"""
Unified avalanche start-zone terrain pipeline.

Usage
-----
    python scripts/pipeline.py <site> [--skip-process] [--skip-analyze]

    python scripts/pipeline.py sisters
    python scripts/pipeline.py star_mtn
    python scripts/pipeline.py sisters --skip-process   # analysis only
    python scripts/pipeline.py star_mtn --skip-analyze  # processing only

Adding a new site
-----------------
Add an entry to SITES below.  Every key is documented inline.
All paths are relative to the project root so the script is
location-independent.
"""
import argparse
import json
import sys
import warnings
from collections import OrderedDict
from glob import glob
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from osgeo import gdal
from osgeo.gdalconst import GA_Update
from scipy import stats
from scipy.ndimage import gaussian_filter
from scipy.stats import kendalltau, kurtosis, pearsonr, spearmanr

import pdal

matplotlib.rcParams['figure.figsize'] = (12, 12)
warnings.filterwarnings('ignore')

ROOT    = Path(__file__).resolve().parent.parent
FIGURES = ROOT / 'figures'

sys.path.insert(0, str(ROOT / 'src'))
from terrain_util import (
    _compute_terrain_attribute,
    crop_tif, get_slope_attributes, load_tif_as_Array,
    plot_3d, plot_attr_vals_probability, show_array,
)
from avl_data import load_avalanches


# =============================================================================
# Site configurations — add new sites here
# =============================================================================

SITES: dict = {
    'sisters': {
        # ── Identity ──────────────────────────────────────────────────────────
        'name':     'Seven Sisters',
        'data_dir': ROOT / 'data' / 'Sisters',

        # ── Processing: LAZ → DEM ────────────────────────────────────────────
        # Glob pattern (relative to data_dir) for raw input tiles
        'laz_glob':   '*.laz',
        # Output file names inside data_dir
        'merged_laz': 'Sisters.laz',
        'ground_laz': 'Sisters_ground.laz',
        # Optional: crop the merged LAZ by EPSG:26913 bounding box before
        # ground filtering.  Set to None to skip.
        'spatial_bounds': dict(x_min=2889242.0, x_max=2891136.0,
                               y_min=1670983.0,  y_max=1672982.0),
        'cropped_laz': 'sisters_crop.laz',   # output name after bounds crop
        # DEM produced from ground-filtered LAZ
        'dem_resolution': 1.2,           # metres per pixel
        'dem_name':       'Sisters_crp.tif',
        # Optional: crop the DEM further by a pixel offset window
        # (offset_x, offset_y, size_x, size_y).  Set to None to skip.
        'pixel_crop_window': None,
        'crp_dem_name':      None,       # only needed when pixel_crop_window set
        # Gaussian gap-fill the DEM in-place after pixel crop?
        'fill_gaps': False,
        # Derive slope/aspect/hillshade/TRI/Roughness rasters from the DEM?
        'compute_terrain_attrs': True,

        # ── Processing: start-zone clipping ──────────────────────────────────
        # KML containing all SZs as one polygon → combined SZ DEM
        # Set to None to skip
        'combined_sz_kml': 'Sisters_SZs.kml',
        'combined_sz_dem': 'Sisters_SZ_DEM.tif',
        # Glob for individual SZ KMLs
        'sz_kml_glob':  'Sister ? SZ.kml',
        # Suffix added to KML stem to produce the output TIF name
        # e.g.  'Sister 1 SZ' + '_DEM' + '.tif' → 'Sister 1 SZ_DEM.tif'
        'sz_dem_suffix': '_DEM',

        # ── Analysis: visualisation ───────────────────────────────────────────
        # Optional satellite image (relative to data_dir)
        'sat_image':  'SatImg.jpg',
        'sat_slice':  (slice(200, None), slice(200, 1500)),
        'sat_labels': [(925, 300, '#1'), (840, 420, '#2'), (770, 420, '#3'),
                       (707, 420, '#4'), (450, 160, '#6')],
        # Labels drawn on the cropped SZ slope map
        'sz_map_labels': [(450, 120, '#1'), (350, 225, '#2'), (280, 225, '#3'),
                          (220, 225, '#4'), (2,   20,  '#6')],
        # richdem attribute name for slope histograms
        'slope_attr': 'slope_degrees',

        # ── Analysis: avalanche data ──────────────────────────────────────────
        # Set avalanche_csv to None to skip avalanche correlation steps
        'avalanche_csv':    ROOT / 'data' / 'avalanches' /
                            'CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv',
        'avalanche_filter': 'Sister',    # substring matched against 'HW Path'
    },

    'star_mtn': {
        # ── Identity ──────────────────────────────────────────────────────────
        'name':     'Star Mountain',
        'data_dir': ROOT / 'data' / 'Star_mtn',

        # ── Processing: LAZ → DEM ────────────────────────────────────────────
        'laz_glob':   'USGS_*.laz',
        'merged_laz': 'Star_mtn.laz',
        'ground_laz': 'Star_mtn_ground.laz',
        'spatial_bounds': None,    # no bounds crop for Star Mtn
        'cropped_laz':    None,
        'dem_resolution': 1.0,
        'dem_name':       'Star_All.tif',
        # Pixel crop applied to the DEM after LAZ conversion
        'pixel_crop_window': (500, 750, 2250, 2050),
        'crp_dem_name':      'Star_mtn_crp.tif',
        'fill_gaps': True,
        'compute_terrain_attrs': False,

        # ── Processing: start-zone clipping ──────────────────────────────────
        'combined_sz_kml': 'Star_All.kml',
        'combined_sz_dem': 'Star_All_SZ.tif',
        'sz_kml_glob':     'Star_?.kml',
        'sz_dem_suffix':   '',   # Star_A.kml → Star_A.tif (no extra suffix)

        # ── Analysis: visualisation ───────────────────────────────────────────
        'sat_image':     None,   # no satellite image for Star Mtn
        'sat_slice':     None,
        'sat_labels':    [],
        'sz_map_labels': [],
        'slope_attr':    'slope_riserun',

        # ── Analysis: avalanche data ──────────────────────────────────────────
        'avalanche_csv':    None,
        'avalanche_filter': None,
    },
}


# =============================================================================
# Shared processing functions
# =============================================================================

def _skip_or(path: Path, action: str) -> bool:
    """Print skip/action message; return True if file already exists."""
    if path.exists():
        print(f"  skip  {path.name}  (already exists)")
        return True
    print(f"  {action:<5} {path.name}")
    return False


def merge_laz_tiles(laz_dir: Path, laz_glob: str, out_file: Path) -> None:
    """Merge LAZ tiles matching laz_glob inside laz_dir into out_file."""
    if _skip_or(out_file, 'merge'):
        return
    pattern  = str(laz_dir / laz_glob)
    pipeline = (
        pdal.Reader.las(filename=pattern)
        | pdal.Filter.merge()
        | pdal.Writer.las(str(out_file))
    )
    pipeline.execute()


def crop_laz_bounds(in_file: Path, out_file: Path, bounds: dict) -> None:
    """Spatially crop a LAZ file to an EPSG:26913 bounding box."""
    if _skip_or(out_file, 'crop'):
        return
    b         = bounds
    bound_str = f"([{b['x_min']},{b['x_max']}],[{b['y_min']},{b['y_max']}])"
    pipeline  = (
        pdal.Reader.las(filename=str(in_file))
        | pdal.Filter.crop(bounds=bound_str)
        | pdal.Writer.las(str(out_file))
    )
    pipeline.execute()


def ground_filter_laz(in_file: Path, out_file: Path) -> None:
    """
    PDAL SMRF ground segmentation:
      reproject (EPSG:26913) → reset classes → ELM noise removal
      → outlier removal → SMRF → extract ground (class 2).
    """
    if _skip_or(out_file, 'smrf'):
        return
    pipeline_def = {
        "pipeline": [
            str(in_file),
            {"type": "filters.reprojection", "out_srs": "EPSG:26913"},
            {"type": "filters.assign",       "assignment": "Classification[:]=0"},
            {"type": "filters.elm"},
            {"type": "filters.outlier"},
            {
                "type":      "filters.smrf",
                "ignore":    "Classification[7:7]",
                "slope":     1.7,
                "window":    16,
                "threshold": 0.45,
                "scalar":    1.2,
            },
            {"type": "filters.range", "limits": "Classification[2:2]"},
            {
                "type":          "writers.las",
                "filename":      str(out_file),
                "minor_version": 1.4,
                "extra_dims":    "all",
            },
        ]
    }
    pdal.Pipeline(json.dumps(pipeline_def)).execute()


def laz_to_dem(in_file: Path, out_file: Path, resolution: float) -> None:
    """Convert ground-filtered LAZ to a mean-elevation DEM GeoTIFF."""
    if _skip_or(out_file, 'dem'):
        return
    pipeline_def = {
        "pipeline": [
            str(in_file),
            {
                "type":        "writers.gdal",
                "filename":    str(out_file),
                "gdaldriver":  "GTiff",
                "output_type": "mean",
                "resolution":  resolution,
            },
        ]
    }
    pdal.Pipeline(json.dumps(pipeline_def)).execute()


def crop_dem_pixels(in_tif: Path, out_tif: Path, window: tuple) -> None:
    """Crop a GeoTIFF by pixel window (offset_x, offset_y, size_x, size_y)."""
    if _skip_or(out_tif, 'crop'):
        return
    gdal.Translate(str(out_tif), str(in_tif), srcWin=list(window))


def fill_dem_gaps(tif_path: Path, sigma: float = 2.0, truncate: float = 4.0) -> None:
    """
    Fill NoData gaps via Gaussian-weighted interpolation (in-place).
    Uses dual convolution to avoid underestimation near gap boundaries.
    """
    print(f"  fill  {tif_path.name}")
    ds  = gdal.Open(str(tif_path), GA_Update)
    dem = np.array(ds.GetRasterBand(1).ReadAsArray(), dtype=float)
    dem[dem == dem.min()] = np.nan
    valid  = np.where(np.isnan(dem), 0.0, dem)
    weight = np.where(np.isnan(dem), 0.0, 1.0)
    filled = (gaussian_filter(valid,  sigma=sigma, truncate=truncate) /
              gaussian_filter(weight, sigma=sigma, truncate=truncate))
    ds.GetRasterBand(1).WriteArray(filled)
    ds.FlushCache()
    ds = None


def compute_terrain_attrs(dem_tif: Path) -> None:
    """Derive slope, aspect, hillshade, TRI, Roughness via GDAL DEMProcessing."""
    ds = gdal.Open(str(dem_tif))
    for attr in ('slope', 'aspect', 'hillshade', 'TRI', 'Roughness'):
        out = dem_tif.parent / f'{dem_tif.stem}_{attr.capitalize()}.tif'
        if _skip_or(out, 'attr'):
            continue
        result = gdal.DEMProcessing(str(out), ds, attr, computeEdges=True)
        result = None
    ds = None


def clip_tif_to_kml(in_tif: Path, out_tif: Path, kml: Path) -> None:
    """Clip a GeoTIFF to the polygon boundary defined in a KML file."""
    if _skip_or(out_tif, 'clip'):
        return
    print(f"        x {kml.name} -> {out_tif.name}")
    ds = gdal.Warp(
        str(out_tif), str(in_tif),
        cutlineDSName=str(kml),
        cropToCutline=True,
        copyMetadata=True,
        dstNodata=np.nan,
    )
    ds = None


# =============================================================================
# Processing pipeline
# =============================================================================

def process_site(cfg: dict) -> None:
    d   = cfg['data_dir']
    print(f"\n{'='*60}")
    print(f"  Processing: {cfg['name']}")
    print(f"{'='*60}\n")

    # 1. Merge raw tiles
    merged_laz = d / cfg['merged_laz']
    print("Step 1 — Merge LAZ tiles")
    merge_laz_tiles(d, cfg['laz_glob'], merged_laz)

    # 2. Optional spatial crop
    laz_for_ground = merged_laz
    if cfg['spatial_bounds']:
        print("\nStep 2 — Spatial bounds crop")
        cropped_laz    = d / cfg['cropped_laz']
        crop_laz_bounds(merged_laz, cropped_laz, cfg['spatial_bounds'])
        laz_for_ground = cropped_laz
    else:
        print("\nStep 2 — Spatial bounds crop  [skipped — not configured]")

    # 3. SMRF ground filter
    ground_laz = d / cfg['ground_laz']
    print("\nStep 3 — SMRF ground segmentation")
    ground_filter_laz(laz_for_ground, ground_laz)

    # 4. LAZ → DEM
    dem_tif = d / cfg['dem_name']
    print("\nStep 4 — LAZ → DEM")
    laz_to_dem(ground_laz, dem_tif, cfg['dem_resolution'])

    # 5. Optional pixel crop
    analysis_dem = dem_tif
    if cfg['pixel_crop_window']:
        print("\nStep 5 — Pixel crop")
        crp_dem     = d / cfg['crp_dem_name']
        crop_dem_pixels(dem_tif, crp_dem, cfg['pixel_crop_window'])
        analysis_dem = crp_dem
    else:
        print("\nStep 5 — Pixel crop  [skipped — not configured]")

    # 6. Optional gap fill
    if cfg['fill_gaps']:
        print("\nStep 6 — Fill DEM gaps")
        fill_dem_gaps(analysis_dem)
    else:
        print("\nStep 6 — Fill DEM gaps  [skipped — not configured]")

    # 7. Terrain attributes
    if cfg['compute_terrain_attrs']:
        print("\nStep 7 — Terrain attributes")
        compute_terrain_attrs(analysis_dem)
    else:
        print("\nStep 7 — Terrain attributes  [skipped — not configured]")

    # 8. Combined SZ DEM
    if cfg['combined_sz_kml']:
        print("\nStep 8 — Clip combined SZ DEM")
        clip_tif_to_kml(analysis_dem,
                         d / cfg['combined_sz_dem'],
                         d / cfg['combined_sz_kml'])
    else:
        print("\nStep 8 — Combined SZ clip  [skipped — not configured]")

    # 9. Individual SZ DEMs
    print("\nStep 9 — Clip individual SZ DEMs")
    for kml in sorted(d.glob(cfg['sz_kml_glob'])):
        out_tif = d / (kml.stem + cfg['sz_dem_suffix'] + '.tif')
        clip_tif_to_kml(analysis_dem, out_tif, kml)

    print(f"\nProcessing done — {d}")


# =============================================================================
# Shared analysis functions
# =============================================================================

def circular_variance(aspects_deg: np.ndarray) -> float:
    r = np.deg2rad(aspects_deg)
    R = np.sqrt(np.nanmean(np.cos(r))**2 + np.nanmean(np.sin(r))**2)
    return float(1 - R)


def circular_std(aspects_deg: np.ndarray) -> float:
    r = np.deg2rad(aspects_deg)
    R = np.sqrt(np.nanmean(np.cos(r))**2 + np.nanmean(np.sin(r))**2)
    return float(np.rad2deg(np.sqrt(-2 * np.log(R)))) if R > 0 else 180.0


def compute_sz_metrics(slope_arr: np.ndarray, aspect_arr: np.ndarray,
                        resolution_m: float = 1.2) -> dict:
    """Terrain variability metrics for one start zone."""
    slopes  = slope_arr.flatten()
    slopes  = slopes[~np.isnan(slopes)]
    slopes  = slopes[(slopes > 0) & (slopes < 90)]
    aspects = aspect_arr.flatten()
    aspects = aspects[~np.isnan(aspects)]
    aspects = aspects[(aspects >= 0) & (aspects <= 360)]
    m = dict(
        slope_mean      = float(np.mean(slopes)),
        slope_std       = float(np.std(slopes)),
        slope_cv        = float(np.std(slopes) / np.mean(slopes)) if np.mean(slopes) else 0.0,
        slope_iqr       = float(np.percentile(slopes, 75) - np.percentile(slopes, 25)),
        slope_kurtosis  = float(stats.kurtosis(slopes)),
        slope_range     = float(np.ptp(slopes)),
        slope_pct_gt_35 = float(np.mean(slopes > 35) * 100),
        slope_pct_gt_40 = float(np.mean(slopes > 40) * 100),
        aspect_circ_var = circular_variance(aspects),
        aspect_circ_std = circular_std(aspects),
        aspect_kurtosis = float(stats.kurtosis(aspects)),
        aspect_n_octants= int(len(np.unique(((aspects + 22.5) % 360 // 45).astype(int)))),
        n_cells         = int(len(slopes)),
        area_m2         = float(len(slopes) * resolution_m**2),
    )
    return m


def correlation_report(x: np.ndarray, y: np.ndarray,
                        x_name: str, y_name: str) -> dict:
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    print(f"\n{'='*55}\n{x_name}  vs  {y_name}  (N={len(x)})\n{'='*55}")
    rho, p_s = spearmanr(x, y)
    tau, p_k = kendalltau(x, y)
    print(f"  Spearman rho = {rho:.3f}  (p = {p_s:.3f})")
    print(f"  Kendall  tau = {tau:.3f}  (p = {p_k:.3f})")
    if len(x) >= 3:
        r, p_p = pearsonr(x, y)
        print(f"  Pearson  r   = {r:.3f}  (p = {p_p:.3f})")
    if len(x) < 10:
        print(f"  WARNING N={len(x)}: p-values unreliable — report effect sizes.")
    return dict(spearman_rho=rho, p_spearman=p_s, kendall_tau=tau, p_kendall=p_k)


def bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                        n_boot: int = 10_000, ci: float = 0.95) -> tuple:
    rng  = np.random.default_rng()
    n    = len(x)
    rhos = np.array([spearmanr(x[i := rng.integers(0, n, n)], y[i])[0]
                     for _ in range(n_boot)])
    lo, hi = np.percentile(rhos, [(1-ci)/2*100, (1-(1-ci)/2)*100])
    obs, _ = spearmanr(x, y)
    print(f"  Spearman rho = {obs:.3f}  [{ci*100:.0f}% CI: {lo:.3f}, {hi:.3f}]")
    return obs, lo, hi


def plot_sz_distributions(sz_paths: list, attrib: str,
                           title: str, out_fig: Path) -> dict:
    """
    Probability histogram + KDE for *attrib* across all SZ DEMs.
    Returns {stem: kurtosis} dict.
    """
    clip_ranges = {
        'slope_riserun': (0.6, 2.0),
        'slope_degrees': (30.0, 65.0),
        'aspect':        (0.0, 360.0),
    }
    lo, hi = clip_ranges.get(attrib, (None, None))

    fig, ax  = plt.subplots()
    kurtosis_d = {}
    for path in sz_paths:
        vals = _compute_terrain_attribute(path, attrib)
        vals = vals[~np.isnan(vals)]

        unique, counts = np.unique(vals, return_counts=True)
        if counts.max() / vals.shape[0] > 0.75:
            vals = vals[vals != unique[counts.argmax()]]

        if lo is not None:
            vals = vals[(vals >= lo) & (vals <= hi)]
        if attrib == 'slope_riserun':
            vals = np.rad2deg(np.arctan(vals))
        elif attrib == 'aspect':
            vals = (vals + 180) % 360

        ku    = round(kurtosis(vals.flatten()), 2)
        stem  = Path(path).stem
        kurtosis_d[stem] = ku
        sns.histplot(vals.flatten(), stat='probability', bins=100, ax=ax,
                     kde=True, line_kws={'lw': 3},
                     label=f"{stem}  (kurtosis: {ku})", alpha=0.4)

    if attrib == 'aspect':
        ax.set_xticklabels(((ax.get_xticks() + 180) % 360).astype(int))

    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(str(out_fig), dpi=150, bbox_inches='tight')
    print(f"  Saved: {out_fig.name}")
    plt.show()
    return kurtosis_d


# =============================================================================
# Analysis pipeline
# =============================================================================

def analyze_site(cfg: dict) -> None:
    d    = cfg['data_dir']
    name = cfg['name']
    print(f"\n{'='*60}")
    print(f"  Analysing: {name}")
    print(f"{'='*60}\n")

    # Determine the primary DEM path
    dem_path = str(d / (cfg['crp_dem_name'] or cfg['dem_name']))

    # ── 1. Satellite image (optional) ────────────────────────────────────────
    if cfg.get('sat_image'):
        print("Step 1 — Satellite image")
        sat = plt.imread(str(d / cfg['sat_image']))
        rs, cs = cfg['sat_slice']
        plt.figure()
        plt.imshow(sat[rs, cs, :])
        for x, y, lbl in cfg['sat_labels']:
            plt.text(x, y, lbl, fontsize=10, color='yellow')
        plt.title(f"{name} — Google Earth view")
        plt.tight_layout()
        plt.show()
    else:
        print("Step 1 — Satellite image  [skipped — not configured]")

    # ── 2. 3-D DEM ───────────────────────────────────────────────────────────
    print("Step 2 — 3-D DEM visualisation")
    import cv2
    dem_arr = load_tif_as_Array(dem_path).astype(float)
    dem_arr[dem_arr == dem_arr.min()] = np.nan
    dem_small = cv2.pyrDown(dem_arr)
    plot_3d(dem_small, title=f"{name} DEM",
            **{'view azi': 130, 'view elev': 10})

    # ── 3. Full-area slope and aspect ─────────────────────────────────────────
    print("Step 3 — Full-area terrain maps")
    full_attrs = get_slope_attributes(dem_path, 'slope_degrees', 'aspect')
    slope_key  = next(k for k in full_attrs if 'slope' in k)
    aspect_key = next(k for k in full_attrs if 'aspect' in k)
    show_array(full_attrs[slope_key][::-1, ::-1],  cmap='magma_r',
               title=f"{name} — Slope (°)")
    show_array(full_attrs[aspect_key][::-1, ::-1], cmap='twilight_shifted',
               title=f"{name} — Aspect")

    # ── 4. Combined SZ DEM slope map ─────────────────────────────────────────
    sz_dem_path = str(d / cfg['combined_sz_dem'])
    if not Path(sz_dem_path).exists() and cfg.get('combined_sz_kml'):
        crop_tif(dem_path, sz_dem_path, str(d / cfg['combined_sz_kml']))

    print("Step 4 — Combined SZ DEM")
    sz_attrs     = get_slope_attributes(sz_dem_path, 'slope_degrees')
    sz_slope_key = list(sz_attrs.keys())[0]
    plt.figure()
    im = plt.imshow(sz_attrs[sz_slope_key][::-1, ::-1], cmap='magma_r')
    plt.colorbar(im, fraction=0.026)
    plt.title(f"{name} — Start zones slope (°)")
    for x, y, lbl in cfg.get('sz_map_labels', []):
        plt.text(x, y, lbl, fontsize=14, color='k')
    plt.tight_layout()
    plt.show()

    # ── 5 & 6. Per-SZ terrain distributions ──────────────────────────────────
    sz_dem_files = sorted(glob(str(d / ('*' + cfg['sz_dem_suffix'] + '.tif')))
                          if cfg['sz_dem_suffix'] else
                          glob(str(d / (cfg['sz_kml_glob'].replace('.kml', '.tif')))))
    # Filter to individual SZ DEMs only (exclude full-area and combined DEMs)
    combined_name = cfg['combined_sz_dem']
    dem_base_name = cfg.get('crp_dem_name') or cfg['dem_name']
    sz_dem_files  = [p for p in sz_dem_files
                     if Path(p).name not in (combined_name, dem_base_name,
                                             cfg['dem_name'])]

    print(f"\nStep 5 — Per-SZ slope distributions  ({len(sz_dem_files)} zones)")
    slope_kurtosis = plot_sz_distributions(
        sz_dem_files,
        attrib=cfg['slope_attr'],
        title=f"{name} — Slope angle distribution",
        out_fig=FIGURES / f"{cfg['name'].lower().replace(' ', '_')}_slope.png",
    )

    print(f"\nStep 6 — Per-SZ aspect distributions")
    aspect_kurtosis = plot_sz_distributions(
        sz_dem_files,
        attrib='aspect',
        title=f"{name} — Aspect distribution",
        out_fig=FIGURES / f"{cfg['name'].lower().replace(' ', '_')}_aspect.png",
    )

    # ── 7–12. Avalanche correlation (Sisters only; skip if no CSV configured) ─
    if not cfg.get('avalanche_csv'):
        print("\nSteps 7-12 — Avalanche analysis  [skipped — no CSV configured]")
        print(f"\nAnalysis done — figures in {FIGURES}")
        return

    DSIZE_ORDER = {f'D{x}': x for x in [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]}

    print("\nStep 7 — CAIC avalanche data")
    avi_df = load_avalanches(str(cfg['avalanche_csv']),
                             cfg['avalanche_filter'],
                             **{'type': ['HS', 'SS'], 'trig': ['N']})
    avi_df['Dsize_numeric'] = avi_df['Dsize'].map(DSIZE_ORDER)
    avi_counts    = avi_df.groupby('HW Path')['Type'].count()
    lg_avi_counts = (avi_df[avi_df['Dsize_numeric'] > 1.0]
                     .groupby('HW Path')['Type'].count())
    total_counts  = pd.concat([avi_counts, lg_avi_counts], axis=1)
    total_counts.columns = ['All HS/SS', '>D1 HS/SS']
    print(total_counts.to_string())

    print("\nStep 8 — Avalanche counts vs kurtosis")
    fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(12, 10))
    total_counts.plot(kind='bar', ax=axes[0], color=['steelblue', 'crimson'])
    axes[0].set_title(f'{name} — Natural HS & SS counts')
    axes[1].bar(avi_counts.index, slope_kurtosis.values(), width=0.5)
    axes[1].axhline(0, color='k')
    axes[1].set_title('Slope kurtosis')
    axes[2].bar(avi_counts.index, aspect_kurtosis.values(), width=0.5)
    axes[2].axhline(0, color='k')
    axes[2].set_title('Aspect kurtosis')
    plt.tight_layout()
    plt.show()

    print("\nStep 9 — Extended terrain metrics")
    sz_metrics = OrderedDict()
    for p in sz_dem_files:
        sid     = Path(p).stem
        s_attrs = get_slope_attributes(p, 'slope_degrees')
        a_attrs = get_slope_attributes(p, 'aspect')
        sz_metrics[sid] = compute_sz_metrics(
            s_attrs[list(s_attrs.keys())[0]],
            a_attrs[list(a_attrs.keys())[0]],
            resolution_m=cfg['dem_resolution'],
        )
    metrics_df = pd.DataFrame(sz_metrics).T
    metrics_df.index.name = 'start_zone'
    print(metrics_df.round(3).to_string())

    analysis_df = metrics_df.copy()
    analysis_df['avi_count']          = total_counts['All HS/SS'].values
    analysis_df['avi_gt_d1']          = total_counts['>D1 HS/SS'].fillna(0).values
    analysis_df['avi_per_area']       = analysis_df['avi_count'] / (analysis_df['area_m2'] / 10_000)
    analysis_df['avi_gt_d1_per_area'] = analysis_df['avi_gt_d1'] / (analysis_df['area_m2'] / 10_000)

    print("\nStep 10 — Nonparametric correlations")
    terrain_vars = [
        ('slope_std',        'Slope SD'),
        ('slope_cv',         'Slope CV'),
        ('slope_iqr',        'Slope IQR'),
        ('slope_kurtosis',   'Slope kurtosis'),
        ('slope_pct_gt_35',  '% cells > 35°'),
        ('aspect_circ_var',  'Aspect circ var'),
        ('aspect_circ_std',  'Aspect circ SD'),
        ('aspect_n_octants', 'N aspect octants'),
    ]
    outcome_vars = [
        ('avi_per_area',       'Avis / 10k m²'),
        ('avi_gt_d1_per_area', '>D1 / 10k m²'),
    ]
    results = []
    for t_col, t_name in terrain_vars:
        for o_col, o_name in outcome_vars:
            r = correlation_report(analysis_df[t_col].values.astype(float),
                                   analysis_df[o_col].values.astype(float),
                                   t_name, o_name)
            r.update(terrain_metric=t_name, outcome=o_name)
            results.append(r)
    corr_df = pd.DataFrame(results)
    print("\nCorrelation summary:")
    print(corr_df[['terrain_metric', 'outcome',
                   'spearman_rho', 'p_spearman']].round(3).to_string(index=False))

    print("\nStep 11 — Scatter plots")
    plot_pairs = [
        ('slope_std',        'avi_per_area',       'Slope SD',        'Avis / 10k m²'),
        ('slope_cv',         'avi_per_area',        'Slope CV',        'Avis / 10k m²'),
        ('slope_pct_gt_35',  'avi_per_area',        '% > 35°',         'Avis / 10k m²'),
        ('aspect_circ_var',  'avi_per_area',        'Aspect circ var', 'Avis / 10k m²'),
        ('aspect_n_octants', 'avi_per_area',        'N octants',       'Avis / 10k m²'),
        ('slope_std',        'avi_gt_d1_per_area',  'Slope SD',        '>D1 / 10k m²'),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f'{name} — Terrain Variability vs Avalanche Frequency', fontsize=13)
    for ax, (xcol, ycol, xlabel, ylabel) in zip(axes.flat, plot_pairs):
        x, y = analysis_df[xcol], analysis_df[ycol]
        ax.scatter(x, y, s=80, zorder=5)
        for idx, row in analysis_df.iterrows():
            ax.annotate(idx, (row[xcol], row[ycol]), fontsize=8, ha='left', va='bottom')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        rho, p = spearmanr(x, y)
        ax.set_title(f'ρ={rho:.2f}, p={p:.2f}', fontsize=10)
    plt.tight_layout()
    scatter_fig = FIGURES / f"{cfg['name'].lower().replace(' ', '_')}_terrain_vs_avis.png"
    plt.savefig(str(scatter_fig), dpi=150, bbox_inches='tight')
    print(f"  Saved: {scatter_fig.name}")
    plt.show()

    print("\nStep 12 — Bootstrap 95% CIs")
    for t_col, t_name in [('slope_std', 'Slope SD'), ('aspect_circ_var', 'Aspect circ var')]:
        for o_col, o_name in [('avi_per_area', 'Avis/area'), ('avi_gt_d1_per_area', '>D1/area')]:
            print(f"\n{t_name}  vs  {o_name}:")
            bootstrap_spearman(analysis_df[t_col].values.astype(float),
                               analysis_df[o_col].values.astype(float))

    print(f"\nAnalysis done — figures in {FIGURES}")


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Avalanche start-zone terrain pipeline.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Available sites: ' + ', '.join(SITES),
    )
    parser.add_argument('--site', choices=list(SITES),
                        help='Site to process/analyse', default='sisters')
    parser.add_argument('--skip-process', action='store_true',
                        help='Skip the LAZ → DEM processing steps')
    parser.add_argument('--skip-analyze', action='store_true',
                        help='Skip the terrain analysis steps')
    args = parser.parse_args()

    FIGURES.mkdir(exist_ok=True)
    cfg = SITES[args.site]

    if not args.skip_process:
        process_site(cfg)
    if not args.skip_analyze:
        analyze_site(cfg)


if __name__ == '__main__':
    main()
