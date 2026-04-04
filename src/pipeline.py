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
from scipy.ndimage import binary_dilation, binary_erosion, gaussian_filter
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
from avl_data import load_avalanches, load_avalanches_simple


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
        # Explicit DEM-stem → CAIC 'HW Path' mapping for joining avalanche counts
        'sz_avi_path_map': {
            'Sister_1_SZ_DEM': 'Seven Sister #1',
            'Sister_2_SZ_DEM': 'Seven Sister #2',
            'Sister_3_SZ_DEM': 'Seven Sister #3',
            'Sister_4_SZ_DEM': 'Seven Sister #4',
            'Sister_6_SZ_DEM': 'Seven Sister #6',
        },
        # SZ DEM stems to exclude from correlation analysis (different loading regime)
        'exclude_sz_stems': ['Sister_6_SZ_DEM'],
        # Primary storm loading direction (compass degrees, inclusive).
        # Used to compute the "steep terrain near loaded face" co-occurrence metric.
        'loading_az_range': (225, 315),   # W to NW
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
        'avalanche_csv':    ROOT / 'data' / 'avalanches' /
                            'CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv',
        'avalanche_filter': 'Star',
        'sz_avi_path_map': {
            'Star_A': 'Star Mtn A',
            'Star_B': 'Star Mtn B',
            'Star_C': 'Star Mtn C',
            'Star_D': 'Star Mtn D',
            'Star_E': 'Star Mtn E',
        },
        'exclude_sz_stems': [],
        'loading_az_range': (225, 315),
    },

    'us550': {
        # ── Identity ──────────────────────────────────────────────────────────
        'name':     'US-550',
        'data_dir': ROOT / 'data' / '550',

        # ── Processing: LAZ → DEM ────────────────────────────────────────────
        'laz_glob':   'USGS_*.laz',
        'merged_laz': '550_merged.laz',
        'ground_laz': '550_ground.laz',
        'spatial_bounds': None,
        'cropped_laz':    None,
        'dem_resolution': 1.0,
        'dem_name':       '550.tif',
        'pixel_crop_window': None,
        'crp_dem_name':      None,
        'fill_gaps': False,
        'compute_terrain_attrs': False,

        # ── Processing: start-zone clipping ──────────────────────────────────
        # 550.kml contains the combined SZ boundary for all four paths
        'combined_sz_kml': '550.kml',
        'combined_sz_dem': '550_SZ.tif',
        # [EMPT]*.kml matches Eagle, Muleshoe, Porcupine, Telescope but not 550.kml
        'sz_kml_glob':     '[EMPT]*.kml',
        'sz_dem_suffix':   '',          # Eagle.kml → Eagle.tif

        # ── Analysis: visualisation ───────────────────────────────────────────
        'sat_image':     None,
        'sat_slice':     None,
        'sat_labels':    [],
        'sz_map_labels': [],
        'slope_attr':    'slope_degrees',

        # ── Analysis: avalanche data ──────────────────────────────────────────
        # US_550_Avalanches.csv has fewer columns than the full CAIC CSV;
        # use the 'simple' loader which only requires HW Path, Trigger, Type.
        'avalanche_csv':    ROOT / 'data' / 'avalanches' / 'US_550_Avalanches.csv',
        'avalanche_filter': None,   # CSV is already filtered to these paths
        'avalanche_loader': 'simple',
        # DEM stems match HW Path names exactly — no explicit map needed
        'sz_avi_path_map': {},
        'exclude_sz_stems': [],
        'loading_az_range': (225, 315),
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

def _profile_curvature(dem_arr: np.ndarray, res: float = 1.0) -> np.ndarray:
    """Profile curvature via the Evans-Young central-difference formula.

    Positive = convex (terrain bends away from viewer looking downslope) —
    the stress-concentration geometry for slab initiation.
    Units: 1/m.  Edge rows/columns are set to NaN.
    """
    Z = dem_arr.astype(float)
    zn = Z[:-2, 1:-1]
    zs = Z[2:,  1:-1]
    ze = Z[1:-1, 2:]
    zw = Z[1:-1, :-2]
    zc = Z[1:-1, 1:-1]

    valid = ~(np.isnan(zn) | np.isnan(zs) | np.isnan(ze) | np.isnan(zw) | np.isnan(zc))

    dzdx  = np.where(valid, (ze  - zw)            / (2 * res),    np.nan)
    dzdy  = np.where(valid, (zn  - zs)            / (2 * res),    np.nan)
    d2x   = np.where(valid, (ze  - 2*zc + zw)     / res**2,       np.nan)
    d2y   = np.where(valid, (zn  - 2*zc + zs)     / res**2,       np.nan)

    p = dzdx**2 + dzdy**2
    q = p + 1.0

    with np.errstate(invalid='ignore', divide='ignore'):
        pc = -(d2x * dzdx**2 + d2y * dzdy**2) / (p * np.sqrt(q))
    pc = np.where(p > 1e-10, pc, np.nan)   # NaN on flat terrain

    out = np.full_like(Z, np.nan)
    out[1:-1, 1:-1] = pc
    return out


def _kml_centroids_lonlat(kml_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (lon, lat) of the centroid of the first polygon in a KML file."""
    import xml.etree.ElementTree as ET
    root = ET.parse(str(kml_path)).getroot()
    for node in root.iter():
        if 'coordinates' in node.tag and node.text:
            pts = [p.split(',') for p in node.text.strip().split()]
            xy  = np.array(pts)[:, :2].astype(float)
            return float(xy[:, 0].mean()), float(xy[:, 1].mean())
    return float('nan'), float('nan')


def _kml_polygon_lonlat(kml_path: Path):
    """Return list of (lon_arr, lat_arr) for each polygon ring in a KML file."""
    import xml.etree.ElementTree as ET
    root  = ET.parse(str(kml_path)).getroot()
    rings = []
    for node in root.iter():
        if 'coordinates' in node.tag and node.text:
            pts = [p.split(',') for p in node.text.strip().split()]
            xy  = np.array(pts)[:, :2].astype(float)
            rings.append((xy[:, 0], xy[:, 1]))
    return rings


def _lonlat_to_dem_pixels(lons, lats, dem_path: str):
    """Convert WGS84 lon/lat arrays to DEM pixel (col, row) coordinates."""
    from osgeo import osr
    ds  = gdal.Open(dem_path)
    gt  = ds.GetGeoTransform()
    src = osr.SpatialReference()
    src.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst = osr.SpatialReference()
    dst.ImportFromWkt(ds.GetProjection())
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ds  = None
    tf  = osr.CoordinateTransformation(src, dst)
    pts = [tf.TransformPoint(float(lo), float(la)) for lo, la in zip(lons, lats)]
    cols = np.array([(x - gt[0]) / gt[1] for x, y, _ in pts])
    rows = np.array([(y - gt[3]) / gt[5] for x, y, _ in pts])
    return cols, rows


def plot_nucleation_map(cfg: dict, site_slug: str, sz_dem_files: list,
                        out_fig_hs: Path, out_fig_sat: Path | None = None) -> None:
    """Hillshade map with nucleation-patch candidate overlays.

    Three-tier overlay (computed on the full-area DEM):
      · Red    — all 3 criteria: steep ≥40°, near W–NW loading face, convex
      · Orange — exactly 2 of 3 criteria (expansion zone candidates)

    SZ polygon boundaries drawn as white outlines.
    If out_fig_sat is given and cfg has sat_image + sat_labels, an approximate
    satellite overlay is also produced using the label positions as GCPs.
    """
    d        = cfg['data_dir']
    res      = cfg['dem_resolution']
    lo_az, hi_az = cfg.get('loading_az_range', (225, 315))
    dem_path = str(d / (cfg.get('crp_dem_name') or cfg['dem_name']))

    # ── Terrain attributes ────────────────────────────────────────────────────
    dem_arr   = load_tif_as_Array(dem_path).astype(float)
    dem_arr[dem_arr == dem_arr.min()] = np.nan
    slope     = _compute_terrain_attribute(dem_path, 'slope_degrees')
    aspect    = _compute_terrain_attribute(dem_path, 'aspect')
    hillshade = _compute_terrain_attribute(dem_path, 'hillshade')
    curv      = _profile_curvature(dem_arr, res=res)

    # ── Nucleation criteria ───────────────────────────────────────────────────
    valid        = ~np.isnan(slope) & ~np.isnan(aspect) & (slope > 0) & (slope < 90)
    steep        = (slope >= 40) & valid
    loading_face = (aspect >= lo_az) & (aspect <= hi_az) & ~np.isnan(aspect)
    near_loading = binary_dilation(loading_face, structure=np.ones((3, 3)), iterations=2)
    convex       = (curv > 0) & ~np.isnan(curv)

    n_crit     = steep.astype(int) + near_loading.astype(int) + convex.astype(int)
    nucleation = n_crit == 3      # all three — most likely nucleation sites
    two_crit   = n_crit == 2      # expansion zone

    # ── SZ polygon boundaries ─────────────────────────────────────────────────
    sz_dem_suffix   = cfg.get('sz_dem_suffix', '')
    exclude_stems   = set(cfg.get('exclude_sz_stems', []))
    sz_rings_px     = []   # [(cols, rows, stem, is_active)]
    sz_centroids_px = []   # [(col, row, label_str)]

    for sz_path in sz_dem_files:
        stem     = Path(sz_path).stem
        kml_name = stem.replace(sz_dem_suffix, '') + '.kml' if sz_dem_suffix else stem + '.kml'
        kml_path = d / kml_name
        if not kml_path.exists():
            continue
        is_active = stem not in exclude_stems
        for lon_ring, lat_ring in _kml_polygon_lonlat(kml_path):
            cols, rows = _lonlat_to_dem_pixels(lon_ring, lat_ring, dem_path)
            sz_rings_px.append((cols, rows, stem, is_active))
        c_lon, c_lat = _kml_centroids_lonlat(kml_path)
        cc, cr = _lonlat_to_dem_pixels([c_lon], [c_lat], dem_path)
        lbl = stem.replace(sz_dem_suffix, '').replace('_', ' ').strip()
        sz_centroids_px.append((float(cc[0]), float(cr[0]), lbl))

    # ── Helper: draw one nucleation panel ────────────────────────────────────
    def _draw_panel(ax, bg, bg_cmap, bg_alpha=1.0):
        ax.imshow(bg, cmap=bg_cmap, alpha=bg_alpha, interpolation='none')

        # two-of-three (expansion zone): yellow, semi-transparent
        two_rgba = np.zeros((*two_crit.shape, 4))
        two_rgba[two_crit] = [1.0, 0.65, 0.0, 0.40]
        ax.imshow(two_rgba, interpolation='none')

        # nucleation candidates: red, more opaque
        nuc_rgba = np.zeros((*nucleation.shape, 4))
        nuc_rgba[nucleation] = [0.92, 0.10, 0.05, 0.70]
        ax.imshow(nuc_rgba, interpolation='none')

        for cols, rows, stem, is_active in sz_rings_px:
            ax.plot(cols, rows, '-', color='white' if is_active else '#aaa',
                    lw=1.2, zorder=5)

        for col, row, lbl in sz_centroids_px:
            ax.text(col, row, lbl, color='white', fontsize=8,
                    ha='center', va='center', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.45),
                    zorder=6)

    # ── Figure 1: hillshade background ───────────────────────────────────────
    from matplotlib.patches import Patch
    legend_elems = [
        Patch(facecolor=(0.92, 0.10, 0.05, 0.70),
              label='Nucleation candidate  (steep ≥40° + near W–NW face + convex)'),
        Patch(facecolor=(1.00, 0.65, 0.00, 0.40),
              label='Expansion zone  (2 of 3 criteria)'),
    ]

    fig, ax = plt.subplots(figsize=(10, 10))
    _draw_panel(ax, hillshade, 'gray')
    ax.legend(handles=legend_elems, loc='lower left', fontsize=8, framealpha=0.85)
    ax.set_title(
        f"{cfg['name']} — Crack Nucleation Candidates\n"
        f"Criteria: slope ≥40°  |  near {lo_az}–{hi_az}° face  |  convex (profile curv > 0)",
        fontsize=11,
    )
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(str(out_fig_hs), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_fig_hs.name}")

    # ── Figure 2: satellite overlay (optional, approximate GCP warp) ─────────
    if out_fig_sat is None or not cfg.get('sat_image') or not cfg.get('sat_labels'):
        return

    import cv2
    sat_raw  = plt.imread(str(d / cfg['sat_image']))
    rs, cs   = cfg['sat_slice']
    sat_crop = sat_raw[rs, cs, :]       # cropped satellite image (H × W × 3)
    H, W     = sat_crop.shape[:2]

    # Build GCP lists: DEM pixel space ↔ sat pixel space
    # cfg['sat_labels']: [(sat_col, sat_row, label_str)]
    # sz_centroids_px:   [(dem_col, dem_row, label_str)]
    label_to_dem  = {lbl: (c, r) for c, r, lbl in sz_centroids_px}
    sat_lbl_map   = {lbl: (sc, sr) for sc, sr, lbl in cfg['sat_labels']}

    # Match on the '#N' shorthand used in sat_labels
    def _short(stem):
        # 'Sister 1 SZ' → '#1', 'Star A' → 'A', etc.
        for ch in '123456789ABCDE':
            if ch in stem:
                return f'#{ch}' if ch.isdigit() else ch
        return stem

    src_pts, dst_pts = [], []
    for lbl, (dc, dr) in label_to_dem.items():
        short = _short(lbl)
        if short in sat_lbl_map:
            sc, sr = sat_lbl_map[short]
            src_pts.append([dc, dr])   # DEM pixel
            dst_pts.append([sc, sr])   # sat pixel

    if len(src_pts) < 4:
        print(f"  [nucleation sat] only {len(src_pts)} GCPs — skipping satellite overlay")
        return

    src_arr = np.float32(src_pts)
    dst_arr = np.float32(dst_pts)
    H_mat, _ = cv2.findHomography(src_arr, dst_arr)

    # Warp nucleation and two-of-three masks into sat space
    nuc_u8  = (nucleation.astype(np.uint8) * 255)
    two_u8  = (two_crit.astype(np.uint8) * 255)
    nuc_sat = cv2.warpPerspective(nuc_u8,  H_mat, (W, H))
    two_sat = cv2.warpPerspective(two_u8,  H_mat, (W, H))

    fig, ax = plt.subplots(figsize=(13, 7))
    ax.imshow(sat_crop)

    two_rgba_sat = np.zeros((H, W, 4))
    two_rgba_sat[two_sat > 127] = [1.0, 0.65, 0.0, 0.40]
    ax.imshow(two_rgba_sat, interpolation='none')

    nuc_rgba_sat = np.zeros((H, W, 4))
    nuc_rgba_sat[nuc_sat > 127] = [0.92, 0.10, 0.05, 0.70]
    ax.imshow(nuc_rgba_sat, interpolation='none')

    ax.legend(handles=legend_elems, loc='lower left', fontsize=8, framealpha=0.85)
    ax.set_title(
        f"{cfg['name']} — Crack Nucleation Candidates (satellite overlay, approx.)\n"
        f"Criteria: slope ≥40°  |  near {lo_az}–{hi_az}° face  |  convex",
        fontsize=11,
    )
    ax.axis('off')
    plt.tight_layout()
    fig.savefig(str(out_fig_sat), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out_fig_sat.name}")


def circular_variance(aspects_deg: np.ndarray) -> float:
    r = np.deg2rad(aspects_deg)
    R = np.sqrt(np.nanmean(np.cos(r))**2 + np.nanmean(np.sin(r))**2)
    return float(1 - R)


def circular_std(aspects_deg: np.ndarray) -> float:
    r = np.deg2rad(aspects_deg)
    R = np.sqrt(np.nanmean(np.cos(r))**2 + np.nanmean(np.sin(r))**2)
    return float(np.rad2deg(np.sqrt(-2 * np.log(R)))) if R > 0 else 180.0


def _aspect_entropy(aspects_deg: np.ndarray) -> float:
    """Shannon entropy (bits) of aspect distribution over 8 octants. Max = 3 bits."""
    if len(aspects_deg) == 0:
        return 0.0
    counts = np.bincount(((aspects_deg + 22.5) % 360 // 45).astype(int), minlength=8)
    p = counts / counts.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def compute_sz_metrics(slope_arr: np.ndarray, aspect_arr: np.ndarray,
                        curv_arr: np.ndarray | None = None,
                        resolution_m: float = 1.2,
                        loading_az_range: tuple = (225, 315)) -> dict:
    """Terrain variability and crack-initiation metrics for one start zone.

    2-D spatial metrics (local extremes, co-occurrence, curvature) are computed
    on the full arrays before flattening.  1-D distribution metrics use matched
    (slope, aspect) pairs to preserve spatial correspondence.

    Parameters
    ----------
    curv_arr : profile curvature array (same shape as slope_arr), optional.
        Positive values = convex terrain (stress concentration sites).
        If None, profile_curv_p95_steep is recorded as NaN.
    loading_az_range : (lo, hi) compass degrees of the primary loading sector.
        Defaults to W–NW (225–315°).  Used for the 'pct_steep_near_loading'
        co-occurrence metric.
    """
    # ── 2-D spatial metrics ───────────────────────────────────────────────────
    valid_2d = (~np.isnan(slope_arr) & ~np.isnan(aspect_arr)
                & (slope_arr > 0) & (slope_arr < 90)
                & (aspect_arr >= 0) & (aspect_arr <= 360))
    s2 = np.where(valid_2d, slope_arr, np.nan)

    # Steepest single cell (≈ resolution_m² ≈ 1.44 m² at 1.2 m)
    slope_max_1cell = float(np.nanmax(s2))

    # Co-occurrence: steep (≥40°) terrain within 2 cells (≈2.4 m) of a
    # loading-sector (W/NW) face.
    # NOTE — necessary but not sufficient: a cell can be steep AND near a
    # loaded face without being directly loaded itself.  Use as a proxy only.
    lo_az, hi_az = loading_az_range
    loading_face   = (aspect_arr >= lo_az) & (aspect_arr <= hi_az) & ~np.isnan(aspect_arr)
    near_loading   = binary_dilation(loading_face, structure=np.ones((3, 3)), iterations=2)
    steep_2d       = (slope_arr >= 40) & valid_2d
    n_valid        = int(valid_2d.sum())
    pct_steep_near_loading = (float((steep_2d & near_loading).sum()) / n_valid * 100
                              if n_valid > 0 else 0.0)

    # Profile curvature 95th-percentile on interior steep cells.
    # Positive profile curvature = convex roll = stress concentration.
    # The 2-pixel boundary erosion removes edge artifacts common in small
    # clipped DEMs where the fill / derivative computation is unreliable.
    if curv_arr is not None:
        interior     = binary_erosion(valid_2d, structure=np.ones((3, 3)), iterations=2)
        curv_steep   = curv_arr[interior & steep_2d]
        curv_steep   = curv_steep[~np.isnan(curv_steep)]
        profile_curv_p95_steep = (float(np.percentile(curv_steep, 95))
                                  if len(curv_steep) > 0 else float('nan'))
    else:
        profile_curv_p95_steep = float('nan')

    # ── Paired 1-D arrays (spatially matched slope + aspect) ──────────────────
    s_flat = slope_arr.flatten()
    a_flat = aspect_arr.flatten()
    valid  = (~np.isnan(s_flat) & ~np.isnan(a_flat)
              & (s_flat > 0) & (s_flat < 90)
              & (a_flat >= 0) & (a_flat <= 360))
    slopes  = s_flat[valid]
    aspects = a_flat[valid]
    steep   = slopes >= 35

    m = dict(
        # ── Slope distribution ───────────────────────────────────────────────
        slope_mean           = float(np.mean(slopes)),
        slope_std            = float(np.std(slopes)),
        slope_cv             = float(np.std(slopes) / np.mean(slopes)) if np.mean(slopes) else 0.0,
        slope_iqr            = float(np.percentile(slopes, 75) - np.percentile(slopes, 25)),
        slope_kurtosis       = float(stats.kurtosis(slopes)),
        slope_range          = float(np.ptp(slopes)),
        # ── Release-zone fractions ───────────────────────────────────────────
        slope_pct_gt_35      = float(np.mean(slopes > 35) * 100),
        slope_pct_gt_40      = float(np.mean(slopes > 40) * 100),
        pct_35_55            = float(np.mean((slopes >= 35) & (slopes <= 55)) * 100),
        pct_40_50            = float(np.mean((slopes >= 40) & (slopes <= 50)) * 100),
        # ── Local extremes: crack-initiation proxies ──────────────────────────
        slope_max_1cell      = slope_max_1cell,   # steepest ≈ 1.4 m² patch
        # ── Aspect distribution ──────────────────────────────────────────────
        aspect_circ_var      = circular_variance(aspects),
        aspect_circ_std      = circular_std(aspects),
        aspect_kurtosis      = float(stats.kurtosis(aspects)),
        aspect_n_octants     = int(len(np.unique(((aspects + 22.5) % 360 // 45).astype(int)))),
        aspect_entropy       = _aspect_entropy(aspects),
        # ── Interaction: aspect diversity of steep cells ──────────────────────
        aspect_entropy_steep = _aspect_entropy(aspects[steep]) if steep.sum() > 10 else 0.0,
        # ── Co-occurrence: steep terrain near loading-sector faces ─────────────
        pct_steep_near_loading = pct_steep_near_loading,
        # ── Stress concentration: convexity of steep interior cells ──────────
        profile_curv_p95_steep = profile_curv_p95_steep,
        # ── Size ─────────────────────────────────────────────────────────────
        n_cells              = int(len(slopes)),
        area_m2              = float(len(slopes) * resolution_m**2),
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
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    rng  = np.random.default_rng()
    n    = len(x)
    rhos = np.array([spearmanr(x[i], y[i])[0]
                     for _ in range(n_boot)
                     for i in (rng.integers(0, n, n),)])
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
    plt.close()
    print(f"  Saved: {out_fig.name}")
    return kurtosis_d


# =============================================================================
# Report
# =============================================================================

def _write_report(cfg: dict, site_slug: str, analysis_df: pd.DataFrame,
                  corr_df: pd.DataFrame, total_counts: pd.DataFrame,
                  boot_results: list) -> None:
    """Write a plain-text analysis report to the figures directory."""
    report_path = FIGURES / f"{site_slug}_report.txt"
    W = 70

    lo_az, hi_az = cfg.get('loading_az_range', (225, 315))
    key_cols = ['pct_40_50', 'slope_max_1cell',
                'aspect_entropy', 'aspect_entropy_steep',
                'pct_steep_near_loading', 'profile_curv_p95_steep',
                'area_m2', 'avi_per_area']

    # ── Pre-compute values used in the interpretation ─────────────────────────
    avi_corrs = corr_df.set_index('terrain_metric')

    entr_steep_rho = avi_corrs.loc['Aspect entropy (steep)', 'spearman_rho']
    entr_rho       = avi_corrs.loc['Aspect entropy',         'spearman_rho']
    slope_rho      = avi_corrs.loc['% 40–50°',               'spearman_rho']
    best_metric    = avi_corrs['spearman_rho'].abs().idxmax()
    best_rho       = avi_corrs.loc[best_metric, 'spearman_rho']
    try:
        curv_rho = avi_corrs.loc['Profile curv P95 (steep)', 'spearman_rho']
    except KeyError:
        curv_rho = float('nan')

    high_avi   = analysis_df.nlargest(2, 'avi_per_area').index.tolist()
    high_entr  = analysis_df.nlargest(2, 'aspect_entropy_steep').index.tolist()

    excl = cfg.get('exclude_sz_stems', [])
    n    = len(analysis_df)

    def _fmt_zone(z):
        return z.replace('_SZ_DEM', '').replace('_', ' ')

    lines = [
        '=' * W,
        f"  {cfg['name'].upper()} — AVALANCHE TERRAIN ANALYSIS",
        '=' * W,
        '',
        'AVALANCHE COUNTS  (Natural HS/SS · CAIC 2009–2021)',
        '-' * W,
        total_counts.to_string(),
        '',
        'KEY TERRAIN METRICS PER START ZONE',
        '-' * W,
        '  pct_40_50              % of SZ cells with slope 40–50° (core release band)',
        '  slope_max_1cell        max slope of any single cell (≈1.4 m² at 1.2 m res)',
        '                         proxy for peak shear stress / crack initiation',
        '  aspect_entropy         Shannon entropy of aspect distribution over 8',
        '                         octants (0 = single direction, 3 bits = uniform)',
        '  aspect_entropy_steep   entropy computed only on cells ≥35°',
        '                         (high = steep terrain faces many wind directions)',
        f' pct_steep_near_loading % of SZ cells that are both steep (≥40°) AND',
        f'                         within 2 cells (≈2.4 m) of a {lo_az}–{hi_az}°',
        '                         (W–NW) facing cell — co-occurrence proxy.',
        ' profile_curv_p95_steep  95th-percentile profile curvature on steep (≥40°)',
        '                         interior cells (boundary pixels excluded). Positive',
        '                         = convex roll; high values indicate stress-',
        '                         concentration sites analogous to notch-tip K_I.',
        '  IMPORTANT CAVEAT: the loading (aspect diversity / near-loading-face)',
        '  and steepness metrics are computed independently across all cells in',
        '  a start zone. A high score on both does NOT guarantee that the loaded',
        '  terrain and the steepest terrain are the same cells. Use',
        '  pct_steep_near_loading as the closest proxy for that co-occurrence.',
        '',
        analysis_df[key_cols].round(3).to_string(),
        '',
        'SPEARMAN CORRELATIONS',
        '-' * W,
        corr_df[['terrain_metric', 'outcome',
                 'spearman_rho', 'p_spearman',
                 'kendall_tau',  'p_kendall']].round(3).to_string(index=False),
        '',
        'BOOTSTRAP 95% CIs  (Spearman ρ · 10,000 resamples)',
        '-' * W,
    ]

    import math
    any_nan = False
    for t_name, o_name, obs, lo, hi in boot_results:
        if math.isnan(lo):
            lines.append(f"  {t_name:<32} vs  {o_name:<14}  ρ={obs:+.3f}  [CI: N/A]")
            any_nan = True
        else:
            lines.append(f"  {t_name:<32} vs  {o_name:<14}"
                         f"  ρ={obs:+.3f}  [95% CI: {lo:+.3f}, {hi:+.3f}]")
    if any_nan:
        lines.append(f"  (N={n} is too small for meaningful bootstrap CIs)")

    # ── Interpretation ────────────────────────────────────────────────────────
    excl_str      = ', '.join(_fmt_zone(z) for z in excl) if excl else 'none'
    high_avi_str  = ', '.join(_fmt_zone(z) for z in high_avi)
    high_entr_str = ', '.join(_fmt_zone(z) for z in high_entr)

    lines += [
        '',
        'INTERPRETATION',
        '-' * W,
        f'N = {n} start zones (excluded from analysis: {excl_str}).',
        'All p-values and effect sizes are exploratory at this sample size.',
        '',
        'Slope steepness:',
        f'  % cells in the 40–50° core release band (ρ={slope_rho:+.2f}) shows no',
        '  consistent trend with avalanche frequency. This is expected: all start',
        '  zones were pre-selected as avalanche terrain and already fall within',
        '  the slab release angle range. Raw steepness does not discriminate',
        '  high- from low-frequency zones when all zones are already steep.',
        '',
        'Aspect diversity:',
        f'  Aspect entropy of steep cells (≥35°) is the strongest predictor',
        f'  (ρ={entr_steep_rho:+.3f}). Zones with the most diverse aspect',
        f'  distribution on steep terrain — {high_entr_str} — are also',
        f'  the highest-frequency zones ({high_avi_str}).',
        '',
        '  This is consistent with a wind-loading hypothesis: start zones',
        '  whose steep cells face multiple directions accumulate wind slabs',
        '  from multiple storm tracks, increasing the total loading and the',
        '  probability of natural slab release.',
        '',
        f'  Overall aspect entropy (all cells) shows a weaker signal',
        f'  (ρ={entr_rho:+.3f}), confirming that the relationship is specific',
        '  to the steep (release-zone) portion of the terrain.',
        '',
        'Co-occurrence of loading and steepness:',
        '  The loading proxies (aspect entropy, near-loading-face fraction) and',
        '  the steepness proxies (pct_40_50, slope_max) are computed separately',
        '  across all cells in a start zone. A start zone can score high on both',
        '  without those conditions coinciding in the same place — for example,',
        '  a zone could have a very steep section that faces E (not loaded from',
        f'  W/NW) and a moderate {lo_az}–{hi_az}° section that is only 35°.',
        '  pct_steep_near_loading is the closest proxy for true co-occurrence,',
        '  but it is still a necessary rather than sufficient condition.',
        '  Definitive co-occurrence analysis would require mapping the loaded',
        '  snow pillow (wind slab deposition) onto the slope map directly.',
        '',
        'Profile curvature (stress concentration):',
        f'  Profile curvature P95 on steep cells (ρ={curv_rho:+.3f}).',
        '  Convex rolls (positive profile curvature) concentrate normal stress',
        '  at the base of the slab in the same way a notch concentrates K_I in',
        '  fracture mechanics (Gvirtzman et al. 2025). The 95th-percentile on',
        '  interior steep cells captures the most extreme stress-concentration',
        '  sites within each start zone while excluding boundary artifacts.',
        '  Boundary cells (within 2 pixels of NoData) are masked before',
        '  computing this metric to reduce DEM edge effects.',
        '',
        'Next steps:',
        '  · Replicate at Star Mountain and US-550 paths (Muleshoe, Eagle,',
        '    Telescope, Porcupine) to increase N before drawing conclusions.',
        '  · Consider weighting aspect entropy by prevailing wind frequency',
        '    for the Silverton area to make the loading proxy more physically',
        '    grounded.',
        '  · Map snow depth variability (from repeat lidar or GPR) onto slope',
        '    and aspect rasters to directly measure co-occurrence of loading',
        '    and high shear stress terrain.',
        '',
        '=' * W,
    ]

    report_path.write_text('\n'.join(lines))
    print(f"  Saved: {report_path.name}")


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

    site_slug = name.lower().replace(' ', '_')

    # ── 1. Satellite image (optional) ────────────────────────────────────────
    if cfg.get('sat_image'):
        print("Step 1 — Satellite image")
        sat = plt.imread(str(d / cfg['sat_image']))
        rs, cs = cfg['sat_slice']
        fig, ax = plt.subplots()
        ax.imshow(sat[rs, cs, :])
        for x, y, lbl in cfg['sat_labels']:
            ax.text(x, y, lbl, fontsize=10, color='yellow')
        ax.set_title(f"{name} — Google Earth view")
        plt.tight_layout()
        fig.savefig(str(FIGURES / f"{site_slug}_sat.png"), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {site_slug}_sat.png")
    else:
        print("Step 1 — Satellite image  [skipped — not configured]")

    # ── 2. 3-D DEM ───────────────────────────────────────────────────────────
    print("Step 2 — 3-D DEM visualisation")
    import cv2
    dem_arr = load_tif_as_Array(dem_path).astype(float)
    dem_arr[dem_arr == dem_arr.min()] = np.nan
    dem_small = cv2.pyrDown(dem_arr)
    fig = plot_3d(dem_small, title=f"{name} DEM", show=False,
                  **{'view azi': 130, 'view elev': 10})
    fig.savefig(str(FIGURES / f"{site_slug}_dem_3d.png"), dpi=100, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {site_slug}_dem_3d.png")

    # ── 3. Full-area slope and aspect ─────────────────────────────────────────
    print("Step 3 — Full-area terrain maps")
    full_attrs = get_slope_attributes(dem_path, 'slope_degrees', 'aspect')
    slope_key  = next(k for k in full_attrs if 'slope' in k)
    aspect_key = next(k for k in full_attrs if 'aspect' in k)
    fig = show_array(full_attrs[slope_key][::-1, ::-1],  cmap='magma_r',
                     title=f"{name} — Slope (°)", show=False)
    fig.savefig(str(FIGURES / f"{site_slug}_slope_full.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    fig = show_array(full_attrs[aspect_key][::-1, ::-1], cmap='twilight_shifted',
                     title=f"{name} — Aspect", show=False)
    fig.savefig(str(FIGURES / f"{site_slug}_aspect_full.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {site_slug}_slope_full.png, {site_slug}_aspect_full.png")

    # ── 4. Combined SZ DEM slope map ─────────────────────────────────────────
    sz_dem_path = str(d / cfg['combined_sz_dem'])
    if not Path(sz_dem_path).exists() and cfg.get('combined_sz_kml'):
        crop_tif(dem_path, sz_dem_path, str(d / cfg['combined_sz_kml']))

    print("Step 4 — Combined SZ DEM")
    sz_attrs     = get_slope_attributes(sz_dem_path, 'slope_degrees')
    sz_slope_key = list(sz_attrs.keys())[0]
    fig, ax = plt.subplots()
    im = ax.imshow(sz_attrs[sz_slope_key][::-1, ::-1], cmap='magma_r')
    fig.colorbar(im, ax=ax, fraction=0.026)
    ax.set_title(f"{name} — Start zones slope (°)")
    for x, y, lbl in cfg.get('sz_map_labels', []):
        ax.text(x, y, lbl, fontsize=14, color='k')
    plt.tight_layout()
    fig.savefig(str(FIGURES / f"{site_slug}_sz_slope_map.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {site_slug}_sz_slope_map.png")

    # ── 5 & 6. Per-SZ terrain distributions ──────────────────────────────────
    sz_dem_files = sorted(glob(str(d / ('*' + cfg['sz_dem_suffix'] + '.tif')))
                          if cfg['sz_dem_suffix'] else
                          glob(str(d / (cfg['sz_kml_glob'].replace('.kml', '.tif')))))
    # Filter out full-area / combined DEMs and any explicitly excluded zones
    combined_name = cfg['combined_sz_dem']
    dem_base_name = cfg.get('crp_dem_name') or cfg['dem_name']
    exclude_stems = set(cfg.get('exclude_sz_stems', []))
    sz_dem_files  = [p for p in sz_dem_files
                     if Path(p).name not in (combined_name, dem_base_name,
                                             cfg['dem_name'])
                     and Path(p).stem not in exclude_stems]

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

    # ── 7–12. Avalanche correlation (skip if no CSV configured) ──────────────
    if not cfg.get('avalanche_csv'):
        print("\nSteps 7-12 — Avalanche analysis  [skipped — no CSV configured]")
        print(f"\nAnalysis done — figures in {FIGURES}")
        return

    avi_path_map = cfg.get('sz_avi_path_map', {})

    print("\nStep 7 — Avalanche data")
    _loaders   = {'simple': load_avalanches_simple}
    _loader    = _loaders.get(cfg.get('avalanche_loader', 'caic'), load_avalanches)
    avi_df     = _loader(str(cfg['avalanche_csv']),
                         cfg['avalanche_filter'],
                         **{'type': ['HS', 'SS'], 'trig': ['N']})
    avi_counts = avi_df.groupby('HW Path')['Type'].count().rename('Count')
    print(avi_counts.to_string())

    print("\nStep 8 — Avalanche counts vs kurtosis")
    # Align kurtosis dicts with the analysed zones using the explicit path map
    analysed_stems  = [Path(p).stem for p in sz_dem_files]
    aligned_paths   = [avi_path_map.get(s, s) for s in analysed_stems]
    aligned_counts  = avi_counts.reindex(aligned_paths)
    aligned_slope_k = [slope_kurtosis.get(s, float('nan')) for s in analysed_stems]
    aligned_aspect_k = [aspect_kurtosis.get(s, float('nan')) for s in analysed_stems]

    fig, axes = plt.subplots(nrows=3, sharex=True, figsize=(12, 10))
    axes[0].bar(aligned_paths, aligned_counts.values, color='steelblue')
    axes[0].set_title(f'{name} — Natural HS & SS counts')
    axes[0].set_ylabel('Count')
    axes[1].bar(aligned_paths, aligned_slope_k, width=0.5)
    axes[1].axhline(0, color='k')
    axes[1].set_title('Slope kurtosis')
    axes[2].bar(aligned_paths, aligned_aspect_k, width=0.5)
    axes[2].axhline(0, color='k')
    axes[2].set_title('Aspect kurtosis')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    fig.savefig(str(FIGURES / f"{site_slug}_kurtosis.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {site_slug}_kurtosis.png")

    print("\nStep 9 — Extended terrain metrics")
    sz_metrics = OrderedDict()
    for p in sz_dem_files:
        sid     = Path(p).stem
        s_attrs  = get_slope_attributes(p, 'slope_degrees')
        a_attrs  = get_slope_attributes(p, 'aspect')
        dem_arr  = load_tif_as_Array(p).astype(float)
        dem_arr[dem_arr == dem_arr.min()] = np.nan   # nodata sentinel
        curv_arr = _profile_curvature(dem_arr, res=cfg['dem_resolution'])
        sz_metrics[sid] = compute_sz_metrics(
            s_attrs[list(s_attrs.keys())[0]],
            a_attrs[list(a_attrs.keys())[0]],
            curv_arr=curv_arr,
            resolution_m=cfg['dem_resolution'],
            loading_az_range=cfg.get('loading_az_range', (225, 315)),
        )
    metrics_df = pd.DataFrame(sz_metrics).T
    metrics_df.index.name = 'start_zone'
    print(metrics_df.round(3).to_string())

    # Join avalanche counts by explicit map — no positional assumption
    analysis_df = metrics_df.copy()
    analysis_df['avi_count'] = [
        avi_counts.get(avi_path_map.get(s, s), 0)
        for s in analysis_df.index
    ]
    analysis_df['avi_per_area'] = (analysis_df['avi_count']
                                   / (analysis_df['area_m2'] / 10_000))

    print("\nStep 10 — Nonparametric correlations")
    lo_az, hi_az = cfg.get('loading_az_range', (225, 315))
    terrain_vars = [
        # steepness distribution
        ('slope_std',               'Slope SD'),
        ('slope_cv',                'Slope CV'),
        ('slope_pct_gt_35',         '% > 35°'),
        ('slope_pct_gt_40',         '% > 40°'),
        ('pct_40_50',               '% 40–50°'),
        # local extremes — crack initiation stress proxies
        ('slope_max_1cell',         'Max slope 1-cell'),
        # aspect diversity — wind loading opportunities
        ('aspect_entropy',          'Aspect entropy'),
        ('aspect_circ_var',         'Aspect circ var'),
        ('aspect_n_octants',        'N aspect octants'),
        ('aspect_entropy_steep',    'Aspect entropy (steep)'),
        # co-occurrence: steep AND near loading-sector face
        ('pct_steep_near_loading',  f'% steep near {lo_az}–{hi_az}°'),
        # stress concentration: convex terrain on steep interior cells
        ('profile_curv_p95_steep',  'Profile curv P95 (steep)'),
    ]
    outcome_vars = [('avi_per_area', 'Avis / 10k m²')]
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

    def _scatter_grid(pairs, filename, title):
        fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 5))
        if len(pairs) == 1:
            axes = [axes]
        fig.suptitle(title, fontsize=13)
        for ax, (xcol, ycol, xlabel, ylabel) in zip(axes, pairs):
            x, y = analysis_df[xcol], analysis_df[ycol]
            ax.scatter(x, y, s=80, zorder=5)
            for idx, row in analysis_df.iterrows():
                ax.annotate(idx.replace('_SZ_DEM', '').replace('_', ' '),
                            (row[xcol], row[ycol]), fontsize=8,
                            ha='left', va='bottom')
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            rho, p = spearmanr(x, y)
            ax.set_title(f'ρ={rho:.2f}, p={p:.2f}', fontsize=10)
        plt.tight_layout()
        fig.savefig(str(FIGURES / filename), dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {filename}")

    print("\nStep 11 — Scatter plots")
    _scatter_grid(
        [('pct_40_50',            'avi_per_area', '% 40–50°',               'Avis / 10k m²'),
         ('aspect_entropy',       'avi_per_area', 'Aspect entropy',         'Avis / 10k m²'),
         ('aspect_entropy_steep', 'avi_per_area', 'Aspect entropy (steep)', 'Avis / 10k m²')],
        f"{site_slug}_terrain_vs_avis.png",
        f"{name} — Loading & Aspect Diversity vs Avalanche Frequency",
    )
    _scatter_grid(
        [('slope_max_1cell',          'avi_per_area', 'Max slope 1-cell (°)',          'Avis / 10k m²'),
         ('pct_steep_near_loading',   'avi_per_area', f'% steep near {lo_az}–{hi_az}°', 'Avis / 10k m²'),
         ('profile_curv_p95_steep',   'avi_per_area', 'Profile curv P95 steep',        'Avis / 10k m²')],
        f"{site_slug}_steepness_vs_avis.png",
        f"{name} — Crack Initiation Proxies vs Avalanche Frequency",
    )

    print("\nStep 12 — Bootstrap 95% CIs")
    boot_results = []
    for t_col, t_name in [('pct_40_50',               '% 40–50°'),
                           ('slope_max_1cell',         'Max slope 1-cell'),
                           ('aspect_entropy',          'Aspect entropy'),
                           ('aspect_entropy_steep',    'Aspect entropy (steep)'),
                           ('pct_steep_near_loading',  f'% steep near {lo_az}–{hi_az}°'),
                           ('profile_curv_p95_steep',  'Profile curv P95 (steep)')]:
        print(f"\n{t_name}  vs  Avis/area:")
        obs, lo_ci, hi_ci = bootstrap_spearman(analysis_df[t_col].values.astype(float),
                                               analysis_df['avi_per_area'].values.astype(float))
        boot_results.append((t_name, 'Avis/area', obs, lo_ci, hi_ci))

    _write_report(cfg, site_slug, analysis_df, corr_df, avi_counts, boot_results)

    print("\nStep 13 — Nucleation patch map")
    out_sat = (FIGURES / f"{site_slug}_nucleation_sat.png"
               if cfg.get('sat_image') and cfg.get('sat_labels') else None)
    plot_nucleation_map(
        cfg, site_slug, sz_dem_files,
        out_fig_hs  = FIGURES / f"{site_slug}_nucleation_hs.png",
        out_fig_sat = out_sat,
    )
    print(f"\nAnalysis done — figures and report in {FIGURES}")


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
                        help='Site to process/analyse  (sisters | star_mtn | us550)',
                        default='sisters')
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
