"""
Unified avalanche start-zone terrain pipeline.

Usage
-----
    python src/pipeline.py [--site SITE] [--skip-process] [--skip-analyze]
                           [--force | --force-from {laz,dem,sz}]

    python src/pipeline.py --site sisters
    python src/pipeline.py --site star_mtn
    python src/pipeline.py --site sisters --skip-process   # analysis only
    python src/pipeline.py --site star_mtn --skip-analyze  # processing only
    python src/pipeline.py --site sisters --force-from sz  # re-clip all SZs
    python src/pipeline.py --site sisters --force          # rebuild everything

Rebuilds
--------
Each processing step is skipped when its output exists AND is newer than
all of its inputs (make-style).  Editing a start-zone KML therefore
re-clips just that SZ; a rebuilt DEM re-clips everything downstream.
Stages, in order:  laz (merge/crop/SMRF) → dem (DEM, pixel crop, gap fill,
terrain attributes) → sz (start-zone clips).  --force-from STAGE rebuilds
that stage and every later one regardless of timestamps; --force is
--force-from laz.

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
from osgeo import gdal, osr
from osgeo.gdalconst import GA_Update
from scipy import stats
from scipy.ndimage import binary_dilation, binary_erosion, gaussian_filter
from scipy.stats import kendalltau, kurtosis, pearsonr, spearmanr

import pdal

matplotlib.rcParams['figure.figsize'] = (12, 12)
warnings.filterwarnings('ignore')
# Raise on GDAL errors instead of returning None — otherwise a failed
# Warp/Translate looks like success and leaves no output behind.
gdal.UseExceptions()

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
        # Optional: pre-made elevation raster (e.g. photogrammetry DTM) to
        # use instead of the LAZ chain.  It is block-averaged to
        # dem_resolution and written as dem_name; LAZ steps are then skipped.
        # Must be BARE EARTH — a DSM with canopy biases slope/aspect/curvature.
        # USGS 3DEP 1 m bare-earth DEM (DRCOG 2020 lidar).  The drone DSM
        # (260922_Sisters_Full_Ground-DSM.tiff) includes forest canopy.
        'source_raster':  'USGS_1M_13_x42y440_CO_DRCOG_2020_B20 (1).tif',
        # Optional crop of source_raster (xmin, ymin, xmax, ymax in its CRS);
        # here the drone-survey footprint.  None = whole raster.
        'source_bounds':  (423279.0, 4391930.0, 424179.0, 4392630.0),
        # Resampling for source_raster → dem_resolution.  'average' suits an
        # integer ratio (5 cm → 1.2 m); for 1.0 → 1.2 m it leaves a 5-pixel
        # ripple in curvature (10× spectral spike), which 'cubic' avoids.
        'source_resample': 'cubic',
        # DEM produced from ground-filtered LAZ (or resampled source_raster)
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
        'sz_kml_glob':  'Sister_?_SZ.kml',
        # Suffix added to KML stem to produce the output TIF name
        # e.g.  'Sister_1_SZ' + '_DEM' + '.tif' → 'Sister_1_SZ_DEM.tif'
        'sz_dem_suffix': '_DEM',

        # ── Analysis: visualisation ───────────────────────────────────────────
        # Optional satellite image (relative to data_dir)
        'sat_image':  'SatImg.jpg',
        # Georeferenced RGB(+NIR) image used as the background of map figures
        # (USGS NAIP 0.6 m, acquired 2019-09-22, from the National Map
        # USGSNAIPImagery service).  None = 3DEP hillshade.
        'basemap':       'NAIP_sisters.tif',
        'basemap_label': 'NAIP 2019-09-22',
        'sat_slice':  (slice(200, None), slice(200, 1500)),
        'sat_labels': [(925, 300, '#1'), (840, 420, '#2'), (770, 420, '#3'),
                       (707, 420, '#4'), (450, 160, '#6')],
        # CAIC database stations {sql_id: [columns]} for src/sql_util.py.
        # Raw CAIC units: tenths of °F and mph; SNOTEL SWE in HUNDREDTHS of
        # an inch, depth in tenths (2000 = flag value 2010–14); pcpac changes
        # scale in 2016 — use SWE for snowfall.  A column may come from
        # only one station.
        #   CALVP — CAIC/CDOT Loveland Pass ridge wind station (wind, temp)
        #   LBAC2 — Loveland Basin SNOTEL (snow depth, SWE, precipitation;
        #           pcpac = accumulated to 2019-09, pcp1 = hourly from 2019-10)
        'weather_stations': {
            'CALVP': ['wspd', 'wdir', 'gust', 'temp', 'rh'],
            'LBAC2': ['depth', 'snowwater', 'pcpac', 'pcp1'],
        },
        # Labels drawn on the cropped SZ slope map: 'auto' places '#N' at each
        # SZ polygon centroid; or a list of (x, y, text) in plot pixels
        'sz_map_labels': 'auto',
        # richdem attribute name for slope histograms
        'slope_attr': 'slope_degrees',

        # ── Analysis: avalanche data ──────────────────────────────────────────
        # Set avalanche_csv to None to skip avalanche correlation steps.
        # Merged file = the named CAIC highway export (2009 – 2021-05-04)
        # + Sisters natural slabs from the 2026 CAIC export (2020 – 2026-03),
        # matched to paths by coordinates (src/assign_avalanche_paths.py,
        # src/merge_avalanche_records.py).  The 2026 export has naturals
        # only, so control-mission (explosive) records end on 2021-05-04.
        'avalanche_csv':    ROOT / 'data' / 'avalanches' /
                            'CAIC_HWY_avalanches_2009-2026_sisters_merged.csv',
        'control_records_end': '2021-05-04',
        'avalanche_filter': 'Sister',    # substring matched against 'HW Path'
        # Explicit DEM-stem → CAIC 'HW Path' mapping for joining avalanche counts
        'sz_avi_path_map': {
            'Sister_1_SZ_DEM': 'Seven Sister #1',
            'Sister_2_SZ_DEM': 'Seven Sister #2',
            'Sister_3_SZ_DEM': 'Seven Sister #3',
            'Sister_4_SZ_DEM': 'Seven Sister #4',
            'Sister_6_SZ_DEM': 'Seven Sister #6',
        },
        # SZ DEM stems to exclude from correlation analysis.  Sister #6 was
        # excluded until its start zone was redrawn (Nov 2026); all five
        # active paths are now analysed.
        'exclude_sz_stems': [],
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
        # Star A & B have a different loading pattern and are not comparable
        # to C, D & E — excluded from correlation analysis.
        'exclude_sz_stems': ['Star_A', 'Star_B'],
        'loading_az_range': (225, 315),
    },

    # 'us550': {
    #     # ── Identity ──────────────────────────────────────────────────────────
    #     'name':     'US-550',
    #     'data_dir': ROOT / 'data' / '550',

    #     # ── Processing: LAZ → DEM ────────────────────────────────────────────
    #     'laz_glob':   'USGS_*.laz',
    #     'merged_laz': '550_merged.laz',
    #     'ground_laz': '550_ground.laz',
    #     'spatial_bounds': None,
    #     'cropped_laz':    None,
    #     'dem_resolution': 1.0,
    #     'dem_name':       '550.tif',
    #     'pixel_crop_window': None,
    #     'crp_dem_name':      None,
    #     'fill_gaps': False,
    #     'compute_terrain_attrs': False,

    #     # ── Processing: start-zone clipping ──────────────────────────────────
    #     # 550.kml contains the combined SZ boundary for all four paths
    #     'combined_sz_kml': '550.kml',
    #     'combined_sz_dem': '550_SZ.tif',
    #     # [EMPT]*.kml matches Eagle, Muleshoe, Porcupine, Telescope but not 550.kml
    #     'sz_kml_glob':     '[EMPT]*.kml',
    #     'sz_dem_suffix':   '',          # Eagle.kml → Eagle.tif

    #     # ── Analysis: visualisation ───────────────────────────────────────────
    #     'sat_image':     None,
    #     'sat_slice':     None,
    #     'sat_labels':    [],
    #     'sz_map_labels': [],
    #     'slope_attr':    'slope_degrees',

    #     # ── Analysis: avalanche data ──────────────────────────────────────────
    #     # US_550_Avalanches.csv has fewer columns than the full CAIC CSV;
    #     # use the 'simple' loader which only requires HW Path, Trigger, Type.
    #     'avalanche_csv':    ROOT / 'data' / 'avalanches' / 'US_550_Avalanches.csv',
    #     'avalanche_filter': None,   # CSV is already filtered to these paths
    #     'avalanche_loader': 'simple',
    #     # DEM stems match HW Path names exactly — no explicit map needed
    #     'sz_avi_path_map': {},
    #     'exclude_sz_stems': [],
    #     'loading_az_range': (225, 315),
    # },
}


# =============================================================================
# Shared processing functions
# =============================================================================

# Processing stages, in pipeline order.  --force-from STAGE rebuilds that
# stage and every later one.
STAGES = ('laz', 'dem', 'sz')
_FORCE_FROM: str | None = None


def _forced(stage: str | None) -> bool:
    return (stage is not None and _FORCE_FROM is not None
            and STAGES.index(stage) >= STAGES.index(_FORCE_FROM))


def _stale_inputs(out: Path, inputs) -> list[Path]:
    """Inputs that exist and were modified after `out`."""
    out_mtime = out.stat().st_mtime
    return [p for p in inputs if p.exists() and p.stat().st_mtime > out_mtime]


def _fmt_names(paths: list[Path], n: int = 3) -> str:
    names = ', '.join(p.name for p in paths[:n])
    return names + (f' +{len(paths) - n} more' if len(paths) > n else '')


def _skip_or(path: Path, action: str, inputs=(), stage: str | None = None) -> bool:
    """
    Return True (and print 'skip') when `path` exists, its stage is not
    forced, and it is newer than every file in `inputs`.  Otherwise print
    why it is being (re)built, delete any stale copy, and return False.
    """
    if not path.exists():
        print(f"  {action:<5} {path.name}")
        return False
    if _forced(stage):
        reason = f"forced from '{_FORCE_FROM}'"
    else:
        newer = _stale_inputs(path, inputs)
        if not newer:
            print(f"  skip  {path.name}  (up to date)")
            return True
        reason = f"newer input: {_fmt_names(newer)}"
    print(f"  {action:<5} {path.name}  (rebuild — {reason})")
    _unlink_with_sidecar(path)
    return False


def _unlink_with_sidecar(path: Path) -> None:
    """Delete a raster and its .aux.xml (cached stats would describe the old data)."""
    path.unlink()
    Path(str(path) + '.aux.xml').unlink(missing_ok=True)


def merge_laz_tiles(laz_dir: Path, laz_glob: str, out_file: Path) -> bool:
    """Merge LAZ tiles matching laz_glob inside laz_dir into out_file."""
    if _skip_or(out_file, 'merge', sorted(laz_dir.glob(laz_glob)), 'laz'):
        return False
    pattern  = str(laz_dir / laz_glob)
    pipeline = (
        pdal.Reader.las(filename=pattern)
        | pdal.Filter.merge()
        | pdal.Writer.las(str(out_file))
    )
    pipeline.execute()
    return True


def crop_laz_bounds(in_file: Path, out_file: Path, bounds: dict) -> bool:
    """Spatially crop a LAZ file to an EPSG:26913 bounding box."""
    if _skip_or(out_file, 'crop', [in_file], 'laz'):
        return False
    b         = bounds
    bound_str = f"([{b['x_min']},{b['x_max']}],[{b['y_min']},{b['y_max']}])"
    pipeline  = (
        pdal.Reader.las(filename=str(in_file))
        | pdal.Filter.crop(bounds=bound_str)
        | pdal.Writer.las(str(out_file))
    )
    pipeline.execute()
    return True


def ground_filter_laz(in_file: Path, out_file: Path) -> bool:
    """
    PDAL SMRF ground segmentation:
      reproject (EPSG:26913) → reset classes → ELM noise removal
      → outlier removal → SMRF → extract ground (class 2).
    """
    if _skip_or(out_file, 'smrf', [in_file], 'laz'):
        return False
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
                "minor_version": 4,
                "extra_dims":    "all",
            },
        ]
    }
    pdal.Pipeline(json.dumps(pipeline_def)).execute()
    return True


def laz_to_dem(in_file: Path, out_file: Path, resolution: float) -> bool:
    """Convert ground-filtered LAZ to a mean-elevation DEM GeoTIFF."""
    if _skip_or(out_file, 'dem', [in_file], 'dem'):
        return False
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
    return True


def resample_source_raster(src: Path, out_tif: Path, resolution: float,
                           bounds: tuple | None = None,
                           alg: str = 'average') -> bool:
    """
    Block-average a (finer) elevation raster to `resolution` metres,
    optionally cropped to `bounds` (xmin, ymin, xmax, ymax in src CRS).
    Reads full-resolution data (overviews ignored — they may be built with a
    different resampling), snaps to a resolution-aligned grid, NoData -9999.
    """
    # A different source, crop or resolution must always rebuild — a
    # swapped-in raster can carry an older mtime than the existing DEM and
    # fool the timestamp check.
    recipe = f"{src.name}|{bounds}|{resolution}|{alg}"
    if out_tif.exists():
        ds = gdal.Open(str(out_tif))
        prev = ds.GetMetadataItem('SOURCE_RECIPE')
        ds = None
        if prev != recipe:
            print(f"  warp  {out_tif.name}  (rebuild — source changed: "
                  f"{prev or 'unknown'} → {recipe})")
            _unlink_with_sidecar(out_tif)
    if _skip_or(out_tif, 'warp', [src], 'dem'):
        return False
    ds = gdal.Warp(str(out_tif), str(src),
                   outputBounds=bounds,
                   xRes=resolution, yRes=resolution, targetAlignedPixels=True,
                   resampleAlg=alg, overviewLevel='NONE', dstNodata=-9999,
                   multithread=True, warpMemoryLimit=4096,
                   creationOptions=['COMPRESS=DEFLATE', 'TILED=YES'])
    ds.SetMetadataItem('SOURCE_RECIPE', recipe)
    ds = None
    return True


def crop_dem_pixels(in_tif: Path, out_tif: Path, window: tuple) -> bool:
    """Crop a GeoTIFF by pixel window (offset_x, offset_y, size_x, size_y)."""
    if _skip_or(out_tif, 'crop', [in_tif], 'dem'):
        return False
    gdal.Translate(str(out_tif), str(in_tif), srcWin=list(window))
    return True


def fill_dem_gaps(tif_path: Path, sigma: float = 2.0, truncate: float = 4.0) -> None:
    """
    Fill NoData gaps via Gaussian-weighted interpolation (in-place).
    Uses dual convolution to avoid underestimation near gap boundaries.

    Only call this on a freshly built DEM: process_site() runs it solely in
    the same run that (re)created the file, so it is never applied twice.
    """
    print(f"  fill  {tif_path.name}")
    ds   = gdal.Open(str(tif_path), GA_Update)
    band = ds.GetRasterBand(1)
    dem  = np.array(band.ReadAsArray(), dtype=float)
    nodata = band.GetNoDataValue()
    dem[dem == (nodata if nodata is not None else dem.min())] = np.nan
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
        if _skip_or(out, 'attr', [dem_tif], 'dem'):
            continue
        result = gdal.DEMProcessing(str(out), ds, attr, computeEdges=True)
        result = None
    ds = None


def clip_tif_to_kml(in_tif: Path, out_tif: Path, kml: Path) -> None:
    """Clip a GeoTIFF to the polygon boundary defined in a KML file."""
    if _skip_or(out_tif, 'clip', [in_tif, kml], 'sz'):
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
# Pre-flight input check
# =============================================================================

def _analysis_dem(cfg: dict) -> Path:
    return cfg['data_dir'] / (cfg['crp_dem_name'] or cfg['dem_name'])


def _sz_pairs(cfg: dict) -> list[tuple[Path, Path]]:
    """(kml, sz_dem) for every individual start-zone KML."""
    d = cfg['data_dir']
    return [(kml, d / (kml.stem + cfg['sz_dem_suffix'] + '.tif'))
            for kml in sorted(d.glob(cfg['sz_kml_glob']))]


def check_inputs(cfg: dict, process: bool, analyze: bool) -> None:
    """
    Verify every input the requested steps need is present, before any
    processing starts.  Exits with a list of problems if not; prints
    warnings for outputs that are stale but will not be rebuilt.
    """
    d = cfg['data_dir']
    errors, warnings_ = [], []

    def rel(p: Path) -> str:
        try:
            return str(p.relative_to(ROOT))
        except ValueError:
            return str(p)

    def need(path: Path, what: str):
        if not path.exists():
            errors.append(f"{what}: {rel(path)} not found")

    # ── Start-zone polygons ──────────────────────────────────────────────────
    if cfg['combined_sz_kml']:
        need(d / cfg['combined_sz_kml'], 'combined SZ KML')
    if not _sz_pairs(cfg):
        errors.append(f"no start-zone KMLs match '{cfg['sz_kml_glob']}' "
                      f"in {rel(d)}")

    # ── Elevation source ─────────────────────────────────────────────────────
    analysis_dem = _analysis_dem(cfg)
    if process:
        # Walk the chain from the finished DEM back to the raw tiles; the
        # first file that exists (and whose stage is not forced) is enough.
        chain = [(analysis_dem, 'dem'), (d / cfg['dem_name'], 'dem'),
                 (d / cfg['ground_laz'], 'laz')]
        if cfg['spatial_bounds']:
            chain.append((d / cfg['cropped_laz'], 'laz'))
        chain.append((d / cfg['merged_laz'], 'laz'))
        if cfg.get('source_raster'):
            need(d / cfg['source_raster'], 'source_raster')
        have_source = (any(p.exists() and not _forced(s) for p, s in chain)
                       or any(d.glob(cfg['laz_glob']))
                       or bool(cfg.get('source_raster')
                               and (d / cfg['source_raster']).exists()))
        if not have_source:
            usable = [p.name for p, s in chain if not _forced(s)]
            errors.append(
                f"no elevation data in {rel(d)} — add LAZ tiles "
                f"matching '{cfg['laz_glob']}'"
                + (f" or one of: {', '.join(dict.fromkeys(usable))}" if usable else ''))
    elif analyze:
        need(analysis_dem, 'DEM (required with --skip-process)')
        for kml, sz_dem in _sz_pairs(cfg):
            if not sz_dem.exists():
                warnings_.append(f"{sz_dem.name} missing — that SZ will be "
                                 f"left out (drop --skip-process to build it)")
            else:
                newer = _stale_inputs(sz_dem, [kml, analysis_dem])
                if newer:
                    warnings_.append(f"{sz_dem.name} is older than "
                                     f"{_fmt_names(newer)} — analysing a stale "
                                     f"clip (drop --skip-process to rebuild)")

    # ── Analysis-only inputs ─────────────────────────────────────────────────
    if analyze:
        if cfg.get('sat_image'):
            need(d / cfg['sat_image'], 'satellite image')
        if cfg.get('avalanche_csv'):
            need(Path(cfg['avalanche_csv']), 'avalanche CSV')

    for w in warnings_:
        print(f"  [warn] {w}")
    if errors:
        print(f"\nCannot run {cfg['name']} — missing inputs:")
        for e in errors:
            print(f"  ✗ {e}")
        sys.exit(1)


# =============================================================================
# Processing pipeline
# =============================================================================

def process_site(cfg: dict) -> None:
    d   = cfg['data_dir']
    print(f"\n{'='*60}")
    print(f"  Processing: {cfg['name']}")
    print(f"{'='*60}\n")

    merged_laz  = d / cfg['merged_laz']
    cropped_laz = d / cfg['cropped_laz'] if cfg['spatial_bounds'] else None
    ground_laz  = d / cfg['ground_laz']
    dem_tif     = d / cfg['dem_name']
    final_dem   = _analysis_dem(cfg)

    def no_input(inp_exists: bool, *downstream: Path | None) -> bool:
        """True when a step's input is absent but a later product exists
        (e.g. only the DEM was copied in) — so the step can be skipped."""
        if inp_exists:
            return False
        have = next((p for p in downstream if p is not None and p.exists()), None)
        if have is not None:
            print(f"  skip  (input absent — using existing {have.name})")
            return True
        return False

    # 0. Optional pre-made elevation raster → DEM (replaces steps 1–4)
    analysis_built = False
    if cfg.get('source_raster'):
        print("Step 0 — Resample source raster")
        analysis_built = resample_source_raster(
            d / cfg['source_raster'], dem_tif, cfg['dem_resolution'],
            cfg.get('source_bounds'), cfg.get('source_resample', 'average'))
        print()

    # 1. Merge raw tiles
    print("Step 1 — Merge LAZ tiles")
    if not no_input(any(d.glob(cfg['laz_glob'])),
                    cropped_laz, ground_laz, dem_tif, final_dem):
        merge_laz_tiles(d, cfg['laz_glob'], merged_laz)

    # 2. Optional spatial crop
    laz_for_ground = merged_laz
    if cropped_laz is not None:
        print("\nStep 2 — Spatial bounds crop")
        if not no_input(merged_laz.exists(), ground_laz, dem_tif, final_dem):
            crop_laz_bounds(merged_laz, cropped_laz, cfg['spatial_bounds'])
        laz_for_ground = cropped_laz
    else:
        print("\nStep 2 — Spatial bounds crop  [skipped — not configured]")

    # 3. SMRF ground filter
    print("\nStep 3 — SMRF ground segmentation")
    if not no_input(laz_for_ground.exists(), dem_tif, final_dem):
        ground_filter_laz(laz_for_ground, ground_laz)

    # 4. LAZ → DEM
    print("\nStep 4 — LAZ → DEM")
    if cfg.get('source_raster'):
        print("  skip  (DEM comes from source_raster — step 0)")
    elif not no_input(ground_laz.exists() or dem_tif.exists(), final_dem):
        analysis_built = laz_to_dem(ground_laz, dem_tif, cfg['dem_resolution'])

    # 5. Optional pixel crop
    analysis_dem = dem_tif
    if cfg['pixel_crop_window']:
        print("\nStep 5 — Pixel crop")
        crp_dem     = d / cfg['crp_dem_name']
        analysis_built = crop_dem_pixels(dem_tif, crp_dem, cfg['pixel_crop_window'])
        analysis_dem = crp_dem
    else:
        print("\nStep 5 — Pixel crop  [skipped — not configured]")

    # 6. Optional gap fill — in-place, so only on a DEM built in this run
    #    (re-filling an already-filled DEM would corrupt it).
    if cfg['fill_gaps']:
        print("\nStep 6 — Fill DEM gaps")
        if analysis_built:
            fill_dem_gaps(analysis_dem)
        else:
            print(f"  skip  {analysis_dem.name}  (filled when it was built)")
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
    for kml, out_tif in _sz_pairs(cfg):
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


def _plot_sz_slope_map(cfg: dict, dem_path: str, out_png: Path,
                       margin_m: float = 90.0, contour_m: float = 10.0) -> None:
    """North-up map: hillshade + contours, slope inside the start zones,
    zone outlines and labels, north arrow and scale bar (map coordinates)."""
    ds = gdal.Open(str(dem_path))
    gt, wkt = ds.GetGeoTransform(), ds.GetProjection()
    z = ds.GetRasterBand(1).ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    if nd is not None:
        z[z == nd] = np.nan
    hs_ds = gdal.DEMProcessing('', ds, 'hillshade', format='MEM', computeEdges=True)
    hs = hs_ds.GetRasterBand(1).ReadAsArray().astype(float); hs_ds = None
    slope, _, _, _, _ = _full_dem_attributes(dem_path, abs(gt[1]))

    zones = [(kml, _kml_mask(kml, z.shape, gt, wkt)) for kml, _ in _sz_pairs(cfg)]
    allz = np.zeros(z.shape, bool)
    for _, m in zones:
        allz |= m
    pad = int(margin_m / abs(gt[1]))
    rr, cc = np.nonzero(allz)
    r0, r1 = max(rr.min() - pad, 0), min(rr.max() + pad, z.shape[0])
    c0, c1 = max(cc.min() - pad, 0), min(cc.max() + pad, z.shape[1])
    win = (slice(r0, r1), slice(c0, c1))
    x0, y0 = gt[0] + c0 * gt[1], gt[3] + r0 * gt[5]
    x1, y1 = gt[0] + c1 * gt[1], gt[3] + r1 * gt[5]
    ext = (x0, x1, y1, y0)                       # left, right, bottom, top

    fig, ax = plt.subplots(figsize=(11, 8))
    basemap = cfg.get('basemap')
    bm_path = cfg['data_dir'] / basemap if basemap else None
    if bm_path is not None and bm_path.exists():
        bm = gdal.Warp('', str(bm_path), format='MEM', outputBounds=(x0, y1, x1, y0),
                       width=c1 - c0, height=r1 - r0, dstSRS=wkt, resampleAlg='bilinear')
        rgb = np.dstack([bm.GetRasterBand(b).ReadAsArray().astype(float) for b in (1, 2, 3)])
        lo, hi = np.percentile(rgb, [1, 99])
        ax.imshow(np.clip((rgb - lo) / (hi - lo), 0, 1), extent=ext)
        base_label, cont_col, slope_alpha = cfg.get('basemap_label', basemap), 'white', .65
    else:
        ax.imshow(hs[win], cmap='gray', extent=ext, vmin=0, vmax=255)
        base_label, cont_col, slope_alpha = '3DEP hillshade', 'saddlebrown', .85
    zw = z[win]
    xs = np.linspace(x0, x1, zw.shape[1]); ys = np.linspace(y0, y1, zw.shape[0])
    lv = np.arange(np.floor(np.nanmin(zw) / contour_m) * contour_m, np.nanmax(zw), contour_m)
    cs = ax.contour(xs, ys, zw, levels=lv, colors=cont_col, linewidths=.4, alpha=.5)
    major = lv[(lv % (5 * contour_m)) == 0]
    ax.contour(xs, ys, zw, levels=major, colors=cont_col, linewidths=.9, alpha=.75)
    ax.clabel(cs, levels=major, fmt='%d m', fontsize=7, inline=True, colors=cont_col)
    im = ax.imshow(np.ma.masked_where(~allz[win], slope[win]), cmap='magma_r',
                   vmin=25, vmax=60, extent=ext, alpha=slope_alpha)
    for kml, m in zones:
        ax.contour(xs, ys, m[win].astype(float), levels=[.5], colors='cyan', linewidths=1.5)
        rr_, cc_ = np.nonzero(m[win])
        num = ''.join(ch for ch in kml.stem if ch.isdigit()) or kml.stem
        ax.text(xs[int(cc_.mean())], ys[rr_.min()] + 12, f'#{num}', fontsize=12, ha='center',
                va='bottom', bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=.85, lw=0))
    fig.colorbar(im, ax=ax, fraction=.03, label='slope (°)')
    # scale bar and north arrow
    sb = 100.0
    bx, by = x0 + .05 * (x1 - x0), y1 + .05 * (y0 - y1)
    ax.plot([bx, bx + sb], [by, by], 'k-', lw=4, solid_capstyle='butt')
    ax.plot([bx, bx + sb], [by, by], 'w-', lw=2, solid_capstyle='butt')
    ax.text(bx + sb / 2, by + .015 * (y0 - y1), f'{sb:.0f} m', ha='center', fontsize=9,
            bbox=dict(fc='white', alpha=.8, lw=0))
    ax.annotate('N', xy=(x1 - .05 * (x1 - x0), y1 + .17 * (y0 - y1)),
                xytext=(x1 - .05 * (x1 - x0), y1 + .06 * (y0 - y1)),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2), ha='center', fontsize=12,
                fontweight='bold', bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=.85, lw=0))
    epsg = osr.SpatialReference(wkt=wkt).GetAuthorityCode(None)
    ax.set_xlabel(f'Easting (m, EPSG:{epsg})'); ax.set_ylabel('Northing (m)')
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.tick_params(labelsize=7)
    ax.set_title(f"{cfg['name']} — start-zone slope on {base_label} "
                 f"({contour_m:.0f} m contours from 3DEP, north up)")
    fig.savefig(str(out_png), dpi=150, bbox_inches='tight')
    plt.close(fig)


def _full_dem_attributes(dem_path: str, res: float):
    """Slope (deg), aspect (deg) and profile curvature on the whole DEM."""
    ds = gdal.Open(str(dem_path))
    z = ds.GetRasterBand(1).ReadAsArray().astype(float)
    nd = ds.GetRasterBand(1).GetNoDataValue()
    z[(z == nd) if nd is not None else (z == z.min())] = np.nan
    # Keep each MEM dataset referenced until read: chaining
    # DEMProcessing(...).GetRasterBand(1) lets Python free it first.
    s_ds = gdal.DEMProcessing('', ds, 'slope', format='MEM', computeEdges=True)
    a_ds = gdal.DEMProcessing('', ds, 'aspect', format='MEM', computeEdges=True)
    slope = s_ds.GetRasterBand(1).ReadAsArray().astype(float)
    aspect = a_ds.GetRasterBand(1).ReadAsArray().astype(float)
    s_ds = a_ds = None
    slope[np.isnan(z)] = np.nan
    aspect[np.isnan(z) | (aspect < 0)] = np.nan     # GDAL flags flat cells −9999
    curv = _profile_curvature(z, res=res)
    return slope, aspect, curv, ds.GetGeoTransform(), ds.GetProjection()


def _kml_mask(kml_path: Path, shape: tuple, gt: tuple, wkt: str) -> np.ndarray:
    """Rasterise a KML polygon onto a grid (cell centres inside)."""
    from osgeo import ogr, osr
    src = ogr.Open(str(kml_path))
    dst_srs = osr.SpatialReference(); dst_srs.ImportFromWkt(wkt)
    dst_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    mem = gdal.GetDriverByName('MEM').Create('', shape[1], shape[0], 1, gdal.GDT_Byte)
    mem.SetGeoTransform(gt); mem.SetProjection(wkt)
    vds = ogr.GetDriverByName('Memory').CreateDataSource('')
    lyr = vds.CreateLayer('z', srs=dst_srs)
    for slyr in src:
        s_srs = slyr.GetSpatialRef()
        if s_srs is not None:
            s_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        tf = osr.CoordinateTransformation(s_srs, dst_srs) if s_srs is not None else None
        for f in slyr:
            g = f.GetGeometryRef().Clone(); g.FlattenTo2D()
            if tf is not None:
                g.Transform(tf)
            nf = ogr.Feature(lyr.GetLayerDefn()); nf.SetGeometry(g); lyr.CreateFeature(nf)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1])
    return mem.GetRasterBand(1).ReadAsArray().astype(bool)


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
        # ── Weak-layer stress under a uniform vertical load HN ────────────────
        # τ = ρ g HN sinθ cosθ ;  σ = ρ g HN cos²θ ;  τ/σ = tanθ.
        # tanθ is the driving (shear/normal) ratio a friction-type failure
        # criterion compares against; it rises monotonically with slope.
        # sinθcosθ is shear stress per unit load; it peaks at 45°.
        tan_mean             = float(np.mean(np.tan(np.radians(slopes)))),
        tan_p90              = float(np.percentile(np.tan(np.radians(slopes)), 90)),
        shear_per_load_mean  = float(np.mean(np.sin(np.radians(slopes))
                                             * np.cos(np.radians(slopes)))),
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

        if len(vals) == 0:
            print(f"  [warn] {Path(path).stem}: no valid values for {attrib}, skipping")
            continue

        unique, counts = np.unique(vals, return_counts=True)
        if counts.max() / vals.shape[0] > 0.75:
            vals = vals[vals != unique[counts.argmax()]]

        if lo is not None:
            vals = vals[(vals >= lo) & (vals <= hi)]

        if len(vals) == 0:
            print(f"  [warn] {Path(path).stem}: no values in range [{lo}, {hi}] for {attrib}, skipping")
            continue
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

def _nucleation_summary(site_slug: str) -> list[str]:
    """Report lines summarising the MicroMet nucleation-potential outputs.

    The nucleation map itself (GeoTIFF, KML of high-potential patches and
    figure) is produced by src/micromet_loading.py, which needs the weather
    record; this only reads its tables.  Returns [] when they are missing.
    """
    patches = FIGURES / f"{site_slug}_nucleation_patches.csv"
    oe_csv  = FIGURES / f"{site_slug}_observed_vs_expected.csv"
    if not patches.exists():
        return []
    p = pd.read_csv(patches).set_index('zone')
    lines = [f'  Patches ≥ P95 of start-zone potential (≥ 2 connected cells):',
             '  ' + p.round(0).astype(int).to_string().replace('\n', '\n  ')]
    if oe_csv.exists():
        oe = pd.read_csv(oe_csv)
        oe = oe[(oe.threshold == 'P99') & (oe.basis == 'hotspot days')]
        oe = oe.assign(zone=oe.zone.str.replace('_SZ', '').str.replace('_', ' #'))
        lines += ['', '  Natural-release days, observed vs expected from hotspot days (P99):',
                  '  ' + oe[['zone', 'observed', 'expected']].to_string(index=False)
                  .replace('\n', '\n  ')]
    lines += ['', f'  Map: {site_slug}_nucleation_potential.png / .tif;'
                  f' patches: {site_slug}_nucleation.kml']
    return lines


def _write_report(cfg: dict, site_slug: str, analysis_df: pd.DataFrame,
                  corr_df: pd.DataFrame, total_counts: pd.DataFrame,
                  boot_results: list, period: str = '',
                  nucleation_lines: list | None = None) -> None:
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

    slope_rho      = avi_corrs.loc['% 40–50°',               'spearman_rho']

    excl = cfg.get('exclude_sz_stems', [])
    n    = len(analysis_df)

    def _fmt_zone(z):
        return z.replace('_SZ_DEM', '').replace('_', ' ')

    lines = [
        '=' * W,
        f"  {cfg['name'].upper()} — AVALANCHE TERRAIN ANALYSIS",
        '=' * W,
        '',
        f'AVALANCHE COUNTS  (Natural HS/SS, CAIC{"; events " + period if period else ""})',
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
        '                         interior cells (boundary pixels excluded); positive',
        '                         = convex roll.  Kept for reference — not used in',
        '                         the current analysis.',
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

    # ── Summary ───────────────────────────────────────────────────────────────
    excl_str = ', '.join(_fmt_zone(z) for z in excl) if excl else 'none'
    ranked   = avi_corrs['spearman_rho'].abs().sort_values(ascending=False)
    top      = [f'    {m}: ρ={avi_corrs.loc[m, "spearman_rho"]:+.2f}'
                for m in ranked.index[:3]]
    lines += [
        '',
        'SUMMARY',
        '-' * W,
        f'N = {n} start zones (excluded: {excl_str}).  At this N, zone-level',
        '  correlations are exploratory: the p-values in the table above are',
        '  asymptotic and too small; exact permutation p (event_stats.py) is',
        '  ≥ 0.017 at N = 5 and ≈ 0.08 for ρ = 0.9.  Outcome is naturals per',
        '  10,000 m², which favours small zones.  Strongest |ρ| this run:',
    ] + top + [
        f'  % in the 40–50° band: ρ={slope_rho:+.2f}.  Per-area rates differ from',
        '  the count-based comparison in the research note (§5.1, mean tanθ vs',
        '  natural count ρ = 0.2), where steepness does not track releases.',
        '  pct_steep_near_loading uses aspect only (no wind model); the loading',
        '  analysis proper is in micromet_loading.py.',
        '',
        'NUCLEATION POTENTIAL  (from src/micromet_loading.py)',
        '-' * W,
    ]
    lines += nucleation_lines or [
        '  Not available — run src/micromet_loading.py --site <site>, then',
        '  rerun this script with --skip-process to include it here.']
    lines += [
        '',
        'See docs/terrain_analysis_report.md for the full analysis.',
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

    # ── 2. Combined SZ DEM slope map ─────────────────────────────────────────
    sz_dem_path = str(d / cfg['combined_sz_dem'])
    if cfg.get('combined_sz_kml'):
        combined_kml = d / cfg['combined_sz_kml']
        if not _skip_or(Path(sz_dem_path), 'clip',
                        [Path(dem_path), combined_kml], 'sz'):
            crop_tif(dem_path, sz_dem_path, str(combined_kml))

    print("Step 4 — Start-zone slope map (north up)")
    _plot_sz_slope_map(cfg, dem_path, FIGURES / f"{site_slug}_sz_slope_map.png")
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
    plot_sz_distributions(
        sz_dem_files,
        attrib=cfg['slope_attr'],
        title=f"{name} — Slope angle distribution",
        out_fig=FIGURES / f"{cfg['name'].lower().replace(' ', '_')}_slope.png",
    )

    print(f"\nStep 6 — Per-SZ aspect distributions")
    plot_sz_distributions(
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

    print("\nStep 8 — Extended terrain metrics")
    # Slope, aspect and curvature come from the FULL DEM and are then masked
    # to each zone.  Computing them on the clipped SZ DEMs biased edge cells:
    # their 3×3 windows reach into NoData and get extrapolated (Sister #6's
    # % > 40° was 31.9 vs 41.9 from the full DEM).
    full_slope, full_aspect, full_curv, full_gt, full_wkt = _full_dem_attributes(
        dem_path, cfg['dem_resolution'])
    sz_metrics = OrderedDict()
    for p in sz_dem_files:
        sid      = Path(p).stem
        kml      = d / ((sid[:-len(cfg['sz_dem_suffix'])] if cfg['sz_dem_suffix'] else sid) + '.kml')
        if not kml.exists():
            print(f"  [warn] {sid}: {kml.name} not found — skipping")
            continue
        m = _kml_mask(kml, full_slope.shape, full_gt, full_wkt)
        if not m.any():
            print(f"  [warn] {sid}: zone does not overlap the DEM — skipping")
            continue
        rows_, cols_ = np.nonzero(m)
        win = (slice(rows_.min(), rows_.max() + 1), slice(cols_.min(), cols_.max() + 1))
        def _masked(a):
            return np.where(m[win], a[win], np.nan)
        sz_metrics[sid] = compute_sz_metrics(
            _masked(full_slope),
            _masked(full_aspect),
            curv_arr=_masked(full_curv),
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
        # weak-layer stress (uniform load): driving ratio and shear per load
        ('tan_mean',                'Mean tanθ (τ/σ)'),
        ('tan_p90',                 'P90 tanθ'),
        ('shear_per_load_mean',     'Mean sinθcosθ (τ/load)'),
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

    print("\nStep 13 — Nucleation potential (from micromet_loading.py)")
    nuc_lines = _nucleation_summary(site_slug)
    for ln in nuc_lines or ['  not available — run src/micromet_loading.py first']:
        print(ln)

    period = ''
    if 'Date' in avi_df.columns:
        yrs = pd.to_datetime(avi_df['Date'], format='mixed', errors='coerce').dt.year
        if yrs.notna().any():
            period = f'{int(yrs.min())}–{int(yrs.max())}'
    _write_report(cfg, site_slug, analysis_df, corr_df, avi_counts, boot_results,
                  period=period, nucleation_lines=nuc_lines)

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
    parser.add_argument('--site', choices=list(SITES), default=None,
                        help='Site to process/analyse (sisters | star_mtn | us550). '
                             'Omit to run all sites.')
    parser.add_argument('--skip-process', action='store_true',
                        help='Skip the LAZ → DEM processing steps')
    parser.add_argument('--skip-analyze', action='store_true',
                        help='Skip the terrain analysis steps')
    force = parser.add_mutually_exclusive_group()
    force.add_argument('--force', action='store_const', const='laz',
                       dest='force_from',
                       help='Rebuild every processing output (same as --force-from laz)')
    force.add_argument('--force-from', choices=STAGES, default=None,
                       help='Rebuild this stage and all later ones regardless of '
                            'timestamps: laz (merge/crop/SMRF), dem (DEM, crop, '
                            'fill, attributes), sz (start-zone clips)')
    args = parser.parse_args()

    global _FORCE_FROM
    _FORCE_FROM = args.force_from

    FIGURES.mkdir(exist_ok=True)
    sites = [args.site] if args.site else list(SITES)

    # Check every site up front so a missing file fails fast, not mid-run.
    for site_key in sites:
        check_inputs(SITES[site_key], process=not args.skip_process,
                     analyze=not args.skip_analyze)

    for site_key in sites:
        cfg = SITES[site_key]
        if not args.skip_process:
            process_site(cfg)
        if not args.skip_analyze:
            analyze_site(cfg)


if __name__ == '__main__':
    main()
