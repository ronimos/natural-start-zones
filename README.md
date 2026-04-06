# Natural Avalanche Start-Zone Terrain Analysis

## What we are trying to achieve

Natural slab avalanches on highway corridors are responsible for road closures and fatalities, yet the terrain features that make one start zone release more frequently than a neighbouring one are not well quantified.  This project asks: **given two start zones on the same mountain, both already steep enough to release an avalanche, what terrain features explain why one releases far more often than the other?**

The answer is not simply "it is steeper" — all start zones in the dataset are pre-selected steep terrain.  We are looking for the finer-scale geometry that controls where a wind slab accumulates and where that slab is most likely to crack.

A secondary goal is to locate the most probable **crack nucleation sites** within each start zone — the small patches where a fracture is most likely to initiate and then propagate into a full slab release.

---

## Study sites

All three sites are along US-550 (the Million Dollar Highway) in the San Juan Mountains of Colorado.  Prevailing storm winds arrive from the W–NW (225–315°).

| Site | Start zones | Data source |
|---|---|---|
| **Seven Sisters** | #1, #2, #3, #4, #6 | USGS 3DEP LiDAR + CAIC |
| **Star Mountain** | A, B, C, D, E | USGS 3DEP LiDAR + CAIC |
| **US-550** | Eagle, Muleshoe, Porcupine, Telescope | USGS 3DEP LiDAR + CAIC |

Avalanche records come from the Colorado Avalanche Information Center (CAIC) highway observer database, 2009–2021.  Only **natural** hard-slab (HS) and soft-slab (SS) avalanches are used.  Human-triggered events are excluded to isolate terrain-controlled release.

---

## Data pipeline

```
Raw LiDAR tiles (LAZ)
        │
        ▼  Step 1  merge tiles
        │
        ▼  Step 2  optional spatial crop (bounding box)
        │
        ▼  Step 3  SMRF ground segmentation
        │
        ▼  Step 4  ground points → DEM (GeoTIFF, 1.0–1.2 m/pixel)
        │
        ▼  Step 5  optional pixel crop
        │
        ▼  Step 6  optional Gaussian gap-fill
        │
        ▼  Step 7  terrain attributes (slope, aspect, hillshade, TRI, roughness)
        │
        ▼  Steps 8–9  clip to start-zone polygons (KML)
        │
        ▼  Steps 10–12  compute metrics, correlate with avalanche frequency
        │
        ▼  Step 13  nucleation patch map + KML export
```

### Why this processing chain?

**LiDAR at 1–1.2 m/pixel.**  Airborne LiDAR gives sub-metre ground elevation even through vegetation.  At this resolution, individual convex rolls and slope breaks — the features most relevant to slab mechanics — are resolved.  Coarser 10 m DEMs smooth away these features entirely.

**SMRF ground filter.**  The Simple Morphological Filter (Pingel et al. 2013) is robust to the dense coniferous canopy and steep rocky ridges found in the San Juans.  We apply ELM (extended local minimum) noise removal and an outlier filter before SMRF, then extract only class-2 (ground) returns.  Parameters: slope 1.7, window 16 m, threshold 0.45 m — tuned for alpine terrain with moderate relief within each filter window.

**Mean-elevation DEM.**  Ground returns within each pixel are averaged.  This is conservative compared to minimum-elevation but reduces sensitivity to isolated rock outcrops and residual noise.

**Per-start-zone clipping.**  Each start zone boundary is drawn in a KML file by hand against ortho imagery.  The full-area DEM is then clipped to each polygon so that metrics are computed only on the terrain actually included in the start zone, not the surrounding runout or terrain traps.

---

## Terrain metrics

All metrics are computed per start zone from the clipped DEM.

### Slope distribution

| Metric | What it measures |
|---|---|
| `slope_mean`, `slope_std`, `slope_cv` | Central tendency and spread of the slope distribution |
| `pct_40_50` | Fraction of cells in the 40–50° core slab release band |
| `slope_max_1cell` | Steepest single cell (≈1–1.4 m²) — proxy for peak driving stress |

**Why slope alone is not enough.**  All start zones were selected as existing avalanche paths, so every zone already has substantial area above 35°.  Raw steepness does not discriminate high- from low-frequency zones; the correlation between `pct_40_50` and avalanche frequency is near zero or negative.  We need metrics that capture *how* steep terrain interacts with wind loading and terrain geometry.

### Aspect diversity

| Metric | What it measures |
|---|---|
| `aspect_entropy` | Shannon entropy (bits) of aspect over 8 compass octants across all cells |
| `aspect_entropy_steep` | Same, restricted to cells ≥35° |
| `aspect_circ_var` | Circular variance of aspect (0 = unimodal, 1 = uniform) |

**Why aspect entropy is the key metric.**  Wind from W–NW loads north- and east-facing lee slopes.  A start zone whose steep cells face only one direction accumulates a slab primarily during storms from that single direction.  A start zone whose steep cells face *multiple* directions accumulates wind slabs from many storm tracks.  Shannon entropy captures this: a high-entropy start zone is exposed to loading from more directions, increasing the probability that some combination of snowpack and wind will produce a natural release.

On the Seven Sisters, `aspect_entropy_steep` is the strongest terrain predictor of avalanche frequency (Spearman ρ = +0.80).  Sisters #4 and #3 — the most complex terrain with aspects ranging from N through E — are also the most active paths.

### Profile curvature (stress concentration)

| Metric | What it measures |
|---|---|
| `profile_curv_p95_steep` | 95th-percentile profile curvature on steep (≥40°) interior cells |

Profile curvature describes how the slope is changing in the downslope direction.  A **positive** (convex) value means the slope is steepening as you look downhill — a convex roll.

This geometry is directly analogous to a notch in fracture mechanics: it concentrates the normal stress at the base of the slab, lowering the energy required to initiate a mode-II (shear) crack.  Gvirtzman et al. (2025, *Nature*) showed that convex terrain features are the preferred nucleation geometry for dry-snow slab fractures.  We use the 95th percentile on interior cells (boundary pixels within 2 cell-lengths of NoData are excluded) to capture the most extreme stress-concentration sites while avoiding DEM edge artifacts.

### Co-occurrence of loading and steepness

`pct_steep_near_loading` — the fraction of start-zone cells that are both steep (≥40°) AND within ~2.4 m of a W–NW-facing cell.

This metric is a proxy for the spatial coincidence of wind-loaded faces and steep terrain.  A high score means that loaded faces and steep release terrain are geographically interleaved, increasing the chance that a wind slab deposited on a loaded face sits directly on a steep slope where driving stress is high.

**Caveat:** this metric is a necessary but not sufficient condition for co-occurrence.  A start zone can score high on both steepness and loading metrics independently without those conditions coinciding in the same cells.  The metric should be interpreted as a spatial proximity proxy, not a direct measurement of wind slab deposition.

---

## Crack nucleation analysis

Beyond per-zone statistics, we map the most likely **crack nucleation sites** — the specific square metres where a fracture is most likely to initiate.

A nucleation site requires all three of:

1. **Steep** — slope ≥ 40° (driving shear stress high enough to overcome weak layer strength)
2. **Wind-loaded lee face** — cell is E-facing (45–135°) AND within ~3 m of N-facing terrain (315–45°)
3. **Convex** — positive profile curvature (stress concentration)

### Why E-facing, near N-facing terrain?

Wind from W–NW (the dominant storm direction here) accelerates over N-facing ridges and deposits snow on the **lee** (east-facing) slope immediately behind the ridgeline.  The critical zone is the **N→E rollover** — the convex break where a N-facing ridge transitions to an E-facing slope.  This is where:

- Snow transport stops and deposition begins (wind shadow)
- The slope steepens from the ridge cap into the face
- Profile curvature peaks (the geometric notch)

The `loaded_zone` mask is computed as E-facing cells within 2 cell-lengths (~2.4 m) of a N-facing cell, ensuring we capture terrain right at the rollover rather than the broader E-facing slope below it.

### Output

- **`figures/{site}_nucleation_hs.png`** — hillshade map with nucleation candidates (red) and two-of-three-criteria cells (yellow) overlaid on start-zone boundaries
- **`figures/{site}_nucleation.kml`** — KML placemarks at the centroid of each connected cluster, with mean slope and aspect on click, for overlay in Google Earth or similar tools

---

## Statistical approach

Avalanche frequency is normalised by start-zone area (avalanches per 10,000 m²) to remove the trivial size effect — a larger zone has more terrain and should produce more avalanches by chance.

Correlations use **Spearman ρ** (rank-based, robust to non-normality and small N) and **Kendall τ** (more conservative, appropriate when N < 10).  Bootstrap 95% confidence intervals are computed with 10,000 resamples.

**Important limitation:** N = 4 at Seven Sisters (Sister #6 is excluded — it has a different loading regime as a convex wind-scoured ridge rather than a lee-loading face).  All p-values are exploratory; effect sizes are reported rather than hypothesis tests.

---

## Running the pipeline

```bash
# Install dependencies (requires Python 3.10+)
uv sync          # or: pip install -r requirements.txt

# Full run for a site
python src/pipeline.py sisters
python src/pipeline.py star_mtn
python src/pipeline.py us550

# Skip the LAZ processing if DEMs already exist
python src/pipeline.py sisters --skip-process

# Skip analysis (processing only)
python src/pipeline.py sisters --skip-analyze
```

Outputs land in `figures/` as PNG figures, a text report, and a KML file.

### Adding a new site

Add an entry to the `SITES` dict in `src/pipeline.py`.  Each key is documented inline.  You need:

- A directory under `data/` with LAZ tiles and KML start-zone polygons
- An avalanche CSV (CAIC format or the simplified two-column format)
- A `sz_avi_path_map` mapping DEM stems to CAIC path names

---

## Project structure

```
data/
  Sisters/          Seven Sisters LiDAR, DEMs, KMLs, satellite image
  Star_mtn/         Star Mountain LiDAR, DEMs, KMLs
  550/              US-550 corridor LiDAR, DEMs, KMLs
  avalanches/       CAIC CSV records
  masks/            Manual snow masks (used in earlier exploratory work)
figures/            All output figures, reports, and KML files
src/
  pipeline.py       Main processing and analysis pipeline
  terrain_util.py   DEM I/O, terrain attribute computation (GDAL + WhiteboxTools)
  avl_data.py       Avalanche record loading and filtering
unused/             Earlier exploratory scripts (kept for reference)
```

---

## Key references

- Gvirtzman, S. et al. (2025). Detection of likely shear initiation location for dry-snow slab avalanches. *Nature*, https://doi.org/10.1038/s41586-024-08287-y
- Pingel, T. J., Clarke, K. C., & McBride, W. A. (2013). An improved simple morphological filter for the terrain classification of airborne LIDAR data. *ISPRS Journal of Photogrammetry and Remote Sensing*, 77, 21–30.
- Colorado Avalanche Information Center (CAIC) highway observer database, 2009–2021.
