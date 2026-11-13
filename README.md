# Natural Avalanche Start Zones

Why do some avalanche paths release naturally far more often than their neighbours, when all of them share the same storms and are controlled together?  This project combines high-resolution terrain, modelled wind transport and station weather to answer that question for the **Seven Sisters** paths above US-6 on Loveland Pass, Colorado, and to locate where within each start zone a slab is most likely to start.

The write-up is in two research notes:

- **[docs/terrain_analysis_report.md](docs/terrain_analysis_report.md)** — *which paths release naturally, and why*: data, terrain, the wind-loading and stress model, hotspots, nucleation potential, the storm-by-storm test and loading days by wind direction.
- **[docs/release_mechanism_note.md](docs/release_mechanism_note.md)** — *when they release, and how a natural slab might start*: the stress threshold, avalanche size, loading rate, hour-by-hour stress before timed releases, and whether natural slabs start in shear rather than as an anticrack.

This README covers what the code does and how to run it.

---

## Approach in brief

The working hypothesis comes from Gvirtzman et al. (2025, *Nature*), a laboratory study of frictional ruptures: a rupture nucleates wherever the **local** shear stress exceeds a threshold, and the nucleation rate grows with the excess stress.  The paper is not about snow; applying it to slab avalanches is this project's analogy.  If it holds, what matters is how often, and by how much, small patches of a start zone exceed a weak-layer stress threshold — not the zone's average steepness.

1. **Terrain.** USGS 3DEP 1 m bare-earth lidar, resampled to 1.2 m (start-zone analysis) and 2 m (wind models).  Start zones are hand-digitised KML polygons.
2. **Wind loading.** MicroMet-style terrain winds (Liston & Elder 2006) driven by the CALVP ridge station, with Pomeroy & Gray (1990) saltation and a downwind erosion/deposition march, run for every hour of the 2009–2026 record.  Forest (from a drone canopy surface) traps snow; cells within 4 m of canopy carry no slab stress.
3. **Stress.** Load = SNOTEL snowfall (LBAC2) + modelled wind deposition; weak-layer stress margin = g·SWE·cos²θ·(tanθ − μ), μ = tan 30°.
4. **Nucleation potential.** Per cell, the excess 72-h stress margin above a threshold, summed over each season.  Patches of high potential are the most likely nucleation sites.
5. **Tests.** Observed vs expected natural-release days per path, a storm-by-storm test (does the model pick the days, and the path, that release?), and the stress threshold implied by the record.

**Headline results** (see the notes for numbers and caveats): steepness alone does not separate the paths; wind landing sites do — Sister #2, the quiet path, has no W–NW landing area.  The stress model predicts natural-slab days (AUC 0.76) and which path releases on those days (0.74).  Days on which W–NW winds load a path's hotspot terrain put the paths in the same order as their natural counts; SW loading days do not.  The record implies a threshold of ≈ 30 Pa of added stress (≈ 14 mm SWE in 72 h on a 45° slope), with the release rate rising about 20-fold well above it.  The last 24 h of loading carry the timing signal (AUC 0.83), and every natural with a recorded time released while snowfall or wind transport was still loading the path — consistent with, but not proof of, natural slabs starting as slow shear failures.

---

## Data

| Data | Source | Location |
|---|---|---|
| Bare-earth DEM | USGS 3DEP 1 m, DRCOG 2020 lidar | `data/Sisters/` (`source_raster` in `SITES`) |
| Canopy | Drone photogrammetry DSM minus 3DEP (2026) | `data/Sisters/260922_Sisters_Full_Ground-DSM.tiff` |
| Imagery | USGS NAIP 0.6 m, 2019-09-22 | `data/Sisters/NAIP_sisters.tif` |
| Start zones | Hand-digitised KML | `data/Sisters/Sister_?_SZ.kml`, `Sisters_SZs.kml` |
| Avalanches | CAIC highway export 2009–2021 + CAIC coordinate export 2020–2026, merged | `data/avalanches/` |
| Weather | CAIC database: CALVP (wind, temperature), LBAC2 SNOTEL (SWE) | `data/weather/sisters.csv` (fetched by `sql_util.py`) |

Only **natural (N) hard and soft slabs (HS, SS)** are analysed.  Five start zones: Sisters #1, #2, #3, #4 and #6 — all the active paths (there is no Sister #5 on the ground, and #7 is inactive).

Star Mountain (`star_mtn`) is configured for the terrain pipeline only; it has no weather or wind-loading analysis yet.  A US-550 site entry is kept, commented out, in `SITES`.

---

## Setup

GDAL and PDAL need full native builds, so the environment comes from conda-forge via [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html):

```bash
micromamba create -p ./.conda -f environment.yml   # once
source .envrc                                      # each new shell (automatic with direnv)
```

In VS Code, select `.conda/bin/python` as the interpreter.  `sql_util.py` needs access to the CAIC MySQL database (localhost).

---

## Running the Seven Sisters analysis

In order (each step reads the previous steps' outputs from `figures/`):

```bash
source .envrc
python src/assign_avalanche_paths.py data/avalanches/caic-avalanche-records-2026-11-12.csv
python src/merge_avalanche_records.py \
    data/avalanches/CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv \
    figures/seven_sisters_new_records_assigned.csv \
    data/avalanches/CAIC_HWY_avalanches_2009-2026_sisters_merged.csv
python src/sql_util.py --site sisters --end 2026-06-01   # weather (slow: ~5 min)
python src/pipeline.py --site sisters                    # DEM, zone clips, zone metrics
python src/wind_exposure.py --site sisters               # Winstral Sx, upwind fetch
python src/storm_analysis.py --site sisters              # weather before each natural
python src/stress_hotspots.py --site sisters             # fetch-index stress hotspots
python src/micromet_loading.py --site sisters            # wind loading, nucleation potential, storm test (~7 min)
python src/event_stats.py --site sisters                 # bootstrap / permutation / period split
python src/threshold_analysis.py --site sisters          # stress threshold, release rate, size
python src/loading_days.py --site sisters                # loading days by wind direction vs natural days
python src/loading_rate.py --site sisters                # 24/48/72-h loading windows; loading rate at fixed load
python src/event_stress_timeline.py --site sisters       # hourly hotspot stress before timed releases (~1 min)
```

`pipeline.py` runs before the wind model, so its text report includes the nucleation-potential summary only from the second run on (`python src/pipeline.py --site sisters --skip-process`).

MicroMet sensitivity runs: `--eta 15|60` (curvature length, m) or `--lf 50|200` (saturation length, m) with `--tag _eta15` etc. to keep outputs apart.  Stress variants for the shear-vs-normal test: `--stress normal|normal30|shear --tag _st_normal` etc., then `python src/stress_variant_test.py --site sisters`.  `--plot-only` redraws figures and re-exports the nucleation GeoTIFF/KML from the saved grids without rerunning the model.

### Scripts

| Script | Does |
|---|---|
| `pipeline.py` | Site config (`SITES`); LAZ or raster → DEM, terrain attributes, start-zone clips; per-zone terrain metrics and correlations with natural counts; text report |
| `assign_avalanche_paths.py` | Match CAIC records that have coordinates but no path name to paths, by tracing each start zone downhill |
| `merge_avalanche_records.py` | Append the matched records to the named highway export (duplicates kept from the highway export) |
| `sql_util.py` | Fetch hourly station data from the CAIC database |
| `wind_exposure.py` | Winstral shelter index and upwind source area per wind direction |
| `storm_analysis.py` | Wind and snowfall in the 72 h before each natural |
| `stress_hotspots.py` | Stress margin and hotspot patches from fetch-index loading |
| `micromet_loading.py` | MicroMet winds + saltation transport; landing sites; zone-day stress; nucleation potential; storm-by-storm test |
| `event_stats.py` | Exact permutation p-values (zone level), cluster bootstrap and date-shift tests (storm level) |
| `threshold_analysis.py` | Stress threshold implied by the record; release rate vs excess stress; avalanche size |
| `loading_days.py` | Days per season each wind direction loads each path's hotspot terrain, against natural days |
| `loading_rate.py` | Whether the last 24 h of loading, or its concentration, predicts naturals beyond the 72-h load |
| `stress_variant_test.py` | Shear margin vs normal stress vs pure shear as the stress per unit load: day skill, path picking, per-path pattern (needs the `--stress` runs of `micromet_loading.py`) |
| `event_stress_timeline.py` | Hour-by-hour hotspot stress (snowfall + MicroMet wind deposition) from 24 h before each natural with a release time, plus station weather around release |
| `terrain_util.py`, `avl_data.py` | DEM I/O and terrain attributes; avalanche CSV loading |

### Main outputs (`figures/`, prefix `seven_sisters_`)

| File | Content |
|---|---|
| `sz_slope_map.png` | Slope inside the start zones on NAIP |
| `loading_patterns.png`, `loading_closeups.png` | Modelled wind accumulation per season, by sector |
| `nucleation_potential.png` | Nucleation-potential map with release points; observed vs expected natural days |
| `nucleation_potential.tif` | The same map as a GeoTIFF (band 1: Pa·day per season; band 2: days per season above threshold) |
| `nucleation.kml` | High-potential patches (≥ P95 of start-zone cells, ≥ 2 cells) as polygons, for Google Earth |
| `stress_hotspots.png` | Fetch-index stress hotspots |
| `threshold_analysis.png` | Skill vs stress threshold; hotspot days vs natural days per path |
| `loading_days.png` | Loading days per season by wind direction and path, with natural days |
| `event_stress_timeline.png`, `event_stress_panels.png` | Hotspot stress from 24 h before each timed natural, aligned at release and per storm |
| `report.txt` | Zone terrain metrics, correlations and nucleation-potential summary |
| `*.csv` | Tables behind every result in the report |

`figures/` is git-ignored; files the report uses are committed with `git add -f`.

---

## Terrain pipeline options

```bash
python src/pipeline.py --site sisters                    # process + analyse
python src/pipeline.py --site sisters --skip-process     # analysis only
python src/pipeline.py --site sisters --skip-analyze     # processing only
python src/pipeline.py --site sisters --force-from sz    # re-clip all start zones
python src/pipeline.py --site sisters --force-from dem   # DEM onward
python src/pipeline.py --site sisters --force            # everything
```

Processing steps are skipped when their output is newer than all of its inputs, so editing a start-zone KML re-clips just that zone.  Missing inputs are reported before anything runs.

A site takes either LAZ tiles (merged, SMRF ground-filtered and gridded) or a pre-made elevation raster via `source_raster`, block-averaged or resampled to `dem_resolution`.  The raster must be **bare earth** — a surface model that includes trees biases slope, aspect and curvature.

### Adding a site

Add an entry to `SITES` in `src/pipeline.py` (every key is documented inline).  The terrain pipeline needs a DEM source, start-zone KMLs, an avalanche CSV and a `sz_avi_path_map` from DEM stems to CAIC path names.  The wind-loading scripts also need `source_raster`, `weather_stations` and a weather CSV from `sql_util.py`; `micromet_loading.py` currently has the CALVP station location hard-coded.

---

## Project structure

```
data/
  Sisters/       3DEP DEM, drone DSM, NAIP, start-zone KMLs, GCPs
  Star_mtn/      Star Mountain LiDAR, DEMs, KMLs
  550/           US-550 corridor LiDAR, DEMs, KMLs
  avalanches/    CAIC exports and the merged Sisters record
  weather/       Hourly station data (from sql_util.py)
  masks/         Manual snow masks (earlier exploratory work)
docs/            Research notes (terrain_analysis_report.md, release_mechanism_note.md)
figures/         Outputs: figures, tables, KML, GeoTIFF
sources/         Reference papers
src/             Analysis scripts (above)
unused/          Earlier exploratory scripts and notebooks
```

---

## Key references

- Gvirtzman, S., Kammer, D. S., Adda-Bedia, M. & Fineberg, J. (2025). How frictional ruptures and earthquakes nucleate and evolve. *Nature* 637, 369–374. https://doi.org/10.1038/s41586-024-08287-y
- Liston, G. E. & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). *Journal of Hydrometeorology* 7, 217–234.
- Pomeroy, J. W. & Gray, D. M. (1990). Saltation of snow. *Water Resources Research* 26, 1583–1594.
- Li, L. & Pomeroy, J. W. (1997). Estimates of threshold wind speeds for snow transport using meteorological data. *Journal of Applied Meteorology* 36, 205–213.
- Winstral, A., Elder, K. & Davis, R. E. (2002). Spatial snow modeling of wind-redistributed snow using terrain-based parameters. *Journal of Hydrometeorology* 3, 524–538.
- Pingel, T. J., Clarke, K. C. & McBride, W. A. (2013). An improved simple morphological filter for the terrain classification of airborne LIDAR data. *ISPRS Journal of Photogrammetry and Remote Sensing* 77, 21–30.
- Colorado Avalanche Information Center (CAIC) highway observer database and avalanche records.
