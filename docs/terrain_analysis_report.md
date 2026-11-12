# Wind Loading, Local Shear Stress and Natural Slab Release
### Seven Sisters, Loveland Pass — internal research note (November 2026)

---

## Summary

- **Storm by storm, a terrain + wind + snowfall stress model predicts natural slab days.** The pre-registered test (72-h load, P99 weak-layer stress margin) gives **AUC 0.73 (95 % CI 0.64–0.83, permutation p = 0.001)** over 32 natural events in 16 seasons. On days with naturals it ranks the released path above the other paths more often than chance (AUC 0.61, p = 0.047).
- **The result holds out of sample.** The model was built on 2009–2021 data; the 17 natural events of 2021–2026 (a 2026 CAIC export, matched to paths by coordinates) score **AUC 0.81 (0.68–0.91, p = 0.001)**.
- **Snowfall sets the timing; terrain and wind pick the path.** Snowfall alone predicts natural days nearly as well (AUC 0.71) but is the same for every path, so it cannot say which path releases.
- **Steepness alone does not explain where naturals happen.** Shear stress per unit load (∝ sinθ cosθ) is almost identical across the five start zones (0.484–0.490, all near 45°), so differences in shear stress must come from differences in load — wind deposition.
- **Wind loading separates the zones, but the transport model is too weak.** Both loading models put Sister #1 and #6 as the most wind-loaded zones (the MicroMet ranking never changes across sensitivity runs), and Sister #6's naturals follow SW-dominated transport (58 % of flux vs 31 % for #1–#4). But MicroMet's saltation-only deposition is small next to direct snowfall, so its zone-level *expected* release counts badly miss Sister #6 (8 observed vs ≈ 1 expected) and Sister #2 (1 vs 8–12).
- **Hotspots match field observations.** On Sister #1, modelled deposition and nucleation potential concentrate on the west (skier's-left) side of the gully, where natural releases have been seen to start in windy periods; the 2020–2026 release points inside start zones sit at above-median potential within their zone (mean 62nd percentile, n = 5, p = 0.18).
- **Zone-level correlations remain non-significant** (five zones; minimum attainable p = 0.017). The storm-level and out-of-sample results are the strongest evidence; the zone-level loading magnitude is the main open problem.

---

## 1. Motivation and framework

Natural slab avalanches on the Seven Sisters paths threaten US-6 over Loveland Pass. All paths share the same storms and are controlled together, yet some release naturally far more often than others. The question is which terrain and loading conditions make a path release on its own.

**Inspiration.** Gvirtzman, Kammer, Adda-Bedia & Fineberg (2025), *How frictional ruptures and earthquakes nucleate and evolve* (Nature 637, 369–374), show in laboratory friction experiments that rupture nucleation is governed by a **local stress threshold**: a slowly expanding patch starts wherever the local shear stress exceeds τ_thresh, and its time to dynamic rupture shortens as the excess stress grows (T ≈ α / (τ − τ_thresh)), so the nucleation rate scales with τ − τ_thresh. The paper does not discuss snow; the application to slabs is this project's analogy. *(An earlier version of this note mis-cited the paper as a snow-avalanche study; that is corrected here.)*

**What carries over to slabs.** In the paper, a nucleating patch transitions to classical fracture when it spans the confining width of the interface. For snow slabs, the critical crack length (tens of centimetres) is far smaller than any start zone's width (tens of metres), so a patch becomes unstable long before it could feel the path boundaries; path width is not expected to matter. What does carry over is that **nucleation is local**: what matters is whether, how often and by how much parts of a start zone exceed the threshold — not the zone's average condition. This note therefore maps stress cell by cell and looks for **hotspots** at or above the critical crack length.

**Stress model.** For a slab of vertical snow load (SWE) on slope θ:

- shear stress τ = g·SWE·sinθ cosθ, normal stress σ = g·SWE·cos²θ, so τ/σ = tanθ;
- for a friction-type weak layer (τ > c + μσ), the stress margin is **g·SWE·cos²θ·(tanθ − μ)**, used throughout with μ = tan 30°.

Load raises the margin only where tanθ > μ, and steepness sets how much each unit of load counts. Because sinθ cosθ is nearly flat between ~37° and 45°, **differences in shear stress between these zones must come mostly from differences in load — that is, from wind deposition.**

---

## 2. Site

Seven Sisters is a set of avalanche paths on the north side of the ridge above US-6, west of the Loveland Pass summit (start zones ≈ 39.677° N, 105.891° W, 3,460–3,600 m). CDOT controls the paths with explosives (hand charges and Avalauncher); every mission includes all paths. Five start zones are analysed: Sisters #1, #2, #3, #4 and #6, which are all the active paths. There is no Sister #5 on the ground (the two CAIC records under that name are presumably misattributions), and Sister #7 is not an active path. Sisters #2–#4 start among spruce forest; Sisters #1 and #6 start in open alpine terrain (Figure 1).

| Zone | Area (m²) | Natural HS/SS 2009–2026 | per 10,000 m² | Slab on control days (2009–2021) |
|---|---|---|---|---|
| Sister #1 | 2,896 | 10 | 34.5 | 67 % |
| Sister #2 | 2,884 | 1 | 3.5 | 45 % |
| Sister #3 | 2,384 | 12 | 50.3 | 61 % |
| Sister #4 | 1,416 | 4 | 28.2 | 55 % |
| Sister #6 | 584 | 9 | 154.1 | 53 % |

*Control-day column: share of the 198 recorded control days on which explosives released a slab on that path.*

Because all paths are shot on every mission, mitigation affects the paths equally, so natural counts are comparable across paths.

![Start-zone slope](../figures/seven_sisters_sz_slope_map.png)

*Figure 1. Slope inside the five start zones on USGS NAIP imagery (0.6 m, 2019-09-22), with 10 m contours from 3DEP (north up; slopes face north, so the top of each start zone is its southern edge).*

---

## 3. Data

**Terrain.** USGS 3DEP 1 m bare-earth lidar DEM (DRCOG 2020; tile `USGS_1M_13_x42y440_CO_DRCOG_2020_B20`), resampled to 1.2 m (cubic) for the start-zone analysis and 2 m for the wind models. A 5 cm drone photogrammetry surface (file dated 2026-09-22; RTK GCPs surveyed 2026-10-20) includes forest canopy (trees up to ~20 m) and is used only to map canopy (drone surface minus 3DEP). After removing a constant ~11.9 m vertical datum offset, the two agree to 0.35 m (NMAD) on open ground with no horizontal shift. Map backgrounds use USGS NAIP imagery (0.6 m, 2019-09-22).

> Averaging 1 m data onto a 1.2 m grid leaves a periodic 5-pixel ripple in curvature (a 10× spectral spike); cubic resampling removes it. Per-zone slope is computed on the full DEM and then masked: computing it on clipped zone DEMs biases edge cells low.

**Start zones.** Hand-digitised KML polygons; Sister #1 was re-digitised in November 2026 (now 2,901 m², centred 40 m further north, i.e. downslope, overlapping the old polygon by 1,126 m²). The combined start-zone polygon is the union of the per-zone polygons.

**Avalanches.** Natural (N) hard and soft slabs (HS, SS), from two sources:

- the named CAIC highway observer export, 2009 – 2021-05-04 (re-exported 2026);
- a 2026 CAIC export of natural slabs with coordinates but no path names (2020 – 2026-03). Records were matched to paths by tracing each start zone downhill on the 3DEP DEM and assigning each record to the nearest traced corridor within 40 m (ties → the path aligned across the slope). In the overlap period, every record found in both exports had the same path, size and type. 24 records matched the Sisters; 19 were new and were added (`src/assign_avalanche_paths.py`, `src/merge_avalanche_records.py`). Two overlapping events are logged as explosive-triggered in the highway export and natural in the 2026 export (2020-11-14 Sister #6, 2021-02-06 Sister #1); the highway export's trigger is kept.

Below-treeline records on Sister #6 are excluded as likely runout misattributions; Sister #6's TL band is alpine and counts with >TL. The 2026 export contains naturals only, so control missions are known only to 2021-05-04.

**Weather** (CAIC database, hourly, Oct 2009 – Jun 2026):

- **CALVP** — CAIC/CDOT ridge station, 3,624 m, on the crest ~380 m SW of the start zones (wind speed, direction, temperature).
- **LBAC2** — Loveland Basin SNOTEL (SWE, in hundredths of an inch). Water year 2011 is missing at both stations and 2019 has half the wind record; snow depth carries a flag value (2000) through 2010–14 and is not used.

---

## 4. Methods

### 4.1 Zone steepness

Mean tanθ, share of cells > 40° and mean sinθ cosθ per zone (Horn slope on the full 1.2 m DEM). Slab-scale sensitivity: the DEM smoothed to 5, 10 and 20 m.

### 4.2 Transport winds

Transport hours at CALVP: wind ≥ 15 mph (6.7 m/s), air temperature below 0 °C and ≥ 0.3 in SWE gain at LBAC2 over the previous 72 h.

### 4.3 Fetch (Winstral)

For each wind direction: Winstral's upwind shelter index Sx (100 m search, ±15° sector) and the upwind source area of each zone — rays marched up to 1 km upwind, stopping at a slope-break snow trap (Sb > 15°, 50 m separation), at forest, or at the data edge (traps within 50 m of the zone ignored: a zone just below a crest is loaded through that crest's own lee break). **Fetch ratio** = source area / zone area.

### 4.4 MicroMet terrain winds and saltation transport

1. **Terrain winds** (Liston & Elder 2006): W = 1 + 0.5·Ωs + 0.5·Ωc (Ωs slope in the wind direction, Ωc curvature over η = 30 m, each scaled to ±0.5); local speed = CALVP speed × W(cell) / W(station).
2. **Threshold** (10-m speed): 5 m/s with ≥ 2.5 mm SWE in the previous 24 h; otherwise Li & Pomeroy (1997) dry snow, 9.43 + 0.18T + 0.0033T² m/s; no transport above 0 °C.
3. **Capacity**: Pomeroy & Gray (1990) saltation flux. **Suspension is omitted, so deposition is a lower bound in strong winds.**
4. **Erosion / deposition**: downwind marching; flux relaxes toward rising capacity over Lf = 100 m and deposits where capacity falls (lees, gully floors, forest).

Every hour of record is run (kg/m² = mm SWE). Sensitivity: η = 15, 60 m and Lf = 50, 200 m.

### 4.5 Stress hotspots and nucleation potential

Stress margin per cell with load = snowfall + wind deposition, reduced to zero between 50° and 60° (snow does not build slabs on very steep rock).

- **Hotspots**: top 5 % of margin among all start-zone cells pooled (1 %, 10 % as checks), in connected patches; two loading variants — the fetch-index geometry (deposition weighted equally with snowfall) and MicroMet season-average loading.
- **Nucleation potential** (following the paper, rate ∝ τ − τ_thresh): for every day, the 72-h stress margin minus a threshold (P99 of start-zone margins; P95 as a check), positive part, summed over all days and expressed per season. **Expected release days** per zone are proportional to the zone's days with a hotspot patch above threshold or to its summed excess, scaled to the observed total. Forest cells are masked in the map.

### 4.6 Storm-by-storm test

For every day November–May, each zone's load is snowfall plus MicroMet wind accumulation in the preceding 72 h (variants: load since the last control mission, evaluated only to 2021-05-04, with and without 5-day settling). Outcome: a natural slab on that path that day. **Primary test, chosen before the results:** 72-h window, P99 stress margin. Statistics: AUC with a cluster bootstrap 95 % CI (event days and seasons resampled), a date-shift permutation test (event days moved ±10–60 days within the season), a within-day permutation test (is the released path ranked above the other paths that day?), and an **out-of-sample split**: 2009–2021 (development) vs 2021–2026 events, which played no part in building the model. Zone-level correlations use exact permutation p-values (minimum two-sided p 0.017 with N = 5, 0.083 with N = 4).

---

## 5. Results

### 5.1 Steepness alone

| Zone | Mean slope | % cells > 40° | Mean tanθ (τ/σ) | Mean sinθ cosθ (τ per load) |
|---|---|---|---|---|
| Sister #1 | 41.4° | 56.6 | 0.895 | 0.490 |
| Sister #2 | 40.6° | 52.2 | 0.865 | 0.490 |
| Sister #3 | 44.7° | 79.9 | 1.010 | 0.490 |
| Sister #4 | 45.1° | 78.9 | 1.035 | 0.487 |
| Sister #6 | 38.9° | 41.9 | 0.813 | 0.484 |

Shear stress per unit load varies by about 1 % between zones. Steepness tracks natural releases poorly: mean tanθ vs natural count gives ρ = 0.4 for Sisters #1–#4 and 0.2 with Sister #6 (exact p ≥ 0.75). Sister #4's steepness is mostly small-scale relief (smoothed to 20 m its share of cells > 40° drops from 79 % to 55 %), and Sister #6 — the least steep zone — has nine naturals.

### 5.2 Transport winds and the weather before naturals

Snow-transporting winds at the ridge are bimodal: **WNW (≈ 290–300°)** and **WSW (≈ 245–250°)**.

![Transport wind rose](../figures/seven_sisters_wind_rose.png)

*Figure 2. CALVP transport winds (direction from).*

Before Sister #6's naturals (10 events), **58 %** of transport flux came from the SW, against **31 %** before Sisters #1–#4's naturals (27 events) and 42 % on all winter days; #1–#4 naturals followed WNW-dominated transport (54 % vs 47 %). Sisters #1–#4 released after larger storms (median 72-h new snow at the 84th percentile of winter days), Sister #6 after ordinary ones (53rd percentile) — consistent with Sister #6 being loaded mainly by wind.

### 5.3 Loading

**Fetch.** Fetch ratios rank #6 (59) ≫ #1 (11) > #3 (3.5) ≈ #2 (3.1) > #4 (1.3), robust to the trap / grace / length settings tested (only the order is interpreted).

**MicroMet.**

| Zone | Wind accumulation, mm SWE / season (mean) | of which SW | of which W–NW | 90th percentile |
|---|---|---|---|---|
| Sister #1 | **+88** | +78 | +12 | +358 |
| Sister #6 | **+23** | +31 | −8 | +212 |
| Sister #2 | −55 | −8 | −44 | −23 |
| Sister #3 | −77 | −20 | −54 | +96 |
| Sister #4 | −94 | −33 | −57 | +1 |

Across sensitivity runs (η 15/30/60 m, Lf 50/100/200 m) the order **#1 > #6 > #2 > #3 > #4 never changes**. In the model, forest upwind of Sisters #2–#4 intercepts drifting snow (the imagery confirms these zones start among trees). On Sister #1, deposition concentrates in a band along the **west side of the gully**, under both SW and W–NW winds; Sister #6 is loaded along its west edge by SW winds.

![Loading patterns](../figures/seven_sisters_loading_patterns.png)

*Figure 3. MicroMet / saltation wind accumulation per season (blue deposition, red erosion, green forest), all transport and by sector. North up.*

![Loading close-ups](../figures/seven_sisters_loading_closeups.png)

*Figure 4. Close-ups of Sister #1 and Sister #6 by wind sector (north up; top of each start zone at the bottom).*

### 5.4 Stress hotspots

| Zone | Fetch-index load: hotspot area / patches | MicroMet load: hotspot area / patches | Naturals |
|---|---|---|---|
| Sister #3 | 124 m² / 4 | 140 m² / 2 | 12 |
| Sister #1 | 295 m² / 11 | 284 m² / 4 | 10 |
| Sister #6 | 26 m² / 2 | 16 m² / 1 | 9 |
| Sister #4 | 17 m² / 3 | 24 m² / 1 | 4 |
| Sister #2 | 4 m² / 1 | 20 m² / 1 | 1 |

With fetch-index loading, hotspot area orders the zones like the natural counts (ρ = 0.90, exact p = 0.083); with MicroMet loading the correlation is weaker (ρ = 0.60, p = 0.35), because MicroMet's deposition is small relative to direct snowfall, so its hotspots mostly follow slope.

**Field check — Sister #1.** In both versions Sister #1's hotspots lie on the west side of the gully (88 % and 99 % of hotspot area in the west half), matching where its naturals have been seen to start. This is consistent but not discriminating: most zones have most of their hotspots on their west side, as expected for prevailing W–WNW winds loading east-facing walls.

![Stress hotspots](../figures/seven_sisters_stress_hotspots.png)

*Figure 5. Fetch-index deposition (left) and stress margin with top-5 % hotspot patches in red (right). North up.*

### 5.5 Nucleation potential and observed vs expected

The nucleation-potential map (Figure 6) sums the excess stress above threshold over every day of record. Within start zones it is highest along the west side of Sister #1, on the west side of Sister #3 and at the top and west edge of Sister #6. Of the 2020–2026 release points, five lie inside start zones; they sit at the 51st–81st percentile of potential within their zone (mean 62nd; random placement would give 50th, p = 0.18; two pairs share coordinates, so they are not fully independent).

At the zone level, **expected counts from MicroMet-based potential do not match observations** (ρ = −0.36 to 0.31): Sister #6 has 8 observed natural days against ≈ 0.6–2 expected, and Sister #2 has 1 against 7–12. Sister #2 is large and moderately steep, so a slope-dominated model expects many releases; Sister #6 is small and less steep, so it expects few. The observations say the opposite — the wind-loading contrast (Sister #2 forest-sheltered, Sister #6 a SW-fed lee) matters more than the MicroMet saltation model credits.

![Nucleation potential](../figures/seven_sisters_nucleation_potential.png)

*Figure 6. Left: nucleation potential (Σ excess stress above the P99 threshold per season; forest masked) on NAIP imagery, with natural release points (stars inside start zones, circles on tracks). Right: observed natural-release days (95 % Poisson CI) vs expected from hotspot days and from summed excess, scaled to the observed total.*

### 5.6 Storm-by-storm skill

| Test | Events | AUC (95 % CI) | Date-shift p | Within-day AUC (p) |
|---|---|---|---|---|
| **Primary: 72-h, P99 stress margin** | 32 | **0.73 (0.64–0.83)** | **0.001** | **0.61 (0.047)** |
| Primary, 2009–2021 (development) | 15 | 0.65 (0.53–0.79) | 0.074 | 0.60 (0.14) |
| **Primary, 2021–2026 (out of sample)** | 17 | **0.81 (0.68–0.91)** | **0.001** | 0.62 (0.11) |
| 72-h snowfall only | 32 | 0.71 (0.62–0.80) | < 0.001 | 0.50 (—) |
| 72-h wind accumulation only | 32 | 0.57 (0.42–0.71) | 0.80 | 0.61 (0.059) |
| Since control, 5-d settling (to 2021) | 15 | 0.65 (0.52–0.79) | 0.019 | 0.58 (0.21) |

*Sensitivity runs: primary AUC 0.730–0.736.*

- **Natural days are predictable from the 72-h stress margin**, and the skill holds on data that played no part in building the model.
- **Timing is mostly snowfall; path selection needs terrain and wind.** On natural days the stress margin and wind accumulation both rank the released path above the others; snowfall alone cannot.
- **Recent loading matters.** Load accumulated since the last control mission without settling does worse; with 5-day settling it is comparable. Natural releases follow loading from the last few days, consistent with older load gaining strength and with the paper's excess-stress / time-to-failure relation.

---

## 6. Limitations

- **Zone count.** Five zones, all the active Sisters paths; no zone-level correlation can be strongly significant and none observed is below 0.08. More zones must come from other sites.
- **Avalanche record.** Many natural dates are estimated. Control missions without results are not recorded, and none are available after 2021-05-04. Records from 2020–2026 are assigned to paths by coordinates (observer-placed points, not exact release locations).
- **Snow data.** SNOTEL SWE resolves 2.5 mm steps; one station represents the whole site; water year 2011 is missing.
- **Wind and transport model.** MicroMet is empirical; saltation only (suspension omitted); unlimited erodible snow; no wind-direction diversion in the transport direction; faint streak artefacts along the wind direction. The zone order is robust to model parameters, but the deposition magnitude is too weak to reproduce zone-level release counts.
- **Deposition not validated against snow.** The canopy mask (drone-derived) affects Sisters #2–#4.
- **Thresholds.** μ = tan 30°, the 50–60° holding cap, the P99 threshold and hotspot percentiles are assumptions; alternatives are reported where tested.

---

## 7. Next steps

1. **Snow-on drone survey** after a W–NW storm and after a SW storm: snow depth = snow-on minus snow-off, processed with the same GCPs and vertical datum (the drone products sit ~11.9 m below 3DEP). Priority targets: Sister #1's west wall and Sister #6's lee. This measures the load directly and calibrates the transport model.
2. **Add suspension transport** (e.g. Pomeroy & Male 1992 or SnowTran-3D); this is what the observed-vs-expected mismatch points to.
3. **WindNinja** runs for a few key storms to check MicroMet's terrain winds near the start zones.
4. **More zones**: extend the analysis to other controlled paths with ridge wind stations.
5. **Control-mission log** from CDOT (including shots with no result, and after 2021) for the since-control analysis.
6. **Pipeline**: replace the legacy nucleation-candidate map with the nucleation-potential model.

---

## Appendix — reproducing

```bash
source .envrc                                                     # conda env (.conda)
python src/assign_avalanche_paths.py data/avalanches/caic-avalanche-records-2026-11-12.csv
python src/merge_avalanche_records.py \
    data/avalanches/CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv \
    figures/seven_sisters_new_records_assigned.csv \
    data/avalanches/CAIC_HWY_avalanches_2009-2026_sisters_merged.csv
python src/sql_util.py --site sisters --end 2026-06-01            # CALVP + LBAC2
python src/pipeline.py --site sisters                             # DEM, zone clips, metrics, Figure 1
python src/wind_exposure.py --site sisters                        # Sx, fetch (§4.3)
python src/storm_analysis.py --site sisters                       # weather before each natural
python src/stress_hotspots.py --site sisters                      # fetch-index hotspots
python src/micromet_loading.py --site sisters                     # MicroMet, potential, storm test
python src/event_stats.py --site sisters                          # bootstrap / permutation / out-of-sample
```

Profile curvature (an earlier "stress concentration" hypothesis) is still computed by the pipeline but is not used here; aspect entropy has been dropped, as explicit wind modelling now addresses the loading question it was a proxy for.

## References

- Gvirtzman, S., Kammer, D. S., Adda-Bedia, M. & Fineberg, J. (2025). How frictional ruptures and earthquakes nucleate and evolve. *Nature* 637, 369–374. https://doi.org/10.1038/s41586-024-08287-y
- Heierli, J., Gumbsch, P. & Zaiser, M. (2008). Anticrack nucleation as triggering mechanism for snow slab avalanches. *Science* 321, 240–243.
- Gaume, J., Gast, T., Teran, J., van Herwijnen, A. & Jiang, C. (2018). Dynamic anticrack propagation in snow. *Nature Communications* 9, 3047.
- Li, L. & Pomeroy, J. W. (1997). Estimates of threshold wind speeds for snow transport using meteorological data. *Journal of Applied Meteorology* 36, 205–213.
- Liston, G. E. & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). *Journal of Hydrometeorology* 7, 217–234.
- Pomeroy, J. W. & Gray, D. M. (1990). Saltation of snow. *Water Resources Research* 26, 1583–1594.
- Schweizer, J., Jamieson, J. B. & Schneebeli, M. (2003). Snow avalanche formation. *Reviews of Geophysics* 41, 1016.
- Winstral, A., Elder, K. & Davis, R. E. (2002). Spatial snow modeling of wind-redistributed snow using terrain-based parameters. *Journal of Hydrometeorology* 3, 524–538.
- Colorado Avalanche Information Center (CAIC) highway observer database, 2009–2021, and CAIC avalanche records export, 2026.
- USGS 3D Elevation Program (3DEP), 1 m DEM, DRCOG 2020 lidar; USGS NAIP imagery (2019).
