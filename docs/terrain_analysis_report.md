# Wind Loading and Natural Slab Release on the Seven Sisters Paths
### Which paths release naturally, and why — internal research note (November 2026)

This is the first of two notes. It covers **which paths release naturally and why**: terrain, wind loading, stress hotspots and the storm-by-storm test. The [companion note](release_mechanism_note.md) covers **when they release and how a natural slab might start**: the stress threshold, loading rate, hour-by-hour stress before timed releases, avalanche size, and whether natural slabs start in shear.

---

## Summary

- **Loading, not terrain, separates the paths.** Shear stress per unit load (∝ sinθ cosθ) is almost identical across the five start zones (0.484–0.490), and zone steepness does not track natural releases (ρ = 0.2). What differs is where wind-transported snow can land.
- **Wind landing sites match the release record and field observation.** Sisters #1, #3, #4 and #6 have steep open cells where W–NW (and, for #1 and #6, SW) winds deposit snow; **Sister #2 has none** (0 m² of W–NW landing area), as observed in the field — and it has one natural slab in 17 seasons.
- **A terrain + wind + snowfall stress model reproduces the zone pattern.** Expected natural-release days per zone (from days with a stress hotspot above threshold) are 6.4 / 1.8 / 9.0 / 4.0 / 10.9 against 8 / 1 / 12 / 3 / 8 observed (Sisters #1 / #2 / #3 / #4 / #6; ρ = 0.82), stable across model settings.
- **Storm by storm, the same model predicts natural slab days and which path releases.** The pre-registered test (72-h load, P99 weak-layer stress margin) gives **AUC 0.76 (95 % CI 0.67–0.84, p = 0.001)** over 32 natural events in 16 seasons, and on natural days it ranks the released path above the others (**AUC 0.74, p < 0.001**). Snowfall sets the timing (AUC 0.71 alone) but cannot pick the path; wind loading does.
- **The 2021–2026 events score as well as the development period** (AUC 0.82). Two model changes — a redrawn Sister #6 start zone and a forest-edge correction — were made after those events were available, though not tuned on them; before the changes the 2021–2026 AUC was already 0.81.
- **W–NW loading of hotspot terrain tracks the paths' natural counts.** Days per season on which W–NW winds add ≥ 10 Pa to a hotspot patch rank the five paths exactly as their natural days (ρ = 0.97, exact p = 0.033, the minimum attainable); SW loading days do not (ρ ≈ 0.7). Sister #2 gets no W–NW loading days at all.
- **When paths release — the stress threshold, loading rate, timed releases, avalanche size, and whether natural slabs start in shear — is in the [companion note](release_mechanism_note.md).**
- **Remaining caveats:** five zones (no zone-level correlation reaches significance), modelled rather than measured loading, and saltation-only transport. A snow-on drone survey is the most valuable next step.

---

## 1. Motivation and framework

Natural slab avalanches on the Seven Sisters paths threaten US-6 over Loveland Pass. All paths share the same storms and are controlled together, yet some release naturally far more often than others. The question is which terrain and loading conditions make a path release on its own.

**Inspiration.** Gvirtzman, Kammer, Adda-Bedia & Fineberg (2025), *How frictional ruptures and earthquakes nucleate and evolve* (Nature 637, 369–374), show in laboratory friction experiments that rupture nucleation is governed by a **local stress threshold**: a slowly expanding patch starts wherever the local shear stress exceeds τ_thresh, and its time to dynamic rupture shortens as the excess stress grows (T ≈ α / (τ − τ_thresh)), so the nucleation rate scales with τ − τ_thresh. The paper does not discuss snow; the application to slabs is this project's analogy. *(An earlier version of this note mis-cited the paper as a snow-avalanche study; that is corrected here.)*

**What carries over to slabs.** In the paper, a nucleating patch transitions to classical fracture when it spans the confining width of the interface. For a slab, the transition comes when the patch reaches the critical crack length. Confinement shapes the nucleation phase up to that point; after it, classical fracture mechanics apply. Path width (tens of metres) is far larger than the critical crack length (tens of centimetres to a few slab thicknesses), so path boundaries are not expected to matter. What does carry over is that **nucleation is local**: what matters is whether, how often and by how much parts of a start zone exceed the threshold — not the zone's average condition. This note therefore maps stress cell by cell and looks for **hotspots** at or above the critical crack length.

**Stress model.** For a slab of vertical snow load (SWE) on slope θ:

- shear stress τ = g·SWE·sinθ cosθ, normal stress σ = g·SWE·cos²θ, so τ/σ = tanθ;
- for a friction-type weak layer (τ > c + μσ), the stress margin is **g·SWE·cos²θ·(tanθ − μ)**, used throughout with μ = tan 30°.

Load raises the margin only where tanθ > μ. Because sinθ cosθ is nearly flat between ~37° and 45°, **differences in shear stress between these zones must come mostly from differences in load — that is, from wind deposition.**
---

## 2. Site

Seven Sisters is a set of avalanche paths on the north side of the ridge above US-6, west of the Loveland Pass summit (start zones ≈ 39.677° N, 105.891° W, 3,460–3,620 m). CDOT controls the paths with explosives (hand charges and Avalauncher); every mission includes all paths. Five start zones are analysed: Sisters #1, #2, #3, #4 and #6, which are all the active paths. There is no Sister #5 on the ground (the two CAIC records under that name are presumably misattributions), and Sister #7 is not an active path. Sisters #2–#4 start among spruce forest; Sisters #1 and #6 start in open alpine terrain (Figure 1).

| Zone | Area (m²) | Natural HS/SS 2009–2026 | Slab on control days (2009–2021) |
|---|---|---|---|
| Sister #1 | 2,896 | 10 | 67 % |
| Sister #2 | 2,884 | 1 | 45 % |
| Sister #3 | 2,384 | 12 | 61 % |
| Sister #4 | 1,416 | 4 | 55 % |
| Sister #6 | 3,840 | 9 | 53 % |

Natural counts are compared per path, not per unit area: what matters is whether a zone has hotspots and whether they get loaded, not what share of the zone they occupy.

*Control-day column: share of the 198 recorded control days on which explosives released a slab on that path.*

Because all paths are shot on every mission, mitigation affects the paths equally, so natural counts are comparable across paths.

**Field observations.** Sister #1's naturals have been seen to start on the west (skier's-left) side of the gully below the top of the start zone, during windy periods. Sister #2 is exposed mainly to W–NW winds with no place for wind-transported snow to land, and sits too low on the slope to be affected by SW winds; Sisters #1, #3 and #4 have locations where W–NW winds can load.

![Start-zone slope](../figures/seven_sisters_sz_slope_map.png)

*Figure 1. Slope inside the five start zones on USGS NAIP imagery (0.6 m, 2019-09-22), with 10 m contours from 3DEP (north up; slopes face north, so the top of each start zone is its southern edge).*
---

## 3. Data

**Terrain.** USGS 3DEP 1 m bare-earth lidar DEM (DRCOG 2020; tile `USGS_1M_13_x42y440_CO_DRCOG_2020_B20`), resampled to 1.2 m (cubic) for the start-zone analysis and 2 m for the wind models. A 5 cm drone photogrammetry surface (file dated 2026-09-22; RTK GCPs surveyed 2026-10-20) includes forest canopy (trees up to ~20 m) and is used only to map canopy (drone surface minus 3DEP). After removing a constant ~11.9 m vertical datum offset, the two agree to 0.35 m (NMAD) on open ground with no horizontal shift. Map backgrounds use USGS NAIP imagery (0.6 m, 2019-09-22).

> Averaging 1 m data onto a 1.2 m grid leaves a periodic 5-pixel ripple in curvature (a 10× spectral spike); cubic resampling removes it. Per-zone slope is computed on the full DEM and then masked: computing it on clipped zone DEMs biases edge cells low.

**Start zones.** Hand-digitised KML polygons. Sister #1 was re-digitised in November 2026 (now 2,901 m², centred 40 m further north, i.e. downslope, overlapping the old polygon by 1,126 m²). **Sister #6 was redrawn in November 2026** (now 3,840 m², up from 579 m²; it contains the old polygon and extends ~31 m north and ~15 m west, taking in the upper west edge where SW loading concentrates). The combined start-zone polygon is the union of the per-zone polygons.

**Avalanches.** Natural (N) hard and soft slabs (HS, SS), from two sources:

- the named CAIC highway observer export, 2009 – 2021-05-04 (re-exported 2026);
- a 2026 CAIC export of natural slabs with coordinates but no path names (2020 – 2026-03). Records were matched to paths by tracing each start zone downhill on the 3DEP DEM and assigning each record to the nearest traced corridor within 40 m (ties → the path aligned across the slope). In the overlap period, every record found in both exports had the same path, size and type. 24 records matched the Sisters; 19 were new and were added (`src/assign_avalanche_paths.py`, `src/merge_avalanche_records.py`). Two overlapping events are logged as explosive-triggered in the highway export and natural in the 2026 export (2020-11-14 Sister #6, 2021-02-06 Sister #1); they are kept as explosive-triggered pending verification, since each would add one natural to a small dataset.

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
### 4.3 MicroMet terrain winds and saltation transport

1. **Terrain winds** (Liston & Elder 2006): W = 1 + 0.5·Ωs + 0.5·Ωc (Ωs slope in the wind direction, Ωc curvature over η = 30 m, each scaled to ±0.5); local speed = CALVP speed × W(cell) / W(station).
2. **Threshold** (10-m speed): 5 m/s with ≥ 2.5 mm SWE in the previous 24 h; otherwise Li & Pomeroy (1997) dry snow, 9.43 + 0.18T + 0.0033T² m/s; no transport above 0 °C.
3. **Capacity**: Pomeroy & Gray (1990) saltation flux (suspension omitted).
4. **Erosion / deposition**: downwind marching; flux relaxes toward rising capacity over Lf = 100 m and deposits where capacity falls (lees, gully floors, forest).

Every hour of record is run (kg/m² = mm SWE). Sensitivity: η = 15, 60 m and Lf = 50, 200 m.
### 4.4 Landing sites, stress hotspots and nucleation potential

- **Landing sites**: steep (35–55°) open cells collecting more than 50 mm SWE per season of wind deposition, by wind sector.
- **Stress margin** per cell with load = snowfall + wind deposition, reduced to zero between 50° and 60° (snow does not build slabs on very steep rock). **Cells within 4 m of forest canopy carry no slab stress**: trees anchor the snowpack, and the transport model dumps the intercepted flux at the tree line. *(This correction was found while investigating Sister #2: its few extreme-deposition cells, holding 70 % of its deposition, all sat at forest edges — not open slab terrain.)*
- **Hotspots**: top 5 % of margin among all start-zone cells pooled, in connected patches (MicroMet season-average loading; a fetch-index geometry is a robustness check).
- **Nucleation potential** (following the paper, rate ∝ τ − τ_thresh): for every day, the 72-h stress margin minus a threshold (P99 of start-zone margins; P95 as a check), positive part, summed over all days per season. **Expected release days** per zone are proportional to the zone's days with a hotspot patch (≥ 2 cells above threshold) or to its summed excess, scaled to the observed total.
- **Wind-deposition multiplier**: load = snowfall + k · wind deposition, k = 0–50, threshold recomputed for each k — does the zone pattern need stronger (or no) wind loading?
### 4.5 Storm-by-storm test

For every day November–May, each zone's load is snowfall plus MicroMet wind accumulation in the preceding 72 h (variants: load since the last control mission, evaluated only to 2021-05-04, with and without 5-day settling; stress on wind-loaded cells only). Outcome: a natural slab on that path that day. **Primary test, chosen before the results:** 72-h window, P99 stress margin. Statistics: AUC with a cluster bootstrap 95 % CI (event days and seasons resampled), a date-shift permutation test (event days moved ±10–60 days within the season), a within-day permutation test (is the released path ranked above the other paths that day?), and a period split: 2009–2021 (development) vs 2021–2026. Zone-level correlations use exact permutation p-values (minimum two-sided p 0.017 with N = 5).

---

## 5. Results

### 5.1 Steepness alone

| Zone | Mean slope | % cells > 40° | Mean tanθ (τ/σ) | Mean sinθ cosθ (τ per load) |
|---|---|---|---|---|
| Sister #1 | 41.4° | 56.6 | 0.895 | 0.490 |
| Sister #2 | 40.6° | 52.2 | 0.865 | 0.490 |
| Sister #3 | 44.7° | 79.9 | 1.010 | 0.490 |
| Sister #4 | 45.1° | 78.9 | 1.035 | 0.487 |
| Sister #6 | 40.5° | 56.4 | 0.862 | 0.488 |

Shear stress per unit load varies by about 1 % between zones, and steepness does not track natural releases (mean tanθ vs count: ρ = 0.2, exact p = 0.78). Sister #2 is as steep as Sisters #1 and #6 yet has one natural against their 10 and 9. Sister #4's steepness is mostly small-scale relief (smoothed to 20 m its share of cells > 40° drops from 79 % to 55 %).
### 5.2 Transport winds and the weather before naturals

Snow-transporting winds at the ridge are bimodal: **WNW (≈ 290–300°)** and **WSW (≈ 245–250°)** (Figure A1).

Before Sister #6's naturals (10 events), **58 %** of transport flux came from the SW, against **31 %** before Sisters #1–#4's naturals (27 events) and 42 % on all winter days; #1–#4 naturals followed WNW-dominated transport (54 % vs 47 %). Sisters #1–#4 released after larger storms (median 72-h new snow at the 84th percentile of winter days), Sister #6 after ordinary ones (53rd percentile) — consistent with Sister #6 being loaded mainly by wind.
### 5.3 Loading

**Fetch** (Winstral source area / zone area; Appendix). Fetch ratios rank #6 (16) > #1 (11) > #3 (3.5) ≈ #2 (3.1) > #4 (1.3) (order robust to the trap / grace / length settings; magnitudes are not).

**MicroMet.**

| Zone | Wind accumulation, mm SWE / season (mean) | of which SW | of which W–NW | 90th percentile |
|---|---|---|---|---|
| Sister #1 | **+85** | +75 | +11 | +354 |
| Sister #6 | **+40** | +53 | −12 | +227 |
| Sister #2 | −52 | −6 | −43 | −19 |
| Sister #3 | −79 | −20 | −56 | +107 |
| Sister #4 | −86 | −25 | −58 | +24 |

Across sensitivity runs (η 15/30/60 m, Lf 50/100/200 m) the order **#1 > #6 > #2 > #3 > #4 never changes**. Zone means hide the pockets that matter: Sisters #3 and #4 are net-eroded but have positive upper percentiles (§5.4).

![Loading patterns](../figures/seven_sisters_loading_patterns.png)

*Figure 2. MicroMet / saltation wind accumulation per season (blue deposition, red erosion, green forest), all transport and by sector. North up.*

![Loading close-ups](../figures/seven_sisters_loading_closeups.png)

*Figure 3. Close-ups of Sister #1 and Sister #6 by wind sector (north up; top of each start zone at the bottom).*
### 5.4 Landing sites — why Sister #2 is quiet

| Zone | W–NW landing area | share of zone | SW landing area | Naturals |
|---|---|---|---|---|
| Sister #1 | 540 m² | 19 % | 1,056 m² | 10 |
| Sister #6 | 312 m² | 8 % | 1,428 m² | 9 |
| Sister #3 | 140 m² | 6 % | 260 m² | 12 |
| Sister #4 | 24 m² | 2 % | 72 m² | 4 |
| **Sister #2** | **0 m²** | **0 %** | 136 m² | **1** |

*Steep (35–55°) open cells (> 4 m from canopy) collecting > 50 mm SWE / season.*

The model independently reproduces the field description: **Sister #2 has no W–NW landing sites** (zero in four of five sensitivity runs), while Sisters #1, #3 and #4 do. The model does give Sister #2 some SW landing area (136 m²), whereas the field observation is that it sits too low on the slope for SW winds — the model may overstate SW exposure at that elevation. Before the forest-edge correction, Sister #2's apparent deposition came from five forest-edge cells; with them excluded its open slab terrain is scoured, not loaded.
### 5.5 Stress hotspots

| Zone | MicroMet load: hotspot area / patches | Fetch-index load: hotspot area / patches | Naturals |
|---|---|---|---|
| Sister #3 | 132 m² / 1 | 105 m² / 3 | 12 |
| Sister #1 | 276 m² / 5 | 251 m² / 9 | 10 |
| Sister #6 | 228 m² / 4 | 266 m² / 8 | 9 |
| Sister #4 | 12 m² / 1 | 9 m² / 1 | 4 |
| Sister #2 | 0 | 0 | 1 |

The fetch-index variant (Winstral source area as the load; Appendix) is a robustness check. Both loading variants put Sister #2 at zero hotspot area and Sister #4 near zero; hotspot area vs natural count gives ρ = 0.70 (exact p = 0.23, N = 5). **Field check — Sister #1:** 99 % of its hotspot area lies in the west half of the zone, matching where its naturals have been seen to start. (Most zones have their hotspots on the west side, as expected for W–WNW winds loading east-facing walls, so this is consistent rather than discriminating.)
### 5.6 Nucleation potential and observed vs expected

| | Sister #1 | Sister #2 | Sister #3 | Sister #4 | Sister #6 |
|---|---|---|---|---|---|
| Observed natural days (with weather data) | 8 | 1 | 12 | 3 | 8 |
| Expected — hotspot days | 6.4 | 1.8 | 9.0 | 4.0 | 10.9 |
| Expected — summed excess | 7.5 | 3.9 | 8.3 | 2.8 | 9.5 |

*P99 threshold; expected values scaled to the observed total of 32. Hotspot-days fit: Poisson deviance 2.8, ρ = 0.82; across sensitivity runs deviance 2.5–5.1, ρ 0.82–0.97, Sister #2 expected 1.5–2.6.*

Every observed count lies within its expected value's Poisson range. This was not true before the two November corrections: with the old Sister #6 polygon and forest-edge cells counted, the model expected about one natural day on Sister #6 (8 observed) and 6–12 on Sister #2 (1 observed).

**Wind multiplier.** With no wind term (k = 0, snowfall × terrain only) the zone pattern fits worse (deviance 7.4–13.3 for the two main bases; restricting stress to wind-loaded cells, which already uses where the wind deposits, fits well even at k = 0, deviance 3.3); with the model's wind deposition as built (k = 1) it fits well (deviance 2.8), and doubling it (k = 2) fits about equally (1.9). Storm-level path picking also peaks at k ≈ 1–2 (within-day AUC 0.64 at k = 0, 0.74 at k = 1–2). The modelled wind loading is therefore roughly the right strength — not an order of magnitude too weak — once deposition is restricted to open slab terrain.

**Release points.** Seven 2020–2026 release points fall inside start zones; they sit at the 0th–81st percentile of potential within their zone (median 58th). One Sister #6 point lies on a zero-potential cell, so this check is inconclusive.

![Nucleation potential](../figures/seven_sisters_nucleation_potential.png)

*Figure 4. Left: nucleation potential (Σ excess stress above the P99 threshold per season; forest and forest-edge cells masked) on NAIP imagery, with natural release points (stars inside start zones, circles on tracks). Right: observed natural-release days on days with weather data — 32 of the 36 naturals; four (2018-04-03 on #1, #4 and #6, and 2018-11-19 on #1) fall in weather-data gaps — with 95 % Poisson CIs, against expected days from hotspot days and from summed excess, scaled to the observed total.*
### 5.7 Storm-by-storm skill

| Test | Events | AUC (95 % CI) | Date-shift p | Within-day AUC (p) |
|---|---|---|---|---|
| **Primary: 72-h, P99 stress margin** | 32 | **0.76 (0.67–0.84)** | **0.001** | **0.74 (< 0.001)** |
| Primary, 2009–2021 | 15 | 0.70 (0.58–0.82) | 0.046 | 0.78 (0.002) |
| Primary, 2021–2026 | 17 | 0.82 (0.69–0.91) | 0.002 | 0.70 (0.017) |
| 72-h snowfall only | 32 | 0.71 (0.62–0.80) | < 0.001 | 0.50 (—) |
| 72-h wind accumulation only | 32 | 0.58 (0.44–0.72) | 0.63 | 0.62 (0.033) |
| Stress on wind-loaded cells only | 32 | 0.73 (0.64–0.82) | 0.002 | 0.70 (0.001) |
| Since control, 5-d settling (to 2021) | 15 | 0.69 (0.58–0.82) | 0.015 | 0.70 (0.018) |

*Sensitivity runs: primary AUC 0.75–0.77, within-day 0.72–0.74.*

- **Natural days are predictable from the 72-h stress margin, and so is the path.** Snowfall alone identifies natural days nearly as well but cannot pick the path; adding terrain and wind loading picks the released path on 74 % of comparisons.
- **The 2021–2026 events score as well as 2009–2021.** They were not used to fit any parameter, but two model changes (the redrawn Sister #6 polygon and the forest-edge rule) were made after they were in hand. Before those changes, the 2021–2026 AUC was 0.81 (within-day 0.62), so the timing skill does not depend on them; the improvement in path picking does.
- **Recent loading matters.** Load accumulated since the last control mission without settling does worse; with 5-day settling it is comparable. Natural releases follow loading from the last few days, consistent with older load gaining strength. The [companion note](release_mechanism_note.md) takes this further: the last 24 h of loading carry the timing signal.
### 5.8 Loading days by wind direction

How often does each wind direction load each path's hotspot terrain? For each path and day, the MicroMet run gives the weak-layer stress that the wind deposition from each 10° direction adds to the path's stress-bearing cells (open, tanθ > μ, < 60°), taken at the 2nd-highest cell (a ≥ 2-cell patch). A **loading day** from a sector is a day on which that sector adds ≥ 10 Pa (≈ 5 mm SWE on a 45° slope) to such a patch (`src/loading_days.py`; 10° bins within a sector summed).

| Path | SW loading days / season | W–NW loading days / season | Natural days |
|---|---|---|---|
| Sister #3 | 51.8 | **46.6** | 12 |
| Sister #6 | 62.4 | **37.7** | 8 |
| Sister #1 | 54.1 | **20.7** | 8 |
| Sister #4 | 42.2 | **7.0** | 3 |
| Sister #2 | 7.7 | **0.0** | 1 |

- **W–NW loading days rank the paths exactly as their natural days do** (ρ = 0.97, exact p = 0.033 — the smallest p attainable with a tie in five paths), at 5, 10 and 20 Pa, and also when the 10° bins are combined by maximum instead of sum (lower bound; ρ = 0.97 at 5–10 Pa, 0.79–0.89 at 20–30 Pa). SW loading days do not (ρ = 0.67–0.82, p ≥ 0.13): Sister #4 gets frequent SW loading but has few naturals.
- By 30° bin, loading from 285–315° tracks natural days best (ρ = 0.97), 255–285° less (0.82), 225–255° least (0.67).
- **Before naturals, the path had been loaded:** 29 of 32 natural days had a loading day on that path in the 72-h window (all zone-days: 57 %). SW loading preceded 22 (69 % vs 40 %, lift 1.7), W–NW loading 14 (44 % vs 24 %, lift 1.8). All 8 Sister #6 natural days followed SW loading (lift 1.8), consistent with its SW exposure; Sister #1's naturals followed W–NW loading twice as often as its other days.
- Per loading day, natural rates are similar across #1, #3 and #6 (0.3–0.5 per 100 days), lower for #4 (0.17) and zero for #2.

This is the zone-level counterpart of §5.4: Sister #2 is the path the W–NW winds never load. The test is exploratory: N = 5, several thresholds and sectors were examined, and the correlation reaches the smallest p the data allow rather than conventional significance by a margin.

![Loading days](../figures/seven_sisters_loading_days.png)

*Figure 5. Left: days per season on which SW and W–NW winds add ≥ 10 Pa to a ≥ 2-cell patch of each path's stress-bearing cells, with natural slab days (diamonds). Right: loading days by 30° wind-from bin.*

---

## 6. Limitations

- **Zone count.** Five zones, all the active Sisters paths; no zone-level correlation reaches significance (best p = 0.08–0.23). More zones must come from other sites.
- **Post-hoc corrections.** The Sister #6 polygon was redrawn and the forest-edge rule added after the zone-level mismatch was seen. Both have physical grounds (the release area as known in the field; trees anchoring snow), but they were not pre-specified.
- **Avalanche record.** Many natural dates are estimated, and two-thirds of 2009–2021 naturals were logged on a control day or the day after, so detection depends on when observers are out. Control missions without results are not recorded, and none are available after 2021-05-04. Records from 2020–2026 are assigned to paths by coordinates (observer-placed points, not exact release locations).
- **Snow data.** SNOTEL SWE resolves 2.5 mm steps; one station represents the whole site; water year 2011 is missing.
- **Wind and transport model.** MicroMet is empirical; saltation only (suspension omitted); unlimited erodible snow; no wind-direction diversion in the transport direction; faint streak artefacts along the wind direction. The zone order is robust to model parameters; absolute amounts are not. The model may overstate SW exposure for the lower zones (Sister #2).
- **Deposition not validated against snow.** The canopy mask (drone-derived) affects Sisters #2–#4.
- **Thresholds.** μ = tan 30°, the 50–60° holding cap, the 4 m forest buffer, the P99 threshold and hotspot percentiles are assumptions; alternatives are reported where tested.
- **Weak-layer state is not modelled.** Loading explains which paths release and on which days, but most loaded days produce no natural ([companion note](release_mechanism_note.md)).

---

## 7. Next steps

1. **Snow-on drone survey** after a W–NW storm and after a SW storm: snow depth = snow-on minus snow-off, processed with the same GCPs and vertical datum (the drone products sit ~11.9 m below 3DEP). Targets: Sister #1's west wall, Sister #6's upper west edge, and Sister #2 (to confirm it has no landing sites). This measures the load directly.
2. **Add suspension transport** (e.g. Pomeroy & Male 1992 or SnowTran-3D) and check whether the zone pattern holds.
3. **WindNinja** runs for a few key storms to check MicroMet's terrain winds near the start zones, especially SW exposure of the lower zones.
4. **More zones**: extend the analysis to other controlled paths with ridge wind stations.
5. **Verify the two trigger conflicts** (2020-11-14 Sister #6, 2021-02-06 Sister #1) and, if possible, obtain the CDOT control-mission log after 2021.

---

## Appendix A — fetch index and transport winds

**Fetch (Winstral).** For each wind direction: Winstral's upwind shelter index Sx (100 m search, ±15° sector) and the upwind source area of each zone — rays marched up to 1 km upwind, stopping at a slope-break snow trap (Sb > 15°, 50 m separation), at forest, or at the data edge (traps within 50 m of the zone ignored). **Fetch ratio** = source area / zone area.

![Transport wind rose](../figures/seven_sisters_wind_rose.png)

*Figure A1. CALVP transport winds (direction from).*

## Appendix B — reproducing

```bash
source .envrc                                                     # conda env (.conda)
python src/assign_avalanche_paths.py data/avalanches/caic-avalanche-records-2026-11-12.csv
python src/merge_avalanche_records.py \
    data/avalanches/CAIC_HWY_avalanches_2009-01-01_2021-05-04.csv \
    figures/seven_sisters_new_records_assigned.csv \
    data/avalanches/CAIC_HWY_avalanches_2009-2026_sisters_merged.csv
python src/sql_util.py --site sisters --end 2026-06-01            # CALVP + LBAC2
python src/pipeline.py --site sisters                             # DEM, zone clips, metrics, Figure 1
python src/wind_exposure.py --site sisters                        # Sx, fetch (Appendix A)
python src/storm_analysis.py --site sisters                       # weather before each natural
python src/stress_hotspots.py --site sisters                      # fetch-index hotspots
python src/micromet_loading.py --site sisters                     # MicroMet, landing sites, potential, storm test
python src/event_stats.py --site sisters                          # bootstrap / permutation / period split
python src/threshold_analysis.py --site sisters                   # stress threshold, release rate, size (companion note)
python src/loading_days.py --site sisters                         # loading days by wind direction (§5.8)
python src/loading_rate.py --site sisters                         # 24/48/72-h windows, loading rate (companion note)
python src/event_stress_timeline.py --site sisters                # hourly stress before timed releases (companion note)
```

The NAIP image (`data/Sisters/NAIP_sisters.tif`) comes from the USGS National Map `USGSNAIPImagery` ImageServer (exportImage, bbox 423100,4391800,424300,4392800, EPSG:26913, 0.6 m). Profile curvature and aspect entropy (earlier hypotheses) are still computed by the pipeline's zone-metric table but are not used here. The earlier rule-based nucleation-candidate map (steep, E-facing near N-facing terrain, convex) has been removed; nucleation sites now come from the nucleation potential, exported by `micromet_loading.py` as `figures/seven_sisters_nucleation_potential.tif` and as a KML of patches at or above the 95th percentile of start-zone potential (`figures/seven_sisters_nucleation.kml`).

## References

- Colorado Avalanche Information Center (CAIC) highway observer database, 2009–2021, and CAIC avalanche records export, 2026.
- Gvirtzman, S., Kammer, D. S., Adda-Bedia, M. & Fineberg, J. (2025). How frictional ruptures and earthquakes nucleate and evolve. *Nature* 637, 369–374. https://doi.org/10.1038/s41586-024-08287-y
- Li, L. & Pomeroy, J. W. (1997). Estimates of threshold wind speeds for snow transport using meteorological data. *Journal of Applied Meteorology* 36, 205–213.
- Liston, G. E. & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). *Journal of Hydrometeorology* 7, 217–234.
- Pomeroy, J. W. & Gray, D. M. (1990). Saltation of snow. *Water Resources Research* 26, 1583–1594.
- USGS 3D Elevation Program (3DEP), 1 m DEM, DRCOG 2020 lidar; USGS NAIP imagery (2019).
- Winstral, A., Elder, K. & Davis, R. E. (2002). Spatial snow modeling of wind-redistributed snow using terrain-based parameters. *Journal of Hydrometeorology* 3, 524–538.
