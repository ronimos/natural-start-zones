"""
micromet_loading.py — terrain-modified winds, snow transport and storm-by-
storm weak-layer stress for each start zone.

1. Terrain winds (MicroMet; Liston & Elder 2006, J. Hydrometeor. 7:217)
     W_d(x) = 1 + γs·Ωs + γc·Ωc
     Ωs = β·cos(ω − ξ)   slope in the wind direction (β slope, ξ aspect,
                          ω wind-from direction); + windward, − lee
     Ωc = curvature over length η (8 neighbours at distance η); + ridges,
          − gullies.  Both scaled to [−0.5, 0.5] over the domain.
   Local speed u(x) = u_station · W_d(x) / W_d(station): the station (CALVP,
   on the ridge crest) measures exposed-ridge wind, not a domain average.
   MicroMet's wind-direction diversion on slopes is not applied to the
   transport direction (simplification).

2. Transport threshold (10-m speed), per hour
     fresh snow (SWE gain ≥ 2.5 mm in the last 24 h): 5 m/s
     otherwise Li & Pomeroy (1997) dry snow: 9.43 + 0.18 T + 0.0033 T² (°C)
     T > 0 °C: no transport (wet snow)

3. Transport capacity — Pomeroy & Gray (1990) saltation
     Q = 0.68 ρ u*t (u*² − u*t²) / (g u*)   [kg m⁻¹ s⁻¹],  u* = κu / ln(z/z0)
   Saltation only (suspension omitted): a lower bound at high wind speeds.

4. Erosion / deposition — march downwind (grid rotated so the wind blows
   along rows).  Where capacity rises, flux relaxes toward it over the
   saturation length Lf (erosion); where it falls (lees, gully floors,
   forest), the excess deposits.  Accumulation A = deposition − erosion
   [kg m⁻² = mm SWE].  Snow supply for erosion is assumed unlimited.

5. Stress per day (Nov–May) from the load accumulated since the last reset
   (explosive control on any Sister path, or a natural release on that path),
   with variants: reset only released paths, 5-day settling decay, and a
   fixed 72-h window.  SWE(x) = h(θ)·max(0, SWE_new + A(x)),
   margin = g·SWE·cos²θ·(tanθ − μ)⁺  [Pa]  (τ − μσ for a friction-type
   weak layer).  h(θ): 1 to 50°, linear to 0 at 60°.

Outputs (figures/): loading maps (all / SW / WNW transport, natural-event
storms), per-zone loading table, zone-day stress table, event test.

Usage:  python src/micromet_loading.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from osgeo import gdal, osr
from scipy import ndimage as ndi

gdal.UseExceptions()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, ROOT, _sz_pairs           # noqa: E402
import wind_exposure as WX                                       # noqa: E402

KAPPA, G, RHO_AIR = 0.4, 9.81, 0.85          # ρ_air at ~3600 m
STATIONS = {'CALVP': (39.6747, -105.894)}    # CAIC stnList lat/lon


# ── 1. Terrain winds ─────────────────────────────────────────────────────────

def curvature(z: np.ndarray, res: float, eta: float) -> np.ndarray:
    n = max(1, int(round(eta / res)))
    def sh(dr, dc):
        return WX._shift(z, dr, dc)
    c = ((z - (sh(0, -n) + sh(0, n)) / 2) / (2 * eta)
         + (z - (sh(-n, 0) + sh(n, 0)) / 2) / (2 * eta)
         + (z - (sh(-n, -n) + sh(n, n)) / 2) / (2 * np.sqrt(2) * eta)
         + (z - (sh(-n, n) + sh(n, -n)) / 2) / (2 * np.sqrt(2) * eta)) / 4
    return c


def scale_half(a: np.ndarray) -> np.ndarray:
    m = np.nanmax(np.abs(a))
    return a / (2 * m) if m > 0 else a


def terrain_weights(z, res, eta, gs=0.5, gc=0.5, dirs=range(0, 360, 10)):
    """{ω: W_ω} for wind-from directions ω."""
    gy, gx = np.gradient(z, res)              # gy: +row = south
    beta = np.arctan(np.hypot(gx, gy))
    # aspect = bearing of steepest DEScent (direction the slope faces)
    xi = np.arctan2(-gx, gy)                  # east = −∂z/∂x, north = +∂z/∂row
    oc = scale_half(curvature(z, res, eta))
    out = {}
    for w in dirs:
        os_ = scale_half(beta * np.cos(np.radians(w) - xi))
        out[w] = 1 + gs * os_ + gc * oc
    return out, np.degrees(beta)


# ── 2–3. Thresholds and capacity ─────────────────────────────────────────────

def u_threshold(temp_c: float, fresh: bool) -> float:
    if temp_c > 0:
        return np.inf
    return 5.0 if fresh else 9.43 + 0.18 * temp_c + 0.0033 * temp_c ** 2


def capacity(u: np.ndarray, ut: float, z=10.0, z0=0.001) -> np.ndarray:
    k = KAPPA / np.log(z / z0)
    us, ust = k * u, k * ut
    q = 0.68 * RHO_AIR * ust * (us ** 2 - ust ** 2) / (G * np.maximum(us, 1e-6))
    return np.where(u > ut, q, 0.0)


# ── 4. Downwind marching on a rotated grid ───────────────────────────────────

def rotate_in(a, ang, cval=np.nan):
    return ndi.rotate(a, ang, reshape=True, order=1, cval=cval, prefilter=False)


def rotate_out(b, ang, shape):
    back = ndi.rotate(b, -ang, reshape=True, order=1, cval=np.nan, prefilter=False)
    r0 = (back.shape[0] - shape[0]) // 2; c0 = (back.shape[1] - shape[1]) // 2
    return back[r0:r0 + shape[0], c0:c0 + shape[1]]


def march(cap: np.ndarray, res: float, lf: float) -> np.ndarray:
    """Accumulation rate (kg m⁻² s⁻¹) for wind blowing along +columns."""
    a = np.zeros_like(cap)
    q = np.full(cap.shape[0], np.nan)
    relax = 1 - np.exp(-res / lf)
    for j in range(cap.shape[1]):
        c = cap[:, j]
        fresh = np.isnan(q) & ~np.isnan(c)
        q = np.where(fresh, c, q)                     # inflow at capacity
        up = c >= q
        qn = np.where(up, q + (c - q) * relax, c)
        a[:, j] = (q - qn) / res                      # + deposition, − erosion
        q = np.where(np.isnan(c), np.nan, qn)
    a[np.isnan(cap)] = np.nan
    return a


# ── Weather ──────────────────────────────────────────────────────────────────

def hourly_weather(csv: Path) -> pd.DataFrame:
    w = pd.read_csv(csv, index_col='time', parse_dates=True)
    h = pd.DataFrame(index=w.index)
    h['u'] = w['wspd'] / 10 * 0.44704                 # m/s
    h['dir'] = w['wdir']
    h['t'] = (w['temp'] / 10 - 32) * 5 / 9            # °C
    swe = (w['snowwater'] / 100 * 25.4).where(lambda x: x > 0)   # mm
    h['swe'] = swe.interpolate(limit=6)
    h['gain24'] = h['swe'] - h['swe'].rolling('24h').min()
    return h.dropna(subset=['u', 'dir', 't'])


def _unflat(vals, idx, shape):
    full = np.full(int(np.prod(shape)), np.nan); full[idx] = vals
    return full.reshape(shape)


def _naip_rgb(cfg, bounds, res, shape, wkt):
    """Basemap RGB (0–1) on the wind grid, or zeros when none is configured."""
    bm = cfg.get('basemap')
    p = cfg['data_dir'] / bm if bm else None
    if p is None or not p.exists():
        return np.zeros(shape + (3,), 'f4')
    ds = gdal.Warp('', str(p), format='MEM', outputBounds=bounds, width=shape[1], height=shape[0],
                   dstSRS=wkt, resampleAlg='bilinear')
    rgb = np.dstack([ds.GetRasterBand(b).ReadAsArray().astype('f4') for b in (1, 2, 3)])
    lo, hi = np.percentile(rgb, [1, 99])
    return np.clip((rgb - lo) / (hi - lo), 0, 1).astype('f4')


def _event_pixels(av, s_srs, d_srs, gt, shape, masks):
    """(row, col, inside_a_zone) of natural records that carry coordinates."""
    out = []
    if not {'Latitude', 'Longitude'} <= set(av.columns):
        return out
    tf = osr.CoordinateTransformation(s_srs, d_srs)
    for _, r in av.dropna(subset=['Latitude', 'Longitude']).iterrows():
        px, py, _ = tf.TransformPoint(float(r['Longitude']), float(r['Latitude']))
        rr, cc = int((py - gt[3]) / gt[5]), int((px - gt[0]) / gt[1])
        if 0 <= rr < shape[0] and 0 <= cc < shape[1]:
            out.append((rr, cc, any(m[rr, cc] for m in masks.values())))
    return out


# ── Figures (from saved grids, so they can be redrawn without rerunning) ─────

def plot_loading(grids: Path, slug: str) -> None:
    g = np.load(grids, allow_pickle=False)
    hs, forest = g['hs'], g['forest']
    names, zones = list(g['zone_names']), g['zones']
    name, eta, lf, res = g['meta']

    def label(n):
        return n.replace('_SZ', '').replace('_', ' #')

    def panel(ax, field, title, vmax, sl=(slice(None), slice(None))):
        ax.imshow(hs[sl], cmap='gray')
        f = np.ma.masked_where(forest[sl], field[sl])
        im = ax.imshow(f, cmap='RdBu', vmin=-vmax, vmax=vmax, alpha=.7)
        ax.imshow(np.ma.masked_where(~forest[sl], forest[sl]), cmap='Greens', vmin=0, vmax=1.5, alpha=.5)
        for n, m in zip(names, zones):
            if m[sl].any():
                ax.contour(m[sl], levels=[.5], colors='k', linewidths=1.2)
                r_, c_ = np.argwhere(m[sl]).mean(axis=0)
                ax.text(c_, r_, label(n), fontsize=8, ha='center', bbox=dict(fc='w', alpha=.7, lw=0))
        ax.set_title(title, fontsize=9); ax.set_axis_off()
        return im

    open_zone = zones.any(axis=0) | ndi.binary_dilation(zones.any(axis=0), iterations=15)
    def vmax_of(f):
        v = np.abs(f[open_zone & ~forest])
        return max(float(np.nanpercentile(v, 98)) if v.size else 1.0, 1e-3)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.6), constrained_layout=True)
    for ax, key, t in ((axs[0], 'all', 'All transport'), (axs[1], 'SW', 'SW winds (200–260°)'),
                       (axs[2], 'WNW', 'W–NW winds (270–330°)')):
        im = panel(ax, g[key], f'{t}: wind accumulation per season (mm SWE)\n'
                   'blue = deposition, red = erosion, green = forest', vmax_of(g[key]))
        fig.colorbar(im, ax=ax, fraction=.035)
    fig.suptitle(f"{name} — MicroMet/saltation loading (north up; η={eta} m, Lf={lf} m; saltation only)")
    fig.savefig(FIGURES / f'{slug}_loading_patterns.png', dpi=140); plt.close(fig)

    # Close-ups: Sister #1 and #6 by sector
    picks = [n for n in names if n.endswith(('1_SZ', '6_SZ'))]
    fig, axs = plt.subplots(len(picks), 3, figsize=(13, 6 * len(picks)), constrained_layout=True)
    axs = np.atleast_2d(axs)
    for row, n in zip(axs, picks):
        m = zones[names.index(n)]
        r, c = np.nonzero(ndi.binary_dilation(m, iterations=int(40 / float(res))))
        sl = (slice(r.min(), r.max()), slice(c.min(), c.max()))
        for ax, key, t in zip(row, ('all', 'SW', 'WNW'), ('All', 'SW', 'W–NW')):
            v = np.abs(g[key][sl][m[sl]]); vm = max(float(np.percentile(v, 98)), 1e-3)
            im = panel(ax, g[key], f'{label(n)} — {t} (mm SWE / season)', vm, sl)
            fig.colorbar(im, ax=ax, fraction=.05)
    fig.suptitle('Close-ups, north up (top of start zone at the bottom)')
    fig.savefig(FIGURES / f'{slug}_loading_closeups.png', dpi=140); plt.close(fig)

    ev, days, paths = g['ev'], list(g['ev_days']), list(g['ev_paths'])
    keep = [i for i in range(len(days)) if np.abs(ev[i]).sum() > 0][:9]
    if keep:
        ncol = 3; nrow = int(np.ceil(len(keep) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(18, 5.2 * nrow), constrained_layout=True)
        axs = np.atleast_1d(axs).ravel()
        for ax, i in zip(axs, keep):
            im = panel(ax, ev[i], f"{days[i]} — natural on {paths[i]}\n72-h wind accumulation (mm SWE)",
                       vmax_of(ev[i]))
            fig.colorbar(im, ax=ax, fraction=.035)
        for ax in axs[len(keep):]:
            ax.set_axis_off()
        fig.savefig(FIGURES / f'{slug}_event_loading.png', dpi=120); plt.close(fig)

    if 'potential' in g.files:
        # Forest is masked: trees anchor the snowpack, and the model's canopy
        # edges collect intercepted snow that would otherwise dominate the map.
        P = np.where(g['region'] & ~forest, g['potential'], np.nan)
        fig, axs = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True,
                                gridspec_kw=dict(width_ratios=[2.2, 1]))
        ax = axs[0]
        if g['naip'].any():
            ax.imshow(g['naip'])
        else:
            ax.imshow(hs, cmap='gray')
        zcells = zones.any(axis=0) & (P > 0)
        vmax = np.nanpercentile(P[zcells], 98) if zcells.any() else 1
        im = ax.imshow(np.ma.masked_where(~(P > 0), P), cmap='inferno', vmin=0, vmax=vmax, alpha=.85)
        for n, m in zip(names, zones):
            ax.contour(m, levels=[.5], colors='cyan', linewidths=1.4)
            r_, c_ = np.argwhere(m).mean(axis=0)
            ax.text(c_, np.argwhere(m)[:, 0].min() - 8, label(n), fontsize=9, ha='center',
                    bbox=dict(fc='w', alpha=.8, lw=0))
        pts = g['pts_xy']
        if len(pts):
            inside = pts[:, 2] == 1
            ax.plot(pts[~inside, 1], pts[~inside, 0], 'o', ms=4, mfc='none', mec='w', label='natural (track / debris point)')
            ax.plot(pts[inside, 1], pts[inside, 0], '*', ms=13, mfc='lime', mec='k', label='natural (inside start zone)')
            ax.legend(loc='lower left', fontsize=8)
        fig.colorbar(im, ax=ax, fraction=.03, extend='max',
                     label='Σ (stress margin − threshold)⁺ per season  [Pa·day]; scale: start-zone P98')
        ax.set_title('Nucleation potential: excess weak-layer stress above threshold, summed over 72-h loads\n'
                     '(terrain + MicroMet winds + SNOTEL snowfall; threshold = P99 of start-zone stress; forest masked; north up)',
                     fontsize=10)
        ax.set_axis_off()
        ax = axs[1]
        zn = [label(str(n)) for n in g['oe_zone']]; o = g['oe_obs']; x = np.arange(len(zn))
        from scipy.stats import chi2
        lo = np.array([chi2.ppf(.025, 2 * k) / 2 if k else 0 for k in o]); hi = chi2.ppf(.975, 2 * (o + 1)) / 2
        ax.bar(x - .2, g['oe_exp_hot'], .4, color='tab:orange', alpha=.8, label='expected: hotspot days')
        ax.bar(x + .2, g['oe_exp_exc'], .4, color='tab:red', alpha=.6, label='expected: summed excess')
        ax.errorbar(x, o, yerr=[o - lo, hi - o], fmt='ko', capsize=4, label='observed (95% Poisson CI)')
        ax.set_xticks(x, zn); ax.set_ylabel('natural slab days (with weather data)')
        ax.set_title('Observed vs expected natural releases\n(expected scaled to the observed total)', fontsize=10)
        ax.legend(fontsize=8)
        fig.savefig(FIGURES / f'{slug}_nucleation_potential.png', dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--res', type=float, default=2.0)
    ap.add_argument('--eta', type=float, default=30.0, help='curvature length scale (m)')
    ap.add_argument('--lf', type=float, default=100.0, help='saturation length (m)')
    ap.add_argument('--mu-deg', type=float, default=30.0)
    ap.add_argument('--window-h', type=int, default=72)
    ap.add_argument('--buffer', type=float, default=1300.0)
    ap.add_argument('--tag', default='', help='suffix for output names (sensitivity runs)')
    ap.add_argument('--map-buffer', type=float, default=100.0,
                    help='ring (m) around the zones kept for the nucleation-potential map')
    ap.add_argument('--plot-only', action='store_true',
                    help='redraw figures from the saved <slug>_loading_grids.npz')
    a = ap.parse_args()
    if a.plot_only:
        slug = SITES[a.site]['name'].lower().replace(' ', '_') + a.tag
        plot_loading(FIGURES / f'{slug}_loading_grids.npz', slug)
        print('figures redrawn'); return
    t0 = time.time()

    cfg = SITES[a.site]; d = cfg['data_dir']
    slug = cfg['name'].lower().replace(' ', '_') + a.tag
    src = d / cfg['source_raster']
    wkt_src = gdal.Open(str(src)).GetProjection()
    proj4 = osr.SpatialReference(wkt=wkt_src).ExportToProj4()
    zones = {k.stem: gp.read_file(k).to_crs(proj4).geometry.union_all() for k, _ in _sz_pairs(cfg)}
    x0, y0, x1, y1 = gp.GeoSeries(list(zones.values())).union_all().buffer(a.buffer).bounds
    z, gt, wkt = WX.load_domain(src, (x0, y0, x1, y1), a.res)
    masks = {k: WX.polygon_mask(g, z.shape, gt, wkt) for k, g in zones.items()}
    dsm = d / '260922_Sisters_Full_Ground-DSM.tiff'
    forest = (WX.canopy_mask(dsm, src, (x0, y0, x1, y1), a.res, z.shape)
              if dsm.exists() else np.zeros(z.shape, bool))

    # Station cell
    s_srs = osr.SpatialReference(); s_srs.ImportFromEPSG(4326)
    s_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    d_srs = osr.SpatialReference(wkt=wkt); d_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    lat, lon = STATIONS['CALVP']
    sx_, sy_, _ = osr.CoordinateTransformation(s_srs, d_srs).TransformPoint(lon, lat)
    srow, scol = int((sy_ - gt[3]) / gt[5]), int((sx_ - gt[0]) / gt[1])
    print(f"domain {z.shape[1]}×{z.shape[0]} @ {a.res} m; CALVP at cell ({srow},{scol}), "
          f"elev {z[srow, scol]:.0f} m (stnList 3624 m)")

    Wd, slope_deg = terrain_weights(z, a.res, a.eta)
    wst = {w: float(np.nanmean(Wd[w][srow - 2:srow + 3, scol - 2:scol + 3])) for w in Wd}
    print("W at station by direction (SW/W/NW): "
          f"{wst[220]:.2f} / {wst[270]:.2f} / {wst[300]:.2f}")

    # Hours → bins (direction 10°, station speed 1 m/s, threshold 0.25 m/s)
    h = hourly_weather(ROOT / 'data' / 'weather' / f'{a.site}.csv')
    h['ut'] = [u_threshold(t, g >= 2.5) for t, g in zip(h['t'], h['gain24'].fillna(0))]
    h['dbin'] = ((h['dir'] / 10).round() * 10 % 360).astype(int)
    h['ubin'] = (h['u']).round().astype(int)
    h['tbin'] = (h['ut'] * 4).round() / 4
    wmax = max(np.nanpercentile(Wd[w], 99.9) / wst[w] for w in Wd)
    h = h[np.isfinite(h['ut']) & (h['u'] * wmax > h['ut'])]   # some cell can transport
    bins = h.groupby(['dbin', 'ubin', 'tbin']).size()
    print(f"{len(h)} candidate transport hours in {len(bins)} bins")

    # Event windows (natural slabs) for maps
    av = pd.read_csv(cfg['avalanche_csv'])
    for c in ('Trigger', 'Type', 'Elev'):
        av[c] = av[c].astype(str).str.strip()
    av = av[av['HW Path'].str.contains(cfg['avalanche_filter'], na=False)
            & (av.Trigger == 'N') & av.Type.isin(['HS', 'SS'])]
    av = av[(av['HW Path'] != 'Seven Sister #6') | av.Elev.isin(['>TL', 'TL'])]
    av['Date'] = pd.to_datetime(av['Date'], format='mixed')
    pmap = {k.replace(cfg['sz_dem_suffix'], ''): v for k, v in cfg['sz_avi_path_map'].items()}
    path2zone = {v: k for k, v in pmap.items()}
    ev_days = sorted(av['Date'].dt.normalize().unique())

    def window_end(day):            # end of local (MST) day, in UTC
        return pd.Timestamp(day) + pd.Timedelta(days=1, hours=7)

    # Zone cells (+ ring) we need per hour; maps we accumulate
    allz = np.zeros(z.shape, bool)
    for m in masks.values():
        allz |= m
    # Cells kept per bin: start zones plus a ring around them (for the
    # nucleation-potential map); zone statistics use the zone cells only.
    region = ndi.binary_dilation(allz, iterations=int(a.map_buffer / a.res))
    zc = np.flatnonzero(region.ravel())
    sector = {'SW': (200, 260), 'WNW': (270, 330)}
    acc_all = np.zeros(z.shape); acc_sec = {k: np.zeros(z.shape) for k in sector}
    acc_ev = {dd: np.zeros(z.shape) for dd in ev_days}
    ev_hours = {dd: h[(h.index > window_end(dd) - pd.Timedelta(hours=a.window_h))
                      & (h.index <= window_end(dd))] for dd in ev_days}
    bin_zone = {}                      # bin → accumulation (kg/m² per hour) at zone cells

    cache_dir = None
    for i, (w, grp) in enumerate(bins.groupby(level=0)):
        flow = (w + 180) % 360
        ang = flow - 90
        Wr = rotate_in(Wd[w] / wst[w], ang)
        Fr = rotate_in(forest.astype(float), ang, cval=0.0) > 0.5
        for (wb, ub, tb), n in grp.items():
            cap = capacity(Wr * ub, tb)
            cap[Fr] = 0.0
            A = rotate_out(march(cap, a.res, a.lf), ang, z.shape) * 3600.0   # kg m⁻² per hour
            A = np.nan_to_num(A)
            acc_all += n * A
            for k, (lo, hi) in sector.items():
                if lo <= w <= hi:
                    acc_sec[k] += n * A
            for dd, eh in ev_hours.items():
                k_ = int(((eh.dbin == wb) & (eh.ubin == ub) & (eh.tbin == tb)).sum())
                if k_:
                    acc_ev[dd] += k_ * A
            bin_zone[(wb, ub, tb)] = A.ravel()[zc].astype('f4')
        print(f"  direction {w:3d}°: {len(grp)} bins  ({time.time() - t0:.0f} s)")

    # ── Per-zone climatological loading ───────────────────────────────────────
    seasons = len(set(h.index.year))
    rows = {}
    for name, m in masks.items():
        rows[name] = dict(
            mean_acc_all=float(acc_all[m].mean()) / seasons,
            p90_acc_all=float(np.percentile(acc_all[m], 90)) / seasons,
            mean_acc_SW=float(acc_sec['SW'][m].mean()) / seasons,
            mean_acc_WNW=float(acc_sec['WNW'][m].mean()) / seasons,
            frac_cells_depositing=float((acc_all[m] > 0).mean()))
    zt = pd.DataFrame(rows).T
    print("\nWind accumulation per season (kg/m² = mm SWE; saltation only):")
    print(zt.round(1).to_string())

    # ── Daily load increments at zone cells ───────────────────────────────────
    sl = np.radians(slope_deg)
    mu = np.tan(np.radians(a.mu_deg))
    shape_term = (np.cos(sl) ** 2 * np.clip(np.tan(sl) - mu, 0, None)
                  * np.clip((60 - slope_deg) / 10, 0, 1)).ravel()[zc]
    zone_of = np.full(len(zc), '', dtype=object)
    for name, m in masks.items():
        zone_of[m.ravel()[zc]] = name
    zsel = {name: zone_of == name for name in masks}

    hw = hourly_weather(ROOT / 'data' / 'weather' / f'{a.site}.csv')
    local_day = lambda idx: (idx - pd.Timedelta(hours=7)).normalize()   # MST day of a UTC hour
    swe_day = hw['swe'].groupby(local_day(hw.index)).mean()
    # New snow: day-to-day gain in mean SWE; gains under one SNOTEL increment
    # (2.54 mm) are treated as sensor jitter.
    gain = swe_day.diff()
    new_day = gain.where(gain >= 2.54, 0.0).clip(lower=0)
    valid_day = hw['swe'].groupby(local_day(hw.index)).count() >= 12
    A_day = {}
    for dd, eh in h.groupby(local_day(h.index)):
        A = np.zeros(len(zc), 'f4')
        for key, n in eh.groupby(['dbin', 'ubin', 'tbin']).size().items():
            A += n * bin_zone[key]
        A_day[dd] = A
    zero = np.zeros(len(zc), 'f4')

    # Resets: explosive control (CAIC AX/AL/AA/AB/AE/AH/AP on any Sister path)
    # and natural releases (that path only).
    raw = pd.read_csv(cfg['avalanche_csv'], usecols=['Date', 'HW Path', 'Trigger'])
    raw['Trigger'] = raw['Trigger'].astype(str).str.strip()
    raw = raw[raw['HW Path'].isin(path2zone)]
    raw['day'] = pd.to_datetime(raw['Date'], format='mixed').dt.normalize()
    expl = raw[raw.Trigger.isin(['AA', 'AB', 'AE', 'AH', 'AL', 'AP', 'AX'])]
    ctrl_days = set(expl['day'])
    ctrl_path_days = {(r.day, path2zone[r['HW Path']]) for _, r in expl.iterrows()}
    nat_path_days = {(pd.Timestamp(r['Date']).normalize(), path2zone[r['HW Path']])
                     for _, r in av.iterrows() if r['HW Path'] in path2zone}
    print(f"\nResets: {len(ctrl_days)} control days, {len(nat_path_days)} natural path-days")

    scenarios = {'72h window': dict(kind='window'),
                 'since control (all paths)': dict(kind='reset', scope='all', T=None),
                 'since control (released paths)': dict(kind='reset', scope='path', T=None),
                 'since control (all), 5-d settling': dict(kind='reset', scope='all', T=5.0),
                 'since control (released), 5-d settling': dict(kind='reset', scope='path', T=5.0)}
    all_days = pd.date_range(swe_day.index.min(), swe_day.index.max(), freq='D')
    # Since-control scenarios are only meaningful while control missions are
    # recorded (the 2026 CAIC export adds naturals only).
    ctrl_end = (pd.Timestamp(cfg['control_records_end'])
                if cfg.get('control_records_end') else None)
    eval_mask = all_days.month.isin([11, 12, 1, 2, 3, 4, 5])

    def run_scenario(sc):
        """{day: {zone: (new_snow_mm, wind_vec, load_vec)}} for evaluation days."""
        out = {}
        if sc['kind'] == 'window':
            for i, dd in enumerate(all_days):
                if not eval_mask[i] or not valid_day.get(dd, False):
                    continue
                win_d = all_days[max(0, i - 2):i + 1]
                nsn = float(sum(new_day.get(x, 0.0) for x in win_d))
                wv = sum((A_day.get(x, zero) for x in win_d), zero.copy())
                out[dd] = {z_: (nsn, wv[zsel[z_]], (nsn + wv)[zsel[z_]]) for z_ in masks}
            return out
        decay = np.exp(-1.0 / sc['T']) if sc['T'] else 1.0
        zz = {z_: np.zeros(int(zsel[z_].sum()), 'f4') for z_ in masks}
        load = {z_: zz[z_].copy() for z_ in masks}; snow = {z_: 0.0 for z_ in masks}
        wind = {z_: zz[z_].copy() for z_ in masks}
        for i, dd in enumerate(all_days):
            if dd.month == 10 and dd.day == 1:                  # new season
                load = {z_: zz[z_].copy() for z_ in masks}; snow = {z_: 0.0 for z_ in masks}
                wind = {z_: zz[z_].copy() for z_ in masks}
            nsn = float(new_day.get(dd, 0.0)); A = A_day.get(dd, zero)
            for z_ in masks:
                Az = A[zsel[z_]]
                snow[z_] = snow[z_] * decay + nsn
                wind[z_] = wind[z_] * decay + Az
                load[z_] = load[z_] * decay + nsn + Az
            if eval_mask[i] and valid_day.get(dd, False) and (ctrl_end is None or dd <= ctrl_end):
                out[dd] = {z_: (snow[z_], wind[z_].copy(), load[z_].copy()) for z_ in masks}
            # resets take effect after the day is evaluated
            for z_ in masks:
                hit = (dd, z_) in nat_path_days or (
                    (dd in ctrl_days) if sc['scope'] == 'all' else (dd, z_) in ctrl_path_days)
                if hit:
                    load[z_] = zz[z_].copy(); snow[z_] = 0.0; wind[z_] = zz[z_].copy()
        return out

    def auc(x, y):
        pos, neg = x[y == 1], x[y == 0]
        if len(pos) == 0 or len(neg) == 0:
            return np.nan
        return float(((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()))

    ev = av.assign(day=av['Date'].dt.normalize().dt.date, zone=av['HW Path'].map(path2zone))[['day', 'zone']]
    ev['natural'] = 1
    tests, zd_main, zd_all = [], None, []
    for sc_name, sc in scenarios.items():
        res = run_scenario(sc)
        # pooled thresholds from all zone-day cells (outcome-blind)
        pool = np.concatenate([G * np.maximum(d_[z_][2], 0) * shape_term[zsel[z_]]
                               for d_ in res.values() for z_ in masks])
        q95, q99 = np.percentile(pool, [95, 99])
        rows = []
        for dd, per in res.items():
            for z_, (nsn, wv, ld) in per.items():
                sel = zsel[z_]
                mz = G * np.maximum(ld, 0) * shape_term[sel]
                rows.append(dict(day=dd.date(), zone=z_, swe_new_mm=nsn,
                                 wind_acc_p90=float(np.percentile(wv, 90)),
                                 margin_p99=float(np.percentile(mz, 99)),
                                 area_gt_q95=float((mz > q95).sum() * a.res ** 2),
                                 area_gt_q99=float((mz > q99).sum() * a.res ** 2)))
        zd = pd.DataFrame(rows).merge(ev.drop_duplicates(), on=['day', 'zone'], how='left').fillna({'natural': 0})
        for col in ('swe_new_mm', 'wind_acc_p90', 'margin_p99', 'area_gt_q95', 'area_gt_q99'):
            x = zd[col].values; y = zd.natural.values
            wins = [auc(g[col].values, g.natural.values)
                    for _, g in zd[zd.day.isin(zd.loc[zd.natural == 1, 'day'])].groupby('day')
                    if g.natural.sum() and (g.natural == 0).sum()]
            tests.append(dict(scenario=sc_name, metric=col, AUC_all=round(auc(x, y), 3),
                              AUC_within_day=round(float(np.nanmean(wins)), 3) if wins else np.nan,
                              n_events=int(y.sum()), thr_q95_Pa=round(float(q95), 1), thr_q99_Pa=round(float(q99), 1)))
        zd_all.append(zd.assign(scenario=sc_name))
        if sc_name == '72h window':
            zd_main = zd
            thr_window = (float(q95), float(q99))
    test = pd.DataFrame(tests)
    print(f"\nZone-days: {len(zd_main)} ({int(zd_main.natural.sum())} with a natural slab on that path)")
    print("AUC (0.5 = no skill; 1 = event zone-days always higher):")
    pd.set_option('display.width', 220)
    print(test.pivot_table(index='metric', columns='scenario', values='AUC_all', sort=False).to_string())
    print("\nWithin event days (released path vs the others):")
    print(test.pivot_table(index='metric', columns='scenario', values='AUC_within_day', sort=False).to_string())

    # ── Nucleation potential (Gvirtzman et al. 2025: rate ∝ τ − τ_thresh) ────
    # Over every evaluated day (72-h load), excess = (margin − thr)⁺ summed per
    # cell; thr = pooled P99 of start-zone cell margins (P95 as a check).
    # Expected natural-release days per zone ∝ days with a hotspot patch
    # (≥ 2 cells above thr) and, alternatively, ∝ the zone's summed excess.
    from scipy.stats import chi2
    inzone = zone_of != ''
    eval_days = [dd for i, dd in enumerate(all_days) if eval_mask[i] and valid_day.get(dd, False)]
    n_seasons = len({season for season in ((d.year + (d.month >= 10)) for d in eval_days)})
    pot = {}; exceed = {}; hot_days = {}; excess_zone = {}
    for qn, thr in (('P99', thr_window[1]), ('P95', thr_window[0])):
        exc = np.zeros(len(zc)); cnt = np.zeros(len(zc)); hd = {z_: 0 for z_ in masks}
        ez = {z_: 0.0 for z_ in masks}
        for dd in eval_days:
            win_d = pd.date_range(dd - pd.Timedelta(days=2), dd, freq='D')
            nsn = float(sum(new_day.get(x, 0.0) for x in win_d))
            wv = sum((A_day.get(x, zero) for x in win_d), zero.copy())
            mg = G * np.maximum(nsn + wv, 0) * shape_term
            e = np.clip(mg - thr, 0, None)
            exc += e; cnt += mg > thr
            for z_ in masks:
                sel = zsel[z_]
                if ((mg[sel] > thr).sum()) >= 2:
                    hd[z_] += 1
                ez[z_] += float(e[sel].sum()) * a.res ** 2
        pot[qn], exceed[qn], hot_days[qn], excess_zone[qn] = exc / n_seasons, cnt / n_seasons, hd, ez

    eval_set = set(d.date() for d in eval_days)
    obs = {z_: int(((ev.zone == z_) & ev.day.isin(eval_set)).sum()) for z_ in masks}
    oe_rows = []
    for qn in ('P99', 'P95'):
        for basis, raw in (('hotspot days', hot_days[qn]), ('summed excess', excess_zone[qn])):
            tot_o = sum(obs.values()); tot_r = sum(raw.values())
            for z_ in masks:
                e_ = tot_o * raw[z_] / tot_r if tot_r else np.nan
                o_ = obs[z_]
                lo_ = chi2.ppf(.025, 2 * o_) / 2 if o_ else 0.0
                hi_ = chi2.ppf(.975, 2 * (o_ + 1)) / 2
                oe_rows.append(dict(threshold=qn, basis=basis, zone=z_, observed=o_, expected=round(e_, 2),
                                    obs_lo95=round(lo_, 2), obs_hi95=round(hi_, 2),
                                    ratio=round(o_ / e_, 2) if e_ else np.nan))
    oe = pd.DataFrame(oe_rows)
    rho_rows = []
    for (qn, basis), g in oe.groupby(['threshold', 'basis'], sort=False):
        from scipy.stats import spearmanr as _sp
        rho_rows.append(dict(threshold=qn, basis=basis, rho=round(float(_sp(g.expected, g.observed)[0]), 2)))
    print(f"\nObserved vs expected natural-release days (evaluated days with weather: {len(eval_days)}, "
          f"{n_seasons} seasons):")
    print(oe.pivot_table(index=['threshold', 'basis', 'zone'], values=['observed', 'expected', 'ratio'],
                         sort=False).to_string())
    print(pd.DataFrame(rho_rows).to_string(index=False))

    # Observed release points inside start zones (2020–2026 coordinates)
    pts_rows = []
    if {'Latitude', 'Longitude'} <= set(av.columns):
        tf = osr.CoordinateTransformation(s_srs, d_srs)
        P = pot['P99']
        full = np.full(z.size, np.nan); full[zc] = P; full = full.reshape(z.shape)
        for _, r in av.dropna(subset=['Latitude', 'Longitude']).iterrows():
            px, py, _ = tf.TransformPoint(float(r['Longitude']), float(r['Latitude']))
            rr_, cc_ = int((py - gt[3]) / gt[5]), int((px - gt[0]) / gt[1])
            if not (0 <= rr_ < z.shape[0] and 0 <= cc_ < z.shape[1]):
                continue
            zone_hit = next((z_ for z_, m in masks.items() if m[rr_, cc_]), None)
            if zone_hit is None:
                continue
            vals = full[masks[zone_hit]]
            pts_rows.append(dict(date=str(r['Date'].date()), zone=zone_hit, potential=round(float(full[rr_, cc_]), 1),
                                 pctile_in_zone=round(float((vals < full[rr_, cc_]).mean() * 100), 0)))
    pts_df = pd.DataFrame(pts_rows)
    if len(pts_df):
        print("\nObserved releases located inside start zones — nucleation potential percentile within zone:")
        print(pts_df.to_string(index=False))
    oe.to_csv(FIGURES / f'{slug}_observed_vs_expected.csv', index=False)
    pts_df.to_csv(FIGURES / f'{slug}_release_points_potential.csv', index=False)

    zt.to_csv(FIGURES / f'{slug}_wind_accumulation.csv')
    pd.concat(zd_all).to_csv(FIGURES / f'{slug}_zone_day_stress.csv', index=False)
    test.to_csv(FIGURES / f'{slug}_stress_event_test.csv', index=False)

    # ── Save grids (zone window, float32) and draw figures ────────────────────
    hs_ds = gdal.Warp('', str(src), format='MEM', outputBounds=(x0, y0, x1, y1),
                      xRes=a.res, yRes=a.res, targetAlignedPixels=True, resampleAlg='cubic')
    hs = gdal.DEMProcessing('', hs_ds, 'hillshade', format='MEM').ReadAsArray()
    rr, cc = np.nonzero(ndi.binary_dilation(allz, iterations=60))
    win = (slice(rr.min(), rr.max()), slice(cc.min(), cc.max()))
    ev_labels = {str(pd.Timestamp(dd).date()): ', '.join(sorted(
        av.loc[av['Date'].dt.normalize() == dd, 'HW Path'].str.replace('Seven Sister ', '')))
        for dd in ev_days}
    grids = FIGURES / f'{slug}_loading_grids.npz'
    np.savez_compressed(
        grids, hs=hs[win].astype('f4'), forest=forest[win],
        slope_deg=slope_deg[win].astype('f4'), res=np.float32(a.res),
        mean_new_snow_mm=np.float32(new_day.groupby((new_day.index + pd.DateOffset(months=3)).year).sum().mean()),
        all=(acc_all / seasons)[win].astype('f4'),
        SW=(acc_sec['SW'] / seasons)[win].astype('f4'),
        WNW=(acc_sec['WNW'] / seasons)[win].astype('f4'),
        zone_names=np.array(list(masks)), zones=np.stack([m[win] for m in masks.values()]),
        ev_days=np.array(list(ev_labels)), ev_paths=np.array(list(ev_labels.values())),
        ev=np.stack([acc_ev[dd][win].astype('f4') for dd in ev_days]),
        meta=np.array([cfg['name'], str(a.eta), str(a.lf), str(a.res)]),
        potential=_unflat(pot['P99'], zc, z.shape)[win].astype('f4'),
        exceed_days=_unflat(exceed['P99'], zc, z.shape)[win].astype('f4'),
        region=region[win],
        naip=_naip_rgb(cfg, (x0, y0, x1, y1), a.res, z.shape, wkt)[win],
        pts_xy=np.array([[r_ - win[0].start, c_ - win[1].start, 1 if inside else 0]
                         for r_, c_, inside in _event_pixels(av, s_srs, d_srs, gt, z.shape, masks)] or
                        np.zeros((0, 3))),
        oe_zone=np.array(list(masks)), oe_obs=np.array([obs[z_] for z_ in masks]),
        oe_exp_hot=oe[(oe.threshold == 'P99') & (oe.basis == 'hotspot days')].set_index('zone').loc[list(masks), 'expected'].values,
        oe_exp_exc=oe[(oe.threshold == 'P99') & (oe.basis == 'summed excess')].set_index('zone').loc[list(masks), 'expected'].values)
    plot_loading(grids, slug)
    print(f"\nSaved figures and tables with prefix {slug}_  ({time.time() - t0:.0f} s)")


if __name__ == '__main__':
    main()
