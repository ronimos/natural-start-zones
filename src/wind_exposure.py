"""
wind_exposure.py — terrain wind-loading metrics per start zone.

For each wind direction d (the direction the wind comes FROM):

  Sx(d)    Winstral et al. (2002) maximum upwind slope: for every cell, the
           largest elevation angle to upwind terrain within `dmax`, averaged
           over d ± 15° in 5° steps.  Positive = sheltered (deposition),
           negative = exposed (erosion).

  fetch(d) Upwind source area of a start zone: rays are marched upwind from
           every zone cell over d ± 15°, up to `lmax`.  A ray stops at the
           first snow trap, at forest (canopy mask), or at the edge of the
           data.  A trap is a slope break (Winstral et al. 2002):
           Sb(d) = Sx(d) − Sx(d) at a point `sep` upwind > `trap_sb`.  Sb is
           ~0 on uniform slopes (where Sx just equals the slope angle) and
           positive in gullies and on abrupt lee drops, where drifts form
           before snow can reach the zone.  fetch = area of the union of
           cells crossed outside the zone.  fetch_ratio = fetch / zone area:
           the snow-source area per unit of receiving area, so extra
           deposition depth scales with it.

Wind weights come from transport hours at the ridge station: wind speed
≥ `u_thresh`, air temperature below freezing, and fresh snow available
(SNOTEL SWE gain ≥ `swe_gain` over the previous 72 h).  Raw CAIC units:
wind tenths of mph, temperature tenths of °F, SWE hundredths of inch.  Directions are
weighted by hours and by Σu³ (snow transport flux scales roughly with u³).

Outputs (figures/):  <site>_wind_exposure.csv  (zone × direction)
                     <site>_wind_loading.csv   (per-zone weighted metrics,
                                                avalanche counts, correlations)
                     <site>_wind_rose.png, <site>_fetch_polar.png,
                     <site>_fetch_map.png

Usage:  python src/wind_exposure.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
from scipy.stats import spearmanr

gdal.UseExceptions()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, ROOT, _sz_pairs   # noqa: E402
from avl_data import load_avalanches                   # noqa: E402

DIRS = np.arange(0, 360, 5)          # wind-from directions evaluated


# ── Terrain ──────────────────────────────────────────────────────────────────

def load_domain(src: Path, bounds, res: float):
    """DEM window (cubic-resampled to res) as array + geotransform + SRS."""
    ds = gdal.Warp('', str(src), format='MEM', outputBounds=bounds,
                   xRes=res, yRes=res, targetAlignedPixels=True,
                   resampleAlg='cubic', overviewLevel='NONE', dstNodata=np.nan)
    z = ds.GetRasterBand(1).ReadAsArray().astype('f8')
    return z, ds.GetGeoTransform(), ds.GetProjection()


def _shift(a: np.ndarray, drow: int, dcol: int) -> np.ndarray:
    """out[r, c] = a[r + drow, c + dcol] (NaN outside)."""
    out = np.full_like(a, np.nan)
    h, w = a.shape
    r0, r1 = max(0, -drow), min(h, h - drow)
    c0, c1 = max(0, -dcol), min(w, w - dcol)
    if r0 < r1 and c0 < c1:
        out[r0:r1, c0:c1] = a[r0 + drow:r1 + drow, c0 + dcol:c1 + dcol]
    return out


def sx_raw(z: np.ndarray, res: float, d_deg: float, dmax: float) -> np.ndarray:
    """Max upwind elevation angle (deg) along bearing d_deg within dmax."""
    th = np.radians(d_deg)
    best = np.full_like(z, -90.0)
    seen = set()
    for s in np.arange(res, dmax + res / 2, res):
        dcol = int(round(s * np.sin(th) / res))
        drow = int(round(-s * np.cos(th) / res))      # north = up = −row
        if (drow, dcol) in seen or (drow, dcol) == (0, 0):
            continue
        seen.add((drow, dcol))
        dist = res * np.hypot(drow, dcol)
        ang = np.degrees(np.arctan((_shift(z, drow, dcol) - z) / dist))
        best = np.fmax(best, ang)
    best[np.isnan(z)] = np.nan
    return best


def sx_sectors(z, res, dmax):
    """{d: Sx averaged over d ± 15° (5° steps)} for every d in DIRS."""
    raw = {d: sx_raw(z, res, d, dmax) for d in DIRS}
    return {d: np.nanmean([raw[(d + k) % 360] for k in range(-15, 20, 5)], axis=0)
            for d in DIRS}


def polygon_mask(geom, shape, gt, wkt) -> np.ndarray:
    mem = gdal.GetDriverByName('MEM').Create('', shape[1], shape[0], 1, gdal.GDT_Byte)
    mem.SetGeoTransform(gt); mem.SetProjection(wkt)
    srs = osr.SpatialReference(); srs.ImportFromWkt(wkt)
    v = ogr.GetDriverByName('Memory').CreateDataSource('')
    lyr = v.CreateLayer('p', srs=srs)
    f = ogr.Feature(lyr.GetLayerDefn()); f.SetGeometry(ogr.CreateGeometryFromWkb(geom.wkb))
    lyr.CreateFeature(f)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1])
    return mem.ReadAsArray().astype(bool)


def canopy_mask(dsm: Path, dem: Path, bounds, res, shape, height=2.0):
    """Forest from (drone DSM − bare-earth DEM), after removing their median
    vertical offset (datum).  Cells outside the DSM are treated as open."""
    opt = dict(format='MEM', outputBounds=bounds, xRes=res, yRes=res,
               targetAlignedPixels=True, resampleAlg='average',
               overviewLevel='NONE', dstNodata=np.nan)
    s = gdal.Warp('', str(dsm), **opt).ReadAsArray().astype('f8')
    g = gdal.Warp('', str(dem), **opt).ReadAsArray().astype('f8')
    s = s[:shape[0], :shape[1]]; g = g[:shape[0], :shape[1]]
    d = s - g
    d -= np.nanmedian(d)
    return np.nan_to_num(d, nan=0.0) > height


def slope_break(sx_d: np.ndarray, res: float, d_deg: float, sep: float) -> np.ndarray:
    """Sb = Sx at the cell minus Sx at the point `sep` metres upwind."""
    th = np.radians(d_deg)
    up = _shift(sx_d, int(round(-sep * np.cos(th) / res)), int(round(sep * np.sin(th) / res)))
    return sx_d - up


def fetch_cells(zone: np.ndarray, sb_d: np.ndarray, blocked: np.ndarray,
                res: float, d_deg: float, lmax: float, trap_sb: float,
                grace: float = 0.0):
    """Union of upwind cells reachable from the zone before a trap/forest/edge.

    Traps within `grace` metres of the zone edge are ignored: a zone just
    below a crest is loaded through that crest's own lee break (cornice /
    rollover), so it is part of the receiving area, not a barrier.
    """
    h, w = zone.shape
    r0, c0 = np.nonzero(zone)
    visited = np.zeros_like(zone)
    for b in range(-15, 20, 5):
        th = np.radians(d_deg + b)
        dr, dc = -np.cos(th), np.sin(th)
        alive = np.ones(len(r0), bool)
        out_dist = np.zeros(len(r0))          # metres travelled outside the zone
        for s in np.arange(res, lmax + res / 2, res):
            r = np.rint(r0 + dr * s / res).astype(int)
            c = np.rint(c0 + dc * s / res).astype(int)
            inside = (r >= 0) & (r < h) & (c >= 0) & (c < w)
            alive &= inside
            if not alive.any():
                break
            rr, cc = r[alive], c[alive]
            in_zone = zone[rr, cc]
            idx = np.flatnonzero(alive)
            out_dist[idx[~in_zone]] += res
            trap = (sb_d[rr, cc] > trap_sb) & (out_dist[idx] > grace)
            stop = blocked[rr, cc] | np.isnan(sb_d[rr, cc]) | trap
            ok = ~stop & ~in_zone
            visited[rr[ok], cc[ok]] = True
            alive[idx[stop & ~in_zone]] = False
    return visited


# ── Wind ─────────────────────────────────────────────────────────────────────

def transport_weights(csv: Path, u_thresh_mph: float, swe_gain_in: float):
    """Per-direction weights from ridge transport hours (raw CAIC units:
    tenths of mph / °F / inches)."""
    w = pd.read_csv(csv, index_col='time', parse_dates=True)
    u = w['wspd'] / 10.0                       # mph
    t = w['temp'] / 10.0                       # °F
    swe = (w['snowwater'] / 100.0).where(lambda x: x > 0)   # hundredths of inch
    swe = swe.resample('1h').mean().interpolate(limit=6).reindex(w.index)
    gain72 = swe - swe.rolling('72h').min()
    ok = (u >= u_thresh_mph) & (t < 32) & (gain72 >= swe_gain_in) & w['wdir'].notna()
    sel = w.loc[ok]
    b = ((sel['wdir'] / 5).round() * 5 % 360).astype(int)
    hours = b.value_counts().reindex(DIRS, fill_value=0)
    flux = ((u[ok] * 0.44704) ** 3).groupby(b).sum().reindex(DIRS, fill_value=0.0)
    return hours, flux, int(ok.sum()), int(((u >= u_thresh_mph) & (t < 32)).sum())


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--res', type=float, default=2.0, help='grid (m), default 2')
    ap.add_argument('--dmax', type=float, default=100.0, help='Sx search distance (m)')
    ap.add_argument('--lmax', type=float, default=1000.0, help='max fetch length (m)')
    ap.add_argument('--trap-sb', type=float, default=15.0,
                    help='slope break Sb (deg) above which a cell traps drifting snow')
    ap.add_argument('--grace', type=float, default=50.0,
                    help='ignore traps within this distance (m) of the zone edge')
    ap.add_argument('--sep', type=float, default=50.0,
                    help='upwind separation (m) for the slope break Sb')
    ap.add_argument('--u-thresh', type=float, default=15.0, help='transport wind (mph)')
    ap.add_argument('--swe-gain', type=float, default=0.3,
                    help='72-h SWE gain (in) for snow to be available')
    ap.add_argument('--buffer', type=float, default=1300.0, help='domain margin (m)')
    a = ap.parse_args()

    cfg = SITES[a.site]; d = cfg['data_dir']
    slug = cfg['name'].lower().replace(' ', '_')
    src = d / cfg['source_raster']
    pairs = _sz_pairs(cfg)

    srs_wkt = gdal.Open(str(src)).GetProjection()
    epsg = osr.SpatialReference(wkt=srs_wkt).GetAuthorityCode(None)
    zones = {kml.stem: gp.read_file(kml).to_crs(f'EPSG:{epsg}').geometry.union_all()
             for kml, _ in pairs}
    x0, y0, x1, y1 = gp.GeoSeries(list(zones.values())).union_all().buffer(a.buffer).bounds
    z, gt, wkt = load_domain(src, (x0, y0, x1, y1), a.res)
    print(f"domain {z.shape[1]}×{z.shape[0]} cells at {a.res} m, valid {np.isfinite(z).mean():.1%}")

    masks = {k: polygon_mask(g, z.shape, gt, wkt) for k, g in zones.items()}
    dsm = d / '260922_Sisters_Full_Ground-DSM.tiff'
    forest = (canopy_mask(dsm, src, (x0, y0, x1, y1), a.res, z.shape)
              if dsm.exists() else np.zeros(z.shape, bool))
    print(f"canopy mask: {forest.mean():.1%} of domain (drone-survey footprint only)")

    print("Sx for 72 directions …")
    sx = sx_sectors(z, a.res, a.dmax)
    sb = {dd: slope_break(sx[dd], a.res, dd, a.sep) for dd in DIRS}

    rows, foot = [], {}
    for name, m in masks.items():
        area = m.sum() * a.res ** 2
        for dd in DIRS:
            v = fetch_cells(m, sb[dd], forest, a.res, dd, a.lmax, a.trap_sb, a.grace)
            foot[(name, dd)] = v
            rows.append(dict(zone=name, dir=dd, zone_area_m2=area,
                             sx_mean=float(np.nanmean(sx[dd][m])),
                             fetch_m2=float(v.sum() * a.res ** 2),
                             fetch_ratio=float(v.sum() * a.res ** 2 / area)))
    per_dir = pd.DataFrame(rows)

    hours, flux, n_tr, n_wind = transport_weights(ROOT / 'data' / 'weather' / f'{a.site}.csv',
                                                  a.u_thresh, a.swe_gain)
    print(f"transport hours: {n_tr} (of {n_wind} cold hours with u ≥ {a.u_thresh} mph)")
    wh = hours / hours.sum(); wf = flux / flux.sum()

    def wavg(col, wt):
        return per_dir.groupby('zone').apply(
            lambda g: float((g.set_index('dir')[col] * wt).sum()), include_groups=False)

    sectors = {'SW (200–240°)': (200, 240), 'W (245–285°)': (245, 285), 'NW (290–330°)': (290, 330)}
    out = pd.DataFrame({
        'zone_area_m2':        per_dir.groupby('zone')['zone_area_m2'].first(),
        'sx_hours_wt':         wavg('sx_mean', wh),
        'fetch_m2_hours_wt':   wavg('fetch_m2', wh),
        'fetch_ratio_hours_wt': wavg('fetch_ratio', wh),
        'fetch_ratio_flux_wt': wavg('fetch_ratio', wf),
    })
    for lab, (lo, hi) in sectors.items():
        sel = per_dir[(per_dir.dir >= lo) & (per_dir.dir <= hi)]
        out[f'fetch_m2 {lab}'] = sel.groupby('zone')['fetch_m2'].mean()

    # Avalanche counts: natural slabs (pipeline filter) and explosive slab days
    csv = cfg['avalanche_csv']
    nat = load_avalanches(str(csv), cfg['avalanche_filter'], type=['HS', 'SS'], trig=['N'])
    nat_n = nat.groupby('HW Path').size()
    raw = pd.read_csv(csv, usecols=['Date', 'HW Path', 'Trigger', 'Type'])
    raw['Trigger'] = raw['Trigger'].str.strip(); raw['Type'] = raw['Type'].str.strip()
    raw = raw[raw['HW Path'].str.contains(cfg['avalanche_filter'], na=False)]
    nat_unf = raw[(raw.Trigger == 'N') & raw.Type.isin(['HS', 'SS'])].groupby('HW Path').size()
    raw['Date'] = pd.to_datetime(raw['Date'], format='mixed')
    ctrl = raw[raw.Trigger.isin(['AX', 'AL'])]
    ctrl_days = ctrl['Date'].dt.normalize().nunique()
    expl = ctrl[ctrl.Type.isin(['HS', 'SS'])].groupby('HW Path')['Date'].apply(
        lambda s: s.dt.normalize().nunique())
    pmap = {k.replace(cfg['sz_dem_suffix'], ''): v for k, v in cfg['sz_avi_path_map'].items()}
    out['nat_slabs'] = [nat_n.get(pmap.get(z_, z_), 0) for z_ in out.index]
    out['nat_slabs_no_elev_filter'] = [nat_unf.get(pmap.get(z_, z_), 0) for z_ in out.index]
    out['nat_per_10k_m2'] = out['nat_slabs'] / (out['zone_area_m2'] / 1e4)
    out['expl_slab_day_frac'] = [expl.get(pmap.get(z_, z_), 0) / ctrl_days for z_ in out.index]

    print("\nPer-zone wind loading (weights: transport hours unless noted)")
    print(out.round(3).to_string())

    # Correlations — with and without excluded zones
    excl = {s.replace(cfg['sz_dem_suffix'], '') for s in cfg.get('exclude_sz_stems', [])}
    metrics = ['sx_hours_wt', 'fetch_m2_hours_wt', 'fetch_ratio_hours_wt', 'fetch_ratio_flux_wt']
    corr = []
    for label, sub in (('all zones', out), (f'excluding {sorted(excl)}', out.drop(index=list(excl & set(out.index))))):
        for outcome in ('nat_per_10k_m2', 'expl_slab_day_frac'):
            for mtr in metrics:
                rho, p = spearmanr(sub[mtr], sub[outcome])
                corr.append(dict(set=label, N=len(sub), outcome=outcome, metric=mtr,
                                 rho=round(float(rho), 3), p=round(float(p), 3)))
    corr = pd.DataFrame(corr)
    print("\nSpearman correlations (exploratory; N ≤ 5)")
    print(corr.to_string(index=False))

    per_dir.to_csv(FIGURES / f'{slug}_wind_exposure.csv', index=False)
    with open(FIGURES / f'{slug}_wind_loading.csv', 'w') as fh:
        fh.write(f"# params: res={a.res} dmax={a.dmax} lmax={a.lmax} trap_sb={a.trap_sb} sep={a.sep} grace={a.grace} "
                 f"u_thresh={a.u_thresh}mph swe_gain={a.swe_gain}in transport_hours={n_tr}\n")
        out.to_csv(fh)
        fh.write('\n'); corr.to_csv(fh, index=False)

    # ── Figures ──────────────────────────────────────────────────────────────
    th = np.radians(DIRS)
    fig = plt.figure(figsize=(6, 6)); ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
    ax.bar(th, wh.values, width=np.radians(5), color='tab:blue', alpha=.8)
    ax.set_title(f"{cfg['name']} — ridge transport winds (direction FROM)\n"
                 f"{n_tr} h: u ≥ {a.u_thresh} mph, T < 32 °F, 72-h SWE gain ≥ {a.swe_gain} in",
                 fontsize=9)
    fig.savefig(FIGURES / f'{slug}_wind_rose.png', dpi=150, bbox_inches='tight'); plt.close(fig)

    fig = plt.figure(figsize=(7, 7)); ax = fig.add_subplot(projection='polar')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
    for name, g in per_dir.groupby('zone'):
        ax.plot(np.r_[th, th[0]], np.r_[g.fetch_ratio.values, g.fetch_ratio.values[0]],
                label=name.replace('_SZ', '').replace('_', ' #'))
    ax.set_title(f"{cfg['name']} — fetch ratio (upwind source area / zone area) by wind direction",
                 fontsize=9)
    ax.legend(loc='lower left', bbox_to_anchor=(-.15, -.12), fontsize=8)
    fig.savefig(FIGURES / f'{slug}_fetch_polar.png', dpi=150, bbox_inches='tight'); plt.close(fig)

    dom = int(DIRS[np.argmax(wh.values)])
    hs = gdal.DEMProcessing('', gdal.Warp('', str(src), format='MEM', outputBounds=(x0, y0, x1, y1),
                            xRes=a.res, yRes=a.res, targetAlignedPixels=True, resampleAlg='cubic'),
                            'hillshade', format='MEM').ReadAsArray()
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(hs, cmap='gray', alpha=.9)
    cols = plt.cm.tab10(np.arange(len(masks)))
    for (name, m), col in zip(masks.items(), cols):
        v = foot[(name, dom)].astype(float)
        ax.imshow(np.ma.masked_where(v == 0, v), cmap=matplotlib.colors.ListedColormap([col]), alpha=.45)
        ax.contour(m, levels=[.5], colors=[col], linewidths=1.5)
        r, c = np.argwhere(m).mean(axis=0)
        ax.text(c, r, name.replace('_SZ', '').replace('_', ' #'), color='w', fontsize=9,
                ha='center', bbox=dict(fc=col, alpha=.8, lw=0))
    ax.imshow(np.ma.masked_where(~forest, forest), cmap='Greens', alpha=.35, vmin=0, vmax=2)
    ax.set_title(f"{cfg['name']} — fetch footprints for the dominant transport wind ({dom}°)\n"
                 f"trap Sb > {a.trap_sb}° (sep {a.sep:.0f} m), max {a.lmax:.0f} m; green = canopy (drone-survey area)",
                 fontsize=10)
    ax.set_axis_off()
    fig.savefig(FIGURES / f'{slug}_fetch_map.png', dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"\nSaved: {slug}_wind_exposure.csv, {slug}_wind_loading.csv, {slug}_wind_rose.png, "
          f"{slug}_fetch_polar.png, {slug}_fetch_map.png")


if __name__ == '__main__':
    main()
