"""
event_stress_timeline.py — hour-by-hour hotspot stress before each timed natural.

For every natural slab with a recorded release time, rebuild the load hour by
hour from 24 h before release to 6 h after: station snowfall (LBAC2 SWE,
3-h median smoothed) plus MicroMet wind deposition (the same terrain winds,
thresholds and saltation march as micromet_loading.py, run only for the
wind/speed/threshold bins that occur in these hours).  Stress margin per cell
= g · load · cos²θ · (tanθ − μ), forest-buffered; the hotspot value is the
2nd-highest cell in each start zone (a patch of ≥ 2 cells), as elsewhere.

Times are local; converted to UTC as MST (UTC−7), so March–May times may be
an hour off.  Records whose date is estimated and time is 12:00 AM are
placeholders and dropped.

Outputs (figures/, prefix <slug>_):
  event_stress_timeline.png   all events aligned at release, one line each
  event_stress_panels.png     one panel per release episode, all five paths
  event_stress_timeline.csv   hourly hotspot stress per event and zone

Usage:  python src/event_stress_timeline.py --site sisters
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
from osgeo import gdal, osr
from scipy import ndimage as ndi

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, ROOT, _sz_pairs           # noqa: E402
import wind_exposure as WX                                       # noqa: E402
from micromet_loading import (G, STATIONS, terrain_weights, u_threshold, capacity,  # noqa: E402
                              rotate_in, rotate_out, march, hourly_weather)

PRE, POST = 24, 6
PATH_COLORS = {'#1': '#2a78d6', '#2': '#eb6834', '#3': '#1baf7a', '#4': '#eda100', '#6': '#e87ba4'}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--res', type=float, default=2.0)
    ap.add_argument('--eta', type=float, default=30.0)
    ap.add_argument('--lf', type=float, default=100.0)
    ap.add_argument('--mu-deg', type=float, default=30.0)
    ap.add_argument('--buffer', type=float, default=1300.0)
    ap.add_argument('--forest-buffer', type=float, default=4.0)
    ap.add_argument('--threshold', type=float, default=30.0, help='stress threshold to draw (Pa)')
    a = ap.parse_args()

    cfg = SITES[a.site]; d = cfg['data_dir']
    slug = cfg['name'].lower().replace(' ', '_')
    src = d / cfg['source_raster']
    proj4 = osr.SpatialReference(wkt=gdal.Open(str(src)).GetProjection()).ExportToProj4()
    zones = {k.stem: gp.read_file(k).to_crs(proj4).geometry.union_all() for k, _ in _sz_pairs(cfg)}
    x0, y0, x1, y1 = gp.GeoSeries(list(zones.values())).union_all().buffer(a.buffer).bounds
    z, gt, wkt = WX.load_domain(src, (x0, y0, x1, y1), a.res)
    masks = {k: WX.polygon_mask(g, z.shape, gt, wkt) for k, g in zones.items()}
    dsm = d / '260922_Sisters_Full_Ground-DSM.tiff'
    forest = (WX.canopy_mask(dsm, src, (x0, y0, x1, y1), a.res, z.shape)
              if dsm.exists() else np.zeros(z.shape, bool))
    s_srs = osr.SpatialReference(); s_srs.ImportFromEPSG(4326)
    s_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    d_srs = osr.SpatialReference(wkt=wkt); d_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    lat, lon = STATIONS['CALVP']
    sx_, sy_, _ = osr.CoordinateTransformation(s_srs, d_srs).TransformPoint(lon, lat)
    srow, scol = int((sy_ - gt[3]) / gt[5]), int((sx_ - gt[0]) / gt[1])
    Wd, slope_deg = terrain_weights(z, a.res, a.eta)
    wst = {w: float(np.nanmean(Wd[w][srow - 2:srow + 3, scol - 2:scol + 3])) for w in Wd}

    # Zone cells, slab-stress shape term (as micromet_loading.py)
    pmap = {k.replace(cfg['sz_dem_suffix'], ''): v for k, v in cfg['sz_avi_path_map'].items()}
    path2zone = {v: k for k, v in pmap.items()}
    allz = np.zeros(z.shape, bool)
    for m in masks.values():
        allz |= m
    zc = np.flatnonzero(allz.ravel())
    sl = np.radians(slope_deg); mu = np.tan(np.radians(a.mu_deg))
    near_forest = ndi.binary_dilation(forest, iterations=max(1, int(round(a.forest_buffer / a.res))))
    shape = ((np.cos(sl) ** 2 * np.clip(np.tan(sl) - mu, 0, None)
              * np.clip((60 - slope_deg) / 10, 0, 1)) * ~near_forest).ravel()[zc]
    zsel = {k: m.ravel()[zc] for k, m in masks.items()}

    # Weather: transport bins (as micromet_loading.py) and smoothed SWE
    h = hourly_weather(ROOT / 'data' / 'weather' / f'{a.site}.csv')
    h = h[~h.index.duplicated()]
    h['ut'] = [u_threshold(t, g >= 2.5) for t, g in zip(h['t'], h['gain24'].fillna(0))]
    h['dbin'] = ((h['dir'] / 10).round() * 10 % 360).astype(int)
    h['ubin'] = h['u'].round().astype(int)
    h['tbin'] = (h['ut'] * 4).round() / 4
    wmax = max(np.nanpercentile(Wd[w], 99.9) / wst[w] for w in Wd)
    ht = h[np.isfinite(h['ut']) & (h['u'] * wmax > h['ut'])]
    swe = hourly_weather(ROOT / 'data' / 'weather' / f'{a.site}.csv')['swe']
    swe = swe[~swe.index.duplicated()].asfreq('h').rolling(3, center=True, min_periods=1).median()

    # Timed naturals
    av = pd.read_csv(cfg['avalanche_csv'])
    for c in ('Trigger', 'Type', 'Elev'):
        av[c] = av[c].astype(str).str.strip()
    av = av[av['HW Path'].str.contains(cfg['avalanche_filter'], na=False)
            & (av.Trigger == 'N') & av.Type.isin(['HS', 'SS']) & av.Time.notna()]
    av = av[(av['HW Path'] != 'Seven Sister #6') | av.Elev.isin(['>TL', 'TL'])]
    av = av[~((av['Date Known'] == 'Estimated') & (av.Time.str.strip() == '12:00 AM'))]
    av['local'] = pd.to_datetime(av.Date + ' ' + av.Time.str.strip(), format='mixed')
    av['t0'] = (av.local + pd.Timedelta(hours=7)).dt.floor('h')
    av['zone'] = av['HW Path'].map(path2zone)
    av = av.dropna(subset=['zone']).sort_values('local')
    # Episode = same day; reference time = the first release that day
    av['episode'] = av.local.dt.date
    ep_t0 = av.groupby('episode').t0.min()
    print(f"{len(av)} timed naturals in {len(ep_t0)} episodes")

    H = pd.Timedelta(hours=1)
    need = pd.concat([ht.loc[t - PRE * H:t + POST * H] for t in ep_t0])
    keys = sorted(set(zip(need.dbin, need.ubin, need.tbin)))
    print(f"{len(keys)} transport bins to compute")
    Fr_cache, A_bin = {}, {}
    for wb, ub, tb in keys:
        ang = (wb + 180) % 360 - 90
        if wb not in Fr_cache:
            Fr_cache[wb] = (rotate_in(Wd[wb] / wst[wb], ang),
                            rotate_in(forest.astype(float), ang, cval=0.0) > 0.5)
        Wr, Fr = Fr_cache[wb]
        cap = capacity(Wr * ub, tb); cap[Fr] = 0.0
        A = np.nan_to_num(rotate_out(march(cap, a.res, a.lf), ang, z.shape) * 3600.0)
        A_bin[(wb, ub, tb)] = A.ravel()[zc]

    rows = []
    for ep, t0 in ep_t0.items():
        hrs = pd.date_range(t0 - PRE * H, t0 + POST * H, freq='h')
        s0 = swe.loc[t0 - (PRE + 3) * H:t0 - PRE * H].min()
        wind = np.zeros(len(zc))
        for t in hrs:
            if t in ht.index:
                r = ht.loc[t]
                wind += A_bin[(r.dbin, r.ubin, r.tbin)]
            snow = max(float(swe.loc[t0 - PRE * H:t].max() - s0), 0.0) if np.isfinite(s0) else np.nan
            load = np.clip(snow + wind, 0, None)
            for zn, sel in zsel.items():
                st = np.sort(G * load[sel] * shape[sel])
                rows.append(dict(episode=ep, t0_utc=t0, hour=(t - t0) / H, zone=zn, snow_mm=snow,
                                 stress_2nd=float(st[-2]),
                                 stress_snow_2nd=float(np.sort(G * snow * shape[sel])[-2])))
    df = pd.DataFrame(rows)
    rel = av.groupby('episode').apply(lambda g: dict(zip(g.zone, (g.t0 - ep_t0[g.name]) / H)),
                                      include_groups=False)
    df['released'] = [z_ in rel[e] for e, z_ in zip(df.episode, df.zone)]
    df['release_hour'] = [rel[e].get(z_, np.nan) for e, z_ in zip(df.episode, df.zone)]
    df.round(2).to_csv(FIGURES / f'{slug}_event_stress_timeline.csv', index=False)

    label = lambda zn: '#' + pmap[zn].split('#')[-1]

    # ── Figure 1: every released path, aligned at its release ────────────────
    fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
    rl = df[df.released].copy(); rl['rh'] = rl.hour - rl.release_hour
    for (ep, zn), g in rl.groupby(['episode', 'zone']):
        c = PATH_COLORS[label(zn)]
        pre = g[g.rh <= 0]; post = g[g.rh >= 0]
        ax.plot(pre.rh, pre.stress_2nd, color=c, lw=1.5, alpha=0.85)
        ax.plot(post.rh, post.stress_2nd, color=c, lw=1.2, alpha=0.5, ls=(0, (2, 2)))
    med = rl[rl.rh.between(-PRE, POST)].groupby(rl.rh.round())['stress_2nd'].median()
    ax.plot(med.index, med.values, color=INK, lw=3, label='median')
    ax.axvline(0, color=INK2, lw=1)
    ax.axhline(a.threshold, color=INK2, lw=1, ls='--')
    ax.text(-PRE + 0.3, a.threshold, f'≈ threshold implied by the record ({a.threshold:.0f} Pa)',
            va='bottom', color=INK2, fontsize=9)
    ax.text(0.2, ax.get_ylim()[1] * 0.97, 'release', color=INK2, fontsize=9, va='top')
    for p, c in PATH_COLORS.items():
        if (rl.zone.map(label) == p).any():
            ax.plot([], [], color=c, lw=2, label=f'Sister {p}')
    ax.plot([], [], color=INK2, lw=1.2, ls=(0, (2, 2)), label='after release')
    ax.legend(frameon=False, loc='upper left', fontsize=9)
    ax.set_xlim(-PRE, POST); ax.set_xticks(range(-PRE, POST + 1, 6))
    ax.set_xlabel('hours from release'); ax.set_ylabel('hotspot stress margin added since −24 h (Pa)')
    ax.set_title(f'Hotspot stress before {len(rl.groupby(["episode", "zone"]))} timed natural slabs '
                 '(2nd-highest cell; station snow + MicroMet wind deposition)', fontsize=11, loc='left')
    ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(FIGURES / f'{slug}_event_stress_timeline.png', dpi=150); plt.close(fig)

    # ── Figure 2: one panel per episode, all paths ────────────────────────────
    eps = list(ep_t0.index); nc = 4; nr = int(np.ceil(len(eps) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4 * nc, 2.7 * nr), sharex=True, facecolor='white')
    for ax, ep in zip(axs.ravel(), eps):
        g = df[df.episode == ep]
        for zn, gz in g.groupby('zone'):
            p = label(zn); r = bool(gz.released.iloc[0])
            ax.plot(gz.hour, gz.stress_2nd, color=PATH_COLORS[p], lw=2.2 if r else 1.0,
                    alpha=1 if r else 0.45, ls='-' if r else (0, (3, 2)))
            if r:
                rh = gz.release_hour.iloc[0]
                ax.plot(rh, np.interp(rh, gz.hour, gz.stress_2nd), 'o', ms=8, mfc=PATH_COLORS[p],
                        mec='white', mew=1.5, zorder=5)
        sn = g.groupby('hour').stress_snow_2nd.max()
        ax.fill_between(sn.index, 0, sn.values, color=GRID, zorder=0, lw=0)
        ax.axhline(a.threshold, color=INK2, lw=0.8, ls='--'); ax.axvline(0, color=INK2, lw=0.6)
        names = ', '.join(sorted(label(z_) for z_ in rel[ep]))
        ax.set_title(f"{pd.Timestamp(ep):%Y-%m-%d} {av[av.episode == ep].local.min():%H:%M}  ({names})",
                     fontsize=9, loc='left')
        ax.grid(color=GRID, lw=0.5); ax.set_axisbelow(True)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.tick_params(labelsize=8, labelbottom=True)
    for ax in axs.ravel()[len(eps):]:
        ax.axis('off')
    for ax in axs[-1]:
        ax.set_xlabel('hours from first release', fontsize=9)
        ax.set_xticks(range(-PRE, POST + 1, 6))
    for ax in axs[:, 0]:
        ax.set_ylabel('stress (Pa)', fontsize=9)
    handles = [plt.Line2D([], [], color=c, lw=2, label=f'Sister {p}') for p, c in PATH_COLORS.items()]
    handles += [plt.Line2D([], [], color=INK2, lw=2.2, label='released (dot = release time)'),
                plt.Line2D([], [], color=INK2, lw=1, ls=(0, (3, 2)), alpha=0.6, label='did not release'),
                plt.Rectangle((0, 0), 1, 1, color=GRID, label='snowfall only (no wind), steepest open cell'),
                plt.Line2D([], [], color=INK2, lw=0.8, ls='--', label=f'{a.threshold:.0f} Pa threshold')]
    fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=9)
    fig.suptitle('Hotspot stress margin, hour by hour, from 24 h before each timed natural release',
                 x=0.01, ha='left', fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(FIGURES / f'{slug}_event_stress_panels.png', dpi=130); plt.close(fig)

    # ── Per-event table: model stress + station weather around release ───────
    h['tr'] = h.u > h.ut
    def gain(t0, t1):
        x = swe.loc[t0:t1].dropna()
        return np.nan if len(x) < 3 else max(float(x.iloc[-1] - x.min()), 0.0)
    ev_rows = []
    for _, r in av.iterrows():
        t = r.t0; ep = r.episode
        g = rl[(rl.episode == ep) & (rl.zone == r.zone)].set_index('rh').stress_2nd
        at = df[(df.episode == ep) & (df.hour == (t - ep_t0[ep]) / H)]
        x6 = h.loc[t - 6 * H:t]
        ds = swe.loc[t - 72 * H:t].diff()
        rise = ds[ds > 0.3]
        ev_rows.append(dict(
            release_local=r.local, path=label(r.zone), type=r.Type, size=r.Dsize,
            time_known=r['Time Known'] if isinstance(r['Time Known'], str) else '',
            stress_at_release=g.get(0.0, np.nan), stress_m6=g.get(-6.0, np.nan), stress_m12=g.get(-12.0, np.nan),
            rise_last3h=g.get(0.0, np.nan) - g.get(-3.0, np.nan),
            rank_among_paths=int((at.stress_2nd > g.get(0.0, np.nan)).sum()) + 1,
            swe_gain_6h=gain(t - 6 * H, t), swe_gain_24h=gain(t - 24 * H, t),
            swe_gain_rest_of_day=gain(t, r.local.normalize() + pd.Timedelta(days=1, hours=7)),
            h_since_snow=(t - rise.index[-1]) / H if len(rise) else np.nan,
            transport_h_last6=int(x6.tr.sum()), temp_c=h.t.get(t, np.nan)))
    ev = pd.DataFrame(ev_rows)
    ev.round(1).to_csv(FIGURES / f'{slug}_timed_naturals.csv', index=False)
    pd.set_option('display.width', 220)
    print(ev.round(1).to_string(index=False))
    n = len(ev)
    print(f"\nmedian stress: −12 h {ev.stress_m12.median():.0f}, −6 h {ev.stress_m6.median():.0f}, "
          f"release {ev.stress_at_release.median():.0f} Pa; ≥ {a.threshold:.0f} Pa at release "
          f"{(ev.stress_at_release >= a.threshold).sum()}/{n}; rising (> 1 Pa in last 3 h) "
          f"{(ev.rise_last3h > 1).sum()}/{n}; released path most stressed {(ev.rank_among_paths == 1).sum()}/{n}")
    print(f"station: transport ≥ 5 of last 6 h {(ev.transport_h_last6 >= 5).sum()}/{n}; snow in last 6 h "
          f"{(ev.swe_gain_6h >= 2.5).sum()}/{n}; snow ≤ 6 h before {(ev.h_since_snow <= 6).sum()}/{n}; "
          f"more snow later that day {(ev.swe_gain_rest_of_day > 0).sum()}/{n} "
          f"(more than the 24 h before: {(ev.swe_gain_rest_of_day > ev.swe_gain_24h).sum()}); "
          f"temperature {ev.temp_c.min():.0f} to {ev.temp_c.max():.0f} °C")
    w = h[h.index.month.isin([11, 12, 1, 2, 3, 4, 5])]
    tr6 = w.tr.astype(float).rolling(6).sum()
    print(f"winter baseline: transport ≥ 5 of 6 h in {(tr6 >= 5).mean():.0%} of hours")


if __name__ == '__main__':
    main()
