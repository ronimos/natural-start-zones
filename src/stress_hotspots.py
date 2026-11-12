"""
stress_hotspots.py — where can slab nucleation patches start?

Nucleation is local: a patch starts where the weak-layer shear stress
exceeds a local threshold.  Zone-average steepness dilutes exactly the spots
that matter, so this maps stress cell by cell and counts high-stress patches.

Per cell x (1.2 m analysis grid):

  Load        SWE(x) ∝ [1 + k · D(x) / D_ref] · h(θ)
    D(x)      wind-deposition index (m of effective upwind fetch), weighted
              over wind directions by ridge transport hours.  For direction d,
              trace upwind from x: every exposed cell (Sx_d ≤ 0) adds its
              length, attenuated by exp(−s_sheltered / λ), where s_sheltered
              is the sheltered distance the snow must cross to reach x (it
              drops out in sheltered terrain).  Only sheltered cells
              (Sx_d(x) > 0) receive; rays stop at forest or after `lmax`.
              This puts deposition just downwind of exposed→sheltered
              transitions: rollovers and lee walls of gullies.
    h(θ)      snow-holding factor: 1 up to 50°, linear to 0 at 60° (snow
              sluffs instead of building a slab on very steep terrain).

  Stress      margin(x) = SWE(x) · cos²θ · (tanθ − μ)⁺
              = (τ − μσ) per unit load for a friction-type weak layer.
              Slope θ from Horn's method on the FULL DEM, then masked per
              zone (clipping first biases edge cells low).

  Hotspots    cells in the top `pct` % of margin across all start-zone cells
              pooled (zones compete on one scale), grouped into 8-connected
              patches of ≥ `min_cells`.

Per zone: hotspot area and fraction, patch count, P99 margin, and the share
of hotspot area in the west / upper half of the zone (a check against field
observations of where releases start).

Parameters are fixed in advance and reported with sensitivity runs.

Usage:  python src/stress_hotspots.py --site sisters
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
from scipy.stats import spearmanr

gdal.UseExceptions()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, ROOT, _sz_pairs           # noqa: E402
from avl_data import load_avalanches                            # noqa: E402
import wind_exposure as W                                       # noqa: E402


def deposition_index(sx: dict, forest: np.ndarray, weights: pd.Series,
                     targets: tuple, res: float, lmax: float, lam: float) -> np.ndarray:
    """D at target cells (rows, cols of the wind grid), Σ_d w_d · dep_d."""
    r0, c0 = targets
    h, w = forest.shape
    D = np.zeros(len(r0))
    for d, wt in weights.items():
        if wt <= 0:
            continue
        sx_d = sx[d]
        th = np.radians(d)
        dr, dc = -np.cos(th), np.sin(th)
        contrib = np.zeros(len(r0)); shel = np.zeros(len(r0))
        alive = np.ones(len(r0), bool)
        for s in np.arange(res, lmax + res / 2, res):
            r = np.rint(r0 + dr * s / res).astype(int)
            c = np.rint(c0 + dc * s / res).astype(int)
            alive &= (r >= 0) & (r < h) & (c >= 0) & (c < w)
            idx = np.flatnonzero(alive)
            if not len(idx):
                break
            rr, cc = r[idx], c[idx]
            v = sx_d[rr, cc]
            dead = forest[rr, cc] | np.isnan(v)
            exposed = (v <= 0) & ~dead
            contrib[idx[exposed]] += res * np.exp(-shel[idx[exposed]] / lam)
            shel[idx[~exposed & ~dead]] += res
            alive[idx[dead]] = False
            alive[idx[shel[idx] > 5 * lam]] = False          # < 1 % would arrive
        receives = sx_d[r0, c0] > 0
        D += wt * contrib * receives
    return D


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--wind-res', type=float, default=2.0)
    ap.add_argument('--lmax', type=float, default=1000.0)
    ap.add_argument('--lam', type=float, default=50.0, help='drop-out length over sheltered terrain (m)')
    ap.add_argument('--k', type=float, default=1.0, help='load gain per D_ref of deposition index')
    ap.add_argument('--d-ref', type=float, default=100.0)
    ap.add_argument('--mu-deg', type=float, default=30.0, help='weak-layer friction angle')
    ap.add_argument('--pct', type=float, default=5.0, help='hotspot = top pct %% of pooled SZ cells')
    ap.add_argument('--min-cells', type=int, default=3)
    a = ap.parse_args()

    cfg = SITES[a.site]; d = cfg['data_dir']
    slug = cfg['name'].lower().replace(' ', '_')
    src = d / cfg['source_raster']
    dem_path = d / (cfg['crp_dem_name'] or cfg['dem_name'])

    # ── Analysis grid (1.2 m) and slope from the full DEM ────────────────────
    ds = gdal.Open(str(dem_path)); gt = ds.GetGeoTransform(); wkt = ds.GetProjection()
    z = ds.ReadAsArray().astype('f8'); z[z == -9999] = np.nan
    res = gt[1]
    slope = np.radians(gdal.DEMProcessing('', ds, 'slope', format='MEM',
                                          computeEdges=True).ReadAsArray().astype('f8'))
    slope[np.isnan(z)] = np.nan
    H, Wd = z.shape
    proj4 = osr.SpatialReference(wkt=wkt).ExportToProj4()
    zones = {k.stem: gp.read_file(k).to_crs(proj4).geometry.union_all() for k, _ in _sz_pairs(cfg)}
    masks = {k: W.polygon_mask(g, z.shape, gt, wkt) for k, g in zones.items()}

    # ── Wind grid (2 m, buffered) — Sx, canopy, transport weights ────────────
    x0, y1 = gt[0], gt[3]; x1, y0 = x0 + Wd * res, y1 - H * res
    wb = (x0 - 1300, y0 - 1300, x1 + 1300, y1 + 1300)
    zw, gtw, wktw = W.load_domain(src, wb, a.wind_res)
    dsm = d / '260922_Sisters_Full_Ground-DSM.tiff'
    forest = (W.canopy_mask(dsm, src, wb, a.wind_res, zw.shape) if dsm.exists()
              else np.zeros(zw.shape, bool))
    print(f"wind grid {zw.shape[1]}×{zw.shape[0]} at {a.wind_res} m; Sx for 72 directions …")
    sx = W.sx_sectors(zw, a.wind_res, 100.0)
    hours, flux, n_tr, _ = W.transport_weights(ROOT / 'data' / 'weather' / f'{a.site}.csv', 15.0, 0.3)
    wh = hours / hours.sum()

    # Target cells: wind-grid cells covering the analysis grid
    rr, cc = np.mgrid[0:zw.shape[0], 0:zw.shape[1]]
    xs = gtw[0] + (cc + .5) * gtw[1]; ys = gtw[3] + (rr + .5) * gtw[5]
    tgt = (xs >= x0) & (xs <= x1) & (ys >= y0) & (ys <= y1)
    print(f"deposition index for {tgt.sum()} cells …")
    Dw = np.full(zw.shape, np.nan)
    Dw[tgt] = deposition_index(sx, forest, wh, (rr[tgt], cc[tgt]), a.wind_res, a.lmax, a.lam)

    # Resample D onto the 1.2 m analysis grid
    mem = gdal.GetDriverByName('MEM').Create('', zw.shape[1], zw.shape[0], 1, gdal.GDT_Float64)
    mem.SetGeoTransform(gtw); mem.SetProjection(wktw); mem.GetRasterBand(1).WriteArray(Dw)
    mem.GetRasterBand(1).SetNoDataValue(np.nan)
    D = gdal.Warp('', mem, format='MEM', outputBounds=(x0, y0, x1, y1), width=Wd, height=H,
                  resampleAlg='bilinear', dstSRS=wkt).ReadAsArray().astype('f8')

    # ── Stress margin ────────────────────────────────────────────────────────
    mu = np.tan(np.radians(a.mu_deg))
    deg = np.degrees(slope)
    hold = np.clip((60 - deg) / 10, 0, 1)
    shape_term = np.cos(slope) ** 2 * np.clip(np.tan(slope) - mu, 0, None)

    def margin(k):
        return (1 + k * np.nan_to_num(D) / a.d_ref) * hold * shape_term

    def zone_metrics(mg, pct):
        allsz = np.zeros(z.shape, bool)
        for m in masks.values():
            allsz |= m
        thr = np.nanpercentile(mg[allsz], 100 - pct)
        hot = (mg >= thr) & allsz
        lab, _ = ndi.label(hot, structure=np.ones((3, 3)))
        sizes = np.bincount(lab.ravel())
        keep = np.isin(lab, np.flatnonzero(sizes >= a.min_cells)[1:] if len(sizes) > 1 else [])
        keep &= lab > 0
        out = {}
        for name, m in masks.items():
            cells = np.argwhere(m); cr, ccn = cells.mean(axis=0)
            hz = keep & m
            nlab = len(np.unique(lab[hz]))
            hr, hc = np.nonzero(hz)
            out[name] = dict(
                hot_area_m2=hz.sum() * res ** 2,
                hot_frac=hz.sum() / m.sum(),
                patches=nlab,
                p99_margin=float(np.nanpercentile(mg[m], 99)),
                hot_west_share=float((hc < ccn).mean()) if len(hc) else np.nan,
                hot_upper_share=float((hr > cr).mean()) if len(hr) else np.nan,  # slope faces N: uphill = south = larger row
            )
        return pd.DataFrame(out).T, keep, thr

    # Avalanche counts (pipeline filter: Sister #6 TL counts as alpine)
    nat = load_avalanches(str(cfg['avalanche_csv']), cfg['avalanche_filter'],
                          type=['HS', 'SS'], trig=['N']).groupby('HW Path').size()
    pmap = {k.replace(cfg['sz_dem_suffix'], ''): v for k, v in cfg['sz_avi_path_map'].items()}
    area = pd.Series({k: m.sum() * res ** 2 for k, m in masks.items()})
    nat_n = pd.Series({k: nat.get(pmap.get(k, k), 0) for k in masks})

    rows, main_tbl, main_keep, main_mg = [], None, None, None
    for k in (0.0, 0.25, 1.0, 4.0):
        mg = margin(k)
        for pct in (1.0, 5.0, 10.0):
            t, keep, thr = zone_metrics(mg, pct)
            t['nat'] = nat_n; t['nat_per_10k'] = nat_n / area * 1e4
            for metric, outcome in (('patches', 'nat'), ('hot_area_m2', 'nat'),
                                    ('hot_frac', 'nat_per_10k'), ('p99_margin', 'nat_per_10k')):
                for lab_, sub in (('N=5', t), ('N=4 (no #6)', t.drop(index=[i for i in t.index if '6' in i]))):
                    rho, p = spearmanr(sub[metric], sub[outcome])
                    rows.append(dict(k=k, pct=pct, metric=metric, outcome=outcome, set=lab_,
                                     rho=round(float(rho), 2), p=round(float(p), 3)))
            if k == a.k and pct == a.pct:
                main_tbl, main_keep, main_mg = t, keep, mg
    corr = pd.DataFrame(rows)

    pd.set_option('display.width', 220)
    print(f"\nMain setting: k={a.k}, top {a.pct}% hotspots, μ=tan{a.mu_deg:.0f}°, λ={a.lam} m, "
          f"patches ≥ {a.min_cells} cells ({a.min_cells * res ** 2:.1f} m²)")
    print(main_tbl.round(3).to_string())
    print("\nSpearman ρ (rows: k = deposition weight, 0 = slope only):")
    print(corr.pivot_table(index=['metric', 'outcome', 'set'], columns=['k', 'pct'],
                           values='rho', sort=False).to_string())

    # ── Outputs ──────────────────────────────────────────────────────────────
    main_tbl.to_csv(FIGURES / f'{slug}_stress_hotspots.csv')
    corr.to_csv(FIGURES / f'{slug}_stress_hotspots_corr.csv', index=False)

    hs = gdal.DEMProcessing('', ds, 'hillshade', format='MEM').ReadAsArray()
    allsz = np.zeros(z.shape, bool)
    for m in masks.values():
        allsz |= m
    ys_, xs_ = np.nonzero(ndi.binary_dilation(allsz, iterations=40))
    win = (slice(ys_.min(), ys_.max()), slice(xs_.min(), xs_.max()))
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    for ax, field, title, cmap in ((axs[0], np.nan_to_num(D), 'Wind-deposition index D (m)', 'Blues'),
                                   (axs[1], main_mg, f'Stress margin (k={a.k}); red = top {a.pct}% patches', 'magma')):
        ax.imshow(hs[win], cmap='gray')
        im = ax.imshow(np.ma.masked_where(field[win] <= 0, field[win]), cmap=cmap, alpha=.6)
        fig.colorbar(im, ax=ax, fraction=.03)
        for name, m in masks.items():
            ax.contour(m[win], levels=[.5], colors='cyan', linewidths=1.2)
            r_, c_ = np.argwhere(m[win]).mean(axis=0)
            ax.text(c_, r_ - 25, name.replace('_SZ', '').replace('_', ' #'), color='w', fontsize=9,
                    ha='center', bbox=dict(fc='k', alpha=.6, lw=0))
        if ax is axs[1]:
            ax.contourf(main_keep[win], levels=[.5, 1.5], colors=['red'], alpha=.9)
        ax.set_title(title); ax.set_axis_off()
    fig.suptitle(f"{cfg['name']} — nucleation hotspots (north up)")
    fig.savefig(FIGURES / f'{slug}_stress_hotspots.png', dpi=150, bbox_inches='tight'); plt.close(fig)

    m1 = masks[[k for k in masks if '1' in k][0]]
    ys_, xs_ = np.nonzero(ndi.binary_dilation(m1, iterations=15))
    w1 = (slice(ys_.min(), ys_.max()), slice(xs_.min(), xs_.max()))
    fig, ax = plt.subplots(figsize=(6, 9))
    ax.imshow(hs[w1], cmap='gray')
    im = ax.imshow(np.ma.masked_where(~m1[w1], main_mg[w1]), cmap='magma', alpha=.7)
    ax.contourf(main_keep[w1] & m1[w1], levels=[.5, 1.5], colors=['red'], alpha=.9)
    ax.contour(m1[w1], levels=[.5], colors='cyan')
    fig.colorbar(im, ax=ax, fraction=.04, label='stress margin')
    ax.set_title('Sister #1 — stress margin and hotspot patches (red)\nnorth up: top of start zone at bottom')
    ax.set_axis_off()
    fig.savefig(FIGURES / f'{slug}_sister1_hotspots.png', dpi=150, bbox_inches='tight'); plt.close(fig)
    print(f"\nSaved: {slug}_stress_hotspots.csv/_corr.csv, {slug}_stress_hotspots.png, "
          f"{slug}_sister1_hotspots.png")


if __name__ == '__main__':
    main()
