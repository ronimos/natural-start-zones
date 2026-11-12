"""
assign_avalanche_paths.py — match avalanche records that have lat/lon but no
path name to the site's paths.

Path corridors: every cell of each start-zone polygon is traced downhill by
steepest descent (D8 on the 3DEP DEM, lightly smoothed) for up to `--trace`
metres, so a corridor = start zone + track.  Each record is assigned to the
nearest corridor when it lies within `--max-dist` metres; the distance to the
second-nearest corridor is reported so ambiguous cases are visible; ties (< 15 m apart) go to the path whose start-zone
centre is nearest across the slope (east–west; these paths run north).  Records
farther than `--max-dist` from every corridor are left unassigned (another
path, or debris that spread outside the traced track).  At Seven Sisters the
five digitised zones are all the active paths (there is no Sister #5, and #7
is inactive).

Usage:
  python src/assign_avalanche_paths.py data/avalanches/caic-avalanche-records-2026-11-12.csv
"""
from __future__ import annotations

import argparse
import html
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

gdal.UseExceptions()
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, _sz_pairs     # noqa: E402
import wind_exposure as WX                          # noqa: E402

OFFS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def d8_receivers(z: np.ndarray, res: float):
    """Index offset (dr, dc) of the steepest-descent neighbour; (0, 0) at pits."""
    best = np.zeros_like(z); rdr = np.zeros(z.shape, int); rdc = np.zeros(z.shape, int)
    for dr, dc in OFFS:
        drop = (z - WX._shift(z, dr, dc)) / (res * np.hypot(dr, dc))
        better = np.nan_to_num(drop, nan=-np.inf) > best
        best = np.where(better, drop, best); rdr[better] = dr; rdc[better] = dc
    return rdr, rdc


def trace(start: np.ndarray, rdr, rdc, steps: int) -> np.ndarray:
    out = start.copy()
    r, c = np.nonzero(start)
    for _ in range(steps):
        dr, dc = rdr[r, c], rdc[r, c]
        moving = (dr != 0) | (dc != 0)
        if not moving.any():
            break
        r, c = r[moving] + dr[moving], c[moving] + dc[moving]
        ok = (r >= 0) & (r < start.shape[0]) & (c >= 0) & (c < start.shape[1])
        r, c = r[ok], c[ok]
        out[r, c] = True
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('csv')
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--res', type=float, default=2.0)
    ap.add_argument('--trace', type=float, default=1200.0, help='max downhill trace (m)')
    ap.add_argument('--max-dist', type=float, default=40.0, help='max distance to assign (m)')
    a = ap.parse_args()

    cfg = SITES[a.site]; d = cfg['data_dir']
    slug = cfg['name'].lower().replace(' ', '_')
    src = d / cfg['source_raster']
    proj4 = osr.SpatialReference(wkt=gdal.Open(str(src)).GetProjection()).ExportToProj4()

    rec = pd.read_csv(a.csv, encoding='utf-8-sig')
    rec = rec.apply(lambda s: s.map(lambda v: html.unescape(v) if isinstance(v, str) else v))
    pts = gp.GeoDataFrame(rec, geometry=gp.points_from_xy(rec['Longitude'], rec['latitude']),
                          crs=4326).to_crs(proj4)
    zones = {kml.stem.replace('_SZ', '').replace('_', ' #'): gp.read_file(kml).to_crs(proj4).geometry.union_all()
             for kml, _ in _sz_pairs(cfg)}
    u = gp.GeoSeries(list(zones.values())).union_all()
    x0, y0, x1, y1 = u.buffer(1500).bounds
    z, gt, wkt = WX.load_domain(src, (x0, y0, x1, y1), a.res)
    zs = ndi.gaussian_filter(np.nan_to_num(z, nan=np.nanmax(z)), 1.0)
    rdr, rdc = d8_receivers(zs, a.res)

    corridors, dist = {}, {}
    for name, g in zones.items():
        m = WX.polygon_mask(g, z.shape, gt, wkt)
        corridors[name] = trace(m, rdr, rdc, int(a.trace / a.res))
        dist[name] = ndi.distance_transform_edt(~corridors[name]) * a.res

    cols = ((pts.geometry.x - gt[0]) / gt[1]).astype(int).values
    rows = ((pts.geometry.y - gt[3]) / gt[5]).astype(int).values
    inside = (rows >= 0) & (rows < z.shape[0]) & (cols >= 0) & (cols < z.shape[1])
    pts = pts[inside].copy(); rows, cols = rows[inside], cols[inside]
    D = pd.DataFrame({n: dist[n][rows, cols] for n in zones}, index=pts.index)
    in_sz = pd.DataFrame({n: zones[n].contains(pts.geometry) for n in zones}, index=pts.index)
    order = np.argsort(D.values, axis=1)
    pts['nearest'] = [D.columns[o[0]] for o in order]
    pts['dist_m'] = D.values[np.arange(len(D)), order[:, 0]].round(0)
    pts['second'] = [D.columns[o[1]] for o in order]
    pts['second_m'] = D.values[np.arange(len(D)), order[:, 1]].round(0)
    pts['in_start_zone'] = [next((n for n in zones if in_sz.loc[i, n]), '') for i in pts.index]
    pts['dem_elev_m'] = z[rows, cols].round(0)
    pts['assigned'] = np.where(pts['dist_m'] <= a.max_dist, pts['nearest'], '')
    pts['ambiguous'] = (pts['assigned'] != '') & (pts['second_m'] - pts['dist_m'] < 15)
    # Tracks merge low on the slope (and run along the highway once they reach
    # it), so corridor distance can tie.  These paths run due north, so ties
    # (< 15 m apart) go to the path whose start-zone centre is nearest in the
    # cross-slope (east–west) direction.  Checked against named CAIC records:
    # 2021-02-06 is Sister #3, which this rule picks.
    cx = {n: g.centroid.x for n, g in zones.items()}
    tied = pts['ambiguous']
    pts.loc[tied, 'assigned'] = [min((pts.loc[i, 'nearest'], pts.loc[i, 'second']),
                                     key=lambda n: abs(pts.loc[i].geometry.x - cx[n]))
                                 for i in pts.index[tied]]
    sz_d = pd.DataFrame({n: pts.distance(g) for n, g in zones.items()}, index=pts.index)
    pts['dist_to_sz_m'] = [round(float(sz_d.loc[i, pts.loc[i, 'assigned']]), 0) if pts.loc[i, 'assigned'] else np.nan
                           for i in pts.index]

    show = pts[pts['dist_m'] <= 300].sort_values('Date')
    cols_out = ['Date', 'Date Known', 'Landmark', 'latitude', 'Longitude', 'Elevation', 'Aspect', 'Type',
                'Trigger', 'Relative Size', 'Destructive Size', 'dem_elev_m', 'in_start_zone',
                'nearest', 'dist_m', 'second', 'second_m', 'assigned', 'ambiguous', 'dist_to_sz_m']
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 30)
    print(f"{len(pts)} records inside the {a.trace:.0f} m-trace domain; "
          f"{(pts.assigned != '').sum()} within {a.max_dist:.0f} m of a Sister corridor\n")
    print(show[cols_out].to_string(index=False))
    out = FIGURES / f'{slug}_new_records_assigned.csv'
    show[cols_out].to_csv(out, index=False)

    hs = gdal.DEMProcessing('', gdal.Warp('', str(src), format='MEM', outputBounds=(x0, y0, x1, y1),
                            xRes=a.res, yRes=a.res, targetAlignedPixels=True, resampleAlg='cubic'),
                            'hillshade', format='MEM').ReadAsArray()
    allc = np.zeros(z.shape, bool)
    for m in corridors.values():
        allc |= m
    rr, cc = np.nonzero(ndi.binary_dilation(allc, iterations=100))
    win = (slice(rr.min(), rr.max()), slice(cc.min(), cc.max()))
    fig, ax = plt.subplots(figsize=(11, 11))
    ax.imshow(hs[win], cmap='gray')
    colors = dict(zip(zones, plt.cm.tab10(np.arange(len(zones)))))
    for n, m in corridors.items():
        ax.contour(m[win], levels=[.5], colors=[colors[n]], linewidths=1.5)
        ax.contour(WX.polygon_mask(zones[n], z.shape, gt, wkt)[win], levels=[.5], colors=[colors[n]], linewidths=3)
        r_, c_ = np.argwhere(WX.polygon_mask(zones[n], z.shape, gt, wkt)[win]).mean(axis=0)
        ax.text(c_, r_ + 40, n, color='w', fontsize=9, ha='center', bbox=dict(fc=colors[n], alpha=.85, lw=0))
    for _, p in show.iterrows():
        c_ = (p.geometry.x - gt[0]) / gt[1] - win[1].start; r_ = (p.geometry.y - gt[3]) / gt[5] - win[0].start
        col = colors.get(p['assigned'], 'k') if p['assigned'] else 'k'
        ax.plot(c_, r_, 'o', ms=8, mfc=col, mec='w', mew=1.5)
        ax.annotate(str(p['Date'])[:10], (c_, r_), xytext=(5, 5), textcoords='offset points', fontsize=7, color='yellow')
    ax.set_title(f"{cfg['name']} — new CAIC records (dots; black = not within {a.max_dist:.0f} m of a "
                 f"Sister corridor)\nthick = start zone, thin = traced track; north up")
    ax.set_axis_off()
    fig.savefig(FIGURES / f'{slug}_new_records_map.png', dpi=140, bbox_inches='tight'); plt.close(fig)
    print(f"\nSaved: {out.name}, {slug}_new_records_map.png")


if __name__ == '__main__':
    main()
