"""
storm_analysis.py — weather before each natural slab avalanche.

For every natural (N) HS/SS avalanche on the site's paths, summarise the
72 h before the end of the avalanche date (local Mountain time; station
times are UTC):

  swe_gain_in      SNOTEL SWE gain (new-snow water, inches)
  transport_h      ridge hours with u ≥ u_thresh and T < 32 °F
  flux_SW / flux_WNW   share of Σu³ transport flux from SW (200–260°) vs
                   W–NW (265–330°) — which sector did the loading
  u3_mean_dir      u³-weighted mean wind direction

and compare with all winter days (Dec–Apr) as the climatological baseline:
the percentile of each event's SWE gain and transport hours.

Output: figures/<site>_storm_events.csv and a printed summary.
Usage:  python src/storm_analysis.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES, ROOT   # noqa: E402

SECTORS = {'SW': (200, 260), 'WNW': (265, 330)}


def load_weather(csv: Path) -> pd.DataFrame:
    w = pd.read_csv(csv, index_col='time', parse_dates=True)
    out = pd.DataFrame(index=w.index)
    out['u'] = w['wspd'] / 10.0                    # mph
    out['dir'] = w['wdir']
    out['t'] = w['temp'] / 10.0                    # °F
    swe = (w['snowwater'] / 100.0).where(lambda x: x > 0)   # inches
    out['swe'] = swe.interpolate(limit=6)
    return out


def window_stats(w: pd.DataFrame, end: pd.Timestamp, hours: int, u_thresh: float) -> dict:
    x = w.loc[end - pd.Timedelta(hours=hours):end]
    if x['swe'].notna().sum() < hours * 0.5 or x['u'].notna().sum() < hours * 0.5:
        return dict(swe_gain_in=np.nan, transport_h=np.nan, u3_mean_dir=np.nan,
                    **{f'flux_{k}': np.nan for k in SECTORS})
    swe = x['swe'].dropna()
    gain = float(max(swe.iloc[-1] - swe.min(), 0.0)) if len(swe) else np.nan
    tr = x[(x.u >= u_thresh) & (x.t < 32) & x.dir.notna()]
    f = (tr.u * 0.44704) ** 3
    tot = f.sum()
    shares = {f'flux_{k}': (float(f[(tr.dir >= lo) & (tr.dir <= hi)].sum() / tot) if tot > 0 else 0.0)
              for k, (lo, hi) in SECTORS.items()}
    if tot > 0:
        a = np.radians(tr.dir)
        mdir = float(np.degrees(np.arctan2((f * np.sin(a)).sum(), (f * np.cos(a)).sum())) % 360)
    else:
        mdir = np.nan
    return dict(swe_gain_in=gain, transport_h=int(len(tr)), u3_mean_dir=mdir, **shares)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--hours', type=int, default=72)
    ap.add_argument('--u-thresh', type=float, default=15.0, help='transport wind (mph)')
    a = ap.parse_args()
    cfg = SITES[a.site]
    slug = cfg['name'].lower().replace(' ', '_')
    w = load_weather(ROOT / 'data' / 'weather' / f'{a.site}.csv')

    av = pd.read_csv(cfg['avalanche_csv'])
    for c in ('Trigger', 'Type', 'Elev', 'Date Known'):
        av[c] = av[c].astype(str).str.strip()
    av = av[av['HW Path'].str.contains(cfg['avalanche_filter'], na=False)
            & (av.Trigger == 'N') & av.Type.isin(['HS', 'SS'])].copy()
    av['Date'] = pd.to_datetime(av['Date'], format='mixed')
    # End of the local (MST, UTC−7) avalanche day, in UTC.
    av['end_utc'] = av['Date'].dt.normalize() + pd.Timedelta(days=1, hours=7)
    # Flag the records the pipeline's Sister #6 elevation filter removes.
    av['pipeline_kept'] = ~((av['HW Path'] == 'Seven Sister #6') & ~av['Elev'].isin(['>TL', 'TL']))

    ev = pd.DataFrame([dict(path=r['HW Path'], date=r['Date'].date(), date_known=r['Date Known'],
                            elev=r['Elev'], dsize=r.get('Dsize'), pipeline_kept=r['pipeline_kept'],
                            **window_stats(w, r['end_utc'], a.hours, a.u_thresh))
                       for _, r in av.iterrows()]).sort_values(['date', 'path'])

    # Climatology: every Dec–Apr day, same window.
    days = pd.date_range(w.index.min().normalize() + pd.Timedelta(days=4),
                         w.index.max().normalize(), freq='D')
    days = days[days.month.isin([12, 1, 2, 3, 4])]
    clim = pd.DataFrame([window_stats(w, d + pd.Timedelta(days=1, hours=7), a.hours, a.u_thresh)
                         for d in days], index=days).dropna(subset=['swe_gain_in'])
    for col in ('swe_gain_in', 'transport_h'):
        ev[f'{col}_pctile'] = [float((clim[col] <= v).mean() * 100) if pd.notna(v) else np.nan
                               for v in ev[col]]

    pd.set_option('display.width', 220)
    cols = ['date', 'path', 'date_known', 'elev', 'pipeline_kept', 'swe_gain_in', 'swe_gain_in_pctile',
            'transport_h', 'transport_h_pctile', 'flux_SW', 'flux_WNW', 'u3_mean_dir']
    print(f"Natural HS/SS avalanches — {a.hours} h before end of avalanche day\n")
    print(ev[cols].round(2).to_string(index=False))

    print(f"\nClimatology (Dec–Apr days, n={len(clim)}): median SWE gain {clim.swe_gain_in.median():.2f} in, "
          f"median transport {clim.transport_h.median():.0f} h, mean SW flux share "
          f"{clim.flux_SW[clim.transport_h > 0].mean():.2f}, WNW {clim.flux_WNW[clim.transport_h > 0].mean():.2f}")
    ev['group'] = np.where(ev.path == 'Seven Sister #6', 'Sister #6', 'Sisters #1–#4')
    g = ev.groupby('group').agg(n=('path', 'size'), swe_gain_median=('swe_gain_in', 'median'),
                                swe_pctile_median=('swe_gain_in_pctile', 'median'),
                                transport_pctile_median=('transport_h_pctile', 'median'),
                                flux_SW_mean=('flux_SW', 'mean'), flux_WNW_mean=('flux_WNW', 'mean'))
    print("\nBy group:"); print(g.round(2).to_string())
    out = FIGURES / f'{slug}_storm_events.csv'
    ev.drop(columns='group').to_csv(out, index=False)
    print(f"\nSaved: {out.name}")


if __name__ == '__main__':
    main()
