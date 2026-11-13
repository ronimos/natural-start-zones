"""
loading_rate.py — does the rate of loading matter, beyond the amount?

Field experience: 30 cm in one day is far less stable than 10 cm a day for
three days.  In the nucleation picture (Gvirtzman et al. 2025), a patch
grows for T ≈ α / (τ − τ_thresh) before dynamic rupture; fast loading raises
τ − τ_thresh quickly, shortening T and leaving the weak layer less time to
strengthen by sintering.  Prediction: at the same 72-h stress, naturals are
more likely when the load arrived in one day than when it was spread out.

Input: the zone-day table from micromet_loading.py — 24-, 48- and 72-h
window scenarios, and for the 72-h window the stress from each single day's
load (margin_2nd_d0/d1/d2: 2nd-highest cell, d0 = the release day) and the
daily new snow (swe_d0/d1/d2).

Measures (per zone-day):
  m72        2nd-highest-cell stress margin from the 72-h load (Pa)
  peak24     largest single-day stress margin of the three days (Pa)
  conc       peak24 / (sum of the three single-day margins): 1 = all load
             in one day, 1/3 = evenly spread; timing-robust (estimated
             release dates make "which day" unreliable, "how concentrated" not)
  conc_snow  the same for station snowfall (SWE)

Tests:
  1. AUC (all zone-days, and within natural days) for 24-, 48-, 72-h
     windows and peak24, with season-block bootstrap CIs.
  2. Release rate by 72-h stress bin × concentration class.
  3. Logistic regression natural ~ log(1 + m72) + conc (and + conc_snow),
     standard errors clustered by season.

Usage:  python src/loading_rate.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES          # noqa: E402

rng = np.random.default_rng(20261113)


def auc(x, y):
    pos, neg = x[y == 1], x[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float(((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean()))


def within_auc(d, col):
    w = [auc(g[col].values, g.natural.values) for _, g in d[d.day.isin(d.loc[d.natural == 1, 'day'])].groupby('day')
         if g.natural.sum() and (g.natural == 0).sum()]
    return float(np.nanmean(w)) if w else np.nan


def season_of(d):
    return d.year + (1 if d.month >= 10 else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--n-boot', type=int, default=1000)
    ap.add_argument('--min-m72', type=float, default=30.0,
                    help='72-h stress (Pa) above which zone-days count as loaded (data-derived threshold)')
    a = ap.parse_args()
    cfg = SITES[a.site]; slug = cfg['name'].lower().replace(' ', '_')
    pd.set_option('display.width', 220)

    zd = pd.read_csv(FIGURES / f'{slug}_zone_day_stress.csv', parse_dates=['day'])
    if 'margin_2nd_d0' not in zd:
        sys.exit('zone-day table has no single-day columns — rerun micromet_loading.py')
    w = {s: zd[zd.scenario == s].set_index(['day', 'zone']) for s in ('24h window', '48h window', '72h window')}
    d = w['72h window'][['natural', 'margin_2nd', 'swe_d0', 'swe_d1', 'swe_d2',
                         'margin_2nd_d0', 'margin_2nd_d1', 'margin_2nd_d2']].rename(columns={'margin_2nd': 'm72'})
    d['m24'] = w['24h window']['margin_2nd']; d['m48'] = w['48h window']['margin_2nd']
    parts = d[['margin_2nd_d0', 'margin_2nd_d1', 'margin_2nd_d2']]
    d['peak24'] = parts.max(axis=1)
    tot = parts.sum(axis=1)
    d['conc'] = np.where(tot > 0, d['peak24'] / tot.where(tot > 0), np.nan)
    sw = d[['swe_d0', 'swe_d1', 'swe_d2']]
    d['conc_snow'] = np.where(sw.sum(axis=1) > 0, sw.max(axis=1) / sw.sum(axis=1).where(sw.sum(axis=1) > 0), np.nan)
    d['peak_day'] = parts.values.argmax(axis=1)            # 0 = release day, 1 = day before, 2 = two days before
    d = d.reset_index()
    d['season'] = d['day'].map(season_of)
    n_nat = int(d.natural.sum())
    print(f"{len(d)} zone-days, {n_nat} natural zone-days, {d.season.nunique()} seasons\n")

    # ── 1. Which window predicts best? ────────────────────────────────────────
    seasons = d.season.unique(); by_s = {s: g for s, g in d.groupby('season')}
    rows = []
    for col, name in (('m24', '24-h window'), ('m48', '48-h window'), ('m72', '72-h window'),
                      ('peak24', 'largest single day of the 72 h')):
        bs = []
        for _ in range(a.n_boot):
            b = pd.concat([by_s[s] for s in rng.choice(seasons, len(seasons))])
            bs.append(auc(b[col].values, b.natural.values))
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        rows.append(dict(load=name, AUC=round(auc(d[col].values, d.natural.values), 3),
                         CI=f'{lo:.2f}–{hi:.2f}', AUC_within_day=round(within_auc(d, col), 3)))
    print("Stress margin (2nd-highest cell) by loading window — AUC for natural zone-days:")
    print(pd.DataFrame(rows).to_string(index=False))

    nd = d[d.natural == 1]
    print(f"\nOn natural days, the heaviest of the three days was: release day {int((nd.peak_day == 0).sum())}, "
          f"day before {int((nd.peak_day == 1).sum())}, two days before {int((nd.peak_day == 2).sum())} "
          f"(of {len(nd)}; dates are partly estimated)")

    # ── 2. Release rate by 72-h stress × concentration ───────────────────────
    L = d[d.m72 >= a.min_m72].copy()
    L['load bin'] = pd.cut(L.m72, [a.min_m72, 60, 110, np.inf],
                           labels=[f'{a.min_m72:.0f}–60 Pa', '60–110 Pa', '> 110 Pa'], right=False)
    L['loading'] = pd.cut(L.conc, [0, 0.5, 0.8, 1.01],
                          labels=['spread (peak day < 50 %)', 'mixed (50–80 %)', 'one day (≥ 80 %)'], right=False)
    t = L.groupby(['load bin', 'loading'], observed=False).agg(zone_days=('natural', 'size'),
                                                               naturals=('natural', 'sum'))
    t['per 1000'] = (1000 * t.naturals / t.zone_days).round(1)
    tt = L.groupby('loading', observed=False).agg(zone_days=('natural', 'size'), naturals=('natural', 'sum'))
    tt['per 1000'] = (1000 * tt.naturals / tt.zone_days).round(1)
    tt['median m72 (Pa)'] = L.groupby('loading', observed=False)['m72'].median().round(0)
    print(f"\nLoaded zone-days (72-h stress ≥ {a.min_m72:.0f} Pa): share of the 3-day stress that came "
          "in the heaviest single day")
    print(tt.to_string())
    print("\nBy 72-h stress bin:")
    print(t.to_string())

    # ── 3. Logistic regression, clustered by season ───────────────────────────
    M = d[d.m72 > 0].dropna(subset=['conc']).copy()
    M['logm72'] = np.log1p(M.m72)
    out = []
    for f in ('natural ~ logm72', 'natural ~ logm72 + conc', 'natural ~ logm72 + conc_snow',
              'natural ~ logm72 + np.log1p(peak24)'):
        Mf = M.dropna(subset=['conc_snow']) if 'conc_snow' in f else M
        r = smf.logit(f, Mf).fit(disp=0, cov_type='cluster', cov_kwds={'groups': Mf['season']})
        for p in r.params.index[1:]:
            ci = r.conf_int().loc[p]
            out.append(dict(model=f, term=p, coef=round(r.params[p], 2),
                            CI=f'{ci[0]:.2f} – {ci[1]:.2f}', p=round(r.pvalues[p], 3),
                            n=int(r.nobs), events=int(Mf.natural.sum()), AIC=round(r.aic, 1)))
    out = pd.DataFrame(out)
    print("\nLogistic regression on zone-days with any 72-h stress (SE clustered by season):")
    print(out.to_string(index=False))
    c = out[(out.model == 'natural ~ logm72 + conc') & (out.term == 'conc')].iloc[0]
    print(f"\n  conc coefficient {c.coef:+.2f}: odds ratio for all-in-one-day (conc 1.0) vs evenly spread "
          f"(1/3) at equal 72-h stress = {np.exp(c.coef * 2 / 3):.2f}")

    d.to_csv(FIGURES / f'{slug}_loading_rate_zone_days.csv', index=False)
    tt.to_csv(FIGURES / f'{slug}_loading_rate_table.csv')
    out.to_csv(FIGURES / f'{slug}_loading_rate_logit.csv', index=False)


if __name__ == '__main__':
    main()
