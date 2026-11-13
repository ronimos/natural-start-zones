"""
stress_variant_test.py — does the shear term carry the skill?

Natural slabs could start as slow shear failures (the stress margin
g·SWE·cos²θ·(tanθ − μ) is what drives them) or by weak-layer collapse
(normal stress g·SWE·cos²θ is what matters).  On a single cell both scale
with the same load, so the only information that separates them is how the
load is weighted across slopes: the shear margin favours steep loaded cells,
normal stress favours flatter ones.  This script compares micromet_loading.py
runs that differ only in the stress per unit load (--stress):

  margin    cos²θ·(tanθ − μ)              shear excess over friction (default)
  normal30  cos²θ, slopes > μ only         same cells, no shear weighting
  normal    cos²θ, all slopes              normal stress / collapse
  shear     sinθ·cosθ                      shear stress, no friction

Everything else (snowfall, wind deposition, forest buffer, 60° cap, P99
thresholds from each variant's own pooled distribution) is identical.

Tests, per variant:
  1. Day skill: AUC of the zone's P99 stress (primary metric) and of the
     hotspot (2nd-highest cell), 72-h and 24-h windows; season-block
     bootstrap CI; paired bootstrap of the difference from 'margin'.
  2. Path picking: within-day AUC on natural days, paired bootstrap over
     event days.
  3. Zone pattern: observed vs expected natural days (P99 hotspot days).
  4. Head to head: logistic natural ~ log(1+margin) + log(1+variant),
     SE clustered by season — which one keeps the signal?

Run first (≈ 7 min each, can run in parallel):
  python src/micromet_loading.py --site sisters --stress normal   --tag _st_normal
  python src/micromet_loading.py --site sisters --stress normal30 --tag _st_normal30
  python src/micromet_loading.py --site sisters --stress shear    --tag _st_shear
Then:  python src/stress_variant_test.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES          # noqa: E402

VARIANTS = {'margin': '', 'normal30': '_st_normal30', 'normal': '_st_normal', 'shear': '_st_shear'}
LABEL = {'margin': 'shear margin cos²θ(tanθ−μ)', 'normal30': 'normal stress cos²θ, > 30° only',
         'normal': 'normal stress cos²θ, all slopes', 'shear': 'shear stress sinθcosθ'}
rng = np.random.default_rng(20261113)


def auc(x, y):
    pos, neg = x[y == 1], x[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float((pos[:, None] > neg[None, :]).mean() + 0.5 * (pos[:, None] == neg[None, :]).mean())


def within_by_day(d, col):
    """{event day: within-day AUC of col}."""
    out = {}
    for day, g in d[d.day.isin(d.loc[d.natural == 1, 'day'])].groupby('day'):
        if g.natural.sum() and (g.natural == 0).sum():
            out[day] = auc(g[col].values, g.natural.values)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--n-boot', type=int, default=2000)
    a = ap.parse_args()
    slug = SITES[a.site]['name'].lower().replace(' ', '_')
    pd.set_option('display.width', 220)

    zd, oe = {}, {}
    for v, tag in VARIANTS.items():
        f = FIGURES / f'{slug}{tag}_zone_day_stress.csv'
        if not f.exists():
            sys.exit(f'missing {f.name} — run micromet_loading.py --stress {v} --tag {tag}')
        t = pd.read_csv(f, parse_dates=['day'])
        zd[v] = t
        oe[v] = pd.read_csv(FIGURES / f'{slug}{tag}_observed_vs_expected.csv')
    key = ['scenario', 'day', 'zone']
    base = zd['margin'][key + ['natural']]
    for v in VARIANTS:                                   # same rows in every run
        assert (zd[v][key].values == base[key].values).all(), v
    base = base.assign(season=base.day.dt.year + (base.day.dt.month >= 10))

    # ── 1–2. Day skill and path picking ───────────────────────────────────────
    rows = []
    for sc, col in (('72h window', 'margin_p99'), ('72h window', 'margin_2nd'), ('24h window', 'margin_2nd')):
        m = (base.scenario == sc).values
        d0 = base[m].reset_index(drop=True)
        cols = {v: zd[v].loc[m, col].values for v in VARIANTS}
        seasons = d0.season.unique(); idx_s = {s: np.flatnonzero(d0.season.values == s) for s in seasons}
        wd = {v: within_by_day(d0.assign(x=cols[v]), 'x') for v in VARIANTS}
        evdays = list(wd['margin'])
        bs = {v: [] for v in VARIANTS}; bw = {v: [] for v in VARIANTS}
        for _ in range(a.n_boot):
            ii = np.concatenate([idx_s[s] for s in rng.choice(seasons, len(seasons))])
            y = d0.natural.values[ii]
            for v in VARIANTS:
                bs[v].append(auc(cols[v][ii], y))
            dd = rng.choice(len(evdays), len(evdays))
            for v in VARIANTS:
                bw[v].append(np.nanmean([wd[v][evdays[k]] for k in dd]))
        for v in VARIANTS:
            diff = np.array(bs[v]) - np.array(bs['margin'])
            wdiff = np.array(bw[v]) - np.array(bw['margin'])
            lo, hi = np.nanpercentile(bs[v], [2.5, 97.5])
            wl, wh = np.nanpercentile(bw[v], [2.5, 97.5])
            rows.append(dict(window=sc[:3], metric=col, variant=v,
                             AUC=round(auc(cols[v], d0.natural.values), 3), AUC_CI=f'{lo:.2f}–{hi:.2f}',
                             dAUC_vs_margin=round(float(np.nanmean(diff)), 3) if v != 'margin' else np.nan,
                             p_worse=round(float((diff >= 0).mean()), 3) if v != 'margin' else np.nan,
                             within=round(float(np.nanmean(list(wd[v].values()))), 3), within_CI=f'{wl:.2f}–{wh:.2f}',
                             dwithin_vs_margin=round(float(np.nanmean(wdiff)), 3) if v != 'margin' else np.nan,
                             p_within_worse=round(float((wdiff >= 0).mean()), 3) if v != 'margin' else np.nan))
    res = pd.DataFrame(rows)
    print("Day skill (AUC, season bootstrap) and path picking (within-day AUC, event-day bootstrap).")
    print("dAUC = variant − shear margin; p_worse = bootstrap share with variant ≥ margin "
          "(small p ⇒ variant is reliably worse):\n")
    print(res.to_string(index=False))

    # ── 3. Zone pattern ──────────────────────────────────────────────────────
    zrows = []
    for v in VARIANTS:
        g = oe[v][(oe[v].threshold == 'P99') & (oe[v].basis == 'hotspot days')]
        o, e = g.observed.values.astype(float), g.expected.values
        dev = 2 * np.sum(np.where(o > 0, o * np.log(np.where(o > 0, o, 1) / e), 0) - (o - e))
        zrows.append(dict(variant=v, **{z.replace('_SZ', '').replace('Sister_', '#'): round(x, 1)
                                        for z, x in zip(g.zone, e)},
                          rho=round(float(spearmanr(e, o)[0]), 2), deviance=round(float(dev), 1)))
    zr = pd.DataFrame(zrows)
    obs = oe['margin'][(oe['margin'].threshold == 'P99') & (oe['margin'].basis == 'hotspot days')]
    print("\nExpected natural days per path (P99 hotspot days, scaled to the observed total); observed: "
          + ', '.join(f"{z.replace('_SZ', '').replace('Sister_', '#')} {o}" for z, o in zip(obs.zone, obs.observed)))
    print(zr.to_string(index=False))

    # ── 4. Head to head ──────────────────────────────────────────────────────
    hrows = []
    for sc, col in (('72h window', 'margin_p99'), ('72h window', 'margin_2nd'), ('24h window', 'margin_2nd')):
        m = (base.scenario == sc).values
        d0 = base[m].reset_index(drop=True)
        for v in ('normal30', 'normal', 'shear'):
            X = d0.assign(lm=np.log1p(zd['margin'].loc[m, col].values), lv=np.log1p(zd[v].loc[m, col].values))
            r = smf.logit('natural ~ lm + lv', X).fit(disp=0, cov_type='cluster', cov_kwds={'groups': X['season']})
            r1 = smf.logit('natural ~ lm', X).fit(disp=0)
            r2 = smf.logit('natural ~ lv', X).fit(disp=0)
            hrows.append(dict(window=sc[:3], metric=col, variant=v,
                              corr=round(float(np.corrcoef(X.lm, X.lv)[0, 1]), 3),
                              coef_margin=round(r.params['lm'], 2), p_margin=round(r.pvalues['lm'], 3),
                              coef_variant=round(r.params['lv'], 2), p_variant=round(r.pvalues['lv'], 3),
                              AIC_margin=round(r1.aic, 1), AIC_variant=round(r2.aic, 1)))
    hh = pd.DataFrame(hrows)
    print("\nHead to head: logit natural ~ log(1+margin) + log(1+variant), SE clustered by season; "
          "single-predictor AICs:")
    print(hh.to_string(index=False))

    res.to_csv(FIGURES / f'{slug}_stress_variant_skill.csv', index=False)
    zr.to_csv(FIGURES / f'{slug}_stress_variant_zones.csv', index=False)
    hh.to_csv(FIGURES / f'{slug}_stress_variant_head2head.csv', index=False)


if __name__ == '__main__':
    main()
