"""
event_stats.py — uncertainty for the zone-level and storm-level results.

Zone level (N = 4–5 start zones)
  Bootstrapping five points is meaningless, so correlations get EXACT
  permutation p-values (all 5! = 120 orderings).  Smallest attainable
  two-sided p: 0.017 (N=5), 0.083 (N=4).
  Hotspots here use the MicroMet season-average load:
      load(x) = h(θ) · max(0, mean seasonal new snow + wind accumulation(x))
      margin  = g · load · cos²θ · (tanθ − μ)⁺
  hotspot = top `pct` % of pooled start-zone cells, patches ≥ `min_cells`.

Storm level (zone-days, ~14 natural events)
  Primary test, named in advance: 72-h window, P99 stress margin.  Everything
  else is exploratory.
  - AUC with a cluster bootstrap 95 % CI: event days (clusters of paths that
    released together) and non-event seasons are resampled with replacement.
  - Date-shift permutation: every event day is moved by a random offset
    (±10–60 days, staying in Nov–May of the same season), paths kept —
    are release days more loaded than other days?
  - Within-day permutation: on each event day, the released path(s) are
    relabelled at random among the zones — does the model pick the path?
  - Out-of-sample split: the primary test re-scored separately on the
    2009–2021 development period and on 2021–2026 events (2026 CAIC
    export), which played no part in building the model.

Usage:  python src/event_stats.py --site sisters
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES                          # noqa: E402
from avl_data import load_avalanches                         # noqa: E402

G = 9.81
rng = np.random.default_rng(20261112)


def exact_perm_spearman(x, y):
    """Spearman ρ and exact two-sided permutation p over all orderings of y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = spearmanr(x, y)[0]
    perms = [spearmanr(x, np.array(p))[0] for p in itertools.permutations(y)]
    perms = np.nan_to_num(np.array(perms))
    return float(rho), float(np.mean(np.abs(perms) >= abs(rho) - 1e-12))


def auc_fast(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    r = rankdata(np.r_[pos, neg])
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def season_of(d):
    d = pd.Timestamp(d)
    return d.year + (1 if d.month >= 10 else 0)


def storm_tests(zd: pd.DataFrame, metric: str, n_boot=2000, n_perm=2000):
    zd = zd.copy()
    zd['day'] = pd.to_datetime(zd['day'])
    zd['season'] = zd['day'].map(season_of)
    val = zd.set_index(['day', 'zone'])[metric]
    pos = zd[zd.natural == 1]; neg = zd[zd.natural == 0]
    obs = auc_fast(pos[metric].values, neg[metric].values)

    # cluster bootstrap
    ev_days = pos['day'].unique(); seasons = neg['season'].unique()
    pos_by_day = {d: g[metric].values for d, g in pos.groupby('day')}
    neg_by_season = {s: g[metric].values for s, g in neg.groupby('season')}
    boots = []
    for _ in range(n_boot):
        p = np.concatenate([pos_by_day[d] for d in rng.choice(ev_days, len(ev_days))])
        n = np.concatenate([neg_by_season[s] for s in rng.choice(seasons, len(seasons))])
        boots.append(auc_fast(p, n))
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])

    # date-shift permutation
    days_avail = set(zd['day'])
    ev = pos[['day', 'zone']]
    null = []
    for _ in range(n_perm):
        shifted = []
        for d, g in ev.groupby('day'):
            for _try in range(50):
                off = int(rng.integers(10, 61)) * (1 if rng.random() < .5 else -1)
                nd = d + pd.Timedelta(days=off)
                if nd in days_avail and nd.month in (11, 12, 1, 2, 3, 4, 5) and season_of(nd) == season_of(d):
                    break
            else:
                continue
            shifted += [(nd, z) for z in g['zone']]
        idx = pd.MultiIndex.from_tuples(shifted, names=['day', 'zone'])
        p = val.reindex(idx).dropna().values
        mask = ~val.index.isin(idx)
        null.append(auc_fast(p, val[mask].values))
    p_shift = float(np.mean(np.array(null) >= obs))

    # within-day permutation
    def within(zz):
        w = []
        for d, g in zz.groupby('day'):
            if 0 < g.natural.sum() < len(g):
                w.append(auc_fast(g.loc[g.natural == 1, metric].values, g.loc[g.natural == 0, metric].values))
        return float(np.mean(w)) if w else np.nan
    evd = zd[zd['day'].isin(ev_days)].copy()
    obs_w = within(evd)
    nullw = []
    for _ in range(n_perm):
        sh = evd.copy()
        sh['natural'] = sh.groupby('day')['natural'].transform(lambda s: rng.permutation(s.values))
        nullw.append(within(sh))
    p_within = float(np.mean(np.array(nullw) >= obs_w - 1e-12))
    return dict(AUC=round(obs, 3), CI95=f"{lo:.2f}–{hi:.2f}", p_date_shift=round(p_shift, 3),
                AUC_within_day=round(obs_w, 3), p_within_day=round(p_within, 3),
                n_events=int(len(pos)), n_event_days=int(len(ev_days)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--mu-deg', type=float, default=30.0)
    ap.add_argument('--min-cells', type=int, default=2, help='patch size at 2 m (2 cells = 8 m²)')
    a = ap.parse_args()
    cfg = SITES[a.site]; slug = cfg['name'].lower().replace(' ', '_')

    # ── Zone level: MicroMet hotspots ─────────────────────────────────────────
    g = np.load(FIGURES / f'{slug}_loading_grids.npz')
    names = [str(n) for n in g['zone_names']]; zones = g['zones']
    slope = g['slope_deg'].astype(float); res = float(g['res'])
    sl = np.radians(slope); mu = np.tan(np.radians(a.mu_deg))
    shape = np.cos(sl) ** 2 * np.clip(np.tan(sl) - mu, 0, None) * np.clip((60 - slope) / 10, 0, 1)
    # same as micromet_loading.py: no slab stress within 4 m of forest canopy
    shape = shape * ~ndi.binary_dilation(g['forest'], iterations=max(1, int(round(4.0 / res))))
    load = np.maximum(float(g['mean_new_snow_mm']) + g['all'].astype(float), 0)
    marg = G * load * shape
    allz = zones.any(axis=0)

    nat = load_avalanches(str(cfg['avalanche_csv']), cfg['avalanche_filter'],
                          type=['HS', 'SS'], trig=['N']).groupby('HW Path').size()
    pmap = {k.replace(cfg['sz_dem_suffix'], ''): v for k, v in cfg['sz_avi_path_map'].items()}
    counts = pd.Series({n: nat.get(pmap.get(n, n), 0) for n in names})
    area = pd.Series({n: zones[i].sum() * res ** 2 for i, n in enumerate(names)})
    print(f"Mean seasonal new snow (SNOTEL): {float(g['mean_new_snow_mm']):.0f} mm SWE")

    tables = {}
    for pct in (1, 5, 10):
        thr = np.percentile(marg[allz], 100 - pct)
        hot = (marg >= thr) & allz
        lab, _ = ndi.label(hot, structure=np.ones((3, 3)))
        sizes = np.bincount(lab.ravel()); big = np.flatnonzero(sizes >= a.min_cells); big = big[big > 0]
        keep = np.isin(lab, big)
        rows = {}
        for i, n in enumerate(names):
            m = zones[i]; hz = keep & m
            r0, c0 = np.argwhere(m).mean(axis=0); hr, hc = np.nonzero(hz)
            rows[n] = dict(hot_area_m2=hz.sum() * res ** 2, patches=len(np.unique(lab[hz])),
                           wind_acc_mean=float(g['all'][m].mean()),
                           west_share=float((hc < c0).mean()) if len(hc) else np.nan,
                           upper_share=float((hr > r0).mean()) if len(hr) else np.nan)
        t = pd.DataFrame(rows).T
        t['nat'] = counts; t['nat_per_10k'] = counts / area * 1e4
        tables[pct] = t
    print("\nMicroMet hotspots (top 5 % of pooled start-zone cells):")
    print(tables[5].round(2).to_string())

    rows = []
    for pct, t in tables.items():
        for metric, outcome in (('hot_area_m2', 'nat'), ('patches', 'nat'), ('wind_acc_mean', 'nat'),
                                ('wind_acc_mean', 'nat_per_10k')):
            for lab_, sub in (('N=5', t), ('N=4 (no #6)', t.drop(index=[n for n in t.index if n.endswith('6_SZ')]))):
                rho, p = exact_perm_spearman(sub[metric], sub[outcome])
                rows.append(dict(pct=pct, metric=metric, outcome=outcome, set=lab_, rho=round(rho, 2), p_exact=round(p, 3)))
    zone_stats = pd.DataFrame(rows)
    print("\nZone-level Spearman ρ with exact permutation p (hotspot metrics vary with pct; wind_acc does not):")
    print(zone_stats.drop_duplicates(subset=['metric', 'outcome', 'set', 'rho', 'p_exact'])
          .to_string(index=False))

    # ── Storm level ───────────────────────────────────────────────────────────
    zd_all = pd.read_csv(FIGURES / f'{slug}_zone_day_stress.csv')
    tests = []
    for sc, metric, role in (('72h window', 'margin_p99', 'PRIMARY'),
                             ('72h window', 'wind_acc_p90', 'exploratory'),
                             ('72h window', 'swe_new_mm', 'exploratory (snowfall only)'),
                             ('72h window', 'margin_p99_loaded', 'exploratory (wind-loaded cells only)'),
                             ('since control (all), 5-d settling', 'margin_p99', 'exploratory')):
        r = storm_tests(zd_all[zd_all.scenario == sc], metric)
        tests.append(dict(role=role, scenario=sc, metric=metric, **r))
    # Out-of-sample: the model was developed on 2009 – 2021-05 data; events
    # after that come from the 2026 CAIC export and were never used to choose
    # anything, so scoring them separately is a genuine test.
    split = pd.Timestamp(cfg.get('control_records_end', '2021-05-04'))
    for period, sel in (('2009–2021 (development)', lambda d: d <= split),
                        ('2021–2026 (out of sample)', lambda d: d > split)):
        zd = zd_all[zd_all.scenario == '72h window'].copy()
        zd = zd[pd.to_datetime(zd['day']).map(sel)]
        if zd.natural.sum() >= 3:
            r = storm_tests(zd, 'margin_p99')
            tests.append(dict(role=f'PRIMARY, {period}', scenario='72h window', metric='margin_p99', **r))
    storm = pd.DataFrame(tests)
    pd.set_option('display.width', 220)
    print("\nStorm-level tests:")
    print(storm.to_string(index=False))

    tables[5].to_csv(FIGURES / f'{slug}_micromet_hotspots.csv')
    zone_stats.to_csv(FIGURES / f'{slug}_zone_stats.csv', index=False)
    storm.to_csv(FIGURES / f'{slug}_storm_stats.csv', index=False)
    print(f"\nSaved: {slug}_micromet_hotspots.csv, {slug}_zone_stats.csv, {slug}_storm_stats.csv")


if __name__ == '__main__':
    main()
