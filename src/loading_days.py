"""
loading_days.py — how often does each wind direction load each path's
hotspot terrain, and do those loading days track natural avalanche days?

Input: <slug>_zone_day_dir_loading.csv from micromet_loading.py — for each
zone, MST day and 10° wind-from bin, the wind deposition on the zone's
stress-bearing cells (open, tanθ > μ, < 60°) and the weak-layer stress it
adds, as the 2nd-highest cell (a patch of ≥ 2 cells, as for hotspots).

A path has a loading day from a sector when the stress added by that day's
wind from the sector reaches T Pa.  Per-direction values within a sector
are summed: an upper bound on the patch stress (patches from neighbouring
10° bins mostly coincide), and the per-bin maximum is a lower bound.

Evaluation days: Nov–May days with weather data (the zone-day table's 72-h
scenario), as in the storm test.  Outputs:
  1. loading days per season by path and sector vs natural days, with
     Spearman ρ and exact permutation p (N = 5 paths) for several T;
  2. naturals per 100 loading days, per path and sector;
  3. event level: share of natural days whose 72-h window held a loading
     day from each sector, against all zone-days (lift);
  4. figure <slug>_loading_days.png and tables <slug>_loading_days*.csv.

Usage:  python src/loading_days.py --site sisters [--t-pa 10]
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES          # noqa: E402

G = 9.81
SECTORS = {'SW': (200, 260), 'W–NW': (270, 330)}
BINS30 = list(range(0, 360, 30))             # 30° bins centred on these


def sector_of(d: int) -> str:
    for k, (lo, hi) in SECTORS.items():
        if lo <= d <= hi:
            return k
    return 'other'


def bin30(d: int) -> int:
    return int(((d + 15) % 360) // 30 * 30)


def swe_equiv(t_pa, slope_deg=45.0, mu_deg=30.0):
    th = np.radians(slope_deg)
    return t_pa / (G * np.cos(th) ** 2 * (np.tan(th) - np.tan(np.radians(mu_deg))))


def exact_perm_spearman(x, y):
    """Spearman ρ and exact two-sided permutation p over all orderings of y."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    r = spearmanr(x, y)[0]
    perms = [spearmanr(x, np.array(p))[0] for p in itertools.permutations(y)]
    return r, float(np.mean(np.abs(perms) >= abs(r) - 1e-12))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--t-pa', type=float, default=10.0,
                    help='added stress (Pa) that makes a loading day (default ≈ 5 mm SWE at 45°)')
    ap.add_argument('--combine', choices=['sum', 'max'], default='sum',
                    help='combine 10° bins within a sector by sum (upper bound) or max (lower bound)')
    a = ap.parse_args()
    cfg = SITES[a.site]; slug = cfg['name'].lower().replace(' ', '_')
    lab = lambda z: z.replace('_SZ', '').replace('_', ' #')
    pd.set_option('display.width', 250)

    dl = pd.read_csv(FIGURES / f'{slug}_zone_day_dir_loading.csv', parse_dates=['day'])
    zd = pd.read_csv(FIGURES / f'{slug}_zone_day_stress.csv', parse_dates=['day'])
    zd = zd[zd['scenario'] == '72h window'][['day', 'zone', 'natural']]
    zones = sorted(zd['zone'].unique())
    eval_days = pd.DatetimeIndex(sorted(zd['day'].unique()))
    seasons = eval_days.map(lambda d: d.year + (d.month >= 10)).nunique()

    dl['sector'] = dl['dir'].map(sector_of)
    dl['b30'] = dl['dir'].map(bin30)
    agg = 'sum' if a.combine == 'sum' else 'max'

    def day_table(key):
        """zone × day × key → added stress (Pa), on evaluation days, zeros filled."""
        t = dl.groupby(['zone', 'day', key])['stress_2nd'].agg(agg).unstack(key, fill_value=0.0)
        idx = pd.MultiIndex.from_product([zones, eval_days], names=['zone', 'day'])
        return t.reindex(idx, fill_value=0.0)

    sec = day_table('sector')
    for k in list(SECTORS) + ['other']:
        if k not in sec:
            sec[k] = 0.0
    sec['any'] = dl.groupby(['zone', 'day'])['stress_2nd'].agg(agg).reindex(sec.index, fill_value=0.0)
    b30 = day_table('b30')
    nat = zd.set_index(['zone', 'day'])['natural'].reindex(sec.index, fill_value=0).astype(int)
    nat_days = nat.groupby(level='zone').sum()

    print(f"{len(eval_days)} evaluation days, {seasons} seasons; loading day = added stress "
          f"≥ T on a ≥ 2-cell patch of stress-bearing cells ({a.combine} over 10° bins in a sector)")
    print("Natural days: " + ', '.join(f"{lab(z)} {nat_days[z]}" for z in zones))

    # ── 1. Loading days per season vs natural days, several thresholds ───────
    cols = list(SECTORS) + ['other', 'any']
    rows, corr = [], []
    for T in sorted({5.0, 10.0, 20.0, 30.0, a.t_pa}):
        cnt = (sec[cols] >= T).groupby(level='zone').sum() / seasons
        for k in cols:
            r, p = exact_perm_spearman(cnt[k].values, nat_days[zones].values)
            corr.append(dict(T_Pa=T, SWE45_mm=round(float(swe_equiv(T)), 1), sector=k,
                             rho=round(r, 2) if r == r else np.nan, p_exact=p))
        if T == a.t_pa:
            rows = cnt.round(1)
    corr = pd.DataFrame(corr)
    tab = rows.copy(); tab.index = [lab(z) for z in tab.index]
    tab['natural days'] = nat_days[zones].values
    print(f"\nLoading days per season at T = {a.t_pa:g} Pa (≈ {swe_equiv(a.t_pa):.1f} mm SWE on 45°):")
    print(tab.to_string())
    print("\nSpearman ρ, loading days per season vs natural days (N = 5 paths; exact p ≥ 0.017):")
    print(corr.pivot_table(index=['T_Pa', 'SWE45_mm'], columns='sector', values=['rho', 'p_exact'])
          .round(3).to_string())

    # ── 2. Naturals per 100 loading days (72-h window: loading day d−2…d) ───
    win = {k: (sec[k] >= a.t_pa).groupby(level='zone')
           .transform(lambda s: s.rolling(3, min_periods=1).max()).astype(bool) for k in cols}
    rate = []
    for z in zones:
        r = {'path': lab(z)}
        for k in cols:
            w = win[k].xs(z, level='zone'); n_ = nat.xs(z, level='zone')
            r[f'{k}: loaded days'] = int(w.sum())
            r[f'{k}: naturals'] = int((w & (n_ == 1)).sum())
            r[f'{k}: per 100'] = round(100 * (w & (n_ == 1)).sum() / w.sum(), 2) if w.sum() else np.nan
        rate.append(r)
    rate = pd.DataFrame(rate).set_index('path')
    print(f"\nNaturals per 100 days with a {a.t_pa:g}-Pa loading day in the 72-h window, by sector:")
    print(rate[[c for c in rate if c.endswith('per 100')]].to_string())

    # ── 3. Event level: was the path loaded from each sector before naturals? ─
    ev = []
    for k in cols:
        w = win[k]
        on_nat = w[nat == 1].mean(); base = w.mean()
        ev.append(dict(sector=k, natural_days_loaded=f"{int(w[nat == 1].sum())}/{int((nat == 1).sum())}",
                       share_natural=round(on_nat, 2), share_all=round(base, 2),
                       lift=round(on_nat / base, 2) if base else np.nan))
        for z in zones:
            wz = w.xs(z, level='zone'); nz = nat.xs(z, level='zone')
            if nz.sum():
                ev.append(dict(sector=k, path=lab(z),
                               natural_days_loaded=f"{int(wz[nz == 1].sum())}/{int(nz.sum())}",
                               share_natural=round(wz[nz == 1].mean(), 2), share_all=round(wz.mean(), 2),
                               lift=round(wz[nz == 1].mean() / wz.mean(), 2) if wz.mean() else np.nan))
    ev = pd.DataFrame(ev)
    print(f"\nShare of natural days with a loading day (≥ {a.t_pa:g} Pa) in the 72-h window, vs all zone-days:")
    print(ev[ev['path'].isna()].drop(columns='path').to_string(index=False))
    print("\nBy path:")
    print(ev[ev['path'].notna()].pivot_table(index='path', columns='sector',
                                             values=['share_natural', 'lift']).round(2).to_string())

    # ── 4. 30° profile and figure ─────────────────────────────────────────────
    prof = (b30 >= a.t_pa).groupby(level='zone').sum() / seasons
    prof = prof.reindex(columns=[b for b in BINS30 if b in prof.columns])
    pc = []
    for b in prof.columns:
        r, p = exact_perm_spearman(prof[b].values, nat_days[zones].values)
        pc.append(dict(dir=b, rho=r, p_exact=p))
    pc = pd.DataFrame(pc).set_index('dir')
    print(f"\nLoading days per season by 30° wind-from bin (≥ {a.t_pa:g} Pa), ρ with natural days:")
    out = prof.copy(); out.index = [lab(z) for z in out.index]
    print(pd.concat([out.round(1), pc.T.round(2).rename(index={'rho': 'ρ', 'p_exact': 'p'})]).to_string())

    tab.to_csv(FIGURES / f'{slug}_loading_days.csv')
    corr.to_csv(FIGURES / f'{slug}_loading_days_corr.csv', index=False)
    rate.to_csv(FIGURES / f'{slug}_loading_days_rate.csv')
    ev.to_csv(FIGURES / f'{slug}_loading_days_events.csv', index=False)
    out.to_csv(FIGURES / f'{slug}_loading_days_30deg.csv')

    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True,
                            gridspec_kw=dict(width_ratios=[1, 1.3]))
    ax = axs[0]
    x = np.arange(len(zones)); w_ = 0.27
    colors = {'SW': 'tab:orange', 'W–NW': 'tab:blue', 'other': 'tab:gray'}
    for i, k in enumerate(['SW', 'W–NW', 'other']):
        ax.bar(x + (i - 1) * w_, rows[k].values, w_, color=colors[k], label=f'{k} loading days')
    ax.set_xticks(x, [lab(z) for z in zones]); ax.set_ylabel('loading days per season')
    ax2 = ax.twinx()
    ax2.plot(x, nat_days[zones].values, 'kD', ms=9, label='natural days (right axis)')
    ax2.set_ylabel('natural slab days'); ax2.set_ylim(0, nat_days.max() * 1.25)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')
    rr = corr[corr.T_Pa == a.t_pa].set_index('sector')
    ax.set_title(f'Days per season on which the wind adds ≥ {a.t_pa:g} Pa to a hotspot patch\n'
                 + ', '.join(f"{k}: ρ = {rr.loc[k, 'rho']:+.2f} (p = {rr.loc[k, 'p_exact']:.2f})"
                             for k in ['SW', 'W–NW', 'any']), fontsize=10)
    ax = axs[1]
    for z in zones:
        ax.plot(prof.columns, prof.loc[z].values, 'o-', label=f'{lab(z)} ({nat_days[z]} naturals)')
    ax.set_xticks(prof.columns, [f'{b}°' for b in prof.columns])
    ax.set_xlabel('wind from (30° bins)'); ax.set_ylabel('loading days per season')
    ax.set_title('Loading days by wind direction', fontsize=10); ax.legend(fontsize=8)
    fig.suptitle(f"{cfg['name']} — directional loading of stress-bearing cells (MicroMet / saltation)")
    fig.savefig(FIGURES / f'{slug}_loading_days.png', dpi=140); plt.close(fig)
    print(f"\nSaved {slug}_loading_days.png and tables")


if __name__ == '__main__':
    main()
