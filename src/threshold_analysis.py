"""
threshold_analysis.py — what stress threshold do the natural avalanches imply,
and how often did each path's hotspots reach it?

Input: the zone-day table written by micromet_loading.py (72-h window), with
margin_2nd = the 2nd-highest cell stress margin in the zone that day.  A zone
has a *hotspot day* at threshold t when margin_2nd > t (a patch of >= 2 cells,
~critical-crack-length scale, exceeds t).

The stress here is the margin added by the last 72 h of loading (snowfall +
wind deposition), τ − μσ in Pa — an increment, not the weak layer's total
strength.  It is also expressed as the equivalent new SWE on a 45° (and 40°)
slope:  SWE = t / (g · cos²θ · (tanθ − μ)).

1. Reference threshold: the model's P99 threshold — per-path hotspot days,
   natural days, hits / misses / false alarms.
2. Threshold scan (5–400 Pa): hit rate H = P(hotspot | natural),
   false-alarm rate F = P(hotspot | no natural), Peirce skill PSS = H − F,
   and the zone-level fit of hotspot days to natural counts.  The threshold
   that maximises PSS is the data's estimate of τ_thresh (as a 72-h
   increment), with a season-block bootstrap 95 % CI.
3. Release rate per 1,000 zone-days in bins of hotspot stress, and a
   logistic trend — the paper predicts onset at a threshold and a rate
   that grows with the excess stress.
4. Avalanche size (D, R) vs hotspot stress, stressed area, new snow and
   wind deposition on the release day, with and without the one D2.5
   outlier (2021-01-24 Sister #6, coordinate-assigned, little new snow).
5. Control missions as censoring (to 2021-05-04): the threshold and stress
   trend recomputed without zone-days on, or the day after, a control
   mission (explosives remove slabs before they can release naturally).

Usage:  python src/threshold_analysis.py --site sisters
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline import SITES, FIGURES          # noqa: E402

G = 9.81
rng = np.random.default_rng(20261112)


def swe_equiv(t_pa, slope_deg, mu_deg=30.0):
    th = np.radians(slope_deg)
    return t_pa / (G * np.cos(th) ** 2 * (np.tan(th) - np.tan(np.radians(mu_deg))))


def season_of(d):
    return d.year + (1 if d.month >= 10 else 0)


def contingency(zd, t):
    hot = zd['margin_2nd'].values > t
    nat = zd['natural'].values == 1
    h = (hot & nat).sum(); m = (~hot & nat).sum(); f = (hot & ~nat).sum(); c = (~hot & ~nat).sum()
    H = h / (h + m) if h + m else np.nan
    F = f / (f + c) if f + c else np.nan
    return h, m, f, c, H, F


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--site', default='sisters', choices=list(SITES))
    ap.add_argument('--n-boot', type=int, default=2000)
    a = ap.parse_args()
    cfg = SITES[a.site]; slug = cfg['name'].lower().replace(' ', '_')

    zd = pd.read_csv(FIGURES / f'{slug}_zone_day_stress.csv')
    zd = zd[zd['scenario'] == '72h window'].copy()
    if 'margin_2nd' not in zd:
        sys.exit('zone-day table has no margin_2nd column — rerun micromet_loading.py')
    zd['day'] = pd.to_datetime(zd['day']); zd['season'] = zd['day'].map(season_of)
    n_seasons = zd['season'].nunique()
    test = pd.read_csv(FIGURES / f'{slug}_stress_event_test.csv')
    t99 = float(test.loc[(test.scenario == '72h window') & (test.metric == 'margin_p99'), 'thr_q99_Pa'].iloc[0])
    zones = sorted(zd['zone'].unique())
    lab = lambda z: z.replace('_SZ', '').replace('_', ' #')

    # ── 1. Reference (P99) threshold, per path ───────────────────────────────
    print(f"Model threshold (pooled P99 of start-zone stress margin, 72-h load): {t99:.1f} Pa")
    print(f"  ≈ {swe_equiv(t99, 45):.0f} mm SWE of new load in 72 h on a 45° slope "
          f"({swe_equiv(t99, 40):.0f} mm at 40°)\n")

    def per_path(t):
        rows = []
        for z in zones:
            g = zd[zd.zone == z]
            h, m, f, c, H, F = contingency(g, t)
            rows.append({'path': lab(z), 'hotspot days': int(h + f),
                         'per season': round((h + f) / n_seasons, 1),
                         'natural days': int(h + m), 'naturals on hotspot days': int(h),
                         'naturals without hotspot': int(m), 'hotspot days without natural': int(f),
                         'P(natural | hotspot day)': round(h / (h + f), 3) if h + f else np.nan,
                         'P(hotspot | natural day)': round(H, 2) if h + m else np.nan})
        return pd.DataFrame(rows)

    pd.set_option('display.width', 250)
    ref = per_path(t99)
    print(f"Per path at the P99 threshold ({t99:.0f} Pa), {n_seasons} seasons:")
    print(ref.to_string(index=False))

    # ── 2. Threshold scan ─────────────────────────────────────────────────────
    grid = np.unique(np.round(np.geomspace(5, 400, 60), 1))
    scan = []
    for t in grid:
        h, m, f, c, H, F = contingency(zd, t)
        hot_z = np.array([(zd[zd.zone == z]['margin_2nd'] > t).sum() for z in zones], float)
        nat_z = np.array([zd[zd.zone == z]['natural'].sum() for z in zones], float)
        rho = pd.Series(hot_z).corr(pd.Series(nat_z), method='spearman') if hot_z.std() > 0 else np.nan
        scan.append(dict(t_Pa=t, SWE45_mm=round(float(swe_equiv(t, 45)), 1), hit_rate=H, false_alarm=F,
                         PSS=H - F, hotspot_zone_days=int(h + f), zone_rho=rho))
    scan = pd.DataFrame(scan)
    best = scan.loc[scan['PSS'].idxmax()]

    # season-block bootstrap of the PSS-optimal threshold
    by_season = {s: g for s, g in zd.groupby('season')}
    seasons = np.array(list(by_season))
    opt = []
    for _ in range(a.n_boot):
        b = pd.concat([by_season[s] for s in rng.choice(seasons, len(seasons))])
        if b['natural'].sum() < 3:
            continue
        pss = []
        for t in grid:
            _, _, _, _, H, F = contingency(b, t)
            pss.append(H - F)
        opt.append(grid[int(np.nanargmax(pss))])
    lo, hi = np.percentile(opt, [2.5, 97.5])
    print(f"\nThreshold that best separates natural days (max Peirce skill H − F): "
          f"{best.t_Pa:.0f} Pa (95 % CI {lo:.0f}–{hi:.0f} Pa)")
    print(f"  ≈ {best.SWE45_mm:.0f} mm SWE of new load in 72 h on a 45° slope "
          f"(CI {swe_equiv(lo, 45):.0f}–{swe_equiv(hi, 45):.0f} mm); "
          f"H = {best.hit_rate:.2f}, F = {best.false_alarm:.2f}, PSS = {best.PSS:.2f}")
    opt_tab = per_path(best.t_Pa)
    print(f"\nPer path at the data-derived threshold ({best.t_Pa:.0f} Pa):")
    print(opt_tab.to_string(index=False))
    print("\nScan (every 5th threshold):")
    print(scan.iloc[::5].round(3).to_string(index=False))

    # ── 3. Release rate vs hotspot stress (paper: rate ∝ τ − τ_thresh) ───────
    from scipy.stats import chi2
    import statsmodels.api as sm
    bins = [0, 15, 30, 45, 60, 80, 110, np.inf]
    zd['stress_bin'] = pd.cut(zd['margin_2nd'], bins, right=False)
    rate = zd.groupby('stress_bin', observed=True).agg(zone_days=('natural', 'size'), naturals=('natural', 'sum'))
    rate['per_1000'] = rate.naturals / rate.zone_days * 1000
    rate['lo95'] = [chi2.ppf(.025, 2 * k) / 2 / n * 1000 if k else 0 for k, n in zip(rate.naturals, rate.zone_days)]
    rate['hi95'] = [chi2.ppf(.975, 2 * (k + 1)) / 2 / n * 1000 for k, n in zip(rate.naturals, rate.zone_days)]
    print("\nNatural releases per 1,000 zone-days by hotspot stress (margin of the 2nd-highest cell):")
    print(rate.round(2).to_string())
    X = sm.add_constant(np.clip(zd['margin_2nd'].values, 0, None) / 10)
    fit = sm.Logit(zd['natural'].values, X).fit(disp=0)
    orr, (olo, ohi) = np.exp(fit.params[1]), np.exp(fit.conf_int()[1])
    print(f"Logistic trend: odds of a natural × {orr:.2f} per +10 Pa (95 % CI {olo:.2f}–{ohi:.2f}, "
          f"p = {fit.pvalues[1]:.1e}; zone-days treated as independent, so the CI is optimistic)")
    rate.to_csv(FIGURES / f'{slug}_release_rate_by_stress.csv')

    # ── 4. Avalanche size vs hotspot stress and stressed area ────────────────
    # Nucleation should follow local excess stress; propagation (size) should
    # depend more on how much slab is loaded above threshold.
    from scipy.stats import spearmanr
    av = pd.read_csv(cfg['avalanche_csv'], dtype=str)
    for c in ('Trigger', 'Type', 'Elev', 'HW Path', 'Dsize', 'Rsize'):
        av[c] = av[c].fillna('').str.strip()
    av = av[av['HW Path'].str.contains(cfg['avalanche_filter']) & (av.Trigger == 'N') & av.Type.isin(['HS', 'SS'])]
    av = av[(av['HW Path'] != 'Seven Sister #6') | av.Elev.isin(['>TL', 'TL'])]
    av['zone'] = av['HW Path'].map({v: k for k, v in
                                    {k.replace(cfg['sz_dem_suffix'], ''): v
                                     for k, v in cfg['sz_avi_path_map'].items()}.items()})
    av['day'] = pd.to_datetime(av['Date'], format='mixed').dt.normalize()
    av['D'] = pd.to_numeric(av['Dsize'].str.replace('D', ''), errors='coerce')
    av['R'] = pd.to_numeric(av['Rsize'].str.replace('R', ''), errors='coerce')
    ev = zd[zd.natural == 1].merge(av[['day', 'zone', 'D', 'R']], on=['day', 'zone'], how='left')
    q95 = float(test.loc[(test.scenario == '72h window') & (test.metric == 'margin_p99'), 'thr_q95_Pa'].iloc[0])

    def perm_rho(x, y, n=20000):
        m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
        r = spearmanr(x, y)[0]
        null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
        return r, float(np.mean(np.abs(null) >= abs(r))), int(m.sum())

    size_rows = []
    preds = (('margin_2nd', 'hotspot stress (Pa)'), ('area_gt_q95', f'area above {q95:.0f} Pa (m²)'),
             ('swe_new_mm', '72-h new snow (mm SWE)'), ('wind_acc_p90', '72-h wind deposition P90 (mm)'))
    for label, sub in (('all naturals', ev), ('without 2021-01-24 (D2.5 outlier)',
                                              ev[ev['day'] != pd.Timestamp('2021-01-24')])):
        for col, name in preds:
            for size in ('D', 'R'):
                r, p, n = perm_rho(sub[col].values.astype(float), sub[size].values.astype(float))
                size_rows.append(dict(set=label, predictor=name, size=size, n=n, rho=round(r, 2), p_perm=round(p, 3)))
    size_tab = pd.DataFrame(size_rows)
    print("\nAvalanche size vs loading on the release day (Spearman, permutation p):")
    print(size_tab.pivot_table(index=['set', 'predictor'], columns='size', values=['rho', 'p_perm'],
                               sort=False).to_string())
    med = ev.groupby('D')[['margin_2nd', 'area_gt_q95', 'swe_new_mm']].median().round(1)
    med['n'] = ev.groupby('D').size()
    print("\nMedian loading by destructive size:"); print(med.to_string())
    size_tab.to_csv(FIGURES / f'{slug}_size_vs_stress.csv', index=False)

    # ── 5. Control missions as censoring (control records end 2021-05-04) ────
    # Explosive control removes slabs on heavily loaded days, so those days
    # cannot produce naturals.  Compare the threshold and stress trend with
    # and without zone-days on, or the day after, a control mission.
    raw = pd.read_csv(cfg['avalanche_csv'], usecols=['Date', 'HW Path', 'Trigger'], dtype=str)
    raw = raw[raw['HW Path'].fillna('').str.contains(cfg['avalanche_filter'])
              & raw['Trigger'].fillna('').str.strip().isin(['AA', 'AB', 'AE', 'AH', 'AL', 'AP', 'AX'])]
    ctrl_days = set(pd.to_datetime(raw['Date'], format='mixed').dt.normalize())
    ctrl = ctrl_days | {d + pd.Timedelta(days=1) for d in ctrl_days}
    end = pd.Timestamp(cfg.get('control_records_end', '2021-05-04'))
    early = zd[zd['day'] <= end]
    cens_rows = []
    for label, sub in (('2009–2021, all zone-days', early),
                       ('2009–2021, control day or day after removed', early[~early['day'].isin(ctrl)])):
        pss = [(t,) + contingency(sub, t)[4:] for t in grid]
        bt, bH, bF = max(pss, key=lambda r: (r[1] - r[2]) if np.isfinite(r[1] - r[2]) else -9)
        Xs = sm.add_constant(np.clip(sub['margin_2nd'].values, 0, None) / 10)
        f_ = sm.Logit(sub['natural'].values, Xs).fit(disp=0)
        h30 = contingency(sub, 30.0)
        cens_rows.append(dict(set=label, zone_days=len(sub), naturals=int(sub['natural'].sum()),
                              best_t_Pa=bt, H=round(bH, 2), F=round(bF, 2), PSS=round(bH - bF, 2),
                              F_at_30Pa=round(h30[5], 2),
                              odds_per_10Pa=round(float(np.exp(f_.params[1])), 2),
                              odds_CI=f"{np.exp(f_.conf_int()[1][0]):.2f}–{np.exp(f_.conf_int()[1][1]):.2f}"))
    cens = pd.DataFrame(cens_rows)
    nat_days = early.loc[early['natural'] == 1, 'day']
    print(f"\nControl missions as censoring ({len(ctrl_days)} control days to {end.date()}; "
          f"{int(nat_days.isin(ctrl).sum())} of {len(nat_days)} natural zone-days fall on a control day "
          f"or the day after):")
    print(cens.to_string(index=False))
    cens.to_csv(FIGURES / f'{slug}_threshold_control_censoring.csv', index=False)

    ref.to_csv(FIGURES / f'{slug}_threshold_per_path_P99.csv', index=False)
    opt_tab.to_csv(FIGURES / f'{slug}_threshold_per_path_optimal.csv', index=False)
    scan.to_csv(FIGURES / f'{slug}_threshold_scan.csv', index=False)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)
    ax = axs[0]
    ax.plot(scan.t_Pa, scan.hit_rate, label='hit rate  P(hotspot | natural day)')
    ax.plot(scan.t_Pa, scan.false_alarm, label='false-alarm rate  P(hotspot | other days)')
    ax.plot(scan.t_Pa, scan.PSS, 'k', lw=2, label='Peirce skill  H − F')
    ax.axvspan(lo, hi, color='0.85', zorder=0, label='95 % CI of best threshold')
    ax.axvline(best.t_Pa, color='k', ls='--')
    ax.axvline(t99, color='tab:red', ls=':', label=f'model P99 threshold ({t99:.0f} Pa)')
    ax.set_xscale('log'); ax.set_xlabel('stress-margin threshold from 72-h loading (Pa)')
    ax.set_ylabel('rate / skill'); ax.legend(fontsize=8, loc='upper right')
    sec = ax.secondary_xaxis('top', functions=(lambda t: swe_equiv(t, 45), lambda s: s * G * 0.5 * (1 - np.tan(np.radians(30)))))
    sec.set_xlabel('≈ new SWE in 72 h on a 45° slope (mm)')
    ax.set_title(f'Which threshold separates natural days?  best ≈ {best.t_Pa:.0f} Pa (CI {lo:.0f}–{hi:.0f})', fontsize=10)

    ax = axs[1]
    x = np.arange(len(zones))
    ax.bar(x - .2, opt_tab['per season'], .4, color='tab:orange', label=f'hotspot days / season (> {best.t_Pa:.0f} Pa)')
    ax2 = ax.twinx()
    ax2.bar(x + .2, opt_tab['natural days'], .4, color='0.3', label='natural slab days (2009–2026, with weather)')
    for i, r in opt_tab.iterrows():
        ax2.text(i + .2, r['natural days'] + .2, f"{r['naturals on hotspot days']}/{r['natural days']}",
                 ha='center', fontsize=8)
    ax.set_xticks(x, opt_tab['path']); ax.set_ylabel('hotspot days per season'); ax2.set_ylabel('natural slab days')
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')
    ax.set_title('Hotspot days vs natural avalanches per path\n(labels: naturals on a hotspot day / all naturals)', fontsize=10)
    fig.savefig(FIGURES / f'{slug}_threshold_analysis.png', dpi=140); plt.close(fig)
    print(f"\nSaved: {slug}_threshold_per_path_*.csv, {slug}_threshold_scan.csv, {slug}_threshold_analysis.png")


if __name__ == '__main__':
    main()
