# =============================================================================
# Improved Statistical Analysis for Start Zone Terrain Variability
# Add these cells to the notebook after the existing terrain computation
# =============================================================================

# %% Cell: Imports for improved analysis
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, kendalltau, pearsonr
from collections import OrderedDict

# %% Cell: Fix D-size ordinal mapping
# BUG FIX: String comparison on D-size is unreliable
# Replace: avalanche_data[avalanche_data['Dsize'] > 'D1']
# With explicit ordinal mapping:

DSIZE_ORDER = {'D1': 1, 'D1.5': 1.5, 'D2': 2, 'D2.5': 2.5,
               'D3': 3, 'D3.5': 3.5, 'D4': 4, 'D4.5': 4.5, 'D5': 5}

avalanche_data['Dsize_numeric'] = avalanche_data['Dsize'].map(DSIZE_ORDER)
larger_avis = avalanche_data[avalanche_data['Dsize_numeric'] > 1.0]


# %% Cell: Circular statistics for aspect
def circular_variance(aspects_deg):
    """
    Compute circular variance (0 = all same direction, 1 = uniform spread).
    Standard kurtosis/std are invalid for circular data like aspect.
    """
    aspects_rad = np.deg2rad(aspects_deg)
    # Mean resultant length (R): 1 = perfectly concentrated, 0 = uniform
    C = np.nanmean(np.cos(aspects_rad))
    S = np.nanmean(np.sin(aspects_rad))
    R = np.sqrt(C**2 + S**2)
    # Circular variance = 1 - R
    return 1 - R


def circular_std(aspects_deg):
    """Circular standard deviation in degrees."""
    aspects_rad = np.deg2rad(aspects_deg)
    C = np.nanmean(np.cos(aspects_rad))
    S = np.nanmean(np.sin(aspects_rad))
    R = np.sqrt(C**2 + S**2)
    # Circular std (in radians), convert to degrees
    if R > 0:
        return np.rad2deg(np.sqrt(-2 * np.log(R)))
    return 180.0  # maximum dispersion


# %% Cell: Compute multiple variability metrics per start zone
def compute_sz_metrics(slope_array, aspect_array):
    """
    Compute terrain variability metrics for a single start zone.
    Returns dict of metrics.
    """
    # Flatten and remove NaN/nodata
    slopes = slope_array.flatten()
    slopes = slopes[~np.isnan(slopes)]
    slopes = slopes[(slopes > 0) & (slopes < 90)]  # valid slope angles

    aspects = aspect_array.flatten()
    aspects = aspects[~np.isnan(aspects)]
    aspects = aspects[(aspects >= 0) & (aspects <= 360)]

    metrics = {}

    # --- Slope variability ---
    metrics['slope_mean'] = np.mean(slopes)
    metrics['slope_std'] = np.std(slopes)
    metrics['slope_cv'] = np.std(slopes) / np.mean(slopes) if np.mean(slopes) > 0 else 0
    metrics['slope_iqr'] = np.percentile(slopes, 75) - np.percentile(slopes, 25)
    metrics['slope_kurtosis'] = stats.kurtosis(slopes)  # keep for comparison
    metrics['slope_range'] = np.ptp(slopes)
    metrics['slope_pct_gt_35'] = np.mean(slopes > 35) * 100  # % of cells > 35 degrees
    metrics['slope_pct_gt_40'] = np.mean(slopes > 40) * 100

    # --- Aspect variability (circular statistics) ---
    metrics['aspect_circ_var'] = circular_variance(aspects)
    metrics['aspect_circ_std'] = circular_std(aspects)
    metrics['aspect_kurtosis'] = stats.kurtosis(aspects)  # keep for comparison

    # Number of distinct aspect octants represented (N/NE/E/SE/S/SW/W/NW)
    octants = ((aspects + 22.5) % 360 // 45).astype(int)
    metrics['aspect_n_octants'] = len(np.unique(octants))

    # --- Start zone size ---
    metrics['n_cells'] = len(slopes)
    # At 1.2m resolution: area in m²
    metrics['area_m2'] = len(slopes) * 1.2 * 1.2

    return metrics


# %% Cell: Build metrics table for all start zones
# Assumes you have individual SZ DEMs clipped and can compute slope/aspect per SZ
# Adapt paths for your expanded dataset (Star Mtn, Muleshoe, etc.)

from glob import glob

sz_metrics = OrderedDict()

for path in sorted(glob(os.path.join(SISTERS_DATA_PATH, 'Sister ? SZ_DEM.tif'))):
    sister_id = path[path.find(' '): path.rfind(' ')].strip()

    slope_data = get_slope_attributes(path, 'slope_degrees')
    aspect_data = get_slope_attributes(path, 'aspect')

    slope_key = [k for k in slope_data.keys() if 'slope' in k][0]
    aspect_key = [k for k in aspect_data.keys() if 'aspect' in k][0]

    sz_metrics[sister_id] = compute_sz_metrics(
        slope_data[slope_key],
        aspect_data[aspect_key]
    )

metrics_df = pd.DataFrame(sz_metrics).T
metrics_df.index.name = 'start_zone'
print(metrics_df.round(3).to_string())


# %% Cell: Merge terrain metrics with avalanche counts
# Normalize counts by area

analysis_df = metrics_df.copy()

# Join avalanche counts — adapt index matching to your naming convention
# This may need adjustment depending on how your SZ names map to HW Path names
analysis_df['avi_count'] = total_avi_counts['All HS/SS avalanches'].values
analysis_df['avi_gt_d1'] = total_avi_counts['>D1 HS/SS avalanches'].fillna(0).values

# Normalize by area (avalanches per 10,000 m²)
analysis_df['avi_per_area'] = analysis_df['avi_count'] / (analysis_df['area_m2'] / 10000)
analysis_df['avi_gt_d1_per_area'] = analysis_df['avi_gt_d1'] / (analysis_df['area_m2'] / 10000)

print(analysis_df[['area_m2', 'avi_count', 'avi_per_area',
                    'avi_gt_d1', 'avi_gt_d1_per_area']].round(3).to_string())


# %% Cell: Correlation analysis — nonparametric for small N
def correlation_report(x, y, x_name, y_name):
    """
    Run Spearman and Kendall correlations (appropriate for small N, 
    ordinal-like data, no normality assumption).
    Pearson included for comparison but interpret cautiously.
    """
    n = len(x)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]

    print(f"\n{'='*60}")
    print(f"{x_name}  vs  {y_name}  (N={len(x)})")
    print(f"{'='*60}")

    # Spearman (rank-based, no distributional assumption)
    rho, p_spearman = spearmanr(x, y)
    print(f"  Spearman rho = {rho:.3f}  (p = {p_spearman:.3f})")

    # Kendall (more robust for small N, handles ties)
    tau, p_kendall = kendalltau(x, y)
    print(f"  Kendall  tau = {tau:.3f}  (p = {p_kendall:.3f})")

    # Pearson (for comparison — assumes linearity + normality)
    if len(x) >= 3:
        r, p_pearson = pearsonr(x, y)
        print(f"  Pearson  r   = {r:.3f}  (p = {p_pearson:.3f})")

    if len(x) < 10:
        print(f"  ⚠ N={len(x)}: p-values unreliable. Report effect sizes, not significance.")

    return {'spearman_rho': rho, 'p_spearman': p_spearman,
            'kendall_tau': tau, 'p_kendall': p_kendall}


# Test key hypotheses
terrain_vars = [
    ('slope_std', 'Slope angle SD'),
    ('slope_cv', 'Slope angle CV'),
    ('slope_iqr', 'Slope angle IQR'),
    ('slope_kurtosis', 'Slope angle kurtosis'),
    ('slope_pct_gt_35', '% cells > 35°'),
    ('aspect_circ_var', 'Aspect circular variance'),
    ('aspect_circ_std', 'Aspect circular SD'),
    ('aspect_n_octants', 'N aspect octants'),
]

outcome_vars = [
    ('avi_per_area', 'Avalanches / 10k m²'),
    ('avi_gt_d1_per_area', '>D1 avalanches / 10k m²'),
]

results = []
for t_col, t_name in terrain_vars:
    for o_col, o_name in outcome_vars:
        r = correlation_report(
            analysis_df[t_col].values,
            analysis_df[o_col].values,
            t_name, o_name
        )
        r['terrain_metric'] = t_name
        r['outcome'] = o_name
        results.append(r)

corr_results = pd.DataFrame(results)


# %% Cell: Visualize the key relationships
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Terrain Variability vs Natural Avalanche Frequency (area-normalized)',
             fontsize=14)

plot_pairs = [
    ('slope_std', 'avi_per_area', 'Slope SD', 'Avis / 10k m²'),
    ('slope_cv', 'avi_per_area', 'Slope CV', 'Avis / 10k m²'),
    ('slope_pct_gt_35', 'avi_per_area', '% cells > 35°', 'Avis / 10k m²'),
    ('aspect_circ_var', 'avi_per_area', 'Aspect circ. var.', 'Avis / 10k m²'),
    ('aspect_n_octants', 'avi_per_area', 'N octants', 'Avis / 10k m²'),
    ('slope_std', 'avi_gt_d1_per_area', 'Slope SD', '>D1 Avis / 10k m²'),
]

for ax, (xcol, ycol, xlabel, ylabel) in zip(axes.flat, plot_pairs):
    x = analysis_df[xcol]
    y = analysis_df[ycol]
    ax.scatter(x, y, s=80, zorder=5)
    for idx, row in analysis_df.iterrows():
        ax.annotate(idx, (row[xcol], row[ycol]),
                    fontsize=8, ha='left', va='bottom')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Spearman rho for annotation
    rho, p = spearmanr(x, y)
    ax.set_title(f'ρ={rho:.2f}, p={p:.2f}', fontsize=10)

plt.tight_layout()
plt.savefig('terrain_variability_vs_avalanches.png', dpi=150, bbox_inches='tight')
plt.show()


# %% Cell: Summary table for the abstract / paper
print("\n" + "="*70)
print("SUMMARY: Terrain Metrics by Start Zone")
print("="*70)

summary_cols = ['area_m2', 'slope_mean', 'slope_std', 'slope_cv',
                'aspect_circ_var', 'aspect_n_octants',
                'avi_count', 'avi_gt_d1', 'avi_per_area']
print(analysis_df[summary_cols].round(3).to_string())

print("\n" + "="*70)
print("CORRELATION SUMMARY (terrain variability → avalanche frequency)")
print("="*70)
print(corr_results[['terrain_metric', 'outcome', 'spearman_rho', 'p_spearman']].round(3).to_string(index=False))


# %% Cell: Bootstrap confidence intervals for correlations (small-sample robust)
def bootstrap_spearman(x, y, n_boot=10000, ci=0.95):
    """
    Bootstrap Spearman correlation for confidence intervals.
    Essential when N is small and p-values from normal theory are unreliable.
    """
    n = len(x)
    rhos = np.zeros(n_boot)
    for i in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        rhos[i], _ = spearmanr(x[idx], y[idx])

    alpha = (1 - ci) / 2
    lo = np.percentile(rhos, alpha * 100)
    hi = np.percentile(rhos, (1 - alpha) * 100)
    observed, _ = spearmanr(x, y)

    print(f"  Spearman rho = {observed:.3f}  [{ci*100:.0f}% CI: {lo:.3f}, {hi:.3f}]")
    return observed, lo, hi


# Run bootstrap on the most promising metric–outcome pair
# Adjust column names based on which metrics look best
print("\nBootstrap CIs for key relationships:")
for t_col, t_name in [('slope_std', 'Slope SD'), ('aspect_circ_var', 'Aspect circ. var')]:
    for o_col, o_name in [('avi_per_area', 'Avis/area'), ('avi_gt_d1_per_area', '>D1/area')]:
        print(f"\n{t_name} vs {o_name}:")
        bootstrap_spearman(
            analysis_df[t_col].values,
            analysis_df[o_col].values
        )


# %% Cell: Poisson regression (if you expand to N >= 10 start zones)
# Count data should be modeled with Poisson or negative binomial, not OLS.
# Uncomment and use when you have enough start zones.

# import statsmodels.api as sm
#
# # Poisson regression: avalanche count ~ terrain variability + offset(log(area))
# # The log(area) offset normalizes for start zone size
#
# X = analysis_df[['slope_std', 'aspect_circ_var']].copy()
# X = sm.add_constant(X)
# y = analysis_df['avi_count'].values
# offset = np.log(analysis_df['area_m2'].values)
#
# poisson_model = sm.GLM(y, X, family=sm.families.Poisson(),
#                         offset=offset).fit()
# print(poisson_model.summary())
#
# # Check for overdispersion
# pearson_chi2 = poisson_model.pearson_chi2 / poisson_model.df_resid
# print(f"\nDispersion parameter: {pearson_chi2:.2f}")
# print("(>1.5 suggests overdispersion → use Negative Binomial instead)")
