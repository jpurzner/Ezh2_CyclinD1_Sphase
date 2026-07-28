"""Distributions for the v44 parameters that are DRAWN from distributions (the population-of-models
heterogeneity), plus the emergent cell-cycle-duration distribution. JP 2026-07-27.

Panels A-C: the three per-lineage input draws (CyclinD1 abundance, birth p27, EZH2 translation).
Panel D: cell-cycle duration -- Nakashima 2015 empirical (mean 15.9h, CV 0.23; lineage-correlated:
siblings CV 0.12 vs non-siblings CV 0.32). Writes fig_param_distributions.pdf + .png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INK, MUT, FAINT, GRID = '#2d2d2d', '#666666', '#999999', '#dddddd'
C_CD, C_P27, C_EZ, C_CYC = '#2c7fb8', '#c0392b', '#7d5ba6', '#4c8c6b'
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans', 'axes.linewidth': 0.8, 'pdf.fonttype': 42})
rng = np.random.default_rng(7)
N = 300000

def lognorm(median, sdlog):
    return median * np.exp(sdlog * rng.standard_normal(N))

def cv_ln(sdlog):
    return np.sqrt(np.exp(sdlog**2) - 1)

def panel(ax, samp, color, title, sub, unit, celltype, extra=None):
    lo, hi = np.percentile(samp, [1, 99])
    ax.hist(samp[(samp >= lo) & (samp <= hi)], bins=90, density=True, color=color, alpha=0.30,
            edgecolor='none')
    from scipy.stats import gaussian_kde
    xs = np.linspace(lo, hi, 400); kde = gaussian_kde(samp[(samp >= lo) & (samp <= hi)])
    ax.plot(xs, kde(xs), color=color, lw=2.2)
    med = np.median(samp); p10, p90 = np.percentile(samp, [10, 90])
    ax.axvline(med, color=color, ls='-', lw=1.4)
    ax.axvspan(p10, p90, color=color, alpha=0.08)
    yl = ax.get_ylim()[1]
    ax.text(med, yl * 1.02, f"median {med:.3g}", color=color, fontsize=8, ha='center', fontweight='bold')
    ax.set_title(title, fontsize=11, fontweight='bold', color=INK, loc='left', pad=16)
    ax.text(0.0, 1.045, sub, transform=ax.transAxes, fontsize=8.4, color=MUT)
    ax.set_xlabel(unit, fontsize=8.5, color=MUT)
    ax.set_yticks([])
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(GRID); ax.tick_params(length=0, labelsize=7.8)
    ax.set_xlim(lo, hi)
    tag = dict(GNP='#2c7fb8', MB='#c0392b', both='#666666')[celltype]
    ax.text(0.985, 0.93, {'GNP': 'P7 GNP', 'MB': 'MB only', 'both': 'both cell types'}[celltype],
            transform=ax.transAxes, ha='right', fontsize=8, color=tag, fontweight='bold')
    if extra: ax.text(0.985, 0.80, extra, transform=ax.transAxes, ha='right', va='top',
                      fontsize=7.6, color=FAINT)

fig, ax = plt.subplots(2, 2, figsize=(13.5, 8.6))
fig.subplots_adjust(hspace=0.5, wspace=0.2, left=0.05, right=0.965, top=0.86, bottom=0.085)
fig.suptitle('v44 distributed parameters — population-of-models heterogeneity', x=0.05, ha='left',
             fontsize=14.5, fontweight='bold', color=INK, y=0.965)
fig.text(0.05, 0.905, 'each proliferating cell draws these independently at birth; the spread is what '
         'produces the population fold-response and cell-to-cell variation',
         fontsize=9.5, color=MUT, ha='left')

# A: CyclinD1 abundance (k_Cd_translation)
sA = lognorm(0.801, 0.70)
panel(ax[0, 0], sA, C_CD, 'CyclinD1 abundance', 'k_Cd_translation · lognormal',
      'protein translation scale (a.u.)', 'both', extra=f"CV ≈ {cv_ln(0.70):.2f}\ngrounds the measured\nCyclinD1 G0/1 CV ≈ 0.70")
# B: birth p27 (P21_div) — MB only
sB = lognorm(0.6, 0.42)
panel(ax[0, 1], sB, C_P27, 'birth p27', 'P21_div · lognormal',
      'inherited p27 at division (a.u.)', 'MB', extra=f"CV ≈ {cv_ln(0.42):.2f}\nGNP fixed at 0\n(no transient G0)")
# C: EZH2 translation (kTlEZ)
sC = lognorm(0.004, 0.51)
panel(ax[1, 0], sC, C_EZ, 'EZH2 translation rate', 'kTlEZ · lognormal',
      'EZH2 translation scale (a.u.)', 'both', extra=f"CV ≈ {cv_ln(0.51):.2f}\nEZH2 buffers CyclinD1\nvariance (governor)")
# D: cell-cycle duration — Nakashima 2015 empirical (lognormal matched to mean 15.9h, CV 0.233)
s = 0.230; mu = np.log(15.9) - s**2 / 2
sD = np.exp(mu + s * rng.standard_normal(N))
panel(ax[1, 1], sD, C_CYC, 'cell-cycle duration', 'Nakashima 2015 · empirical (P7 GNP)',
      'cycle length (h)', 'GNP', extra="mean 15.9 h · SD 3.7 · CV 0.23\nlineage-correlated:\nsiblings CV 0.12 ≪ non-sib CV 0.32")
ax[1, 1].text(0.02, 0.80, "⚠ model output here is near-delta\n(cycle is growth-timed, ~invariant to\nthe CyclinD1 draw) — a known gap",
              transform=ax[1, 1].transAxes, fontsize=7.4, color='#c0392b', va='top')

out = os.path.join(ROOT, 'simulations', 'fig_param_distributions')
fig.savefig(out + '.pdf', bbox_inches='tight'); fig.savefig(out + '.png', dpi=150, bbox_inches='tight')
print('wrote', out + '.pdf / .png')
