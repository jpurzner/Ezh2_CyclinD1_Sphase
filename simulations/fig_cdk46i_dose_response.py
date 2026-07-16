"""Figure: CDK4/6i DOSE-RESPONSE, GNP vs MB -- the pRb-Ser807/811⁺ residual is dose-limited (not a resistant clone).
Reads cdk46i_dose_response_results.json (sim_cdk46i_dose_response.py). Every dividing cell is pRb⁺; a drug-
independent reversible reserve sits out of cycle at every dose. Run: ./venv/bin/python simulations/fig_cdk46i_dose_response.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cdk46i_dose_response_results.json')))
resids = np.array(R['resids']); inhib = 100 * (1 - resids)
COL = dict(GNP='#2f7d54', MB='#c0392b')

fig, ax = plt.subplots(figsize=(10.4, 6.0))
fig.subplots_adjust(left=0.10, right=0.97, top=0.87, bottom=0.20)

for coh in ('GNP', 'MB'):
    y = np.array(R['cohorts'][coh]['pct_dividing_prb_pos'])
    ax.plot(inhib, y, '-o', color=COL[coh], lw=2.4, ms=6, zorder=3, label=f'{coh}: % dividing (all pRb-Ser807/811⁺)')
    res = 100 * R['reserve_frac'][coh]
    ax.axhline(res, color=COL[coh], ls=':', lw=1.2, alpha=0.7, zorder=1)
ax.fill_between(inhib, 0, np.array(R['cohorts']['MB']['pct_dividing_prb_pos']), color=COL['MB'], alpha=0.06, zorder=0)

# JP's palbo data points
ax.axhspan(16, 18, color='#7f8c8d', alpha=0.16, zorder=0)
ax.text(1.5, 34, "1 µM: 16–18% pRb⁺\n(sub-saturating, reversible)", fontsize=9.5, color='#333', va='center')
ax.annotate("5 µM: ~2% pRb⁺\n(saturating; slower to\nwash out, still reversible)", xy=(88, 2), xytext=(58, 32),
            fontsize=9.5, color='#333', ha='left', arrowprops=dict(arrowstyle='->', color='#666', lw=1.1))
ax.text(99, 20.5, "MB reserve ~20%", fontsize=8.2, color=COL['MB'], alpha=0.85, va='bottom', ha='right')
ax.text(99, 2.5, "GNP reserve ~2%", fontsize=8.2, color=COL['GNP'], alpha=0.85, va='bottom', ha='right')

ax.set_xlabel('CDK4/6 inhibition  (dose →)', fontsize=12)
ax.set_ylabel('% of cells dividing  (= % pRb⁺)', fontsize=12)
ax.set_xlim(-2, 102); ax.set_ylim(-3, 103)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('CDK4/6i dose-response: a dose-limited, pRb-positive, reversible residual',
             fontsize=12.5, fontweight='bold', pad=12)
ax.legend(loc='upper right', fontsize=10, frameon=False)

foot = ('v44 populations (CyclinD1 CV 0.70 + CKI CV 0.42); palbo = graded residual CDK4/6 activity. Every dividing '
        'cell keeps Rb hyperphosphorylated (pRb-Ser807/811⁺) — the residual\nis the high-CyclinD tail beating partial '
        'inhibition, reversible and gone at saturating dose. Dotted lines = drug-independent reversible reserve. No '
        'resistant clone.')
fig.text(0.10, 0.015, foot, fontsize=7.6, color='#666', va='bottom', linespacing=1.5)
out = os.path.join(HERE, 'fig_cdk46i_dose_response.png')
fig.savefig(out, dpi=150)
print('MB dividing (no drug -> saturating):', [round(v) for v in R['cohorts']['MB']['pct_dividing_prb_pos']])
print('->', out)
