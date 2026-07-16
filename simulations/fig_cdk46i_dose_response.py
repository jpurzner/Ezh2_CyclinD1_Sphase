"""Figure: CDK4/6i DOSE-RESPONSE on a real concentration axis, fit to JP's two palbo points (1 µM -> 16-18% pRb⁺,
5 µM -> ~2%). Reads cdk46i_dose_response_results.json (resid-based curve from sim_cdk46i_dose_response.py) and maps
residual CDK4/6 activity to concentration via a Langmuir fit resid = K/(K+C), K≈1.5 µM (EMPIRICAL cell-culture
constant; the in-vitro IC50 ~11 nM is irrelevant to cell dosing). Every dividing cell is pRb-Ser807/811⁺; a drug-
independent reversible reserve sits out of cycle. Run: ./venv/bin/python simulations/fig_cdk46i_dose_response.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cdk46i_dose_response_results.json')))
K_UM = 1.5                                        # empirical concentration->activity fit (µM); NOT the in-vitro IC50
resid = np.array(R['resids'])
# resid = K/(K+C)  ->  C = K*(1/resid - 1); drop resid=0 (C->inf)
keep = resid > 0.02
conc = K_UM * (1.0 / resid[keep] - 1.0)
COL = dict(GNP='#2f7d54', MB='#c0392b')

fig, ax = plt.subplots(figsize=(10.2, 6.0))
fig.subplots_adjust(left=0.10, right=0.97, top=0.87, bottom=0.20)

for coh in ('GNP', 'MB'):
    y = np.array(R['cohorts'][coh]['pct_dividing_prb_pos'])[keep]
    order = np.argsort(conc)
    ax.plot(conc[order], y[order], '-o', color=COL[coh], lw=2.4, ms=5, zorder=3,
            label=f'{coh}: model % dividing (all pRb-Ser807/811⁺)')
    ax.axhline(100 * R['reserve_frac'][coh], color=COL[coh], ls=':', lw=1.1, alpha=0.6, zorder=1)

# JP's two MEASURED points
ax.errorbar([1.0], [17], yerr=[1], fmt='D', ms=9, color='#111', zorder=5, capsize=4)
ax.annotate('1 µM: 16–18% pRb⁺\n(measured, reversible)', xy=(1.0, 17), xytext=(1.35, 40),
            fontsize=9.5, color='#111', arrowprops=dict(arrowstyle='->', color='#111', lw=1.1))
ax.errorbar([5.0], [2], yerr=[1], fmt='D', ms=9, color='#111', zorder=5, capsize=4)
ax.annotate('5 µM: ~2% pRb⁺\n(measured, saturating,\nslower washout)', xy=(5.0, 2), xytext=(3.4, 30),
            fontsize=9.5, color='#111', arrowprops=dict(arrowstyle='->', color='#111', lw=1.1))
ax.text(5.9, 100 * R['reserve_frac']['MB'] + 0.5, 'MB reserve ~20%', fontsize=8, color=COL['MB'], ha='right', va='bottom')

ax.set_xlabel('palbociclib (µM)   —  empirical fit K ≈ 1.5 µM (NOT the in-vitro IC50)', fontsize=12)
ax.set_ylabel('% of cells dividing  (= % pRb⁺)', fontsize=12)
ax.set_xlim(0, 6); ax.set_ylim(-3, 103)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('CDK4/6i dose-response fit to the two measured doses (pRb-positive, reversible)',
             fontsize=12.5, fontweight='bold', pad=12)
ax.legend(loc='upper right', fontsize=10, frameon=False)
foot = ('Model = v44 MB/GNP populations (CyclinD1 CV 0.70 + CKI CV 0.42), palbo as graded residual CDK4/6 activity, '
        'resid=K/(K+C). Fits both measured doses; every divider is pRb⁺.\nCAVEAT: the bistable R-point makes the drop '
        'STEEP (IC50 ~0.7 µM, ~saturated by 2 µM) — a testable prediction; widening heterogeneity does not soften it. '
        'Dotted = drug-independent reserve.')
fig.text(0.10, 0.015, foot, fontsize=7.6, color='#666', va='bottom', linespacing=1.5)
out = os.path.join(HERE, 'fig_cdk46i_dose_response.png')
fig.savefig(out, dpi=150)
fig.savefig(out.replace('.png', '.pdf'))
print('conc(µM) at resid grid:', [round(c, 2) for c in sorted(conc)][:8], '...')
print('->', out)
