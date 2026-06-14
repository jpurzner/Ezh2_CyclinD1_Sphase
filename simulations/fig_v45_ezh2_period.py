"""v45 Figure -- EZH2 cell-cycle dynamics + the period mixture + phase count-fractions.

Panels:
  (A) EZH2 vs cell-cycle phase (G0/G1/S/G2M), GNP & MB -- EZH2 is cycle-driven and peaks in S/G2
      (Fig 4A/B, G2/G0 ~2x). EZH2 integrates over the cycle and is halved at division.
  (B) EZH2 represses CyclinD1: cycle-mean effective CyclinD1 vs EZH2 inhibition (-/+ EZH2i) for GNP
      and MB -- removing EZH2 de-represses CyclinD1 (the rescue mechanism, ~Fig 3C).
  (C) Period distribution (GNP & MB): a right-skewed MIXTURE of fast committed cyclers + a G0 tail,
      from birth noise alone. Median ~22h; the fast mode ~15h.
  (D) MB cell-cycle COUNT-fractions vs the flow target (2N / S / G2M).

Run:  ./venv/bin/python simulations/fig_v45_ezh2_period.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

N = 400
print('running GNP / MB ensembles ...')
ENS = {}
for nm in ('GNP', 'MB', 'GNP+EZH2i', 'MB+EZH2i'):
    drive, ezh2i = v45.condition_drive(nm)
    ENS[nm] = v45.ensemble(drive, EZH2i=ezh2i, N=N)

fig, axes = plt.subplots(2, 2, figsize=(12, 8.6))

# ---- (A) EZH2 by phase ----
axA = axes[0, 0]
phases = ['G0', 'G1', 'S', 'G2M']
for nm, c in [('GNP', '#2166AC'), ('MB', '#B2182B')]:
    means = []
    for ph in phases:
        vals = []
        for cell in ENS[nm]:
            sel = cell['ph'][ph]
            if sel.any():
                vals.append(float(cell['traj']['EZH2'][sel].mean()))
        means.append(np.mean(vals) if vals else np.nan)
    means = np.array(means) / means[0]                # normalize to G0
    axA.plot(phases, means, '-o', color=c, lw=2, ms=8, label=nm)
axA.axhspan(1.8, 2.5, color='#DDDDDD', alpha=0.6)
axA.text(3, 2.15, 'data G1/S, G2/M\n(Fig 4A/B ~2x)', fontsize=7.5, ha='right', color='#555')
axA.set_ylabel('EZH2 (relative to G0)'); axA.set_ylim(0.8, 2.6)
axA.set_title('A  EZH2 is cell-cycle driven (peaks S/G2)', loc='left', fontweight='bold', fontsize=10.5)
axA.legend(fontsize=9)

# ---- (B) EZH2 represses CyclinD1 (de-repression by EZH2i) ----
axB = axes[0, 1]
pairs = [('GNP', 'GNP+EZH2i', '#2166AC'), ('MB', 'MB+EZH2i', '#B2182B')]
xb = np.arange(2)
for i, (base, ezi, c) in enumerate(pairs):
    cd0 = v45.ezh2_stats(ENS[base])['cd_mean']
    cd1 = v45.ezh2_stats(ENS[ezi])['cd_mean']
    axB.plot(xb + i * 0.02, [cd0, cd1], '-o', color=c, lw=2, ms=9, label=base)
    axB.annotate(f'{cd1/cd0:.1f}x', (1, cd1), textcoords='offset points', xytext=(8, 0),
                 fontsize=9, color=c, fontweight='bold')
axB.set_xticks(xb); axB.set_xticklabels(['EZH2 on\n(DMSO)', 'EZH2i\n(de-repressed)'])
axB.set_yscale('log'); axB.set_ylabel('effective CyclinD1 (cycle-mean)')
axB.set_title('B  EZH2 $\\dashv$ CyclinD1 (EZH2i de-represses)', loc='left', fontweight='bold', fontsize=10.5)
axB.legend(fontsize=9, loc='center left')

# ---- (C) period distribution mixture ----
axC = axes[1, 0]
for nm, c in [('GNP', '#2166AC'), ('MB', '#B2182B')]:
    pers = np.array([cell['per'] for cell in ENS[nm]])
    axC.hist(pers, bins=np.arange(10, 60, 2.5), alpha=0.55, color=c, label=f'{nm} (median {np.median(pers):.0f}h)')
    axC.axvline(np.median(pers), color=c, ls='--', lw=1.4)
axC.set_xlabel('cell-cycle period (h)'); axC.set_ylabel('cells')
axC.set_title('C  Period = mixture (fast cyclers + G0 tail)', loc='left', fontweight='bold', fontsize=10.5)
axC.legend(fontsize=8.5)
axC.text(0.97, 0.6, 'right-skew from\nbirth noise alone', transform=axC.transAxes, ha='right',
         fontsize=7.8, color='#444')

# ---- (D) MB count-fractions vs flow target ----
axD = axes[1, 1]
cf, _ = v45.count_fractions(ENS['MB'])
model = [cf['G0'] + cf['G1'], cf['S'], cf['G2M']]
target = [68.2, 15.7, 16.1]                          # 2N / S / G2M (flow; G2M 4N-gate, soft)
labels = ['2N (G0+G1)', 'S', 'G2M']
xd = np.arange(3)
axD.bar(xd - 0.2, model, width=0.4, color='#762A83', label='v45 model')
axD.bar(xd + 0.2, target, width=0.4, color='#BBBBBB', label='flow data (MB DMSO)')
axD.text(2 + 0.2, target[2] + 1, '4N gate\n(soft; direct\nG2+M~2.5h)', fontsize=6.6, ha='center', color='#777')
axD.set_xticks(xd); axD.set_xticklabels(labels); axD.set_ylabel('count-fraction (%)')
axD.set_title('D  MB phase count-fractions', loc='left', fontweight='bold', fontsize=10.5)
axD.legend(fontsize=8.5)

fig.suptitle('v45: EZH2 cell-cycle feedback, the period mixture, and phase fractions', fontsize=12, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig('simulations/fig_v45_ezh2_period.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v45_ezh2_period.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_ezh2_period.{png,pdf}')
