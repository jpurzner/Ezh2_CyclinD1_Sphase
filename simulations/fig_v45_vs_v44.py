"""v44 vs v45 -- a behavioral comparison of the two models on the same conditions.

v44 (committed working model): Heldt G1/S core + explicit replication + CyclinB/CDK1 mitosis +
  growth-gated restriction point + Skp2-p27 feedforward. A SINGLE DETERMINISTIC cell simulated
  continuously; commitment is a phenomenological GROWTH/SIZE gate; "arrest" is a population/fractional
  effect (parameter heterogeneity), not single-cell.
v45 (stochastic-commitment successor): a reduced kernel with a BISTABLE CDK2-p27 toggle + Rb-dilution
  G1 timer + mitogen-gated CyclinD bootstrap + mass-capped quiescence + dynamic EZH2. A STOCHASTIC
  ENSEMBLE of single cells (noisy birth); commitment is an emergent bifurcation; arrest is real
  (mass-capped) at the single-cell level.

Panels: (A) MB cell-cycle COUNT-fractions, both models vs flow data; (B) cell-cycle period -- v44 single
deterministic value vs the v45 distribution (mixture); (C) the rescue, each model's native read-out.

Run:  ./venv/bin/python simulations/fig_v45_vs_v44.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from simulations.validate_v44 import run as run44, count_divisions, classify, P27_THR, P16_MB, P18_MB, KSYP21_MB
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

_MB = dict(p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB)
print('running v44 conditions (each a ~200h deterministic sim) ...')
v44 = {
    'MB':           run44(shh=0.5, ptch1_cn=0.1, gdc=0.0, ezh2i=0.0, mycn_amp=2.8, **_MB),
    'MB+HHi':       run44(shh=0.5, ptch1_cn=0.1, gdc=1.0, ezh2i=0.0, mycn_amp=2.8, **_MB),
    'MB+HHi+EZH2i': run44(shh=0.5, ptch1_cn=0.1, gdc=1.0, ezh2i=1.0, mycn_amp=2.8, **_MB),
    'GNP':          run44(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0),
}
print('running v45 ensembles ...')
v45_cells = {}
for nm in ('MB', 'MB+HHi', 'MB+HHi+EZH2i', 'GNP'):
    drive, ezh2i = v45.condition_drive(nm)
    v45_cells[nm] = v45.ensemble(drive, EZH2i=ezh2i, N=400, return_stats=True)

C44, C45, CD = '#2166AC', '#762A83', '#BBBBBB'
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# (A) MB count-fractions: 2N / S / G2M, both models vs data
axA = axes[0]
cf44 = classify(v44['MB'], P27_THR)[0]
m44 = [cf44['G0'] + cf44['G1'], cf44['S'], cf44['G2']]
cells_mb = v45_cells['MB'][0]; cf45, _ = v45.count_fractions(cells_mb)
m45 = [cf45['G0'] + cf45['G1'], cf45['S'], cf45['G2M']]
data = [68.2, 15.7, 16.1]
x = np.arange(3); w = 0.27
axA.bar(x - w, m44, w, color=C44, label='v44')
axA.bar(x, m45, w, color=C45, label='v45')
axA.bar(x + w, data, w, color=CD, label='flow data (MB)')
axA.set_xticks(x); axA.set_xticklabels(['2N (G0+G1)', 'S', 'G2M']); axA.set_ylabel('count-fraction (%)')
axA.set_title('A  MB phase fractions (vs data)', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=8.5)
axA.text(0.02, 0.92, 'v44 over-counts 2N / under-counts S;\nv45 closer (mixture + count-frac)',
         transform=axA.transAxes, fontsize=6.8, color='#444', va='top')

# (B) period: v44 deterministic vs v45 distribution
axB = axes[1]
pers45 = np.array([c['per'] for c in v45_cells['GNP'][0]])
axB.hist(pers45, bins=np.arange(8, 56, 2.5), color=C45, alpha=0.65, label=f'v45 GNP (median {np.median(pers45):.0f}h)')
_, _, per44_arr = count_divisions(v44['GNP'])
per44 = float(np.mean(per44_arr)) if len(per44_arr) else np.nan
axB.axvline(per44, color=C44, lw=2.5, label=f'v44 GNP (deterministic {per44:.0f}h)')
axB.set_xlabel('cell-cycle period (h)'); axB.set_ylabel('cells (v45)')
axB.set_title('B  Period: deterministic vs distribution', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=8)
axB.text(0.97, 0.6, 'v44: one growth-clock period\nv45: a mixture (fast cyclers\n+ a G0 tail) from birth noise',
         transform=axB.transAxes, ha='right', fontsize=6.8, color='#444')

# (C) the rescue, each model's native read-out (normalized to MB)
axC = axes[2]
conds = ['MB', 'MB+HHi', 'MB+HHi+EZH2i']
div44 = np.array([count_divisions(v44[c])[0] for c in conds], float)
# v45: cycling fraction = 1 - quiescent fraction
cyc45 = []
for c in conds:
    cells, st = v45_cells[c]
    cf, _ = v45.count_fractions(cells)
    quiesc = 100 * st['arrest_frac'] + (1 - st['arrest_frac']) * cf['G0']
    cyc45.append(100 - quiesc)
cyc45 = np.array(cyc45)
xx = np.arange(3)
axC.plot(xx, div44 / div44[0], '-o', color=C44, lw=2.2, ms=8, label='v44: divisions/168h (norm.)')
axC.plot(xx, cyc45 / cyc45[0], '-s', color=C45, lw=2.2, ms=8, label='v45: cycling fraction (norm.)')
axC.set_xticks(xx); axC.set_xticklabels(['MB', '+HHi', '+HHi\n+EZH2i'])
axC.set_ylabel('proliferation (norm. to MB)')
axC.set_title('C  Rescue: both reproduce it', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=8, loc='lower left')
axC.text(0.5, 0.92, 'HHi $\\downarrow$ then EZH2i rescues', transform=axC.transAxes, ha='center', fontsize=7.5, color='#444')

fig.suptitle('v44 (deterministic growth-gated cell-cycle engine) vs v45 (stochastic bistable-commitment ensemble): same biology, different mechanism',
             fontsize=10.5, y=1.01)
fig.tight_layout()
fig.savefig('simulations/fig_v45_vs_v44.png', dpi=160, bbox_inches='tight')
fig.savefig('simulations/fig_v45_vs_v44.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_vs_v44.{png,pdf}')
print(f'  MB fractions 2N/S/G2: v44 {[round(v) for v in m44]}  v45 {[round(v) for v in m45]}  data {data}')
print(f'  GNP period: v44 {per44:.0f}h (deterministic), v45 median {np.median(pers45):.0f}h (mixture)')
print(f'  v44 divisions MB/+HHi/+rescue: {div44}; v45 cycling%: {[round(v) for v in cyc45]}')
