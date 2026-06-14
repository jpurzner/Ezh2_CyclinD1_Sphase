"""v45 Figure -- the population rescue + the MYCN floor (the headline result).

Drives the v45 stochastic-commitment ensemble with the v44 Gli/MYCN CyclinD1 drive per condition;
EZH2 represses CyclinD1 INSIDE v45. Panels:
  (A) Quiescent fraction by condition, split into PERMANENT arrest (pRb-/Ki67- exit) vs transient G0.
      Shows Hh-withdrawal raises quiescence, EZH2i rescues, and the MYCN floor (MB+HHi does NOT
      permanently arrest while GNP+HHi does).
  (B) Effective CyclinD1 (cycle-mean) and EZH2 (relative) by condition -- EZH2 represses CyclinD1,
      EZH2i de-represses (the rescue mechanism).
  (C) The rescue arc: MB -> MB+HHi -> MB+HHi+EZH2i quiescent%.

Run:  ./venv/bin/python simulations/fig_v45_rescue.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

N = 400
ORDER = ['GNP', 'GNP-SHH', 'GNP+HHi', 'GNP+EZH2i', 'MB', 'MB+HHi', 'MB+HHi+EZH2i', 'MB+EZH2i']


def summarize(name):
    drive, ezh2i = v45.condition_drive(name)
    cells, st = v45.ensemble(drive, EZH2i=ezh2i, N=N, return_stats=True)
    arr = st['arrest_frac']
    cf, _ = v45.count_fractions(cells) if cells else ({}, None)
    ez = v45.ezh2_stats(cells) if cells else {'mean': np.nan, 'cd_mean': np.nan}
    g0 = cf.get('G0', 0.0)
    return dict(drive=drive, arrest=100 * arr, transG0=(1 - arr) * g0,
                quiesc=100 * arr + (1 - arr) * g0, cd=ez['cd_mean'], ezh2=ez['mean'])


print("running v45 ensembles per condition ...")
S = {n: summarize(n) for n in ORDER}
gE = S['GNP']['ezh2']

fig = plt.figure(figsize=(13.5, 4.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.2, 1.0], wspace=0.34)

# ---- (A) quiescent fraction: permanent arrest + transient G0 ----
axA = fig.add_subplot(gs[0])
x = np.arange(len(ORDER))
arr = np.array([S[n]['arrest'] for n in ORDER])
tg0 = np.array([S[n]['transG0'] for n in ORDER])
axA.bar(x, arr, color='#762A83', label='permanent arrest (pRb$^-$/Ki67$^-$)')
axA.bar(x, tg0, bottom=arr, color='#C2A5CF', label='transient G0 (cycling pool)')
axA.set_xticks(x); axA.set_xticklabels(ORDER, rotation=40, ha='right', fontsize=8)
axA.set_ylabel('quiescent fraction (%)'); axA.set_ylim(0, 100)
axA.set_title('A  Quiescence by condition', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=7.5, loc='upper left', framealpha=0.9)
# MYCN-floor annotation: GNP+HHi arrests, MB+HHi does not
gi, mi = ORDER.index('GNP+HHi'), ORDER.index('MB+HHi')
axA.annotate('MYCN floor:\nMB+HHi does NOT\npermanently arrest',
             xy=(mi, tg0[mi] + 2), xytext=(mi + 0.1, 88), fontsize=7.3, color='#B2182B',
             ha='center', va='top', arrowprops=dict(arrowstyle='->', color='#B2182B', lw=1.1))

# ---- (B) effective CyclinD1 + EZH2 (relative) ----
axB = fig.add_subplot(gs[1])
cd = np.array([S[n]['cd'] for n in ORDER])
ezrel = np.array([S[n]['ezh2'] / gE for n in ORDER])
axB.bar(x - 0.2, cd, width=0.4, color='#1B7837', label='CyclinD1 (cycle-mean, eff.)')
axB.set_xticks(x); axB.set_xticklabels(ORDER, rotation=40, ha='right', fontsize=8)
axB.set_ylabel('effective CyclinD1 (GNP=1)', color='#1B7837')
axB.set_yscale('log'); axB.tick_params(axis='y', colors='#1B7837')
axBt = axB.twinx()
axBt.bar(x + 0.2, ezrel, width=0.4, color='#E08214', label='EZH2 (rel. to GNP)')
axBt.set_ylabel('EZH2 (GNP=1)', color='#E08214'); axBt.tick_params(axis='y', colors='#E08214')
axBt.set_ylim(0, 1.6)
axB.set_title('B  CyclinD1 $\\dashv$ EZH2 (EZH2i de-represses)', loc='left', fontweight='bold', fontsize=11)

# ---- (C) the rescue arc ----
axC = fig.add_subplot(gs[2])
arc = ['MB', 'MB+HHi', 'MB+HHi+EZH2i']
yq = [S[n]['quiesc'] for n in arc]
axC.plot(range(3), yq, '-o', color='#762A83', lw=2.2, ms=9)
for i, (n, y) in enumerate(zip(arc, yq)):
    axC.annotate(f'{y:.0f}%', (i, y), textcoords='offset points', xytext=(0, 9),
                 ha='center', fontsize=9, fontweight='bold')
axC.set_xticks(range(3)); axC.set_xticklabels(['MB', '+HHi\n(vismo)', '+HHi\n+EZH2i'], fontsize=8.5)
axC.set_ylabel('quiescent fraction (%)'); axC.set_ylim(0, max(yq) * 1.4)
axC.set_title('C  EZH2i rescues', loc='left', fontweight='bold', fontsize=11)
axC.annotate('Hh withdrawal\n$\\to$ quiescence', xy=(1, yq[1]), xytext=(0.15, yq[1] + 6),
             fontsize=7.5, color='#444')
axC.annotate('EZH2i\nrescue', xy=(2, yq[2]), xytext=(1.55, yq[2] + 5), fontsize=7.5, color='#1B7837')

fig.suptitle('v45 stochastic-commitment model: mitogen-withdrawal arrest, the MYCN floor, and the EZH2i rescue',
             fontsize=11.5, y=1.02)
fig.savefig('simulations/fig_v45_rescue.png', dpi=170, bbox_inches='tight')
fig.savefig('simulations/fig_v45_rescue.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_rescue.{png,pdf}')
for n in ORDER:
    s = S[n]
    print(f'  {n:16s} arrest {s["arrest"]:4.0f}%  quiesc {s["quiesc"]:4.0f}%  Cd_eff {s["cd"]:5.2f}  EZH2rel {s["ezh2"]/gE:.2f}')
