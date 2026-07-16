"""Figure: kinetics of EZH2/PRC2 inhibition on CyclinD1, the role of H3K27me3, and historical vs current models.
Reads ezh2_prc2_kinetics_results.json (sim_ezh2_prc2_kinetics.py). 2x2:
  A  EZH2i (catalytic, complex stays) vs EED KO (complex removed) -> CyclinD1 fold(t)
  B  the H3K27me3 mark decays after EZH2i -> that decay IS the de-repression timescale
  C  more H3K27me3 (a_rw sweep) -> more baseline repression -> bigger, still-gradual de-repression
  D  historical -> current: instant(direct-EZH2) / explicit-mark memory / PRC2-occupancy / serial me-chain
Run: ./venv/bin/python simulations/fig_ezh2_prc2_kinetics.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'ezh2_prc2_kinetics_results.json')))

fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.0))
fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.09, hspace=0.34, wspace=0.22)
(axA, axB), (axC, axD) = axes

def plot(ax, block, key, colors, lw=2.4, ls='-'):
    d = R[block][key]; t = np.array(d['t'])
    ax.plot(t, np.array(d['cd']), ls, color=colors, lw=lw, label=key)
    return t

# A: EZH2i vs EED KO (CyclinD1 fold)
plot(axA, 'A', 'EZH2i (catalytic)', '#2f7d54')
plot(axA, 'A', 'EED KO (complex removed)', '#7b241c')
axA.axhline(1.0, color='#999', ls=':', lw=1)
axA.set_title('A. EZH2 inhibition vs EED knockout', fontsize=12, fontweight='bold')
axA.set_ylabel('CyclinD1 fold (vs baseline)', fontsize=11)
axA.legend(fontsize=9.5, frameon=False, loc='upper right')
axA.text(20, 1.05, 'EZH2i = catalytic (SET-domain) block: complex STAYS\nbound → PARTIAL de-repression (mark decays but the\nPRC2 occupancy floor remains). EED KO = complex\nremoved → FULLER de-repression.',
         fontsize=8.3, color='#444', va='bottom')

# B: normalized H3K27me3 mark decay for the historical variants that carry a mark (explains D's timescales)
vmark = {'explicit-mark memory': ('#e0a13c', '-'), 'PRC2-occupancy (lumped)': ('#2e86c1', '-'),
         'serial me-chain (current)': ('#1e5638', '-')}
for key, (c, ls) in vmark.items():
    if key in R['C']:
        d = R['C'][key]; mk = np.array(d['mk'])
        if mk[0] > 1e-6:
            axB.plot(np.array(d['t']), mk / mk[0], ls, color=c, lw=2.4, label=key)
axB.set_title('B. The mark decay sets the timescale', fontsize=12, fontweight='bold')
axB.set_ylabel('H3K27me3 mark (fraction of baseline)', fontsize=11)
axB.legend(fontsize=9, frameon=False, loc='upper right')
axB.text(20, 0.5, 'after EZH2i, writing stops and the mark decays\n(turnover + Gli eraser + replication dilution).\nA SLOWER-decaying mark (memory) → slower\nde-repression; the chain adds a me2→me3 lag.',
         fontsize=8.3, color='#444', va='center')

# C: mark amount sweep (a_rw)
cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len([k for k in R['B']])))
for (key, d), c in zip(R['B'].items(), cmap):
    axC.plot(np.array(d['t']), np.array(d['cd']), '-', color=c, lw=2.3, label=key.replace('a_rw ', 'H3K27me3 read-write '))
axC.axhline(1.0, color='#999', ls=':', lw=1)
axC.set_title('C. Changing the amount of H3K27me3', fontsize=12, fontweight='bold')
axC.set_xlabel('hours after EZH2i', fontsize=11); axC.set_ylabel('CyclinD1 fold (vs baseline)', fontsize=11)
axC.legend(fontsize=8.8, frameon=False, loc='upper left', title='mark level', title_fontsize=9)

# D: historical vs current
styles = {'instant (direct EZH2)': ('#c0392b', '--'), 'explicit-mark memory': ('#e0a13c', '-'),
          'PRC2-occupancy (lumped)': ('#2e86c1', '-'), 'serial me-chain (current)': ('#1e5638', '-')}
for key, (c, ls) in styles.items():
    if key in R['C']:
        d = R['C'][key]; axD.plot(np.array(d['t']), np.array(d['cd']), ls, color=c, lw=2.5 if 'current' in key else 2.1, label=key)
axD.axhline(1.0, color='#999', ls=':', lw=1)
axD.set_title('D. Historical models vs current', fontsize=12, fontweight='bold')
axD.set_xlabel('hours after EZH2i', fontsize=11); axD.set_ylabel('CyclinD1 fold (vs baseline)', fontsize=11)
axD.legend(fontsize=9, frameon=False, loc='center right')
axD.text(20, 5.6, 'instant = old wrong model (step, no mark timescale)\nmemory = slow mark-decay (over-large fold)\ncurrent chain = gradual, partial (complex stays)',
         fontsize=8.3, color='#444', va='top')

for ax in (axA, axB, axC, axD):
    ax.set_xlabel(ax.get_xlabel() or 'hours after EZH2i', fontsize=11)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
    ax.axvspan(-1, 0, color='none')

fig.suptitle('Kinetics of EZH2 / PRC2 inhibition on CyclinD1 — the role and amount of H3K27me3, historical vs current',
             fontsize=13.5, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_ezh2_prc2_kinetics.png')
fig.savefig(out, dpi=140)
fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
