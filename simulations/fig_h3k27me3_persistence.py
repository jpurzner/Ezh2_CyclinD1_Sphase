"""Figure: why H3K27me3 decays too fast, the eraser<->ChIP trade-off, and the bistable resolution.
A: reducing the Gli eraser lengthens the mark (intrinsic me3 half-life) but loses MB<GNP (ChIP) -- no regime does
   both in the continuous model. B: a bistable read-write (cooperative + suppressed nucleation) latches the mark ON
   (GNP) vs OFF (MB, tipped by higher transcription) -> MB<GNP with persistence and NO eraser (reduced 1-locus proto).
Run: ./venv/bin/python simulations/fig_h3k27me3_persistence.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

UPH = 62.8; GLI_MB = 0.19; DEL = 0.0006
erasers = np.array([0.22, 0.11, 0.05, 0.02, 0.0])
halflife_h = (np.log(2) / (DEL + erasers * GLI_MB)) / UPH        # intrinsic me3 half-life (MB)
mb_gnp = np.array([0.56, 0.79, 0.92, 0.98, 1.02])               # from build_model sweep

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.4, 5.4))
fig.subplots_adjust(left=0.07, right=0.93, top=0.86, bottom=0.14, wspace=0.42)

# A: trade-off (dual axis)
x = np.arange(len(erasers))
c1, c2 = '#1e5638', '#c0392b'
axA.plot(x, halflife_h, '-o', color=c1, lw=2.6, ms=7, label='mark half-life')
axA.set_yscale('log'); axA.set_ylabel('intrinsic me3 half-life (h, log)', color=c1, fontsize=11)
axA.tick_params(axis='y', labelcolor=c1)
axA.axhspan(5, 30, color=c1, alpha=0.07)
axA.text(0.1, 8, 'genuine memory\n(~1 cycle)', fontsize=8.5, color=c1)
axB2 = axA.twinx()
axB2.plot(x, mb_gnp, '-s', color=c2, lw=2.6, ms=7, label='MB/GNP mark')
axB2.axhline(1.0, color=c2, ls=':', lw=1, alpha=0.6)
axB2.axhspan(0.4, 0.7, color=c2, alpha=0.07)
axB2.set_ylabel('MB/GNP mark ratio (ChIP)', color=c2, fontsize=11); axB2.tick_params(axis='y', labelcolor=c2)
axB2.set_ylim(0.3, 1.15); axB2.text(3.0, 0.5, 'ChIP target\nMB<GNP', fontsize=8.5, color=c2, ha='center')
axA.set_xticks(x); axA.set_xticklabels([f'{e:g}' for e in erasers])
axA.set_xlabel('Gli eraser strength  k_jmjd3_gli  (← current 0.22 ...... 0 →)', fontsize=11)
axA.set_title('A. The eraser ↔ ChIP trade-off (continuous model)', fontsize=12, fontweight='bold')
for s in ('top',): axA.spines[s].set_visible(False)
axA.text(1.0, 0.35, 'reducing the eraser lengthens the mark\nbut sends MB/GNP → 1 (ChIP lost).\nNO regime does both.',
         fontsize=8.6, color='#444', ha='left', transform=axA.get_xaxis_transform())

# B: bistable resolution (prototype: n=2, g=0.15, ke=0.03, continuous)
states = ['GNP\n(low transcription)', 'MB\n(high transcription)']
vals = [0.67, 0.01]; cols = ['#2f7d54', '#c0392b']
axB.bar([0, 1], vals, width=0.55, color=cols, edgecolor='white')
axB.set_xticks([0, 1]); axB.set_xticklabels(states, fontsize=10)
axB.set_ylabel('H3K27me3 mark (latched)', fontsize=11); axB.set_ylim(0, 0.85)
for i, v in enumerate(vals):
    axB.text(i, v + 0.02, f'{v:.2f}\n{"latched ON" if v>0.3 else "tipped OFF"}', ha='center', fontsize=9, fontweight='bold')
axB.set_title('B. Bistable read-write: MB<GNP with NO eraser', fontsize=12, fontweight='bold')
for s in ('top', 'right'): axB.spines[s].set_visible(False)
axB.text(0.5, 0.76, 'cooperative read-write (M²) + suppressed nucleation\n→ latched memory; MB tipped OFF by its higher\ntranscription-eviction. MB/GNP = 0.01, no eraser.',
         fontsize=8.4, color='#444', ha='center')

fig.suptitle('H3K27me3 persistence: the mark decays too fast because the Gli eraser is strong; the tension and its resolution',
             fontsize=12.8, fontweight='bold', y=0.965)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_h3k27me3_persistence.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
