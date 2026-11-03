"""Kymographs: what the nucleosome lattice actually does over time.

  python simulations/lattice_kymograph.py

(A) one cell establishing the domain from naive chromatin
(B) one cell maintaining it through divisions (per-nucleosome partitioning)
(C) ensemble-mean me3, showing the sawtooth of replicative dilution
(D) core vs flank time courses, with and without division
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'lattice'))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from promoter_lattice import PromoterLattice, Geometry, Params, Inputs, U, AC, ME1, ME2, ME3

OUT = os.path.dirname(os.path.abspath(__file__))
GEOM = Geometry(up_bp=3000, dn_bp=8000, core_up_bp=400, core_dn_bp=700)
PAR = Params(a_rw_prc2=3.5)
INP = Inputs(Gli1=0.3, activity=0.4, tx=0.4)
HOURS, DT, REC = 300.0, 0.05, 1.0
CYCLE = 24.0

STATE_COLS = ['#e8e8e8', '#2e9e5b', '#bcd8ef', '#6aa6d6', '#8e2b3f']
STATE_NAMES = ['unmodified', 'K27ac', 'me1', 'me2', 'me3']


def trace(seed, init_me3, cycle_h=None, n_cells=1):
    """Return (times, states[t,x]) for one cell, or (times, mean me3[t,x])."""
    frames, times = [], []
    lats = [PromoterLattice(GEOM, PAR, seed=seed + i, init_me3=init_me3)
            for i in range(n_cells)]
    nst = int(round(HOURS / DT)); rec_n = int(round(REC / DT))
    t = 0.0; nxt = cycle_h if cycle_h else np.inf
    div_times = []
    for i in range(nst):
        for lat in lats:
            lat.step(DT, INP)
        t += DT
        if t >= nxt:
            for lat in lats:
                lat.replicate()
            div_times.append(t); nxt += cycle_h
        if i % rec_n == 0:
            times.append(t)
            frames.append(lats[0].state.copy() if n_cells == 1
                          else np.mean([l.me3() for l in lats], 0))
    return np.array(times), np.array(frames), np.array(div_times), lats[0].x


fig = plt.figure(figsize=(14.5, 9.2))
gs = fig.add_gridspec(2, 2, hspace=.32, wspace=.22)
cmap = ListedColormap(STATE_COLS)
norm = BoundaryNorm(np.arange(-0.5, 5.1, 1), cmap.N)

# ---- (A) establishment, single cell ----
tA, SA, _, x = trace(1, init_me3=0.0)
axA = fig.add_subplot(gs[0, 0])
axA.imshow(SA, aspect='auto', cmap=cmap, norm=norm, origin='lower',
           extent=[x[0] / 1000, x[-1] / 1000, tA[0], tA[-1]], interpolation='nearest')
axA.set_title('(A) One cell establishing the domain from naive chromatin',
              fontweight='bold', fontsize=10.5)

# ---- (B) maintenance through division, single cell ----
tB, SB, divB, _ = trace(1, init_me3=0.0, cycle_h=CYCLE)
axB = fig.add_subplot(gs[0, 1])
axB.imshow(SB, aspect='auto', cmap=cmap, norm=norm, origin='lower',
           extent=[x[0] / 1000, x[-1] / 1000, tB[0], tB[-1]], interpolation='nearest')
for d in divB[::2]:
    axB.axhline(d, color='k', lw=.4, alpha=.45)
axB.set_title(f'(B) Same cell dividing every {CYCLE:g} h\n'
              '(thin lines = S phase; each nucleosome independently replaced, p=1/2)',
              fontweight='bold', fontsize=10.5)

for ax in (axA, axB):
    ax.axvspan(-GEOM.core_up_bp / 1000, GEOM.core_dn_bp / 1000, color='k', alpha=.08, zorder=3)
    ax.axvline(0, color='k', lw=.8, ls=':')
    ax.set_xlabel('position relative to TSS (kb)'); ax.set_ylabel('time (h)')
axA.legend(handles=[Patch(facecolor=c, label=n, edgecolor='#999', lw=.4)
                    for c, n in zip(STATE_COLS, STATE_NAMES)],
           loc='lower right', fontsize=7.5, framealpha=.95, title='H3K27 state',
           title_fontsize=7.5)

# ---- (C) ensemble mean me3 with divisions ----
tC, MC, divC, _ = trace(500, init_me3=0.0, cycle_h=CYCLE, n_cells=20)
axC = fig.add_subplot(gs[1, 0])
im = axC.imshow(MC, aspect='auto', cmap='magma', origin='lower', vmin=0, vmax=1,
                extent=[x[0] / 1000, x[-1] / 1000, tC[0], tC[-1]], interpolation='bilinear')
axC.axvspan(-GEOM.core_up_bp / 1000, GEOM.core_dn_bp / 1000, color='w', alpha=.18)
axC.axvline(0, color='w', lw=.8, ls=':')
axC.set_xlabel('position relative to TSS (kb)'); axC.set_ylabel('time (h)')
axC.set_title('(C) Ensemble mean me3 (20 cells) — the ChIP-comparable view',
              fontweight='bold', fontsize=10.5)
fig.colorbar(im, ax=axC, label='me3 occupancy', fraction=.046)

# ---- (D) time courses ----
axD = fig.add_subplot(gs[1, 1])
core_m = (x >= -GEOM.core_up_bp) & (x <= GEOM.core_dn_bp)
flank_m = x > 2000
for cyc, ls, lab in [(None, '-', 'no division'), (CYCLE, '-', f'{CYCLE:g} h cycle')]:
    tt, MM, dd, _ = trace(900, init_me3=0.0, cycle_h=cyc, n_cells=20)
    axD.plot(tt, MM[:, flank_m].mean(1), ls, lw=1.9,
             color='#8e2b3f' if cyc is None else '#d98a98', label=f'flank — {lab}')
    axD.plot(tt, MM[:, core_m].mean(1), ls, lw=1.9,
             color='#2f6fb0' if cyc is None else '#9ec4e4', label=f'core — {lab}')
axD.set_xlabel('time (h)'); axD.set_ylabel('me3 occupancy')
axD.set_title('(D) Core stays depleted; flank fills, then is held down by dilution',
              fontweight='bold', fontsize=10.5)
axD.legend(fontsize=8); axD.set_ylim(bottom=0)

fig.suptitle('H3K27me3 nucleosome lattice over time  '
             f'({len(x)} nucleosomes, 190 bp pitch, Gli1={INP.Gli1:g})',
             fontweight='bold', fontsize=13)
p = os.path.join(OUT, 'lattice_kymograph')
fig.savefig(p + '.png', dpi=145, bbox_inches='tight')
fig.savefig(p + '.pdf', bbox_inches='tight')
print('wrote', p + '.png/.pdf')
