"""Consolidated: H3K27me3 mark dynamics in G0 arrest — plateaus, a race, and reversible (monostable).

Three-beat story, one figure:
  (A) The mark does NOT run away in arrest — it overshoots then settles to a modest plateau, because
      the writer EZH2 is E2f-gated and collapses when the cell stops cycling (arrest time-course).
  (B) Whether the arrested mark ends up above or below the cycling mark is a RACE between EZH2-writer
      persistence and mark turnover; the measured EZH2 half-life places GNPs on the RELAX side, near the
      fence (phase diagram; reads the g0_mark_regime cache).
  (C) The mark does NOT make G0 bistable: from a high-mark (arrested) or low-mark (cycling) start the
      cell settles to the SAME steady state at every drive -> monostable, reversible (a sharp but
      single-valued transition). "Sticky/metastable", not a bistable lock. (reads the bistability cache)

Run:  ./venv/bin/python simulations/fig_v44_g0_mark_dynamics.py
Regenerates panel A (fast); B and C read the cached grids (run g0_mark_regime.py / g0_bistability_confirm.py first).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import TwoSlopeNorm
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

LN2 = np.log(2.0); hl = lambda r: LN2 / r / 60.0
C_MARK = '#7a5aa6'; C_TX = '#b0402f'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})

fig = plt.figure(figsize=(14.5, 8.6))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 0.92], hspace=0.42, wspace=0.28)

# ---------- (A) arrest time-course: the mark plateaus ----------
axA = fig.add_subplot(gs[0, 0])
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
try: r.integrator.setValue('maximum_num_steps', 200000)
except Exception: pass
for shh, col, lab in [(1.0, '#117a65', 'cycling (SHH 1.0)'), (0.10, '#2980b9', 'arrested (SHH 0.10)')]:
    r.reset()
    for k, v in GNP.items(): r[k] = v
    r['SHH'] = shh
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            d = r.simulate(0, 42000, 14000, selections=['time', 'Mk', 'EZH2']); break
        except Exception:
            r.reset()
            for k, v in GNP.items(): r[k] = v
            r['SHH'] = shh
    axA.plot(d['time'] / 60 / 24, d['Mk'], color=col, lw=1.3, label=lab)
    if shh == 0.10:
        axA.annotate('overshoot\n(dilution stops)', xy=(2.5, 0.37), xytext=(6, 0.45), fontsize=8,
                     color='#2980b9', arrowprops=dict(arrowstyle='->', color='#2980b9'))
        axA.annotate('settles ~0.31\n(EZH2 collapses to basal)', xy=(22, 0.31), xytext=(14, 0.19), fontsize=8,
                     color='#2980b9', arrowprops=dict(arrowstyle='->', color='#2980b9'))
axA.set_xlabel('time in arrest (days)'); axA.set_ylabel('H3K27me3  Mk')
axA.set_title('(A) The mark plateaus — it does not run away\n(writer EZH2 is E2f-gated, collapses in G0)')
axA.legend(fontsize=8.5, loc='upper right'); axA.grid(alpha=0.15); axA.set_ylim(0.1, 0.58)

# ---------- (B) phase diagram: the race ----------
axB = fig.add_subplot(gs[0, 1])
zr = np.load('simulations/fig_v44_g0_mark_regime_cache.npz')
KDEEZ, DELMK, CYC, ARR = zr['KDEEZ'], zr['DELMK'], zr['CYC'], zr['ARR']
DIFF = ARR - CYC; HLE, HLM = hl(KDEEZ), hl(DELMK)


def _edges(c):
    lc = np.log10(c); e = np.empty(len(c) + 1)
    e[1:-1] = 10 ** (0.5 * (lc[:-1] + lc[1:])); e[0] = 10 ** (lc[0] - 0.5 * (lc[1] - lc[0]))
    e[-1] = 10 ** (lc[-1] + 0.5 * (lc[-1] - lc[-2])); return e


XE, YE = _edges(HLE), _edges(HLM); vlim = np.nanmax(np.abs(DIFF))
pm = axB.pcolormesh(XE, YE, DIFF, cmap='RdBu_r', norm=TwoSlopeNorm(0, -vlim, vlim), shading='flat')
try:
    axB.contour(HLE, HLM, DIFF, levels=[0], colors='k', linewidths=2)
except Exception:
    pass
axB.plot(hl(1.5e-4), hl(7e-4), 'o', ms=12, mfc='yellow', mec='k', mew=1.6, zorder=6)
axB.text(hl(1.5e-4) * 1.08, hl(7e-4) * 1.05, ' measured\n (77 h)', fontsize=8)
axB.set_xscale('log'); axB.set_yscale('log')
axB.text(0.05, 0.92, 'RELAX\n(reversible)', transform=axB.transAxes, fontsize=9, va='top', color='#1a4a7a', fontweight='bold')
axB.text(0.62, 0.16, 'ACCUMULATE', transform=axB.transAxes, fontsize=9, va='top', color='#8a2020', fontweight='bold')
axB.set_xlabel('EZH2 half-life (h) →'); axB.set_ylabel('mark turnover half-life (h) →')
axB.set_title('(B) Which way? A race: EZH2 persistence\nvs mark turnover (measured → RELAX, near fence)')
fig.colorbar(pm, ax=axB, label='arrested − cycling mark', fraction=0.046, pad=0.03)

# ---------- (C) two-IC test: monostable / reversible ----------
axC = fig.add_subplot(gs[1, :])
zb = np.load('simulations/fig_v44_g0_bistability_confirm_cache.npz', allow_pickle=True)
SHHG = zb['SHHG']; names = list(zb['names']); CD_A, CD_C = zb['CD_ARR'], zb['CD_CYC']
jacc = [i for i, n in enumerate(names) if 'ACCUM' in n][0]
axC.plot(SHHG, CD_C[jacc], '-o', color='#2980b9', ms=6, lw=2.2, label='from CYCLING start (low mark)')
axC.plot(SHHG, CD_A[jacc], '-s', color='#c0392b', ms=6, lw=2.2, label='from ARRESTED start (high mark)')
axC.axvspan(0.46, 0.51, color='0.85', alpha=0.6)
axC.text(0.485, axC.get_ylim()[1] * 0.35 if False else 1.55, 'sharp but\nsingle-valued\ntransition', ha='center', fontsize=8.5, color='#555')
axC.set_xlabel('mitogen drive  SHH  (fixed; each point equilibrated from BOTH starts to steady state)')
axC.set_ylabel('mean CyclinD1  Cd')
axC.set_title('(C) The mark does NOT make G0 bistable — both starts converge to one state at every drive → MONOSTABLE, reversible (metastable "sticky" arrest, not a bistable lock)  [ACCUMULATE regime]')
axC.legend(fontsize=9, loc='upper left'); axC.grid(alpha=0.15)

fig.suptitle('H3K27me3 in arrested G0: it plateaus (A), which way it goes is a writer-vs-turnover race (B), and it stays reversible — no bistable lock (C)',
             fontsize=12.5, fontweight='bold', y=0.99)
plt.savefig('simulations/fig_v44_g0_mark_dynamics.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_mark_dynamics.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_g0_mark_dynamics.png')
