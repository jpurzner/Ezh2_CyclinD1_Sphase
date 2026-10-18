"""How replicative dilution and Ezh2 replenishment set the H3K27me3 level (JP 2026-09-19, conceptual figure).
Ezh2 is an S-phase-associated writer: at each round of DNA replication the mark is diluted across newly
synthesized histones, and Ezh2 rebuilds it between/through S-phase. As mitogen drive speeds up division,
the mark is maintained only when replenishment keeps pace. Three regimes (Ezh2 output scaled via kTlEZ):
FAILING (replenishment < dilution -> mark erodes, Cyclin D1 de-repressed), SUFFICIENT (balanced -> mark
held), EXCESSIVE (replenishment >> dilution -> mark saturated, beyond what maintenance requires).

Generated from the baked model. All three cells start from the same high mark, then diverge under their
Ezh2 capacity while continuing to cycle.

Run:  ./venv/bin/python simulations/fig_v44_dilution_replenishment.py   (~8 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

KTLEZ0 = 0.015098862630302894
C_FAIL, C_SUFF, C_EXC = '#c0392b', '#e0a800', '#1f6f4a'
C_MK, C_EZ, C_CD = '#6a2c8f', '#2472a4', '#b0402f'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})

M = build_model_v44()
RR = te.loada(M); RR.integrator.setValue('relative_tolerance', 1e-7); RR.integrator.setValue('absolute_tolerance', 1e-9)
RR.integrator.setValue('maximum_num_steps', 12000000)


def mb_conditions():
    RR['SHH'] = 0.5; RR['Ptch1_copy_number'] = 0.3
    try: RR['Cd2_expr'] = RR['CD2_EXPR_MB']
    except Exception: pass
    RR['p18'] = 1.73; RR['p19'] = 0.58


def run_regime(scale, T_eq=6000, T_show=13000, npts=13000):
    """Equilibrate at a high (excess) writer to establish a high mark, then switch to `scale` and record."""
    mb_conditions()
    RR['kTlEZ'] = KTLEZ0 * 2.0            # common high-mark start
    RR.reset()
    RR.simulate(0, T_eq, 200, selections=['time'])
    RR['kTlEZ'] = KTLEZ0 * scale          # switch to this regime's Ezh2 output
    r = RR.simulate(T_eq, T_eq + T_show, npts,
                    selections=['time', 'Mk', 'EZH2', 'Cd', 'Dna', 'aRc'])
    t = (r['time'] - T_eq) / 60.0
    return dict(t=t, Mk=r['Mk'], EZH2=r['EZH2'], Cd=r['Cd'], Dna=r['Dna'], aRc=r['aRc'])


fail = run_regime(0.05)
suff = run_regime(0.25)
exc = run_regime(2.0)


def divisions(d):
    dna = d['Dna']; idx = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    return d['t'][idx]


fig = plt.figure(figsize=(14.5, 9))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 0.95], hspace=0.42, wspace=0.28)

# ---------- (A) three regimes: mark sawtooth, common high start, diverging ----------
axA = fig.add_subplot(gs[0, :])
for d, col, lab in [(fail, C_FAIL, 'FAILING — replenishment < dilution'),
                    (suff, C_SUFF, 'SUFFICIENT — replenishment ≈ dilution'),
                    (exc, C_EXC, 'EXCESSIVE — replenishment ≫ dilution')]:
    axA.plot(d['t'], d['Mk'], color=col, lw=1.6, label=lab)
for x in divisions(suff):
    axA.axvline(x, color='#bbb', lw=0.6, alpha=0.5, zorder=0)
axA.set_xlabel('time (h)   — cells dividing under mitogen drive'); axA.set_ylabel('H3K27me3 at the locus (occupancy)')
axA.set_ylim(0, 1.03); axA.set_xlim(0, suff['t'][-1])
axA.set_title('(A) All cells start with a full mark; its fate depends on whether Ezh2 replenishment keeps pace with replicative dilution')
axA.legend(fontsize=8.5, loc='center right', framealpha=0.95)
axA.text(0.01, 0.06, 'grey lines = cell divisions (each dilutes the mark)', transform=axA.transAxes, fontsize=8, color='#666', style='italic')
axA.annotate('mark eroded → gene de-repressed', (fail['t'][-1]*0.6, fail['Mk'][-3000:].mean()),
             (fail['t'][-1]*0.35, 0.28), fontsize=8, color=C_FAIL,
             arrowprops=dict(arrowstyle='->', color=C_FAIL))

# ---------- (B) mechanism (schematic): Ezh2 (S-phase writer) vs replication dilution over one cycle ----------
axB = fig.add_subplot(gs[1, 0])
tc = np.linspace(0, 24, 700)
S0, S1 = 8.0, 16.0            # S phase window
axB.axvspan(S0, S1, color='#cfe3f5', alpha=0.6, zorder=0)
# idealized H3K27me3: high in G1, diluted across S as new histones are laid down, re-methylated by Ezh2 into G2
mark = np.where(tc < S0, 0.92,
        np.where(tc < S1, 0.92 - 0.34 * (tc - S0) / (S1 - S0),      # dilution during replication
                 0.58 + 0.34 * (tc - S1) / (24 - S1)))              # Ezh2 re-methylation
# idealized Ezh2 writer: low in G1, elevated through S phase (S-phase-associated), declines late
ez = 0.25 + 0.70 / (1 + np.exp(-(tc - 9) / 1.6)) * np.exp(-np.maximum(tc - 18, 0) / 4)
axB.plot(tc, mark, color=C_MK, lw=2.2, label='H3K27me3 (mark)', zorder=3)
axr = axB.twinx()
axr.plot(tc, ez, color=C_EZ, lw=1.8, ls='--', label='Ezh2 (writer)')
axr.set_ylim(0, 1.15); axr.set_ylabel('Ezh2 activity (a.u.)', color=C_EZ); axr.tick_params(axis='y', labelcolor=C_EZ)
axB.set_ylim(0.5, 1.0); axB.set_xlim(0, 24)
axB.set_xticks([4, 12, 20]); axB.set_xticklabels(['G1', 'S', 'G2/M'])
axB.set_xlabel('cell-cycle phase (schematic)'); axB.set_ylabel('H3K27me3 (a.u.)', color=C_MK); axB.tick_params(axis='y', labelcolor=C_MK)
axB.set_title('(B) Mechanism: replication dilutes; Ezh2 (S-phase writer) re-methylates')
axB.annotate('replication\ndilutes', (12, 0.70), (12, 0.55), fontsize=8.5, color='#3a6ea5', ha='center',
             arrowprops=dict(arrowstyle='->', color='#3a6ea5'))
axB.annotate('Ezh2\nre-methylates', (20, 0.80), (21.5, 0.62), fontsize=8.5, color=C_EZ, ha='center',
             arrowprops=dict(arrowstyle='->', color=C_EZ))
h1, l1 = axB.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
axB.legend(h1 + h2 + [plt.Rectangle((0, 0), 1, 1, fc='#cfe3f5', alpha=0.6)], l1 + l2 + ['S phase'],
           fontsize=8, loc='lower center', framealpha=0.9)

# ---------- (C) consequence: Cyclin D1 (the repressed gene) in each regime ----------
axC = fig.add_subplot(gs[1, 1])
labels = ['FAILING', 'SUFFICIENT', 'EXCESSIVE']
mk_ss = [fail['Mk'][-4000:].mean(), suff['Mk'][-4000:].mean(), exc['Mk'][-4000:].mean()]
cd_ss = [fail['Cd'][-4000:].mean(), suff['Cd'][-4000:].mean(), exc['Cd'][-4000:].mean()]
x = np.arange(3); w = 0.36
b1 = axC.bar(x - w/2, mk_ss, w, color=C_MK, label='H3K27me3 (mark)')
axr2 = axC.twinx()
b2 = axr2.bar(x + w/2, cd_ss, w, color=C_CD, label='Cyclin D1')
axC.set_xticks(x); axC.set_xticklabels(labels, fontsize=9)
axC.set_ylabel('H3K27me3 (steady)', color=C_MK); axC.tick_params(axis='y', labelcolor=C_MK); axC.set_ylim(0, 1.03)
axr2.set_ylabel('Cyclin D1 (steady)', color=C_CD); axr2.tick_params(axis='y', labelcolor=C_CD)
axC.set_title('(C) Consequence: low mark → Cyclin D1 up; excess → floored')
axC.legend(handles=[b1, b2], labels=['H3K27me3', 'Cyclin D1'], fontsize=8.5, loc='upper center')

fig.suptitle('Replicative dilution vs Ezh2 replenishment of H3K27me3 — the mark is held only when '
             'replenishment keeps pace with division', fontsize=12.5, fontweight='bold', y=0.995)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_dilution_replenishment')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
print('steady Mk  fail/suff/exc:', [round(x,3) for x in mk_ss])
print('steady Cd  fail/suff/exc:', [round(x,2) for x in cd_ss])
print('divisions in window (suff):', len(divisions(suff)))
