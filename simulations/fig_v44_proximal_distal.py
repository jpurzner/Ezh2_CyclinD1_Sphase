"""Two-compartment PRC2 / H3K27me3 model — the reservoir-sustained transient-G0 amplifier (JP 2026-08-09).

JP's molecular mechanism: CDK6 & CyclinD1 are heavily H3K27me3-marked yet highly EXPRESSED because the
PROXIMAL promoter is kept depleted of the PRC2 complex (Pol2 elongation clears it), while the broad DISTAL
region holds H3K27me3 as a RESERVOIR. Repression is the PRC2 COMPLEX (not the mark); the reservoir recruits
it via read-write. When a cell stops cycling (Pol2 stops elongating) the reservoir rapidly re-loads the
proximal promoter -> reservoir-SUSTAINED repression that deepens & prolongs the arrest, and RELEASES when
mitogen/cycling returns (reversible = TRANSIENT G0).

Four beats:
  (A) MARKED-BUT-EXPRESSED: cycling cell — distal reservoir Mk HIGH, proximal occupancy P_prox ~0
      (Pol2 clears it) -> CDK6/CyclinD1 highly expressed. The paradox resolved.
  (B) AMPLIFIER ENGAGES IN ARREST: drive a cell to arrest (high CDKI) -> E2f collapses -> elongation stops
      -> the persistent reservoir re-loads P_prox -> CDK6 repressed (the reservoir SUSTAINS PRC2 as EZH2 falls).
  (C) SCENARIO MAP: amplifier strength -> transient dwell length & validation. Reversible TRANSIENT window
      (30/32, dwell ~90-150 h) between a leaky limit and a PERMANENT-LOCK cliff; the deep permanent latch is
      MORE data-excluded (<=29/32; loses Palbo-EZH2 + MB+HHi) = the known chromatin-latch bound, new direction.
  (D) TRANSIENT-G0 HYSTERESIS: arrest, then drop CDKI -> reservoir-sustained dwell delays re-entry
      (transient C) vs immediate (amplifier OFF) vs never (permanent A). The reversible multi-cycle dwell.

Run:  ./venv/bin/python simulations/fig_v44_proximal_distal.py   (~3-4 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)
C_RES = '#7a5aa6'   # reservoir (distal Mk)
C_PROX = '#b0402f'  # proximal PRC2 occupancy
C_EXPR = '#117a65'  # expressed gene (CDK6/CyclinD1)
C_E2F = '#2472a4'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})


def make_rr(amp=None):
    P = dict(BP)
    if amp:
        P.update(amp)
    M = build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_proximal_distal=True, params=P)
    rr = te.loada(M)
    rr.integrator.setValue('relative_tolerance', 1e-7)
    rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr, P


def setp(rr, P, mult):
    for k, v in P.items():
        try: rr[k] = v
        except Exception: pass
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult; rr['p19'] = 0.58 * mult; rr['kSyP21'] = BP['kSyP21'] * mult


fig = plt.figure(figsize=(15, 9.2))
gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.30)

# ---------- (A) marked-but-expressed: cycling cell ----------
axA = fig.add_subplot(gs[0, 0])
rr, P = make_rr()
setp(rr, P, 1.0); rr.reset()
d = rr.simulate(0, 12000, 6000, selections=['time', 'Mk', 'P_prox', 'cdk6'])
th = d['time'] / 60
axA.plot(th, d['Mk'], color=C_RES, lw=1.6, label='distal reservoir  H3K27me3 (Mk)')
axA.plot(th, d['P_prox'], color=C_PROX, lw=1.6, label='proximal PRC2 occupancy (P$_{prox}$)')
axA.set_ylim(-0.03, 1.02); axA.set_xlabel('time (h)'); axA.set_ylabel('occupancy [0–1]')
axA.set_title('(A) Marked but expressed — cycling cell')
axr = axA.twinx()
axr.plot(th, d['cdk6'], color=C_EXPR, lw=1.4, ls='--', label='CDK6 (expressed)')
axr.set_ylabel('CDK6 level', color=C_EXPR); axr.tick_params(axis='y', labelcolor=C_EXPR)
axr.set_ylim(0, max(6, d['cdk6'].max() * 1.15))
h1, l1 = axA.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
axA.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='center right', framealpha=0.9)
axA.text(0.03, 0.06, 'reservoir HIGH, proximal ~0\n→ Pol2 elongates → gene ON',
         transform=axA.transAxes, fontsize=8, style='italic',
         bbox=dict(boxstyle='round', fc='#f4f1fa', ec=C_RES, alpha=0.9))

# ---------- (B) amplifier engages in sustained arrest ----------
axB = fig.add_subplot(gs[0, 1])
rr, P = make_rr()
setp(rr, P, 6.0); rr.reset()
d = rr.simulate(0, 30000, 12000, selections=['time', 'E2f', 'elong_gate', 'P_prox', 'cdk6'])
th = d['time'] / 60
axB.plot(th, d['E2f'], color=C_E2F, lw=1.5, label='E2f (proliferation)')
axB.plot(th, d['elong_gate'], color='#888', lw=1.3, ls=':', label='elongation gate')
axB.plot(th, d['P_prox'], color=C_PROX, lw=1.8, label='proximal PRC2 (P$_{prox}$) — LOADS')
axB.set_ylim(-0.03, 1.05); axB.set_xlabel('time (h)'); axB.set_ylabel('level [0–1]')
axB.set_title('(B) Sustained arrest — reservoir re-loads the proximal promoter')
axr = axB.twinx()
axr.plot(th, d['cdk6'], color=C_EXPR, lw=1.5, ls='--', label='CDK6 — repressed')
axr.set_ylabel('CDK6 level', color=C_EXPR); axr.tick_params(axis='y', labelcolor=C_EXPR)
axr.set_ylim(0, 6)
h1, l1 = axB.get_legend_handles_labels(); h2, l2 = axr.get_legend_handles_labels()
axB.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='center right', framealpha=0.9)
axB.text(0.03, 0.55, 'E2f collapses → Pol2 stops →\nreservoir sustains PRC2\n(does NOT need EZH2)',
         transform=axB.transAxes, fontsize=8, style='italic',
         bbox=dict(boxstyle='round', fc='#fbeeea', ec=C_PROX, alpha=0.9))

# ---------- (C) validation cost vs amplifier strength (the strong regime is data-excluded) ----------
axC = fig.add_subplot(gs[1, 0])
# amp_scenario_wf results: (a_P, passed/32, verdict)
scen = [(0.05, 30, 'T'), (0.08, 30, 'T'), (0.10, 30, 'T'), (0.12, 30, 'P'), (0.15, 30, 'P'),
        (0.30, 29, 'P'), (0.50, 28, 'P')]
ax_ap = [s[0] for s in scen]; ax_ps = [s[1] for s in scen]
axC.plot(ax_ap, ax_ps, 'o-', color=C_PROX, lw=1.7, ms=7)
axC.axhline(31, color='#117a65', ls='--', lw=1); axC.text(0.42, 31.06, 'default (flag OFF) 31/32', fontsize=7.5, color='#117a65')
axC.axvspan(0.11, 0.52, color='#d9534f', alpha=0.09)
axC.axvline(0.11, color='#a02c24', ls='--', lw=1)
axC.text(0.30, 28.4, 'PERMANENT LOCK\n(bistable latch — loses\nPalbo-EZH2 + MB+HHi;\ndata-excluded)', fontsize=8,
         ha='center', color='#a02c24', style='italic')
axC.text(0.075, 29.4, 'reversible\ntransient regime', fontsize=8, ha='center', color=C_PROX, style='italic')
axC.set_xlabel('amplifier strength  a$_P$  (reservoir→proximal load)')
axC.set_ylabel('validation targets passed / 32')
axC.set_xlim(0.02, 0.52); axC.set_ylim(27.5, 31.4)
axC.set_title('(C) Validation cost — safe transient (30/32) vs excluded permanent lock')

# ---------- (D) dwell vs CDKI excursion: safe amplifier ~= OFF; only permanent locks (population_amplifier.py) ----------
axD = fig.add_subplot(gs[1, 1])
mags = [2.0, 3.0, 4.0, 5.0, 6.0]
dwell_off = [1.2, 10.3, 25.3, 147.8, 145.2]      # amplifier OFF  (intrinsic CDKI-recovery)
dwell_tr = [1.2, 10.3, 25.3, 147.8, 145.2]       # transient (C)  -- IDENTICAL to OFF
dwell_pm = [1.2, 10.3, 25.3, None, None]         # permanent (A)  -- locks deep excursions
axD.plot(mags, dwell_off, 'o-', color='#888', lw=3.2, ms=9, alpha=0.7, label='amplifier OFF')
axD.plot(mags, dwell_tr, 's--', color=C_PROX, lw=1.6, ms=6, label='transient (C) — ≡ OFF')
pm_x = [m for m, d in zip(mags, dwell_pm) if d is not None]; pm_y = [d for d in dwell_pm if d is not None]
axD.plot(pm_x, pm_y, '^', color='#a02c24', ms=8, label='permanent (A)')
for m, d in zip(mags, dwell_pm):
    if d is None:
        axD.annotate('LOCK', (m, 6), color='#a02c24', fontsize=8, ha='center', va='bottom', fontweight='bold')
axD.axhline(24, color='k', ls=':', lw=0.8); axD.text(2.05, 27, 'one cycle (~24 h)', fontsize=7.5)
axD.set_xlabel('excursion CDKI (× baseline)'); axD.set_ylabel('transient-G0 dwell after re-entry (h)')
axD.set_ylim(0, 160)
axD.set_title('(D) Dwell = intrinsic CDKI-recovery — safe amplifier adds nothing; only permanent locks')
axD.legend(fontsize=8, loc='upper left')

fig.suptitle('Two-compartment PRC2 / H3K27me3 — reservoir-sustained transient-G0 amplifier '
             '(with_proximal_distal, flag; default OFF)', fontsize=12.5, fontweight='bold', y=0.995)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_proximal_distal')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight')
fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
