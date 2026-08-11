"""MB transient G0 via the CDK4/6-alone Rb-hyper-P escape (Yang 2020) — high CyclinD1/CDK6 pulls MB cells out
of the two-step-Rb permanent lock (JP 2026-08-09). Flag with_cd_hyper_escape (default OFF).

The two-step Rb makes a deep MB arrest PERMANENT (CyclinD-CDK4/6 only mono-Ps Rb; hyper-P needs CyclinE/A
which needs E2f = deadlock). Letting CyclinD-CDK4/6 WEAKLY hyper-P Rb (Yang 2020: CDK4/6 alone commits)
lets the huge MB CyclinD1/CDK6 pull the cell back out -> reversible TRANSIENT G0.

  (A) proliferation vs CDKI: OFF collapses to permanent arrest at 6x; escape keeps deep-CDKI cells cycling.
  (B) time-course at fixed 6x CDKI: OFF = permanent (no division); escape = MULTI-PERIODIC transient G0 --
      CyclinD1 accumulates in the dwell, pulls the cell out, it divides, re-arrests (the reversible dwell).
  (C) G0 dwell (period) grows with arrest depth -- deeper CDKI = longer transient G0.
  (D) validation vs escape strength: 31/32 up to w_cd_hyper~0.2; >0.2 breaks the GNP-SHH arrest (self-bound).

Run:  ./venv/bin/python simulations/fig_v44_cdk46_escape.py   (~4 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)
C_OFF = '#888888'; C_ESC = '#1f6f4a'; C_CD = '#b0402f'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})


def make(wch):
    P = dict(BP)
    if wch is not None:
        P['w_cd_hyper'] = wch
    rr = te.loada(build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_cd_hyper_escape=(wch is not None), params=P))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
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


def prolif(wch, mults, T=60000):
    rr, P = make(wch); out = []
    for m in mults:
        setp(rr, P, m); rr.reset()
        r = rr.simulate(0, T, T // 10, selections=['time', 'Dna'])
        d = np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]
        per = float(np.mean(np.diff(r['time'][d])) / 60) if len(d) >= 2 else None
        out.append((len(d), per))
    return out

fig = plt.figure(figsize=(15, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.28)
mults = [1, 2, 3, 4, 5, 6, 7, 8]

# (A) proliferation (divisions/60000min) vs CDKI
axA = fig.add_subplot(gs[0, 0])
for wch, col, lab in [(None, C_OFF, 'escape OFF'), (0.05, '#7fb069', 'escape 0.05'), (0.2, C_ESC, 'escape 0.20')]:
    res = prolif(wch, mults)
    axA.plot(mults, [r[0] for r in res], 'o-', color=col, lw=1.7, ms=6, label=lab)
axA.set_xlabel('CDKI (× baseline)'); axA.set_ylabel('divisions / 60000 min')
axA.set_title('(A) Deep-CDKI arrest: OFF → permanent (0); escape keeps cycling')
axA.legend(fontsize=8.5); axA.axhline(0, color='k', lw=0.6)
axA.annotate('OFF: permanent\narrest at 6×', (6, 2), fontsize=8, color=C_OFF, ha='center')

# (B) time-course at 6x CDKI: OFF permanent vs escape multi-periodic transient G0
axB = fig.add_subplot(gs[0, 1])
for wch, col, lab in [(None, C_OFF, 'escape OFF (permanent)'), (0.2, C_ESC, 'escape 0.20 (transient G0)')]:
    rr, P = make(wch); setp(rr, P, 6.0); rr.reset()
    r = rr.simulate(0, 55000, 27500, selections=['time', 'Dna', 'Cd'])
    axB.plot(r['time'] / 60, r['Cd'], color=col, lw=1.4, label=lab)
    dd = np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]
    for x in r['time'][dd] / 60:
        axB.axvline(x, color=col, ls=':', lw=0.8, alpha=0.5)
axB.set_xlabel('time (h)'); axB.set_ylabel('CyclinD1 (Cd)')
axB.set_title('(B) At fixed 6× CDKI — escape → CyclinD1 pulls cell out, divides, re-arrests')
axB.legend(fontsize=8.5, loc='center right')
axB.annotate('G0 dwell\n(CyclinD1 rebuilds)', (18, 10), fontsize=8, color=C_ESC, ha='center')

# (C) G0 dwell (period) vs CDKI depth (escape 0.2)
axC = fig.add_subplot(gs[1, 0])
res = prolif(0.2, mults)
xs = [m for m, r in zip(mults, res) if r[1] is not None]
ys = [r[1] for r in res if r[1] is not None]
axC.plot(xs, ys, 'o-', color=C_ESC, lw=1.8, ms=7)
axC.axhline(24, color='k', ls=':', lw=0.8); axC.text(1.2, 26, 'normal cycle ~24 h', fontsize=8)
axC.set_xlabel('CDKI (× baseline)'); axC.set_ylabel('cycle period = G0 dwell + cycle (h)')
axC.set_title('(C) Deeper arrest → longer transient-G0 dwell (drive-graded)')

# (D) validation vs escape strength (from escape_eval workflow)
axD = fig.add_subplot(gs[1, 1])
w = [0.0, 0.03, 0.05, 0.10, 0.20, 0.40]; passed = [31, 31, 31, 31, 31, 30]
axD.plot(w, passed, 'o-', color=C_ESC, lw=1.8, ms=7)
axD.axhline(31, color='#117a65', ls='--', lw=1); axD.text(0.01, 31.05, 'full baseline 31/32', fontsize=8, color='#117a65')
axD.axvspan(0.25, 0.42, color='#d9534f', alpha=0.09)
axD.text(0.33, 30.4, 'too strong:\nbreaks GNP-SHH\narrest (self-bound)', fontsize=8, ha='center', color='#a02c24', style='italic')
axD.set_xlabel('escape strength  w_cd_hyper'); axD.set_ylabel('validation passed / 32')
axD.set_ylim(29.6, 31.4); axD.set_xlim(-0.02, 0.42)
axD.set_title('(D) Validation preserved (31/32) up to w_cd_hyper ~0.2')

fig.suptitle('MB transient G0 via CDK4/6-alone escape — high CyclinD1/CDK6 pulls cells out of the '
             'two-step-Rb lock (with_cd_hyper_escape, flag; default OFF)', fontsize=12, fontweight='bold', y=0.995)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_cdk46_escape')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
