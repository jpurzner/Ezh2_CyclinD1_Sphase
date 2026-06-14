"""v45 -- does the RAPIDITY of mitogen withdrawal change the CyclinD1<->EZH2 response?

EZH2 is a SLOW, stable integrator; the mitogen change has its own timescale. Whether the feedback
buffers (helps) or overshoots (hurts) depends on the ratio. We cut mitogen from M_hi to an intermediate
M_lo (a MYCN-floor level) over a withdrawal time T_cut, sweeping T_cut from an instantaneous STEP to a
slow ramp, and integrate the CyclinD1<->EZH2 module (v45 forms), feedback ON vs OFF.

Prediction: a FAST cut outruns EZH2 -- EZH2 lags high -> transiently OVER-represses CyclinD1 (an
undershoot below its new steady level) -> CyclinD1 dips LOWER than no-feedback for a while. A SLOW cut
lets EZH2 track -> quasi-static buffering, no undershoot. So the feedback's sign of effect depends on
the cut rate.

Run:  ./venv/bin/python simulations/fig_v45_cutoff_rate.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

Kcd, K = 2.0, 1.0
kEZbas, kEZsyn, kdEZ = 0.08, 0.49, 0.16          # k* = data-pinned feedback strength (fig_v45_feedback_prediction)
ksyn, kdeg = 0.30, 0.30
M_HI, M_LO, T0 = 5.0, 1.0, 10.0                   # cut from 5 to an intermediate floor of 1


def cd_ss(M, ezh2i=0):
    cd = M
    for _ in range(300):
        E = (kEZbas + kEZsyn * cd / (Kcd + cd)) / kdEZ
        cd2 = (ksyn / kdeg) * M * K / (K + E * (1 - ezh2i))
        if abs(cd2 - cd) < 1e-10:
            break
        cd = 0.5 * cd + 0.5 * cd2
    return cd, E


def M_of_t(t, Tcut):
    if t < T0:
        return M_HI
    if Tcut <= 0:
        return M_LO
    f = min(1.0, (t - T0) / Tcut)
    return M_HI * (1 - f) + M_LO * f


def run(Tcut, ezh2i):
    cd0, E0 = cd_ss(M_HI, ezh2i)

    def dydt(y, t):
        E, cd = y
        M = M_of_t(t, Tcut)
        return [kEZbas + kEZsyn * cd / (Kcd + cd) - kdEZ * E, ksyn * M * K / (K + E * (1 - ezh2i)) - kdeg * cd]
    t = np.linspace(0, T0 + 120, 3000)
    sol = odeint(dydt, [E0, cd0], t)
    return t, sol[:, 1], sol[:, 0]              # time, CyclinD1, EZH2


CFB, CNF = '#762A83', '#1B7837'
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# (A) CyclinD1 time-courses for several cut rates (feedback ON)
axA = axes[0]
cd_lo_fb = cd_ss(M_LO, 0)[0]
for Tcut, c in [(0.0, '#40004B'), (4.0, '#762A83'), (20.0, '#9970AB'), (80.0, '#C2A5CF')]:
    t, cd, _ = run(Tcut, 0)
    axA.plot(t - T0, cd, color=c, lw=2, label=f'cut over {Tcut:.0f} h' if Tcut else 'instant step')
axA.axhline(cd_lo_fb, color='#888', ls=':', lw=1.2); axA.text(110, cd_lo_fb, ' new steady', va='bottom', ha='right', fontsize=7, color='#555')
axA.set_xlabel('time after withdrawal starts (h)'); axA.set_ylabel('CyclinD1 (feedback ON)')
axA.set_title('A  Faster cut $\\to$ transient UNDERSHOOT', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=8, loc='upper right'); axA.set_xlim(-8, 112)
# zoom inset on the undershoot + recovery (small, because the data-pinned feedback is weak)
axin = axA.inset_axes([0.42, 0.42, 0.5, 0.4])
for Tcut, c in [(0.0, '#40004B'), (4.0, '#762A83'), (20.0, '#9970AB')]:
    t, cd, _ = run(Tcut, 0); axin.plot(t - T0, cd, color=c, lw=1.6)
axin.axhline(cd_lo_fb, color='#888', ls=':', lw=1)
axin.set_xlim(-2, 45); axin.set_ylim(cd_lo_fb - 0.06, cd_lo_fb + 0.10)
axin.set_title('zoom: undershoot $\\sim$7%\n(recovers in $\\sim$6h = 1/kdEZ)', fontsize=6.5)
axin.tick_params(labelsize=6); axA.indicate_inset_zoom(axin, edgecolor='#999')

# (B) undershoot depth vs cut time, feedback vs no-feedback
axB = axes[1]
Tcuts = np.concatenate([[0], np.geomspace(0.5, 120, 26)])
for ezh2i, c, lab in [(0, CFB, 'EZH2 feedback'), (1, CNF, 'no feedback')]:
    cdf = cd_ss(M_LO, ezh2i)[0]
    und = []
    for Tc in Tcuts:
        _, cd, _ = run(Tc, ezh2i)
        und.append(100 * max(0.0, (cdf - cd.min())) / cd_ss(M_HI, ezh2i)[0])
    axB.plot(Tcuts, und, '-o', color=c, lw=2, ms=3, label=lab)
axB.set_xscale('symlog', linthresh=0.5)
axB.set_xlabel('mitogen cut duration T_cut (h)  [0 = instant]'); axB.set_ylabel('CyclinD1 undershoot below new steady (%)')
axB.set_title('B  Undershoot only with feedback + fast cut', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=8.5, loc='upper right')

# (C) consequence: integrated CyclinD1 above a cycling threshold (~ divisions) vs cut time
axC = axes[2]
thr = 0.45 * cd_ss(M_HI, 0)[0]                    # a nominal cycling threshold
for ezh2i, c, lab in [(0, CFB, 'EZH2 feedback'), (1, CNF, 'no feedback')]:
    integ = []
    for Tc in Tcuts:
        t, cd, _ = run(Tc, ezh2i)
        m = t >= T0
        integ.append(np.trapz(np.clip(cd[m] - thr, 0, None), t[m]))
    axC.plot(Tcuts, integ, '-o', color=c, lw=2, ms=3, label=lab)
axC.set_xscale('symlog', linthresh=0.5)
axC.set_xlabel('mitogen cut duration T_cut (h)  [0 = instant]'); axC.set_ylabel('integrated CyclinD1 over threshold ($\\sim$ divisions)')
axC.set_title('C  Slower cut $\\to$ more proliferation (both)', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=8.5, loc='upper left')
axC.annotate('feedback < no-feedback at EVERY rate\n(EZH2 only represses); the gap is\nLARGEST for a fast cut (lag adds extra arrest)',
             xy=(2, 0), xytext=(1.2, max(integ) * 0.18), fontsize=6.8, color='#555')

fig.suptitle('v45: rapidity of mitogen withdrawal -- a fast cut exposes the slow-EZH2 lag (transient over-repression / undershoot); a slow cut is quasi-static. Feedback never adds proliferation.',
             fontsize=9.8, y=1.01)
fig.tight_layout()
fig.savefig('simulations/fig_v45_cutoff_rate.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v45_cutoff_rate.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_cutoff_rate.{png,pdf}')
_, cd_step, _ = run(0.0, 0)
print(f'  instant step: CyclinD1 undershoot to {cd_step.min():.3f} (new steady {cd_ss(M_LO,0)[0]:.3f}); recovers as EZH2 falls')
print(f'  EZH2 timescale 1/kdEZ = {1/kdEZ:.1f} h sets the undershoot/recovery time')
