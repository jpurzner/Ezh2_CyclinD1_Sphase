"""v45 -- PARAMETER-SPACE analysis of the CyclinD1<->EZH2 negative feedback.

The question is not "do we see buffering / extra divisions / reduced mitogen sensitivity" but UNDER WHAT
PARAMETERS, HOW BROADLY, and WITH WHAT CONSEQUENCE. We isolate the CyclinD1<->EZH2 module (the v45
functional forms) and sweep its two governing dimensions:
  * LOOP STRENGTH  = EZH2 synthesis gain kEZsyn (how hard EZH2 tracks/represses CyclinD1).
  * LOOP DELAY     = EZH2 turnover kdEZ (EZH2 stability; slow EZH2 = long delay relative to mitogen change).
(K, Kcd, ksyn/kdeg fixed at the v45-matched module values.)

Module (steady state solved by fixed-point; dynamics by integration):
  EZH2_ss   = (kEZbas + kEZsyn*Cd_p/(Kcd+Cd_p)) / kdEZ
  CyclinD1_mRNA = M * K/(K + EZH2*(1-EZH2i)) ;  dCd_p/dt = ksyn*mRNA - kdeg*Cd_p

Outputs mapped over (kEZsyn, kdEZ):
  (A) steady-state mitogen DOSE-RESPONSE for a few loop strengths -> mitogen SENSITIVITY (the slope).
  (B) heatmap of mitogen sensitivity S = mean d ln(CyclinD1)/d ln(mitogen): 1 = linear, <1 = buffered/homeostatic.
  (C) heatmap of OVERSHOOT after a mitogen step (the delayed-negative-feedback consequence: ringing).
  (D) representative time courses at 3 regimes (v45 default; strong+slow = oscillatory; strong+fast = buffered).

Run:  ./venv/bin/python simulations/fig_v45_feedback_paramspace.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.integrate import odeint

Kcd, K = 2.0, 1.0
kEZbas = 0.08
ksyn, kdeg = 0.30, 0.30
DEF_kEZsyn, DEF_kdEZ = 0.24, 0.16          # the v45-matched default (EZH2 swing ~2.5x, repression ~2.5x)


def cd_ss(M, kEZsyn, kdEZ, ezh2i):
    """Steady-state CyclinD1 protein at mitogen M (fixed-point of the feedback)."""
    cd = M
    for _ in range(200):
        E = (kEZbas + kEZsyn * cd / (Kcd + cd)) / kdEZ
        cd_new = (ksyn / kdeg) * M * K / (K + E * (1 - ezh2i))
        if abs(cd_new - cd) < 1e-9:
            break
        cd = 0.5 * cd + 0.5 * cd_new
    return cd, E


def sensitivity(kEZsyn, kdEZ, Mgrid):
    """mean d ln(CyclinD1)/d ln(mitogen) over Mgrid (1=linear, <1=buffered)."""
    cd = np.array([cd_ss(M, kEZsyn, kdEZ, 0)[0] for M in Mgrid])
    lnM, lncd = np.log(Mgrid), np.log(cd)
    return np.mean(np.diff(lncd) / np.diff(lnM))


def overshoot(kEZsyn, kdEZ, Mhi=5.0, Mlo=2.0):
    """ring/overshoot after a step DOWN in mitogen: max |Cd_p - Cd_final| beyond a monotone approach,
    as a fraction of the step. Slow strong feedback -> Cd_p undershoots then recovers (damped oscillation)."""
    cd0, E0 = cd_ss(Mhi, kEZsyn, kdEZ, 0)
    cdf, _ = cd_ss(Mlo, kEZsyn, kdEZ, 0)

    def dydt(y, t):
        E, cd = y
        mRNA = Mlo * K / (K + E)
        return [kEZbas + kEZsyn * cd / (Kcd + cd) - kdEZ * E, ksyn * mRNA - kdeg * cd]
    t = np.linspace(0, 80, 1600)
    sol = odeint(dydt, [E0, cd0], t)
    cd = sol[:, 1]
    under = (cdf - cd.min()) / (cd0 - cdf + 1e-9)        # how far it dips BELOW the final value
    return max(0.0, under)


# ---- grids ----
Mgrid = np.linspace(0.6, 5.0, 30)
KS = np.linspace(0.0, 0.7, 36)                            # loop strength (0 = no feedback)
KD = np.linspace(0.05, 0.6, 34)                           # EZH2 turnover (small = slow = long delay)
Sens = np.array([[sensitivity(ks, kd, Mgrid) for ks in KS] for kd in KD])
Over = np.array([[overshoot(ks, kd) for ks in KS] for kd in KD])
print('computing grids done')

fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# (A) dose-response for a few loop strengths
axA = axes[0, 0]
for ks, c, lab in [(0.0, '#888888', 'no feedback'), (0.12, '#9970AB', 'weak'),
                   (DEF_kEZsyn, '#762A83', 'v45 default'), (0.5, '#40004B', 'strong')]:
    cd = [cd_ss(M, ks, DEF_kdEZ, 0)[0] for M in Mgrid]
    axA.plot(Mgrid, cd, color=c, lw=2.2, label=f'{lab} (kEZsyn={ks})')
axA.set_xlabel('mitogen (CyclinD1 drive)'); axA.set_ylabel('steady-state CyclinD1 (a.u.)')
axA.set_title('A  Mitogen dose-response vs loop strength', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=8, loc='upper left')
axA.annotate('stronger feedback\n$\\to$ flatter $\\to$ LOWER\nmitogen sensitivity\n(homeostasis)',
             xy=(4.2, cd_ss(4.2, 0.5, DEF_kdEZ, 0)[0]), xytext=(2.4, 1.0), fontsize=7.5, color='#40004B',
             arrowprops=dict(arrowstyle='->', color='#40004B', lw=1))

# (B) mitogen sensitivity heatmap
axB = axes[0, 1]
im = axB.pcolormesh(KS, 1 / KD, Sens, shading='auto', cmap='viridis', vmin=0.3, vmax=1.0)
cs = axB.contour(KS, 1 / KD, Sens, levels=[0.5, 0.7, 0.85], colors='white', linewidths=0.8)
axB.clabel(cs, fmt='%.2f', fontsize=7)
# DATA CONSTRAINT: CyclinD1 MB/GNP = 5.07x at drive ratio ~7.4 -> sensitivity S ~ ln(5.07)/ln(7.4) ~ 0.81.
# Strong feedback (S far below this) would FLATTEN the mitogen response and violate the measured fold.
S_data = np.log(5.07) / np.log(7.4)
cd_line = axB.contour(KS, 1 / KD, Sens, levels=[S_data], colors='#E41A1C', linewidths=2.0, linestyles='--')
axB.clabel(cd_line, fmt=lambda v: 'data fold 5.07x', fontsize=7.5, colors='#E41A1C')
axB.plot(DEF_kEZsyn, 1 / DEF_kdEZ, 'o', color='#E41A1C', ms=11, mec='white', mew=1.5)
axB.text(DEF_kEZsyn, 1 / DEF_kdEZ, '  v45', color='#E41A1C', fontsize=8, fontweight='bold', va='center')
axB.set_xlabel('loop strength (EZH2 synthesis gain kEZsyn)'); axB.set_ylabel('EZH2 timescale 1/kdEZ (delay)')
axB.set_title('B  Mitogen sensitivity  (1=linear, <1=buffered)', loc='left', fontweight='bold', fontsize=11)
fig.colorbar(im, ax=axB, label='d ln CyclinD1 / d ln mitogen')

# (C) overshoot (delayed-feedback ringing) heatmap
axC = axes[1, 0]
im2 = axC.pcolormesh(KS, 1 / KD, Over, shading='auto', cmap='magma')
cs2 = axC.contour(KS, 1 / KD, Over, levels=[0.02, 0.08, 0.15], colors='cyan', linewidths=0.8)
axC.clabel(cs2, fmt='%.2f', fontsize=7)
axC.plot(DEF_kEZsyn, 1 / DEF_kdEZ, 'o', color='#4DAF4A', ms=11, mec='white', mew=1.5)
axC.text(DEF_kEZsyn, 1 / DEF_kdEZ, '  v45', color='#4DAF4A', fontsize=8, fontweight='bold', va='center')
axC.set_xlabel('loop strength (EZH2 synthesis gain kEZsyn)'); axC.set_ylabel('EZH2 timescale 1/kdEZ (delay)')
axC.set_title('C  Overshoot after mitogen step (ringing)', loc='left', fontweight='bold', fontsize=11)
fig.colorbar(im2, ax=axC, label='undershoot fraction')

# (D) representative time courses
axD = axes[1, 1]
Mhi, Mlo = 5.0, 2.0
regimes = [('v45 default', DEF_kEZsyn, DEF_kdEZ, '#762A83'),
           ('strong + slow EZH2\n(oscillatory)', 0.6, 0.06, '#B2182B'),
           ('strong + fast EZH2\n(buffered, no ring)', 0.6, 0.5, '#1B7837')]
for lab, ks, kd, c in regimes:
    cd0, E0 = cd_ss(Mhi, ks, kd, 0); cdf, _ = cd_ss(Mlo, ks, kd, 0)

    def dydt(y, t, ks=ks, kd=kd):
        E, cd = y
        return [kEZbas + ks * cd / (Kcd + cd) - kd * E, ksyn * Mlo * K / (K + E) - kdeg * cd]
    tt = np.linspace(0, 80, 1600); sol = odeint(dydt, [E0, cd0], tt)
    axD.plot(tt, sol[:, 1] / sol[0, 1], color=c, lw=2, label=lab)
axD.axhline(cd_ss(Mlo, DEF_kEZsyn, DEF_kdEZ, 0)[0] / cd_ss(Mhi, DEF_kEZsyn, DEF_kdEZ, 0)[0],
            color='#888', ls=':', lw=1)
axD.set_xlabel('time after mitogen step-down (a.u.)'); axD.set_ylabel('CyclinD1 (norm. to pre-step)')
axD.set_title('D  Consequence: dynamics by regime', loc='left', fontweight='bold', fontsize=11)
axD.legend(fontsize=7.5, loc='upper right')

fig.suptitle('v45 CyclinD1$\\dashv$EZH2 feedback across parameter space: mitogen sensitivity (homeostasis) and dynamic ringing as functions of loop STRENGTH and DELAY',
             fontsize=10.5, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig('simulations/fig_v45_feedback_paramspace.png', dpi=160, bbox_inches='tight')
fig.savefig('simulations/fig_v45_feedback_paramspace.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_feedback_paramspace.{png,pdf}')
print(f'  v45 default: sensitivity={sensitivity(DEF_kEZsyn,DEF_kdEZ,Mgrid):.2f} (1=linear), overshoot={overshoot(DEF_kEZsyn,DEF_kdEZ):.3f}')
print(f'  no feedback: sensitivity={sensitivity(0.0,DEF_kdEZ,Mgrid):.2f}')
print(f'  fraction of swept space buffered (S<0.85): {100*np.mean(Sens<0.85):.0f}%; with ringing (overshoot>0.05): {100*np.mean(Over>0.05):.0f}%')
