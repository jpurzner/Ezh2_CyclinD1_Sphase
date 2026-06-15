"""v45 Figure -- CyclinD1 transcript & protein under a LINEAR mitogen decline, with vs without the
EZH2 feedback.

The CyclinD1<->EZH2 negative-feedback module from v45 (Fig 4K), the same functional forms used in the
validated cell-cycle model:
  * EZH2 is a Rb-E2f target -> synthesis tracks the CyclinD1 PROTEIN (mitogen-dose, saturating):
        dEZH2/dt = kEZbas + kEZsyn*(CyclinD1_prot/(Kcd+CyclinD1_prot)) - kdEZ*EZH2     (STABLE, slow kdEZ)
  * EZH2 REPRESSES CyclinD1 TRANSCRIPTION (H3K27me3 at the promoter):
        CyclinD1_mRNA = mitogen * K/(K + EZH2*(1-EZH2i))
  * CyclinD1 PROTEIN follows the transcript with a synthesis/degradation lag:
        dCyclinD1_prot/dt = ksyn*CyclinD1_mRNA - kdeg*CyclinD1_prot

We ramp the mitogen (Gli/MYCN drive) DOWN linearly over time and integrate, for EZH2i=0 (feedback ON)
and EZH2i=1 (feedback OFF, EZH2 cannot repress -> CyclinD1_mRNA = mitogen). With the feedback, EZH2
tracks the falling mitogen and FALLS, relieving repression -> CyclinD1 transcript & protein are BUFFERED
(decline sub-linearly, maintained longer) -- homeostasis. Without it, CyclinD1 follows mitogen linearly.

Run:  ./venv/bin/python simulations/fig_v45_cyclind1_decline.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# module parameters (v45 forms; a.u.). Kcd, K mirror the v45 kernel (Kcd=2.0, K_EZH2_Cd=1.0).
Kcd, K = 2.0, 1.0
kEZbas, kEZsyn, kdEZ = 0.08, 0.24, 0.16           # EZH2: stable (slow kdEZ) -> integrates; synth saturates in CyclinD1;
                                                  # basal floor so EZH2 swings ~3x (1.5 high -> 0.5 low), matching ~2x MB/GNP
ksyn, kdeg = 0.30, 0.30                            # CyclinD1 protein from transcript (lags)
M0, T_RAMP = 5.0, 60.0                             # mitogen 5 (MB-like) -> 0 over 60 h
T_SETTLE = 24.0                                    # hold at M0 first to reach steady cycling


def mitogen(t):
    if t < T_SETTLE:
        return M0
    return M0 * max(0.0, 1.0 - (t - T_SETTLE) / T_RAMP)


def deriv(y, t, EZH2i):
    EZH2, Cd_p = y
    mRNA = mitogen(t) * K / (K + EZH2 * (1 - EZH2i))
    dEZH2 = kEZbas + kEZsyn * (Cd_p / (Kcd + Cd_p)) - kdEZ * EZH2
    dCd_p = ksyn * mRNA - kdeg * Cd_p
    return [dEZH2, dCd_p]


t = np.linspace(0, T_SETTLE + T_RAMP + 25, 2000)
out = {}
for fb, ezh2i in [('feedback', 0), ('no_feedback', 1)]:
    # initialize at the steady state for M0
    y0 = [2.0, 4.0] if fb == 'feedback' else [2.0, 5.0]
    for _ in range(3):                             # settle the initial condition
        y0 = odeint(deriv, y0, np.linspace(0, T_SETTLE, 200), args=(ezh2i,))[-1]
    sol = odeint(deriv, y0, t, args=(ezh2i,))
    EZH2 = sol[:, 0]; Cd_p = sol[:, 1]
    mRNA = np.array([mitogen(tt) for tt in t]) * K / (K + EZH2 * (1 - ezh2i))
    out[fb] = dict(EZH2=EZH2, mRNA=mRNA, prot=Cd_p)

M = np.array([mitogen(tt) for tt in t])
CFB, CNF, CM, CE = '#762A83', '#1B7837', '#888888', '#E08214'
m0 = (t >= T_SETTLE - 1) & (t < T_SETTLE)         # window for the "high mitogen" normalization

fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4))

# (A) the inputs: mitogen (linear) + EZH2 (tracks it)
axA = axes[0, 0]
axA.plot(t, M, color=CM, lw=2, ls=':', label='mitogen (Gli/MYCN drive)')
axA.plot(t, out['feedback']['EZH2'], color=CE, lw=2.2, label='EZH2 (feedback ON, tracks mitogen)')
axA.plot(t, out['no_feedback']['EZH2'], color=CE, lw=1.4, ls='--', alpha=0.7, label='EZH2 (feedback OFF)')
axA.axvspan(T_SETTLE, T_SETTLE + T_RAMP, color='#F2EEF6', alpha=0.8)
axA.text(T_SETTLE + T_RAMP / 2, M0 * 0.95, 'mitogen withdrawal', ha='center', fontsize=8, color='#555')
axA.set_xlabel('time (h)'); axA.set_ylabel('level (a.u.)')
axA.set_title('A  Inputs: mitogen $\\downarrow$, EZH2 falls with it', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=8, loc='upper right')

# (B) CyclinD1 TRANSCRIPT (absolute)
axB = axes[0, 1]
axB.plot(t, M, color=CM, lw=1.3, ls=':', label='mitogen (no-repression reference)')
axB.plot(t, out['no_feedback']['mRNA'], color=CNF, lw=2.2, label='no feedback (= mitogen)')
axB.plot(t, out['feedback']['mRNA'], color=CFB, lw=2.2, label='EZH2 feedback')
axB.axvspan(T_SETTLE, T_SETTLE + T_RAMP, color='#F2EEF6', alpha=0.8)
axB.set_xlabel('time (h)'); axB.set_ylabel('CyclinD1 transcript (a.u.)')
axB.set_title('B  CyclinD1 transcript (absolute)', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=7.5, loc='upper right')
axB.annotate('feedback REPRESSES\n(lower absolute level)', xy=(T_SETTLE - 2, out['feedback']['mRNA'][m0].mean()),
             xytext=(2, M0 * 0.30), fontsize=7.2, color=CFB, arrowprops=dict(arrowstyle='->', color=CFB, lw=1))

# (C) CyclinD1 PROTEIN (absolute, lags transcript)
axC = axes[1, 0]
axC.plot(t, out['no_feedback']['prot'], color=CNF, lw=2.2, label='no feedback')
axC.plot(t, out['feedback']['prot'], color=CFB, lw=2.2, label='EZH2 feedback')
axC.axvspan(T_SETTLE, T_SETTLE + T_RAMP, color='#F2EEF6', alpha=0.8)
axC.set_xlabel('time (h)'); axC.set_ylabel('CyclinD1 protein (a.u.)')
axC.set_title('C  CyclinD1 protein (lags transcript)', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=7.5, loc='upper right')

# (D) NORMALIZED to the high-mitogen value -> the BUFFERING (rate of decline)
axD = axes[1, 1]
mref = M[m0].mean()
axD.plot(t, M / mref, color=CM, lw=1.6, ls=':', label='mitogen (linear)')
for sp, c, lab in [('no_feedback', CNF, 'no feedback'), ('feedback', CFB, 'EZH2 feedback')]:
    axD.plot(t, out[sp]['mRNA'] / out[sp]['mRNA'][m0].mean(), color=c, lw=2.2, label=f'CyclinD1 mRNA, {lab}')
axD.axvspan(T_SETTLE, T_SETTLE + T_RAMP, color='#F2EEF6', alpha=0.8)
axD.set_xlabel('time (h)'); axD.set_ylabel('CyclinD1 transcript (norm. to high mitogen)')
axD.set_title('D  Normalized: feedback BUFFERS the decline', loc='left', fontweight='bold', fontsize=11)
axD.legend(fontsize=7.5, loc='upper right'); axD.set_ylim(-0.03, 1.08)
axD.annotate('EZH2 falls $\\to$ de-represses\n$\\to$ transcript declines SUB-linearly',
             xy=(T_SETTLE + T_RAMP * 0.6, np.interp(T_SETTLE + T_RAMP * 0.6, t, out['feedback']['mRNA'] / out['feedback']['mRNA'][m0].mean())),
             xytext=(2, 0.30), fontsize=7.2, color=CFB, arrowprops=dict(arrowstyle='->', color=CFB, lw=1))

fig.suptitle('v45 CyclinD1$\\dashv$EZH2 feedback under linear mitogen decline: the feedback REPRESSES CyclinD1 (lower) AND, because EZH2 tracks the falling mitogen, BUFFERS its decline (sub-linear)',
             fontsize=10, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig('simulations/fig_v45_cyclind1_decline.png', dpi=170, bbox_inches='tight')
fig.savefig('simulations/fig_v45_cyclind1_decline.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_cyclind1_decline.{png,pdf}')
ez = out['feedback']['EZH2']
print(f'  EZH2 swing (feedback): {ez[m0].mean():.2f} (high mitogen) -> {ez[-1]:.2f} (low)  ~{ez[m0].mean()/max(ez[-1],1e-6):.1f}x')
# half-decline times (normalized) -> the buffering
for sp, lab in [('no_feedback', 'no feedback'), ('feedback', 'feedback')]:
    nm = out[sp]['mRNA'] / out[sp]['mRNA'][m0].mean()
    thalf = t[np.argmax(nm < 0.5)] - T_SETTLE if np.any(nm < 0.5) else np.nan
    print(f'  CyclinD1 transcript ({lab}): falls to 50% at t-T_settle = {thalf:.0f}h into withdrawal')
