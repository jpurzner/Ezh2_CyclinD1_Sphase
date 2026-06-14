"""v45 -- the FLIPPED analysis: from the measured folds to a falsifiable prediction.

Instead of assuming the feedback strength and reading off the dynamics, we ask: given the measured
steady-state folds, what feedback strength is implied, what does it PREDICT, and what experiment would
distinguish a strong- vs weak-feedback EZH2?

Isolated CyclinD1<->EZH2 module (v45 forms), mitogen drive ratio MB/GNP ~7.4 (condition_drive: GNP 2.56,
MB 19.0). Sweep loop strength kEZsyn and compute, at steady state:
  * CyclinD1 MB/GNP fold = cd_ss(M_MB)/cd_ss(M_GNP)   (data: 5.07x, Fig 4I)
  * EZH2 MB/GNP fold     = E_ss(M_MB)/E_ss(M_GNP)      (data: 2.05x, Fig 4J)
  * mitogen-ENTRY threshold shift = how much MORE mitogen is needed to reach a fixed commitment level
    CyclinD1=T_commit WITH the feedback vs WITHOUT (EZH2i) = the rShh-dose-to-cycle shift that Taz would reveal.

Key structural fact: EZH2 can only LOWER CyclinD1 (cd_fb = M*K/(K+EZH2) <= M = cd_noFB always), so EZH2
inhibition can NEVER give fewer divisions than the feedback -- i.e. "EZH2-high -> extra divisions" is
impossible; the only sign is EZH2i -> MORE proliferation (the paper's rescue). The falsifiable content is
therefore the THRESHOLD SHIFT (and the dose-response slope), pinned by the CyclinD1 fold.

Run:  ./venv/bin/python simulations/fig_v45_feedback_prediction.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt

Kcd, K = 2.0, 1.0
kEZbas, kdEZ = 0.08, 0.16
ksyn, kdeg = 0.30, 0.30
M_GNP, M_MB = 2.56, 19.0                 # condition_drive normalized drives (ratio ~7.4)
FOLD_CD_DATA, FOLD_EZ_DATA = 5.07, 2.05  # Fig 4I, Fig 4J


def cd_ss(M, kEZsyn, ezh2i=0):
    cd = M
    for _ in range(300):
        E = (kEZbas + kEZsyn * cd / (Kcd + cd)) / kdEZ
        cd_new = (ksyn / kdeg) * M * K / (K + E * (1 - ezh2i))
        if abs(cd_new - cd) < 1e-10:
            break
        cd = 0.5 * cd + 0.5 * cd_new
    return cd, E


def entry_mitogen(kEZsyn, T_commit, ezh2i=0):
    """smallest mitogen whose steady CyclinD1 reaches T_commit (the dose-to-cycle)."""
    Ms = np.linspace(0.05, 12, 1200)
    cd = np.array([cd_ss(M, kEZsyn, ezh2i)[0] for M in Ms])
    i = np.argmax(cd >= T_commit)
    return Ms[i] if cd[i] >= T_commit else np.nan


KS = np.linspace(0.0, 0.9, 60)
foldCd = np.array([cd_ss(M_MB, k)[0] / cd_ss(M_GNP, k)[0] for k in KS])
foldEz = np.array([cd_ss(M_MB, k)[1] / cd_ss(M_GNP, k)[1] for k in KS])

# pin the strength from the CyclinD1 fold (monotone decreasing)
k_star = np.interp(-FOLD_CD_DATA, -foldCd, KS)            # interp on decreasing foldCd
foldEz_at_star = np.interp(k_star, KS, foldEz)

# threshold shift: T_commit set so that WITHOUT feedback the dose-to-cycle is M_GNP (GNP just cycles)
T_commit = cd_ss(M_GNP, 0.0)[0] * 0.85
shift = np.array([entry_mitogen(k, T_commit) / entry_mitogen(0.0, T_commit) for k in KS])
shift_at_star = np.interp(k_star, KS, shift)

print(f'k* (pinned by CyclinD1 fold {FOLD_CD_DATA}) = {k_star:.3f}')
print(f'  predicted EZH2 fold there = {foldEz_at_star:.2f} (data {FOLD_EZ_DATA}) -> {"CONSISTENT" if abs(foldEz_at_star-FOLD_EZ_DATA)<0.4 else "tension"}')
print(f'  predicted mitogen-entry threshold shift (EZH2 vs EZH2i) = {shift_at_star:.2f}x')

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# (A) the two folds pin the same strength (self-consistency)
axA = axes[0]
axA.plot(KS, foldCd, color='#1B7837', lw=2.4, label='CyclinD1 MB/GNP fold')
axA.plot(KS, foldEz, color='#E08214', lw=2.4, label='EZH2 MB/GNP fold')
axA.axhline(FOLD_CD_DATA, color='#1B7837', ls=':', lw=1.4); axA.axhline(FOLD_EZ_DATA, color='#E08214', ls=':', lw=1.4)
axA.axvline(k_star, color='#E41A1C', ls='--', lw=1.8)
axA.text(k_star, 7.2, ' k* (data-pinned)', color='#E41A1C', fontsize=8, fontweight='bold')
axA.scatter([k_star, k_star], [FOLD_CD_DATA, foldEz_at_star], color='#E41A1C', zorder=5, s=35)
axA.set_xlabel('feedback strength (EZH2 synthesis gain kEZsyn)'); axA.set_ylabel('MB/GNP fold')
axA.set_title('A  The two measured folds pin one strength', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=8.5, loc='center right'); axA.set_ylim(0.8, 7.8)
axA.text(0.02, 6.6, 'weaker $\\to$', fontsize=7.5, color='#555'); axA.text(0.62, 1.2, '$\\leftarrow$ stronger (buffered)', fontsize=7.5, color='#555')

# (B) the PREDICTION: mitogen-entry threshold shift
axB = axes[1]
axB.plot(KS, shift, color='#762A83', lw=2.6)
axB.axvline(k_star, color='#E41A1C', ls='--', lw=1.8)
axB.scatter([k_star], [shift_at_star], color='#E41A1C', zorder=5, s=45)
axB.annotate(f'PREDICTION at data-pinned k*:\nEZH2 raises the mitogen\nthreshold to cycle by ~{shift_at_star:.1f}x\n(test: rShh dose-to-cycle $\\pm$ Taz)',
             xy=(k_star, shift_at_star), xytext=(k_star + 0.08, shift_at_star - 0.9),
             fontsize=8, color='#762A83', arrowprops=dict(arrowstyle='->', color='#762A83', lw=1.2))
axB.set_xlabel('feedback strength (EZH2 synthesis gain kEZsyn)'); axB.set_ylabel('mitogen-entry threshold shift (EZH2 / EZH2i)')
axB.set_title('B  Falsifiable PREDICTION: the threshold shift', loc='left', fontweight='bold', fontsize=11)
axB.set_ylim(0.9, shift.max() * 1.05)

# (C) the discriminating experiment: dose-response slope; + the structural impossibility
axC = axes[2]
Ms = np.linspace(0.3, 12, 200)
for k, c, lab in [(0.0, '#888888', 'no feedback (EZH2i)'),
                  (k_star, '#1B7837', f'data-pinned k* (weak; fold {FOLD_CD_DATA})'),
                  (0.8, '#B2182B', 'strong feedback (fold ~%.1f)' % (cd_ss(M_MB, 0.8)[0] / cd_ss(M_GNP, 0.8)[0]))]:
    cd = [cd_ss(M, k)[0] for M in Ms]
    axC.plot(Ms, cd, color=c, lw=2.2, label=lab)
axC.scatter([M_GNP, M_MB], [cd_ss(M_GNP, k_star)[0], cd_ss(M_MB, k_star)[0]], color='#1B7837', zorder=5, s=30)
axC.text(M_GNP, cd_ss(M_GNP, k_star)[0], ' GNP', fontsize=7.5); axC.text(M_MB, cd_ss(M_MB, k_star)[0], ' MB', fontsize=7.5)
axC.set_xlabel('mitogen (rShh / Gli-MYCN drive)'); axC.set_ylabel('steady-state CyclinD1')
axC.set_title('C  Discriminator: CyclinD1 dose-response slope', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=7.3, loc='upper left')
axC.text(0.97, 0.04, 'EZH2 only LOWERS CyclinD1\n$\\Rightarrow$ EZH2i always $\\geq$ proliferation\n(the rescue); no "EZH2 $\\to$ extra divisions"',
         transform=axC.transAxes, ha='right', fontsize=7, color='#444',
         bbox=dict(boxstyle='round', fc='#F5F5F5', ec='#bbb'))

fig.suptitle('v45 flipped analysis: the measured CyclinD1 (5.07x) + EZH2 (2.05x) folds pin a WEAK feedback, which PREDICTS the mitogen-entry threshold shift (a separable experiment)',
             fontsize=10, y=1.01)
fig.tight_layout()
fig.savefig('simulations/fig_v45_feedback_prediction.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v45_feedback_prediction.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_feedback_prediction.{png,pdf}')
