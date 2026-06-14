"""v44 -- mitogen sensitivity: how CyclinD1, EZH2 and proliferation respond to mitogen (SHH/Hh) dose,
with vs without the EZH2 -> CyclinD1 feedback. Generated on the committed Heldt-engine model.

Sweep SHH (the physiological mitogen, ~ rShh in Fig 4H) in a GNP context (Ptch1_copy_number=1), with the
EZH2->CyclinD1 feedback ON (EZH2i=0) vs OFF (EZH2i=1). Read the settled CyclinD1 transcript (Cd_mRNA) and
protein (Cd), EZH2, and the number of divisions in 168 h.

  (A) CyclinD1 transcript & protein vs mitogen, +/- feedback -- the dose-response. The feedback REPRESSES
      CyclinD1 and lowers its mitogen SENSITIVITY (flatter / buffered), consistent with the v45 module analysis.
  (B) EZH2 vs mitogen -- EZH2 tracks Hh/CyclinD1 dose (Rb-E2f target; Fig 4H), the source of the feedback.
  (C) proliferation (divisions / 168 h) vs mitogen, +/- feedback -- the feedback RAISES the mitogen
      threshold to enter cycle (the mitogen sensitivity of proliferation).

Run:  ./venv/bin/python simulations/fig_v44_mitogen_sensitivity.py    (~10 min; 2 x SHH-sweep of v44 runs)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from simulations.validate_v44 import run as run44, count_divisions, mean_settled

SHH = np.array([0.0, 0.08, 0.15, 0.22, 0.3, 0.4, 0.5, 0.65, 0.85, 1.1, 1.4])
T_END, N_PTS = 10080, 16000          # 168 h
WINDOW_H = (T_END - 3000) / 60.0      # settled window for divisions/period


def sweep(ezh2i):
    cd, mrna, ez, div = [], [], [], []
    for shh in SHH:
        r = run44(shh=float(shh), ptch1_cn=1.0, gdc=0.0, ezh2i=ezh2i, mycn_amp=1.0,
                  t_end=T_END, n_pts=N_PTS)
        cd.append(mean_settled(r, 'Cd')); mrna.append(mean_settled(r, 'Cd_mRNA'))
        ez.append(mean_settled(r, 'EZH2'))
        n, _, _ = count_divisions(r); div.append(n / WINDOW_H * 168.0)
        print(f'  EZH2i={ezh2i} SHH={shh:.2f}: Cd={cd[-1]:.2f} mRNA={mrna[-1]:.2f} EZH2={ez[-1]:.2f} div/168h={div[-1]:.1f}')
    return dict(cd=np.array(cd), mrna=np.array(mrna), ez=np.array(ez), div=np.array(div))


_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mitogen_sensitivity_cache.npz')
if os.path.exists(_CACHE) and '--fresh' not in sys.argv:
    z = np.load(_CACHE)
    FB = {k[3:]: z[k] for k in z.files if k.startswith('fb_')}
    NF = {k[3:]: z[k] for k in z.files if k.startswith('nf_')}
    print(f'loaded cached sweep ({_CACHE}); pass --fresh to recompute')
else:
    print('sweeping SHH WITH feedback (EZH2i=0) ...');  FB = sweep(0.0)
    print('sweeping SHH WITHOUT feedback (EZH2i=1) ...'); NF = sweep(1.0)
    np.savez(_CACHE, **{f'fb_{k}': v for k, v in FB.items()}, **{f'nf_{k}': v for k, v in NF.items()})

CFB, CNF = '#762A83', '#1B7837'
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# (A) CyclinD1 transcript & protein vs mitogen
axA = axes[0]
axA.plot(SHH, NF['mrna'], '--', color=CNF, lw=1.8, label='transcript, no feedback')
axA.plot(SHH, FB['mrna'], '-', color=CFB, lw=2.2, label='transcript, EZH2 feedback')
axA.plot(SHH, NF['cd'], '--o', color=CNF, lw=1.4, ms=3, alpha=0.7, label='protein, no feedback')
axA.plot(SHH, FB['cd'], '-o', color=CFB, lw=1.4, ms=3, alpha=0.7, label='protein, EZH2 feedback')
axA.set_xlabel('mitogen (SHH / Hh dose)'); axA.set_ylabel('CyclinD1 (transcript / protein, a.u.)')
axA.set_title('A  CyclinD1 dose-response', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=7.3, loc='upper left')
axA.annotate('feedback represses &\nflattens (lower mitogen\nsensitivity)', xy=(1.1, FB['mrna'][-2]),
             xytext=(0.45, FB['mrna'].max() * 0.45), fontsize=7.3, color=CFB,
             arrowprops=dict(arrowstyle='->', color=CFB, lw=1))

# (B) EZH2 vs mitogen
axB = axes[1]
axB.plot(SHH, FB['ez'], '-o', color='#E08214', lw=2.2, ms=4, label='EZH2 (feedback ON)')
axB.plot(SHH, NF['ez'], '--o', color='#E08214', lw=1.4, ms=4, alpha=0.6, label='EZH2 (feedback OFF)')
axB.set_xlabel('mitogen (SHH / Hh dose)'); axB.set_ylabel('EZH2 (a.u.)')
axB.set_title('B  EZH2 tracks mitogen (Fig 4H)', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=8, loc='lower right')

# (C) proliferation vs mitogen -> threshold shift
axC = axes[2]
axC.plot(SHH, NF['div'], '--s', color=CNF, lw=2, ms=5, label='no feedback (EZH2i)')
axC.plot(SHH, FB['div'], '-o', color=CFB, lw=2.2, ms=5, label='EZH2 feedback')
# 50%-of-max thresholds
def thr(d):
    half = 0.5 * np.nanmax(d)
    above = np.where(d >= half)[0]
    return np.interp(half, d[:above[0] + 1], SHH[:above[0] + 1]) if len(above) and above[0] > 0 else np.nan
tFB, tNF = thr(FB['div']), thr(NF['div'])
for t, c in [(tNF, CNF), (tFB, CFB)]:
    if not np.isnan(t):
        axC.axvline(t, color=c, ls=':', lw=1.2, alpha=0.7)
axC.set_xlabel('mitogen (SHH / Hh dose)'); axC.set_ylabel('proliferation (divisions / 168 h)')
axC.set_title('C  Mitogen threshold to cycle', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=8, loc='center right'); axC.set_ylim(-0.6, np.nanmax(NF['div']) * 1.15)
# no-feedback cycles even at SHH=0 (basal CyclinD1 suffices) -> threshold <= 0; feedback installs a finite threshold
nf_cycles_at_0 = NF['div'][0] > 0.5 * np.nanmax(NF['div'])
note = (f'WITHOUT feedback: GNP cycles even at SHH=0\n(basal CyclinD1 suffices). The feedback\nINSTALLS the Hh threshold (SHH~{tFB:.2f}).'
        if nf_cycles_at_0 else f'feedback raises threshold {tNF:.2f}$\\to${tFB:.2f}')
axC.axvline(tFB, color=CFB, ls=':', lw=1.2, alpha=0.7)
axC.annotate(note, xy=(tFB, 0.5 * np.nanmax(FB['div'])), xytext=(tFB + 0.12, 0.30 * np.nanmax(FB['div'])),
             fontsize=7.0, color=CFB, arrowprops=dict(arrowstyle='->', color=CFB, lw=1))

fig.suptitle('v44 (Heldt-engine model): mitogen sensitivity -- the EZH2$\\dashv$CyclinD1 feedback represses/buffers CyclinD1 and raises the Hh threshold to enter cycle',
             fontsize=10.5, y=1.01)
fig.tight_layout()
fig.savefig('simulations/fig_v44_mitogen_sensitivity.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v44_mitogen_sensitivity.pdf', bbox_inches='tight')
print('wrote simulations/fig_v44_mitogen_sensitivity.{png,pdf}')
print(f'  mitogen threshold to cycle: no-feedback SHH~{tNF:.2f}, feedback SHH~{tFB:.2f}')
