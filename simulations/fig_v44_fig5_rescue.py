"""Figure 5A-C (main text) + Supp Fig 8 I,J  from the v44 model.

Fig 5: predicted cell-cycle oscillations (Cyclin B) in
  (A) Ptch+/- MB                -> sustained cycling
  (B) MB + HHi (vismodegib)     -> arrest
  (C) MB + HHi + EZH2i          -> rescue (cycling restored by CyclinD1 de-repression)
Supp Fig 8 I,J: MB + CDK4/6i (arrest) and MB + CDK4/6i + EZH2i (NO rescue) -- the
  downstream-block control showing EZH2i rescues upstream (HHi) but not downstream (CDK4/6i).

Uses the calibrated parameter set from validate_v44.PARAMS. Real-time minutes -> hours.
Run:  ./venv/bin/python simulations/fig_v44_fig5_rescue.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from simulations.validate_v44 import run, count_divisions, P16_MB, P18_MB, KSYP21_MB

# MB CDK-inhibitor brake (INK4 p16+p18 + CIP/KIP p21/p27), applied to every MB condition below
_MB = dict(p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB)

T_END = 10080      # 168 h in minutes
N_PTS = 20160

RESCUE = {
    'Ptch+/- MB':            dict(shh=0.5, ptch1_cn=0.1, hhi=0.0, ezh2i=0.0, mycn_amp=2.8, **_MB),
    'MB + HHi':              dict(shh=0.5, ptch1_cn=0.1, hhi=1.0, ezh2i=0.0, mycn_amp=2.8, **_MB),
    'MB + HHi + EZH2i':      dict(shh=0.5, ptch1_cn=0.1, hhi=1.0, ezh2i=1.0, mycn_amp=2.8, **_MB),
    'MB + CDK4/6i':          dict(shh=0.5, ptch1_cn=0.1, hhi=0.0, ezh2i=0.0, mycn_amp=2.8, cdk46i=True, **_MB),
    'MB + CDK4/6i + EZH2i':  dict(shh=0.5, ptch1_cn=0.1, hhi=0.0, ezh2i=1.0, mycn_amp=2.8, cdk46i=True, **_MB),
}
COLORS = ['#762A83', '#E08214', '#4A1486', '#C2185B', '#D81B60']

sims = {n: run(**kw, t_end=T_END, n_pts=N_PTS) for n, kw in RESCUE.items()}
names = list(RESCUE)
CB_YMAX = max(sims[n]['Cb'].max() for n in names) * 1.08    # shared y-scale (arrest reads as flat 0)

# ---- Fig 5 A-C: the three main-text rescue panels ----
fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))
for i, name in enumerate(names[:3]):
    ax = axes[i]; r = sims[name]
    n, _, per = count_divisions(r)
    ax.plot(r['time'] / 60.0, r['Cb'], color=COLORS[i], linewidth=1.0)
    p = f", {np.mean(per):.1f} h" if len(per) else ""
    ax.set_title(f"{name}\n({n} div{p})", fontsize=10, fontweight='bold', color=COLORS[i])
    ax.set_xlabel('Time (h)', fontsize=9)
    ax.set_ylabel('Cyclin B (a.u.)' if i == 0 else '', fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 168); ax.set_ylim(-0.05, CB_YMAX)
    if 'EZH2i' in name:
        ax.patch.set_facecolor('#e8f5e9'); ax.patch.set_alpha(0.3)
fig.suptitle('Figure 5A-C: EZH2 inhibition rescues HHi-arrested Ptch+/- MB (v44 model)',
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('simulations/fig_v44_fig5_rescue.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_v44_fig5_rescue.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_fig5_rescue.png/pdf")

# ---- Supp Fig 8 I,J: CDK4/6i (no rescue) control ----
fig2, axes2 = plt.subplots(1, 2, figsize=(9, 3.6))
for j, name in enumerate(names[3:]):
    ax = axes2[j]; r = sims[name]
    n, _, per = count_divisions(r)
    ax.plot(r['time'] / 60.0, r['Cb'], color=COLORS[3 + j], linewidth=1.0)
    p = f", {np.mean(per):.1f} h" if len(per) else ""
    ax.set_title(f"{name}\n({n} div{p})", fontsize=10, fontweight='bold', color=COLORS[3 + j])
    ax.set_xlabel('Time (h)', fontsize=9)
    ax.set_ylabel('Cyclin B (a.u.)' if j == 0 else '', fontsize=9)
    ax.grid(True, alpha=0.2); ax.set_xlim(0, 168); ax.set_ylim(-0.05, CB_YMAX)
    ax.patch.set_facecolor('#ffebee'); ax.patch.set_alpha(0.3)
fig2.suptitle('Supp Fig 8 I,J: CDK4/6 inhibition cannot be rescued by EZH2i (v44 model)',
              fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig('simulations/fig_v44_cdk46i_norescue.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_v44_cdk46i_norescue.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_cdk46i_norescue.png")

print("\nDivision summary:")
for name in names:
    n, _, per = count_divisions(sims[name])
    print(f"  {name:24s}: {n:2d} div" + (f", period {np.mean(per):.1f} h" if len(per) else ", arrest"))
