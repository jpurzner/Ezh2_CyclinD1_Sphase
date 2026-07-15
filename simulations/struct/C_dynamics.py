"""STRUCTURE ANALYSIS — step C: dynamical portraits of the core loops.

Reads the actual settled limit-cycle trajectory and draws phase portraits that expose each core loop's dynamical
TYPE: the (Cdc20,MPF) mitotic relaxation oscillator, the (E2f,pRb) R-point excursion, the (p27,CyclinE)
commitment, and the (Cd,EZH2)/(Cd,Mk) feedback. Plus a fine birth-p27 sweep showing the SHARP commitment
threshold (the fingerprint of the underlying bistable toggle inside the hybrid oscillator). A rigorous
saddle-node/Hopf continuation on reduced subsystems is the deep-dive; this is the tractable evidence of type.

Run: ./venv/bin/python simulations/struct/C_dynamics.py
"""
import os, sys, json
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import tellurium as te
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import (run, classify, PTCH1_MB, MYCN_AMP_MB, P16_MB, P18_MB, KSYP21_MB)

OUT = os.path.dirname(os.path.abspath(__file__))
SEL = ['time', 'MPF', 'Cdc20', 'E2f', 'pRb', 'Rbm', 'Rb', 'P21', 'Ce', 'Ca', 'Cd', 'EZH2', 'Mk', 'mass', 'Dna', 'Skp2']
rr = te.loada(build_model_v44())
rr['SHH'] = 0.5
rr.reset(); rr.simulate(0, 5000, 2)         # settle
r = rr.simulate(5000, 5000 + 3 * 1400, 6000, selections=SEL)   # ~3 cycles, fine
tr = {s: np.array(r[s]) for s in SEL}

# ---- birth-p27 sweep: commitment sharpness (MB condition) ----
p27s = np.linspace(0.6, 2.2, 17)
g0c = []
for p in p27s:
    res = run(shh=0.5, ptch1_cn=PTCH1_MB, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=float(p))
    cnt, _, _ = classify(res, 1.0)
    g0c.append(cnt.get('G0', np.nan))

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
fig.suptitle('v44 core loops — dynamical portraits from the limit cycle (+ the commitment threshold)', fontsize=14, fontweight='bold')

def portrait(a, x, y, xl, yl, title, c='#39598a'):
    a.plot(tr[x], tr[y], '-', color=c, lw=.9, alpha=.85)
    a.plot(tr[x][0], tr[y][0], 'o', color='#b34a43', ms=5)
    a.set_xlabel(xl); a.set_ylabel(yl); a.set_title(title, fontsize=11, fontweight='bold'); a.grid(alpha=.2)

portrait(ax[0, 0], 'Cdc20', 'MPF', 'Cdc20', 'MPF', '(N1) Mitotic relaxation oscillator\nMPF-Cdc20 — sharp spike + reset = relaxation', '#7a3b9a')
portrait(ax[0, 1], 'E2f', 'pRb', 'E2f (free)', 'pRb (hyper)', '(L1) Rb-E2f R-point excursion\nlow-E2f/low-pRb  ->  jump to high (commitment)', '#b34a43')
portrait(ax[0, 2], 'P21', 'Ce', 'p27 (P21)', 'CyclinE-CDK2 (Ce)', '(L2) Skp2-p27 commitment toggle\nhigh-p27/low-CDK2  <->  low-p27/high-CDK2', '#12756c')
portrait(ax[1, 0], 'Cd', 'EZH2', 'CyclinD1 (Cd)', 'EZH2', '(N2) EZH2 -| CyclinD1 feedback\ntight loop = monostable buffering, no hysteresis', '#9a6a12')
portrait(ax[1, 1], 'mass', 'pRb', 'mass (size)', 'pRb (commitment)', '(N5) growth-gated cycle\nmass ramps, resets at division (event)', '#2e7d32')

axc = ax[1, 2]
axc.plot(p27s, g0c, '-o', color='#b34a43', lw=2, ms=5)
axc.axvspan(1.3, 1.45, color='#ccc', alpha=.4)
axc.set_xlabel('birth-p27'); axc.set_ylabel('MB G0 count-fraction (%)')
axc.set_title('Commitment THRESHOLD is sharp (switch)\nfingerprint of the bistable Skp2-p27 toggle', fontsize=11, fontweight='bold')
axc.grid(alpha=.25)
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'C_dynamics.png'), dpi=140, bbox_inches='tight')
print("wrote C_dynamics.png")

# quantify the threshold sharpness (max slope) + whether EZH2-Cd loop shows hysteresis (loop area)
g0 = np.array(g0c)
dslope = np.max(np.abs(np.diff(g0)) / np.diff(p27s))
# EZH2-Cd loop "thickness": how much the (Cd,EZH2) trajectory is a thin curve (monostable) vs fat loop
cd = tr['Cd']; ez = tr['EZH2']
cdn = (cd - cd.min()) / (cd.ptp() + 1e-9); ezn = (ez - ez.min()) / (ez.ptp() + 1e-9)
loop_area = 0.5 * abs(np.sum(cdn[:-1] * ezn[1:] - cdn[1:] * ezn[:-1]))  # normalized enclosed area
print(f"commitment threshold max slope: {dslope:.1f} %/unit birth-p27 (sharp = switch)")
print(f"(Cd,EZH2) normalized loop area: {loop_area:.3f}  (small = thin curve = monostable buffering)")
json.dump(dict(birth_p27=list(p27s), mb_g0_count=g0c, threshold_max_slope=float(dslope),
               ezh2_cd_loop_area=float(loop_area)),
          open(os.path.join(OUT, 'C_dynamics.json'), 'w'), indent=2)
print("wrote C_dynamics.json")
