"""v44 phase-proportion + EZH2-gradient metrics vs experimental targets (MB context).

Phase classification from a single cycling trajectory (time-fraction = asynchronous
population proportion). Per the manuscript, G0 vs G1 is delineated by Rb phosphorylation:
  G0  : pre-replication, Rb HYPO-phosphorylated (pRb < thr)   -> quiescent-like
  G1  : pre-replication, Rb phosphorylated (pRb >= thr)        -> committed, pre-S
  S   : origins fired & replicating (aRc high, Dna < 0.98)
  G2M : replication complete (Dna >= 0.98) until division
The pRb threshold is anchored so DMSO G0 matches the (renormalized) experimental value;
the HU-induced G0/G1/S/G2 redistribution is then a prediction.

Targets (Section F, MB; renormalized to the classified fraction so phases sum to 100%):
  DMSO  G0/G1/S/G2 = 24.8 / 43.4 / 15.7 / 16.1 %   (from 21.7/37.9/13.7/14.1, sum 87.4)
  HU folds vs DMSO: S x1.36 (up), G2 x0.23 (down), G0 x1.19 (up), G1 x1.16
EZH2 gradient (Sections D/E): cycling/G0 ~2.0 (transcript); G2/G0 protein 1.48x;
  HU 10uM -> EZH2 in S / DMSO = 1.31x.

MB = MYCN_amplification ~2.8. HU=1.0 := 10 uM (calibratable via KmHU_fork).
Run:  ./venv/bin/python simulations/v44_phase_metrics.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

M = build_model_v44(with_ezh2=True, with_hh=True)
SEL = ["time", "pRb", "E2f", "aRc", "Dna", "MPF", "EZH2", "Ce", "Ca"]


def run(hu, mycn_amp=2.8, hhi=0.0, t_end=14000, n_pts=56000):
    rr = te.loada(M)
    rr['HU'] = hu; rr['MYCN_amplification'] = mycn_amp; rr['HHi'] = hhi
    return rr.simulate(0, t_end, n_pts, selections=SEL)


def classify(res, pRb_thr, settle=5000):
    t = res["time"]; m = t >= settle
    pRb = res["pRb"][m]; aRc = res["aRc"][m]; Dna = res["Dna"][m]
    in_S = (aRc > 0.05) & (Dna < 0.98)
    in_G2M = (Dna >= 0.98)
    preS = ~in_S & ~in_G2M
    G0 = preS & (pRb < pRb_thr)
    G1 = preS & (pRb >= pRb_thr)
    frac = dict(G0=G0.mean(), G1=G1.mean(), S=in_S.mean(), G2=in_G2M.mean())
    s = sum(frac.values())
    frac = {k: 100 * v / s for k, v in frac.items()}   # renormalize to 100%
    ez = res["EZH2"][m]
    ezph = {ph: float(ez[sel].mean()) if sel.any() else np.nan
            for ph, sel in [("G0", G0), ("G1", G1), ("S", in_S), ("G2", in_G2M)]}
    return frac, ezph


def period(res, settle=5000):
    t = res["time"]; mpf = res["MPF"]; m = t >= settle
    pk, _ = find_peaks(mpf[m], prominence=0.2, distance=50)
    return float(np.mean(np.diff(t[m][pk]))) / 60 if len(pk) > 1 else np.nan


# Anchor pRb threshold so DMSO G0 ~ 24.8% (target)
dmso = run(0.0)
prb = dmso["pRb"][dmso["time"] >= 5000]
# search threshold
best = None
for thr in np.linspace(prb.min() + 0.05, prb.max() - 0.05, 60):
    f, _ = classify(dmso, thr)
    d = abs(f["G0"] - 24.8)
    if best is None or d < best[0]:
        best = (d, thr)
pRb_thr = best[1]

TARGET_DMSO = dict(G0=24.8, G1=43.4, S=15.7, G2=16.1)
TARGET_HU_FOLD = dict(G0=1.19, G1=1.16, S=1.36, G2=0.23)

f0, ez0 = classify(dmso, pRb_thr)
hu = run(1.0)
f1, ez1 = classify(hu, pRb_thr)

print("=" * 76)
print(f"v44 PHASE METRICS vs targets (MB; pRb G0-thr={pRb_thr:.2f}; "
      f"period DMSO={period(dmso):.1f}h)")
print("=" * 76)
print(f"\n{'phase':>6} | {'model DMSO':>11}{'target':>8} | {'model HU':>10}"
      f"{'mdl fold':>9}{'tgt fold':>9}")
for ph in ("G0", "G1", "S", "G2"):
    mfold = f1[ph] / f0[ph] if f0[ph] else np.nan
    print(f"{ph:>6} | {f0[ph]:>10.1f}%{TARGET_DMSO[ph]:>7.1f}% | {f1[ph]:>9.1f}%"
          f"{mfold:>9.2f}{TARGET_HU_FOLD[ph]:>9.2f}")

print(f"\nEZH2 by phase (ratio vs G0):  target cycling/G0~2.0, G2/G0 protein 1.48x")
print(f"  DMSO:  G1/G0={ez0['G1']/ez0['G0']:.2f}  S/G0={ez0['S']/ez0['G0']:.2f}  "
      f"G2/G0={ez0['G2']/ez0['G0']:.2f}")
print(f"  HU 10uM: EZH2 in S / DMSO EZH2 in S = {ez1['S']/ez0['S']:.2f}  (target 1.31x)")
