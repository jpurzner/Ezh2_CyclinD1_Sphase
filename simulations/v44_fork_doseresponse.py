"""v44 dose-response: does reducing fork speed (kSyDna = HU target) lengthen S and the
cycle CONCENTRATION-DEPENDENTLY in the cycling model, with mitosis held until Dna=1?

This is the property the whole rebuild exists to get (GG could not). Vary kSyDna,
measure per-cycle G1 / S / total-period from a stable cycle, confirm: (a) S lengthens
as forks slow, (b) mitosis (MPF spike) occurs only after Dna~1 (checkpoint holds),
(c) cells keep dividing (no arrest).

Run:  ./venv/bin/python simulations/v44_fork_doseresponse.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
import tellurium as te
from v44_build_test import build  # reuse the augmented model

BASE_KSYDNA = 0.0093


def run(mult, t_end=8000, n_pts=32000):
    rr = te.loada(build())
    rr.kSyDna = BASE_KSYDNA * mult
    res = rr.simulate(0, t_end, n_pts, selections=["time", "Dna", "MPF", "E2f", "Ca", "Cb"])
    return res


def analyze(res, settle=2500):
    t = res["time"]; dna = res["Dna"]; mpf = res["MPF"]
    mit, _ = find_peaks(mpf, prominence=0.2, distance=50)
    mit = mit[t[mit] > settle]
    if len(mit) < 3:
        return None
    per = np.mean(np.diff(t[mit]))
    # division = Dna reset edge (>0.9 then <0.1); cycle starts at the low (post-reset) sample
    div = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0] + 1
    div = div[t[div] > settle]
    g1s, ss, g2s = [], [], []
    for k in range(len(div) - 1):
        d0, d1 = div[k], div[k + 1]
        seg_d = dna[d0:d1]; seg_t = t[d0:d1]
        a05 = np.where(seg_d > 0.05)[0]; a95 = np.where(seg_d > 0.95)[0]
        if len(a05) and len(a95):
            ts0 = seg_t[a05[0]]; ts1 = seg_t[a95[0]]
            g1s.append(ts0 - seg_t[0]); ss.append(ts1 - ts0); g2s.append(seg_t[-1] - ts1)
    return dict(period=per, n_mit=len(mit),
                g1=np.mean(g1s) if g1s else np.nan,
                s=np.mean(ss) if ss else np.nan,
                g2=np.mean(g2s) if g2s else np.nan,
                dna_at_mit=float(np.mean(dna[mit])))


if __name__ == "__main__":
    print("=" * 78)
    print("v44 FORK-SPEED DOSE RESPONSE (kSyDna = replication fork speed = HU target)")
    print("=" * 78)
    print(f"{'fork x':>8}{'kSyDna':>9}{'period_h':>10}{'G1_h':>8}{'S_h':>8}{'G2M_h':>8}"
          f"{'mitoses':>9}{'Dna@mit':>9}")
    for mult in (2.0, 1.0, 0.6, 0.4, 0.25):
        res = run(mult)
        a = analyze(res)
        if a is None:
            print(f"{mult:>8.2f}{BASE_KSYDNA*mult:>9.5f}   no stable cycle / arrest")
            continue
        print(f"{mult:>8.2f}{BASE_KSYDNA*mult:>9.5f}{a['period']/60:>10.2f}"
              f"{a['g1']/60:>8.2f}{a['s']/60:>8.2f}{a['g2']/60:>8.2f}{a['n_mit']:>9d}"
              f"{a['dna_at_mit']:>9.2f}")
    print("\nExpect: S_h grows as forks slow; G1_h ~flat; mitoses>0 (no arrest); Dna@mit~1.0")
