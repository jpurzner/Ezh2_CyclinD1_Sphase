"""v44 phase-balance calibration loop. Edit PARAMS, re-run, read the table.

Goal: a G1-dominated ~22 h MB cycle matching Section F proportions.
Coarse targets (renormalized to 100%):  pre-S(G0+G1)=68%, S=16%, G2=16%, period~22h.
HU 10uM redistribution: S x1.36 up, G2 x0.23 down. EZH2: G2/G0~1.5-2.0, HU S-boost 1.31x.

Run:  ./venv/bin/python simulations/v44_calibrate_phases.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

# ---- calibration candidate (edit me) ----
PARAMS = {
    # G1 length: slow Rb phosphorylation (dominant phase)
    "kPhRbCd": 0.2, "kPhRbCe": 0.3, "kPhRbCa": 0.3,
    # cyclin synthesis (G1->S drive)
    "kSyCe": 0.01, "kSyCa": 0.02,
    # S length: replication fork speed
    "kSyDna": 0.0093,
    # G2 length: CycB synthesis / mitotic entry
    "kSyCb": 0.01,
}

SEL = ["time", "pRb", "aRc", "Dna", "MPF", "EZH2"]


def run(params, hu, t_end=20000, n_pts=60000):
    m = build_model_v44(with_ezh2=True, with_hh=True, params=params)
    rr = te.loada(m)
    rr['HU'] = hu; rr['MYCN_amplification'] = 2.8
    return rr.simulate(0, t_end, n_pts, selections=SEL)


def metrics(res, settle=7000, prb_frac=0.5):
    t = res["time"]; m = t >= settle
    pRb = res["pRb"][m]; aRc = res["aRc"][m]; Dna = res["Dna"][m]; ez = res["EZH2"][m]
    in_S = (aRc > 0.05) & (Dna < 0.98)
    in_G2 = (Dna >= 0.98)
    preS = ~in_S & ~in_G2
    thr = prb_frac * pRb.max()
    G0 = preS & (pRb < thr); G1 = preS & (pRb >= thr)
    tot = len(pRb)
    fr = dict(preS=100*preS.mean(), S=100*in_S.mean(), G2=100*in_G2.mean(),
              G0=100*G0.mean(), G1=100*G1.mean())
    def emean(sel): return float(ez[sel].mean()) if sel.any() else np.nan
    ez_ph = dict(G0=emean(G0), G1=emean(G1), S=emean(in_S), G2=emean(in_G2))
    mpf = res["MPF"][m]
    pk, _ = find_peaks(mpf, prominence=0.2, distance=50)
    per = float(np.mean(np.diff(t[m][pk]))) / 60 if len(pk) > 1 else np.nan
    return fr, ez_ph, per, len(pk)


if __name__ == "__main__":
    print("=" * 78)
    print("v44 PHASE-BALANCE CALIBRATION (MB)   targets: preS 68% | S 16% | G2 16% | ~22h")
    print("=" * 78)
    print("PARAMS:", PARAMS)
    d, ezd, perd, nd = metrics(run(PARAMS, 0.0))
    h, ezh, perh, nh = metrics(run(PARAMS, 1.0))
    print(f"\nperiod DMSO={perd:.1f}h ({nd} mitoses)   HU={perh:.1f}h ({nh} mitoses)   target ~22h")
    print(f"\n{'':8}{'preS':>8}{'S':>8}{'G2':>8}  |{'(G0':>7}{'G1)':>7}")
    print(f"{'DMSO':8}{d['preS']:>7.1f}%{d['S']:>7.1f}%{d['G2']:>7.1f}%  |"
          f"{d['G0']:>6.1f}%{d['G1']:>6.1f}%")
    print(f"{'target':8}{68:>7.0f}%{16:>7.0f}%{16:>7.0f}%  |{24.8:>6.1f}%{43.4:>6.1f}%")
    print(f"{'HU':8}{h['preS']:>7.1f}%{h['S']:>7.1f}%{h['G2']:>7.1f}%")
    print(f"\nHU folds:  S x{h['S']/d['S']:.2f} (tgt 1.36)   G2 x{h['G2']/d['G2']:.2f} (tgt 0.23)"
          f"   preS x{h['preS']/d['preS']:.2f}")
    print(f"\nEZH2 gradient (DMSO):  G1/G0={ezd['G1']/ezd['G0']:.2f}  S/G0={ezd['S']/ezd['G0']:.2f}"
          f"  G2/G0={ezd['G2']/ezd['G0']:.2f}  (tgt cycling/G0~2, G2/G0 1.48-2.0)")
    print(f"EZH2 HU S-boost: {ezh['S']/ezd['S']:.2f}  (tgt 1.31)")
