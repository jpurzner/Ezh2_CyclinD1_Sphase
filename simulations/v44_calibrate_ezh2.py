"""v44 EZH2 calibration: within-cycle gradient (transcript + protein) and HU S-boost.

The model carries BOTH EZH2 transcript (EZH2m) and protein (EZH2), so we score:
  - EZH2m (transcript) gradient vs G0  -> Section D (scRNA-seq, MB): G1/S ~1.8-2.5, G2/M ~1.9-2.5
  - EZH2  (protein)    gradient vs G0  -> Section E (IF, Supp.6D):   G2/G0 = 1.48x
  - HU 10uM -> EZH2(protein) in S / DMSO -> Section E (Fig.4G):       1.31x

Phases by Rb-phosphorylation (G0/G1) + replication state (S/G2), MB context.
Edit EZH2 params via PARAMS and re-run.  ./venv/bin/python simulations/v44_calibrate_ezh2.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

PARAMS = {
    "kEZbas": 0.0003, "kEZE2f": 0.010, "K_E2f_EZ": 0.3,
    "K_Ce_EZ": 0.5, "K_Ca_EZ": 0.8, "wCe": 0.5,
    "kDeEZm": 0.02, "kTlEZ": 0.004, "kDeEZ": 0.0003,
}
SEL = ["time", "pRb", "aRc", "Dna", "EZH2", "EZH2m"]


def run(params, hu, t_end=16000, n_pts=64000):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=params))
    rr['HU'] = hu; rr['MYCN_amplification'] = 2.8
    return rr.simulate(0, t_end, n_pts, selections=SEL)


def by_phase(res, var, settle=6000, prb_frac=0.5):
    t = res["time"]; m = t >= settle
    pRb = res["pRb"][m]; aRc = res["aRc"][m]; Dna = res["Dna"][m]; v = res[var][m]
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = (Dna >= 0.98); preS = ~in_S & ~in_G2
    thr = prb_frac * pRb.max()
    G0 = preS & (pRb < thr); G1 = preS & (pRb >= thr)
    e = lambda s: float(v[s].mean()) if s.any() else np.nan
    return dict(G0=e(G0), G1=e(G1), S=e(in_S), G2=e(in_G2))


if __name__ == "__main__":
    print("=" * 72)
    print("v44 EZH2 CALIBRATION (MB)")
    print("PARAMS:", PARAMS)
    d = run(PARAMS, 0.0); h = run(PARAMS, 1.0)
    for var, lab, tgt in [("EZH2m", "transcript", "D: G1/S~1.8-2.5, G2/M~1.9-2.5"),
                          ("EZH2", "protein", "E: G2/G0=1.48")]:
        p = by_phase(d, var)
        print(f"\n{lab} ({var}) gradient vs G0  [{tgt}]")
        print(f"   G1/G0={p['G1']/p['G0']:.2f}  S/G0={p['S']/p['G0']:.2f}  "
              f"G2/G0={p['G2']/p['G0']:.2f}")
    ps_d = by_phase(d, "EZH2"); ps_h = by_phase(h, "EZH2")
    print(f"\nHU 10uM -> EZH2(protein) in S / DMSO = {ps_h['S']/ps_d['S']:.2f}  (target 1.31x)")
