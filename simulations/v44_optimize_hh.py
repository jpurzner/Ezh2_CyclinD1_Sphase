"""Automated fit of the HH/MYCN -> CyclinD1 parameters to the Section A between-condition
ratios. Random search (robust to CVODE crashes / quiescence) + Nelder-Mead refinement.

Fits CyclinD1 transcript (Cd_mRNA, matches RNA-seq), MYCN, Gli1 ratios across:
  GNP (Ptch1=1, MYCN_amp=1), GNP+HHi, MB (Ptch1=0.1, MYCN_amp=2.8), MB+HHi.

Run (background):  ./venv/bin/python simulations/v44_optimize_hh.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from scipy.optimize import minimize
from src.build_model_v44_heldt import build_model_v44

COND = {
    "GNP":     dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=0.0),
    "GNP+HHi": dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=1.0),
    "MB":      dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=0.0),
    "MB+HHi":  dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=1.0),
}
SPEC = ["Cd_mRNA", "MYCN", "Gli1"]

# parameter name -> (low, high, default).  Fit the HH/MYCN -> CyclinD1 block.
PSPEC = [
    ("k_Cd_tx_basal",   0.05, 0.6,  0.3),
    ("k_Cd_tx_Gli_max", 2.0,  25.0, 4.0),
    ("K_Gli_act_CycD",  0.2,  3.0,  0.4),
    ("k_Cd_tx_MYCN",    4.0,  40.0, 15.0),
    ("K_MYCN_Cd",       0.5,  2.0,  1.5),
    ("n_MYCN_Cd",       2.0,  6.0,  3.0),
    ("K_Ptch_Smo",      0.05, 0.6,  0.4),
]
NAMES = [p[0] for p in PSPEC]
LO = np.array([p[1] for p in PSPEC]); HI = np.array([p[2] for p in PSPEC])
X0 = np.array([p[3] for p in PSPEC])

# (species, num_cond, den_cond, target, weight)
TARGETS = [
    ("Cd_mRNA", "GNP+HHi", "GNP", 0.14, 2.0),
    ("Cd_mRNA", "MB+HHi",  "MB",  0.40, 2.0),
    ("Cd_mRNA", "MB",      "GNP", 5.07, 3.0),
    ("MYCN",    "GNP+HHi", "GNP", 0.78, 1.0),
    ("MYCN",    "MB+HHi",  "MB",  0.86, 1.0),
    ("MYCN",    "MB",      "GNP", 2.80, 1.0),
    ("Gli1",    "GNP+HHi", "GNP", 0.01, 1.0),
]


def measure(x):
    params = dict(zip(NAMES, x))
    try:
        m = build_model_v44(with_ezh2=True, with_hh=True, params=params)
    except Exception:
        return None
    out = {}
    for name, inp in COND.items():
        try:
            rr = te.loada(m)
            for k, v in inp.items():
                rr[k] = v
            rr.integrator.setValue("absolute_tolerance", 1e-9)
            r = rr.simulate(0, 2500, 5000, selections=["time"] + SPEC)
            sel = r["time"] >= 1500
            out[name] = {s: max(float(np.mean(r[s][sel])), 1e-6) for s in SPEC}
        except Exception:
            return None
    return out


def objective(x):
    if np.any(x < LO) or np.any(x > HI):
        return 1e6
    o = measure(x)
    if o is None:
        return 1e5
    err = 0.0
    for sp, a, b, tgt, w in TARGETS:
        ratio = o[a][sp] / o[b][sp]
        if "Gli1" in sp and tgt < 0.02:                 # >99% reduction: penalize if not small
            err += w * (max(ratio - 0.02, 0.0) * 20) ** 2
        else:
            err += w * (np.log(ratio / tgt)) ** 2
    return err


def report(x, label=""):
    o = measure(x)
    print(f"\n--- {label}  obj={objective(x):.3f} ---")
    print("params: " + ", ".join(f"{n}={v:.3g}" for n, v in zip(NAMES, x)))
    if o is None:
        print("  (failed)"); return
    for sp, a, b, tgt, w in TARGETS:
        print(f"  {sp:8} {a:>7}/{b:<4} = {o[a][sp]/o[b][sp]:6.3f}   target {tgt}")


if __name__ == "__main__":
    print("=" * 64)
    print("v44 HH/MYCN -> CyclinD1 OPTIMIZER (random search + Nelder-Mead)")
    report(X0, "baseline (v42 defaults)")

    rng = np.random.default_rng(0)
    best_x, best_f = X0.copy(), objective(X0)
    N = 400
    print(f"\nrandom search ({N} samples)...")
    for i in range(N):
        x = LO + (HI - LO) * rng.random(len(LO))
        f = objective(x)
        if f < best_f:
            best_f, best_x = f, x.copy()
            print(f"  [{i}] obj={f:.3f}")
    report(best_x, "best random")

    print("\nNelder-Mead refine...")
    res = minimize(objective, best_x, method="Nelder-Mead",
                   options=dict(maxiter=600, xatol=1e-3, fatol=1e-3))
    report(res.x, "REFINED BEST")
    print("\nFINAL PARAMS:")
    for n, v in zip(NAMES, res.x):
        print(f'    "{n}": {v:.4g},')
