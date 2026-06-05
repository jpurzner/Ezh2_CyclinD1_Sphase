"""Re-fit the HH/MYCN->CyclinD1 module to the bulk RNA-seq folds (docs/PARAMETERIZATION.md A).

Fixes the Hedgehog SATURATION: Gli1 must rise ~6.9x GNP->MB (model previously ~1.1x), and CyclinD1
must track it (~7.6x), with MYCN as the HHi-resistant residual. Adds the Gli-pathway saturation
knobs (K_Ptch_Smo, K_Smo_Gli_switch, K_Gli_act_Gli1, Vmax_Gli1_tx) to the search so Gli1 can
de-saturate. Random search + Nelder-Mead on steady-state levels (the HH module is quasi-static).

MYCN_amplification is a context INPUT (set per condition from the data Mycn fold), not fit.
Run (background):  ./venv/bin/python simulations/v44_recalibrate_gli.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.optimize import minimize
from src.build_model_v44_heldt import build_model_v44

# context = Ptch1 dosage + MYCN amplification (from the data Mycn fold), +/- HHi
COND = {
    "GNP":     dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0,  GDC0449=0.0),
    "GNP+HHi": dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0,  GDC0449=1.0),
    "P7Ptch":  dict(SHH=0.5, Ptch1_copy_number=0.5, MYCN_amplification=1.38, GDC0449=0.0),
    "MB":      dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.86, GDC0449=0.0),
    "MB+HHi":  dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.86, GDC0449=1.0),
}
SPEC = ["Cd_mRNA", "MYCN", "Gli1"]

PSPEC = [  # name, lo, hi, default     (Gli-pathway de-saturation + CyclinD1 transcription)
    ("K_Ptch_Smo",       0.02, 0.6,  0.4),
    ("K_Smo_Gli_switch", 0.2,  1.5,  0.6),
    ("K_Gli_act_Gli1",   0.1,  2.0,  0.3),
    ("Vmax_Gli1_tx",     0.5,  4.0,  1.2),
    ("k_Cd_tx_basal",    0.05, 0.5,  0.3),
    ("k_Cd_tx_Gli_max",  2.0,  35.0, 4.0),
    ("K_Gli_act_CycD",   0.2,  4.0,  0.4),
    ("k_Cd_tx_MYCN",     2.0,  40.0, 15.0),
    ("K_MYCN_Cd",        0.5,  2.5,  1.5),
    ("n_MYCN_Cd",        2.0,  6.0,  3.0),
]
NAMES = [p[0] for p in PSPEC]
LO = np.array([p[1] for p in PSPEC]); HI = np.array([p[2] for p in PSPEC]); X0 = np.array([p[3] for p in PSPEC])

# (species, num, den, target, weight)
TARGETS = [
    ("Gli1",    "MB",     "GNP", 6.90, 3.0),
    ("Gli1",    "P7Ptch", "GNP", 1.50, 1.0),
    ("Gli1",    "GNP+HHi","GNP", 0.01, 1.5),
    ("Cd_mRNA", "MB",     "GNP", 7.58, 2.0),
    ("Cd_mRNA", "P7Ptch", "GNP", 1.90, 1.0),
    ("Cd_mRNA", "GNP+HHi","GNP", 0.14, 2.0),
    ("Cd_mRNA", "MB+HHi", "MB",  0.40, 1.5),
]


# load the model ONCE; all fitted params + context inputs are runtime-settable (5-10x faster)
_RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)


def measure(x):
    params = dict(zip(NAMES, x))
    out = {}
    for name, inp in COND.items():
        try:
            _RR.reset()
            for k, v in params.items():
                _RR[k] = v
            for k, v in inp.items():
                _RR[k] = v
            r = _RR.simulate(0, 2500, 5000, selections=["time"] + SPEC)
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
        if tgt < 0.02 and "Gli1" in sp:                 # >99% reduction
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
        print(f"  {sp:8} {a:>7}/{b:<4} = {o[a][sp]/o[b][sp]:7.3f}   target {tgt}")


if __name__ == "__main__":
    print("=" * 66)
    print("v44 HH DE-SATURATION recalibration (Gli1 6.9x, CyclinD1 7.6x MB/GNP)")
    report(X0, "baseline (current)")
    rng = np.random.default_rng(0)
    best_x, best_f = X0.copy(), objective(X0)
    N = 700
    print(f"\nrandom search ({N})...")
    for i in range(N):
        x = LO + (HI - LO) * rng.random(len(LO))
        f = objective(x)
        if f < best_f:
            best_f, best_x = f, x.copy()
            print(f"  [{i}] obj={f:.3f}")
    print("\nNelder-Mead refine...")
    res = minimize(objective, best_x, method="Nelder-Mead", options=dict(maxiter=900, xatol=1e-3, fatol=1e-3))
    if objective(res.x) < best_f:
        best_x = res.x
    report(best_x, "REFINED BEST")
    print("\nFINAL PARAMS:")
    for n, v in zip(NAMES, best_x):
        print(f'    "{n}": {v:.4g},')
