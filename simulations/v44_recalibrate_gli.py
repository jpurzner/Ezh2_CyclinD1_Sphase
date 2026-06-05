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
SPEC = ["Cd_mRNA", "MYCN", "Gli1", "Ptch1_mRNA"]

PSPEC = [  # name, lo, hi, default(=currently-baked HH params; seed NM here, not from scratch)
    ("K_Ptch_Smo",       0.02, 0.6,  0.1935),
    ("K_Smo_Gli_switch", 0.2,  1.5,  1.448),
    ("K_Gli_act_Gli1",   0.1,  2.0,  0.5395),
    ("Vmax_Gli1_tx",     0.5,  4.0,  0.5837),
    ("k_Cd_tx_basal",    0.02, 0.5,  0.137),
    ("k_Cd_tx_Gli_max",  2.0,  70.0, 70.0),
    ("K_Gli_act_CycD",   0.1,  4.0,  0.6246),
    ("k_Cd_tx_MYCN",     1.0,  40.0, 14.48),
    ("K_MYCN_Cd",        0.5,  2.5,  1.968),
    ("n_MYCN_Cd",        2.0,  6.0,  2.431),
    ("k_Smo_act",        0.5,  4.0,  1.5),     # Smo activation gain
    ("k_Ptch1_basal",    0.05, 1.0,  0.3),     # Gli->Ptch1 negative feedback: basal floor
    ("k_Ptch1_Gli",      0.1,  3.0,  0.6),     #   Gli-induced Ptch1 (closes the loop)
    ("K_Gli_Ptch",       0.1,  2.0,  0.5),     #   half-max Gli_act for Ptch1 induction
]
NAMES = [p[0] for p in PSPEC]
LO = np.array([p[1] for p in PSPEC]); HI = np.array([p[2] for p in PSPEC]); X0 = np.array([p[3] for p in PSPEC])

# (species, num, den, target, weight)   -- folds from docs/PARAMETERIZATION.md bulk RNA-seq
# KEY UPDATE: MB_GDC0449 row shows vismo crashes MB CyclinD1 by 85.6% (0.144), not 0.40.
# So Gli must DOMINATE MB CyclinD1 (~86%); Mycn/basal is only the ~14% HHi-resistant residual.
# Ptch1 is a Gli TARGET (negative feedback): elevated in MB (co-regulated with Gli1), dropping under HHi.
TARGETS = [
    ("Gli1",      "MB",     "GNP", 6.90, 3.0),
    ("Gli1",      "P7Ptch", "GNP", 1.50, 1.0),
    ("Gli1",      "GNP+HHi","GNP", 0.01, 1.5),
    ("Gli1",      "MB+HHi", "MB",  0.023, 2.0),   # vismo crashes Gli1 to ~2% of MB
    ("Cd_mRNA",   "MB",     "GNP", 7.58, 2.0),
    ("Cd_mRNA",   "P7Ptch", "GNP", 1.90, 1.0),
    ("Cd_mRNA",   "GNP+HHi","GNP", 0.157, 2.0),
    ("Cd_mRNA",   "MB+HHi", "MB",  0.144, 3.5),   # <-- the new vismo-on-MB CyclinD1 drop
    ("Ptch1_mRNA","MB",     "GNP", 5.0,  0.5),    # soft: Ptch1 co-regulated w/ Gli1 (no exact bulk value)
    ("Ptch1_mRNA","MB+HHi", "MB",  0.10, 0.5),    # soft: vismo collapses Gli -> Ptch1 transcription too
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
        if tgt < 0.05 and "Gli1" in sp:                 # near-total Gli block: one-sided
            err += w * (max(ratio - tgt, 0.0) * 20) ** 2  # only penalize if ABOVE target
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
    BEST_PATH = os.path.join(os.path.dirname(__file__), "v44_recalibrate_gli_best.txt")

    def save_best(x):  # write incrementally so a kill still leaves the current best
        lines = [f'    "{n}": {v:.4g},' for n, v in zip(NAMES, x)]
        with open(BEST_PATH, "w") as fh:
            fh.write(f"obj={objective(x):.4f}\n" + "\n".join(lines) + "\n")

    report(X0, "baseline (currently-baked params)")
    rng = np.random.default_rng(0)
    best_x, best_f = X0.copy(), objective(X0)
    save_best(best_x)

    # (1) refine straight from the baked X0 -- the known-good basin
    print("\nNelder-Mead from baked X0...")
    res0 = minimize(objective, X0, method="Nelder-Mead", options=dict(maxiter=1500, xatol=1e-4, fatol=1e-4))
    if objective(res0.x) < best_f:
        best_f, best_x = objective(res0.x), res0.x.copy()
        save_best(best_x); print(f"  -> obj={best_f:.3f}")

    # (2a) GLOBAL random search (full bounds) -- the Gli->Ptch1 feedback needs a different basin
    #      (stronger feedback: lower GNP Gli, bigger MB/GNP contrast, higher MB Ptch1 mRNA)
    print("\nglobal random search (full bounds)...")
    for i in range(400):
        x = LO + (HI - LO) * rng.random(len(LO))
        f = objective(x)
        if f < best_f:
            best_f, best_x = f, x.copy()
            save_best(best_x); print(f"  [g{i}] obj={f:.3f}")
    print("  refining best global hit..."); r = minimize(objective, best_x, method="Nelder-Mead",
                  options=dict(maxiter=1200, xatol=1e-4, fatol=1e-4))
    if objective(r.x) < best_f:
        best_f, best_x = objective(r.x), r.x.copy(); save_best(best_x); print(f"  -> obj={best_f:.3f}")

    # (2b) local random restarts AROUND the current best (wider jitter), refine each
    print("\nlocal random restarts around best...")
    for i in range(16):
        jitter = np.exp(rng.normal(0.0, 0.5, len(X0)))
        x = np.clip(best_x * jitter, LO, HI)
        r = minimize(objective, x, method="Nelder-Mead", options=dict(maxiter=600, xatol=1e-3, fatol=1e-3))
        f = objective(r.x)
        if f < best_f:
            best_f, best_x = f, r.x.copy()
            save_best(best_x); print(f"  [{i}] obj={f:.3f}")
    report(best_x, "REFINED BEST")
    print("\nFINAL PARAMS:")
    lines = [f'    "{n}": {v:.4g},' for n, v in zip(NAMES, best_x)]
    print("\n".join(lines))
    with open(os.path.join(os.path.dirname(__file__), "v44_recalibrate_gli_best.txt"), "w") as fh:
        fh.write(f"obj={objective(best_x):.4f}\n" + "\n".join(lines) + "\n")
