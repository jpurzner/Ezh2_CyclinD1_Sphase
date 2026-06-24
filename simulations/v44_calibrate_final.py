"""Constrained whole-model calibration for v44 figure generation.

Random search over the HH/MYCN->CyclinD1 + EZH2-repression parameters with a HARD GUARD:
any candidate that (a) crashes CVODE on any key condition, or (b) breaks the figure-critical
cycling behavior (GNP+SHH & MB must cycle; GNP+HHi, MB+HHi must arrest; MB+HHi+EZH2i must
rescue) is rejected. Among survivors, minimize a weighted error on the figure-relevant ratios.

This deliberately bounds k_Cd_tx_MYCN / n_MYCN_Cd away from the stiff regime that crashes
MB+HHi (the raw transcript-only optimizer hit n=5.6, k=38 and crashed the rescue conditions).

Run (background):  ./venv/bin/python simulations/v44_calibrate_final.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

COND = {
    "GNP":     dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=0.0, EZH2i=0.0),
    "GNP+HHi": dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=1.0, EZH2i=0.0),
    "GNP+EZH2i":dict(SHH=0.5,Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=0.0, EZH2i=1.0),
    "MB":      dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=0.0, EZH2i=0.0),
    "MB+HHi":  dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=1.0, EZH2i=0.0),
    "MB+HHi+EZH2i": dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=1.0, EZH2i=1.0),
}
SPEC = ["time", "Cd_mRNA", "MYCN", "Gli1", "EZH2", "MPF"]

PSPEC = [   # name, low, high
    ("k_Cd_tx_basal",   0.10, 0.45),
    ("k_Cd_tx_Gli_max", 4.0,  18.0),
    ("K_Gli_act_CycD",  0.30, 0.80),
    ("k_Cd_tx_MYCN",    10.0, 26.0),
    ("K_MYCN_Cd",       1.0,  1.8),
    ("n_MYCN_Cd",       3.0,  4.0),
    ("K_Ptch_Smo",      0.10, 0.40),
    ("K_EZH2_repression", 0.5, 1.0),
]
NAMES = [p[0] for p in PSPEC]
LO = np.array([p[1] for p in PSPEC]); HI = np.array([p[2] for p in PSPEC])

# figure-relevant ratio targets: (species, num, den, target, weight)
TARGETS = [
    ("Cd_mRNA", "GNP+HHi",       "GNP", 0.14, 2.0),
    ("Cd_mRNA", "MB+HHi",        "MB",  0.40, 2.0),
    ("Cd_mRNA", "MB",            "GNP", 5.07, 0.7),   # soft (architectural tension)
    ("Cd_mRNA", "GNP+EZH2i",     "GNP", 2.20, 1.5),
    ("EZH2",    "MB",            "GNP", 2.05, 1.0),
    ("MYCN",    "MB",            "GNP", 2.80, 0.5),
]


def ndiv(r):
    t = r["time"]; m = t >= 3000; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r["MPF"][m], prominence=0.15, distance=int(200 / dt))
    return len(pk)


def measure(x):
    params = dict(zip(NAMES, x))
    try:
        model = build_model_v44(with_ezh2=True, with_hh=True, params=params)
    except Exception:
        return None
    out = {}
    for name, inp in COND.items():
        try:
            rr = te.loada(model)
            rr.integrator.setValue("absolute_tolerance", 1e-9)
            rr.integrator.setValue("relative_tolerance", 1e-6)
            for k, v in inp.items():
                rr[k] = v
            r = rr.simulate(0, 8000, 16000, selections=SPEC)
            sel = r["time"] >= 4000
            out[name] = {s: max(float(np.mean(r[s][sel])), 1e-6) for s in SPEC[1:]}
            out[name]["div"] = ndiv(r)
        except Exception:
            return None                       # any crash -> reject candidate
    return out


def guard_ok(o):
    return (o["GNP"]["div"] > 0 and o["MB"]["div"] > 0 and
            o["GNP+HHi"]["div"] == 0 and o["MB+HHi"]["div"] == 0 and
            o["MB+HHi+EZH2i"]["div"] > 0)


def objective(o):
    err = 0.0
    for sp, a, b, tgt, w in TARGETS:
        ratio = o[a][sp] / o[b][sp]
        err += w * (np.log(ratio / tgt)) ** 2
    return err


def report(x, o, label=""):
    print(f"\n--- {label}  obj={objective(o):.3f}  guard={'OK' if guard_ok(o) else 'FAIL'} ---")
    print("  params: " + ", ".join(f"{n}={v:.3g}" for n, v in zip(NAMES, x)))
    print("  div: " + ", ".join(f"{c}={o[c]['div']}" for c in COND))
    for sp, a, b, tgt, w in TARGETS:
        print(f"    {sp:8} {a:>12}/{b:<4} = {o[a][sp]/o[b][sp]:6.3f}  target {tgt}")


if __name__ == "__main__":
    print("=" * 70)
    print("v44 CONSTRAINED CALIBRATION (figure-safe; guard on arrest/rescue/stability)")
    rng = np.random.default_rng(1)
    best = None
    survivors = 0
    N = 300
    for i in range(N):
        x = LO + (HI - LO) * rng.random(len(LO))
        o = measure(x)
        if o is None or not guard_ok(o):
            continue
        survivors += 1
        f = objective(o)
        if best is None or f < best[0]:
            best = (f, x.copy(), o)
            print(f"  [{i}] NEW BEST obj={f:.3f}  "
                  f"Cd MB/GNP={o['MB']['Cd_mRNA']/o['GNP']['Cd_mRNA']:.2f}  "
                  f"MB+HHi/MB={o['MB+HHi']['Cd_mRNA']/o['MB']['Cd_mRNA']:.2f}")
    print(f"\nsurvivors: {survivors}/{N}")
    if best:
        report(best[1], best[2], "BEST")
        print("\nFINAL PARAMS:")
        for n, v in zip(NAMES, best[1]):
            print(f'    "{n}": {v:.4g},')
    else:
        print("NO SURVIVORS — relax guard or ranges.")
