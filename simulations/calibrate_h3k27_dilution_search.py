"""Automated fold-preserving search for the replicative-dilution H3K27me3 params {k_w,k0,del,K,n}.
Hybrid local-Gaussian + wide-restart search; 6 hard guards (GNP/MB cycle, GNP+HHi/MB+HHi arrest,
MB+HHi+EZH2i rescue, GNP+EZH2i de-repress) + soft targets (MB/GNP fold 5.07, EZH2i fold 2.2,
GNP/MB Cd ~1.9/8.3, periods ~22h). JSONL + _best.json incremental.
Run: ./venv/bin/python simulations/calibrate_h3k27_dilution_search.py [seed] [hours]
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
HOURS = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
PREFIX = f'simulations/h3k27_dilution_search_s{SEED}'
OUT, BEST = f'{PREFIX}.jsonl', f'{PREFIX}_best.json'
rng = np.random.default_rng(1000 + SEED)

_RR = te.loada(build_model_v44(with_hh=True, with_h3k27_dilution=True))
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)
PMK = ['k_w_mk', 'k0_mk', 'del_mk', 'K_mk', 'n_mk']
LO = np.array([0.0012, 0.0001, 0.0010, 0.22, 3.0])
HI = np.array([0.0035, 0.0005, 0.0025, 0.34, 6.5])
SEED_P = np.array([0.0023, 0.0002, 0.0015, 0.27, 4.0])

COND = {
    'GNP':          dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, GDC0449=0, EZH2i=0),
    'GNP+HHi':      dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, GDC0449=0.95, EZH2i=0),
    'GNP+EZH2i':    dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, GDC0449=0, EZH2i=1),
    'MB':           dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, GDC0449=0, EZH2i=0),
    'MB+HHi':       dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, GDC0449=0.95, EZH2i=0),
    'MB+HHi+EZH2i': dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, GDC0449=0.95, EZH2i=1),
}


def _run(p, cond):
    for hz, npts, settle in [(11000, 22000, 7000), (8000, 16000, 5000), (5500, 11000, 3500)]:
        try:
            _RR.reset()
            for k, v in zip(PMK, p):
                _RR[k] = v
            for k, v in cond.items():
                _RR[k] = v
            _RR['Mk'] = 0.2
            r = _RR.simulate(0, hz, npts, selections=['time', 'MPF', 'Cd'])
            m = r['time'] >= settle
            pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=400)
            per = float(np.mean(np.diff(r['time'][m][pk])) / 60.0) if len(pk) > 1 else float('nan')
            return len(pk), float(r['Cd'][m].mean()), per
        except Exception:
            continue
    return -1, float('nan'), float('nan')


def evaluate(p):
    R = {c: _run(p, COND[c]) for c in COND}
    div = {c: R[c][0] for c in COND}
    cd = {c: R[c][1] for c in COND}
    per = {c: R[c][2] for c in COND}
    guards = [div['GNP'] > 0, div['MB'] > 0, div['GNP+HHi'] == 0, div['MB+HHi'] == 0,
              div['MB+HHi+EZH2i'] > 0, div['GNP+EZH2i'] > 0]
    nguard = sum(guards)
    fold = cd['MB'] / cd['GNP'] if cd['GNP'] else 0
    ezf = cd['GNP+EZH2i'] / cd['GNP'] if cd['GNP'] else 0
    soft = [('fold', fold, 5.07, 0.15), ('ezf', ezf, 2.2, 0.15), ('gnpCd', cd['GNP'], 1.9, 0.25),
            ('mbCd', cd['MB'], 8.3, 0.25), ('gnpPer', per['GNP'], 22, 0.2), ('mbPer', per['MB'], 22, 0.2)]
    relerr, npass = 0.0, 0
    for _, val, tgt, tol in soft:
        if not np.isfinite(val):
            relerr += 1.0; continue
        re = abs(val - tgt) / tgt
        relerr += min(re, 1.0)
        if re <= tol:
            npass += 1
    score = nguard * 100 + npass * 10 - relerr
    return score, dict(nguard=nguard, npass=npass, fold=round(fold, 2), ezf=round(ezf, 2),
                       div=div, cd={k: round(v, 2) for k, v in cd.items()}, per={k: round(v, 1) for k, v in per.items()})


def sample(best_p):
    if rng.random() < 0.2:
        return LO + rng.random(5) * (HI - LO)
    step = (HI - LO) * 0.08 * rng.standard_normal(5)
    return np.clip(best_p + step, LO, HI)


def main():
    t0 = time.time()
    best_p = SEED_P.copy()
    best_score, best_det = evaluate(best_p)
    open(OUT, 'w').close()
    with open(BEST, 'w') as f:
        json.dump(dict(score=best_score, params=dict(zip(PMK, best_p.tolist())), detail=best_det), f, indent=2)
    print(f'seed score={best_score:.1f} (guards {best_det["nguard"]}/6, soft {best_det["npass"]}/6) {best_det}', flush=True)
    n = 0
    while time.time() - t0 < HOURS * 3600:
        n += 1
        p = sample(best_p)
        sc, det = evaluate(p)
        with open(OUT, 'a') as f:
            f.write(json.dumps(dict(n=n, score=round(sc, 2), params=dict(zip(PMK, p.tolist())), det=det)) + '\n')
        if sc > best_score:
            best_score, best_p, best_det = sc, p, det
            with open(BEST, 'w') as f:
                json.dump(dict(score=best_score, params=dict(zip(PMK, best_p.tolist())), detail=best_det, eval=n), f, indent=2)
            print(f'  [{n}] NEW BEST {sc:.1f} g{det["nguard"]}/6 s{det["npass"]}/6 fold={det["fold"]} ezf={det["ezf"]}', flush=True)
        elif n % 20 == 0:
            print(f'  [{n}] {(time.time()-t0)/3600:.1f}h best={best_score:.1f}', flush=True)
    print(f'DONE: {n} evals, best {best_score:.1f} -> {BEST}', flush=True)


main()
