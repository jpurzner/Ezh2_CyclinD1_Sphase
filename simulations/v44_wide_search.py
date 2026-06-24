"""WIDE unattended parameter optimization for v44 (overnight run).

Searches the key tunable parameters against the FULL validation objective (the same targets as
validate_v44, plus the new wide-rShh EZH2 dose-response JP flagged), maximizing the number of soft
targets met subject to hard figure-critical guards. Hybrid search: mostly local Gaussian moves around
the best-so-far, with periodic wide random restarts for exploration. Seeds from the current (good)
parameter set.

Writes incrementally to a JSONL log (every eval) and a `_best.json` (best-so-far) so progress survives
and is reviewable any time. Runs until a wall-clock limit.

Run:  ./venv/bin/python simulations/v44_wide_search.py [seed] [hours] [out_prefix]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 1
HOURS = float(sys.argv[2]) if len(sys.argv) > 2 else 9.0
PREFIX = sys.argv[3] if len(sys.argv) > 3 else 'simulations/v44_wide_search'
OUT, BEST = f'{PREFIX}.jsonl', f'{PREFIX}_best.json'
T_END, N_PTS, SETTLE = 9000, 14000, 4000

_RR = te.loada(build_model_v44(with_hh=True))   # builder default (single-step Rb) -- matches validate_v44; exposes p16
_RR.integrator.setValue("relative_tolerance", 1e-6)
rng = np.random.default_rng(1000 + SEED)

# searched params: name -> (current/seed value, lo, hi, log?)
P = {
    'kPhRbCd':         (0.5,   0.35, 0.7,  False),   # commitment threshold (lower = higher threshold)
    'K_CdRb':          (0.357, 0.25, 0.55, False),
    'p18':             (0.4,   0.3,  1.0,  False),   # GNP baseline INK4 brake (MB overrides per-cond)
    'k_Cd_tx_basal':   (0.3086,0.12, 0.5,  False),
    'k_Cd_tx_Gli_max': (46.31, 28.0, 75.0, False),
    'k_Cd_tx_MYCN':    (35.22, 18.0, 55.0, False),
    'K_EZH2_repression':(0.75, 0.4,  1.2,  False),
    'Kez_cd':          (2.0,   1.0,  4.5,  False),   # EZH2 mitogen-dose saturation
    'kEZE2f':          (0.026, 0.012,0.045, True),
    'kEZbas':          (0.0003,0.0001,0.0012, True),
    'kDeEZ':           (5e-5,  2e-5, 2.5e-4, True),
    'k_Cd_translation':(0.75,  0.55, 0.95, False),
    'P16_MB':          (0.15,  0.05, 0.6,  False),   # MB CKI tones (applied only to MB conditions)
    'P18_MB':          (1.5,   0.9,  2.2,  False),
    'KSYP21_MB':       (0.004, 0.002,0.009, False),
    'KmHU_fork':       (0.35,  0.2,  0.6,  False),
    'wCe':             (0.5,   0.2,  0.8,  False),   # EZH2 cyclE/cyclA gate weight (S/G0 gradient)
}
KEYS = list(P)

# conditions: (SHH, ptch1, hhi, ezh2i, mycn, is_mb, cdk46i, hu)
CONDS = {
    'GNP':      (0.5, 1.0, 0, 0, 1.0, 0, 0, 0),
    'GNP-SHH':  (0.0, 1.0, 0, 0, 1.0, 0, 0, 0),
    'GNP+HHi':  (0.5, 1.0, 1, 0, 1.0, 0, 0, 0),
    'GNP+EZH2i':(0.5, 1.0, 0, 1, 1.0, 0, 0, 0),
    'MB':       (0.5, 0.1, 0, 0, 2.8, 1, 0, 0),
    'MB+HHi':   (0.5, 0.1, 1, 0, 2.8, 1, 0, 0),
    'MB+EZH2i': (0.5, 0.1, 0, 1, 2.8, 1, 0, 0),
    'MB+HHi+EZH2i':(0.5,0.1, 1, 1, 2.8, 1, 0, 0),
    'MB+CDK46i':(0.5, 0.1, 0, 0, 2.8, 1, 1, 0),
    'MB+CDK46i+EZH2i':(0.5,0.1,0,1, 2.8, 1, 1, 0),
}
DEF_kPhRbCd_key = 'kPhRbCd'
SEL = ["time", "MPF", "Cd", "Cd_mRNA", "MYCN", "Gli1", "EZH2"]


def _set(p, cond):
    shh, ptch, hhi, ezh2i, mycn, is_mb, cdk46i, hu = cond
    _RR['SHH'] = shh; _RR['Ptch1_copy_number'] = ptch; _RR['HHi'] = hhi
    _RR['EZH2i'] = ezh2i; _RR['MYCN_amplification'] = mycn; _RR['HU'] = hu
    for k in KEYS:
        if k in ('P16_MB', 'P18_MB', 'KSYP21_MB'):
            continue
        try: _RR[k] = p[k]
        except Exception: pass
    # CKI tones: GNP uses p16=0 + searched p18; MB uses searched MB tones
    if is_mb:
        _RR['p16'] = p['P16_MB']; _RR['p18'] = p['P18_MB']; _RR['kSyP21'] = p['KSYP21_MB']
    else:
        _RR['p16'] = 0.0; _RR['p18'] = p['p18']
    _RR['kPhRbCd'] = 0.0 if cdk46i else p['kPhRbCd']


def _run(p, cond):
    for te_end, te_pts in ((T_END, N_PTS), (7000, 11000), (5500, 8500)):   # horizon retry for stiff (MB+HHi+EZH2i) cond
        settle = min(SETTLE, te_end - 2500)
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            _RR.reset(); _set(p, cond)
            try: _RR.integrator.setValue("maximum_num_steps", 300000)
            except Exception: pass
            _RR.integrator.setValue("absolute_tolerance", atol)
            try:
                r = _RR.simulate(0, te_end, te_pts, selections=SEL)
                m = r['time'] >= settle; dt = r['time'][1] - r['time'][0]
                pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
                return dict(div=len(pk), Cd=float(np.mean(r['Cd'][m])), EZH2=float(np.mean(r['EZH2'][m])),
                            MYCN=float(np.mean(r['MYCN'][m])), Gli1=float(np.mean(r['Gli1'][m])))
            except Exception:
                continue
    return None


def _ez_dose(p):
    """EZH2 at GNP SHH=0.5,1,2,4 -> dose-response slope (should be wide/increasing)."""
    vals = []
    for shh in (0.5, 1.0, 2.0, 4.0):
        r = _run(p, (shh, 1.0, 0, 0, 1.0, 0, 0, 0))
        vals.append(r['EZH2'] if r else np.nan)
    return vals


def evaluate(p):
    S = {n: _run(p, c) for n, c in CONDS.items()}
    if any(S[n] is None for n in ('GNP', 'GNP-SHH', 'GNP+HHi', 'MB', 'MB+HHi+EZH2i')):
        return -1e3, {'fail': 'integration'}
    d = lambda n: S[n]['div']
    # ---- hard guards ----
    guards = {
        'GNP_cycles': d('GNP') > 0, 'GNP-SHH_arrest': d('GNP-SHH') == 0, 'GNP+HHi_arrest': d('GNP+HHi') == 0,
        'MB_cycles': d('MB') > 0, 'MB+HHi_arrest': d('MB+HHi') == 0, 'rescue': d('MB+HHi+EZH2i') > 0,
        'CDK46i_arrest': (S['MB+CDK46i'] is None) or d('MB+CDK46i') == 0,
        'CDK46i_norescue': (S['MB+CDK46i+EZH2i'] is None) or d('MB+CDK46i+EZH2i') == 0,
    }
    nguard = sum(guards.values())
    cd = lambda n: S[n]['Cd'] if S[n] else np.nan
    ez = lambda n: S[n]['EZH2'] if S[n] else np.nan
    my = lambda n: S[n]['MYCN'] if S[n] else np.nan
    gl = lambda n: S[n]['Gli1'] if S[n] else np.nan
    edose = _ez_dose(p)
    # ---- soft targets: (value, target, rel_tol) ----
    T = {
        'CyclinD1 MB/GNP': (cd('MB') / cd('GNP'), 5.07, 0.30),
        'CyclinD1 MB+HHi/MB': (cd('MB+HHi') / cd('MB'), 0.144, 0.30),
        'CyclinD1 GNP+HHi/GNP': (cd('GNP+HHi') / cd('GNP'), 0.157, 0.30),
        'EZH2 MB/GNP': (ez('MB') / ez('GNP'), 2.05, 0.25),
        'EZH2i CycD1 fold': (cd('GNP+EZH2i') / cd('GNP'), 2.2, 0.30),
        'MYCN MB/GNP': (my('MB') / my('GNP'), 2.8, 0.20),
        'MYCN GNP+HHi/GNP': (my('GNP+HHi') / my('GNP'), 0.78, 0.18),
        'MYCN MB+HHi/MB': (my('MB+HHi') / my('MB'), 0.86, 0.15),
        'Gli1 MB/GNP': (gl('MB') / gl('GNP'), 6.9, 0.30),
        'EZH2 dose wide (SHH4/0.5)': (edose[3] / edose[0] if edose[0] else np.nan, 1.9, 0.45),
        'EZH2 dose monotonic': (1.0 if (np.all(np.diff([e for e in edose if e == e]) > -0.02)) else 0.0, 1.0, 0.01),
    }
    npass = sum(1 for v, t, tol in T.values() if (v == v) and abs(v - t) <= tol * abs(t))
    # score: guards dominate; then soft passes; tie-break by negative total relative error
    relerr = sum(min(2.0, abs(v - t) / abs(t)) for v, t, tol in T.values() if v == v)
    score = nguard * 100 + npass * 10 - relerr
    return score, dict(guards=guards, nguard=nguard, npass=npass,
                       targets={k: round(v[0], 4) if v[0] == v[0] else None for k, v in T.items()},
                       edose=[round(e, 3) if e == e else None for e in edose])


def sample(best_p):
    """local Gaussian move around best, with 20% wide random restart."""
    p = {}
    wide = rng.random() < 0.20
    for k, (cur, lo, hi, islog) in P.items():
        b = best_p[k]
        if wide:
            v = (10 ** rng.uniform(np.log10(lo), np.log10(hi))) if islog else rng.uniform(lo, hi)
        else:
            s = 0.12 * (np.log10(hi) - np.log10(lo)) if islog else 0.12 * (hi - lo)
            v = 10 ** (np.log10(b) + rng.normal(0, s)) if islog else b + rng.normal(0, s)
        p[k] = float(np.clip(v, lo, hi))
    return p


def main():
    t0 = time.time()
    best_p = {k: P[k][0] for k in KEYS}
    best_score, best_det = evaluate(best_p)
    n = 0
    with open(BEST, 'w') as f:
        json.dump(dict(score=best_score, params=best_p, detail=best_det), f, indent=2)
    open(OUT, 'w').close()
    print(f'seed {SEED}: seed score={best_score:.1f} (guards {best_det.get("nguard")}/8, soft {best_det.get("npass")}/11)', flush=True)
    while time.time() - t0 < HOURS * 3600:
        p = sample(best_p)
        try:
            sc, det = evaluate(p)
        except Exception as e:
            sc, det = -1e4, {'err': str(e)[:80]}
        n += 1
        with open(OUT, 'a') as f:
            f.write(json.dumps(dict(n=n, t=round(time.time() - t0, 1), score=round(sc, 2),
                                    nguard=det.get('nguard'), npass=det.get('npass'), params=p)) + '\n')
        if sc > best_score:
            best_score, best_p, best_det = sc, p, det
            with open(BEST, 'w') as f:
                json.dump(dict(score=best_score, params=best_p, detail=best_det, eval=n), f, indent=2)
            print(f'  [{n}] NEW BEST score={sc:.1f} guards={det["nguard"]}/8 soft={det["npass"]}/11', flush=True)
        if n % 25 == 0:
            print(f'  [{n}] {(time.time()-t0)/3600:.1f}h elapsed; best={best_score:.1f} (g{best_det.get("nguard")}/s{best_det.get("npass")})', flush=True)
    print(f'DONE seed {SEED}: {n} evals, best score {best_score:.1f} -> {BEST}', flush=True)


if __name__ == '__main__':
    main()
