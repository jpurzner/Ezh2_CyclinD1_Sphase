"""Broad UNATTENDED parameter search for the vismodegib / EZH2i population rescue.

Target (JP's pRb+ data): MB(DMSO) ~100%  ->  MB+HHi ~25%  ->  MB+HHi+EZH2i ~75%  ;  CDK4/6i 0%.
Current model (p16 competitive brake) gives ~100/18/58 -- this search tunes the rescue-determining
parameters to best hit 100/25/75 while keeping the figure-critical guards (GNP cycles, GNP+HHi arrests,
MB cycles, CDK4/6i not rescued).

Searched params (the rest of the model -- HH/expression calibration -- is held fixed; it already
passes its targets and varying it needs the expensive multi-condition expression fit):
  p16_MB    : MB competitive CDK4/6 brake (eff K_CdRb = K_CdRb*(1+p16))   [1.0 .. 6.0]
  sigma_cd  : cell-to-cell CyclinD1 setpoint heterogeneity (lognormal)    [0.5 .. 1.2]
  K_EZ      : K_EZH2_repression (EZH2i de-repression strength)            [0.25 .. 0.9]
  ktl0      : k_Cd_translation (CyclinD1 protein scale; GNP must cycle)   [0.70 .. 0.95]

Robust integration: atol retries 1e-9..1e-5 + max steps; cells that fail ALL tolerances are EXCLUDED
from the fraction denominator (so integration crashes don't bias the rescue downward).

Run (one of several in parallel):
  ./venv/bin/python simulations/v44_rescue_search.py <seed> <out.jsonl> <time_limit_h>
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEED   = int(sys.argv[1]) if len(sys.argv) > 1 else 1
OUT    = sys.argv[2] if len(sys.argv) > 2 else f"simulations/v44_rescue_search_{SEED}.jsonl"
TLIMIT = float(sys.argv[3]) if len(sys.argv) > 3 else 6.5     # hours

T_END, N_PTS, SETTLE = 12000, 24000, 4000
N_MB, N_HHI, N_RES, N_CDK = 60, 90, 90, 15      # cells per condition
RANGES = dict(p16=(1.0, 6.0), sigma=(0.5, 1.2), kez=(0.25, 0.9), ktl0=(0.70, 0.95))
TGT = dict(MB=100.0, HHi=25.0, RES=75.0, CDK=0.0)

_RR = te.loada(build_model_v44())
_RR.integrator.setValue("relative_tolerance", 1e-6)
DEF_kPhRbCd = _RR['kPhRbCd']
rng = np.random.default_rng(1000 + SEED)
Z = rng.normal(0.0, 1.0, max(N_MB, N_HHI, N_RES))   # fixed standard-normal draws -> sigma scales them


def _cycles(ktl, p16, kez, gdc, ezh2i, cdk46i):
    """True/False cycling, or None if integration fails at all tolerances."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        _RR.reset()
        _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 0.1; _RR['MYCN_amplification'] = 2.86
        _RR['GDC0449'] = gdc; _RR['EZH2i'] = ezh2i; _RR['p16'] = p16
        _RR['k_Cd_translation'] = ktl; _RR['K_EZH2_repression'] = kez
        _RR['kPhRbCd'] = 0.0 if cdk46i else DEF_kPhRbCd
        _RR.integrator.setValue("absolute_tolerance", atol)
        try:
            _RR.integrator.setValue("maximum_num_steps", 200000)
        except Exception:
            pass
        try:
            r = _RR.simulate(0, T_END, N_PTS, selections=["time", "MPF"])
        except Exception:
            continue
        t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
        pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
        return len(pk) >= 2
    return None


def _frac(ktl_arr, p16, kez, gdc, ezh2i, cdk46i):
    """Cycling % over cells that integrated (crashes excluded). Returns (pct, n_ok, n_crash)."""
    cyc = ok = crash = 0
    for ktl in ktl_arr:
        c = _cycles(ktl, p16, kez, gdc, ezh2i, cdk46i)
        if c is None:
            crash += 1
        else:
            ok += 1; cyc += int(c)
    return (100.0 * cyc / ok if ok else 0.0), ok, crash


def _gnp_cycles(ktl, kez, gdc):
    _RR.reset()
    _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 1.0; _RR['MYCN_amplification'] = 1.0
    _RR['GDC0449'] = gdc; _RR['EZH2i'] = 0; _RR['p16'] = 0.0
    _RR['k_Cd_translation'] = ktl; _RR['K_EZH2_repression'] = kez; _RR['kPhRbCd'] = DEF_kPhRbCd
    for atol in (1e-9, 1e-8, 1e-7):
        _RR.reset(); _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 1.0; _RR['MYCN_amplification'] = 1.0
        _RR['GDC0449'] = gdc; _RR['p16'] = 0.0; _RR['k_Cd_translation'] = ktl; _RR['K_EZH2_repression'] = kez
        _RR.integrator.setValue("absolute_tolerance", atol)
        try:
            r = _RR.simulate(0, T_END, N_PTS, selections=["time", "MPF"]); break
        except Exception:
            r = None
    if r is None:
        return -1
    t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt)); return len(pk)


def evaluate(p16, sigma, kez, ktl0):
    cd = np.clip(np.exp(Z * sigma), 0.2, 5.0)
    ktl = ktl0 * cd
    # guards (cheap, deterministic)
    gnp = _gnp_cycles(ktl0, kez, 0)
    gnp_hhi = _gnp_cycles(ktl0, kez, 1)
    guard_ok = (gnp >= 2) and (gnp_hhi == 0)
    mb, mb_ok, mb_cr = _frac(ktl[:N_MB], p16, kez, 0, 0, False)
    hhi, hhi_ok, hhi_cr = _frac(ktl[:N_HHI], p16, kez, 1, 0, False)
    res, res_ok, res_cr = _frac(ktl[:N_RES], p16, kez, 1, 1, False)
    cdk, cdk_ok, cdk_cr = _frac(ktl[:N_CDK], p16, kez, 0, 0, True)
    obj = (0.02 * (mb - TGT['MB'])**2 + 0.04 * (hhi - TGT['HHi'])**2
           + 0.04 * (res - TGT['RES'])**2 + 0.05 * (cdk - TGT['CDK'])**2)
    if not guard_ok:
        obj += 1000.0
    return dict(p16=p16, sigma=sigma, kez=kez, ktl0=ktl0, obj=round(obj, 3),
                MB=round(mb, 1), HHi=round(hhi, 1), RES=round(res, 1), CDK=round(cdk, 1),
                gnp_div=gnp, gnp_hhi_div=gnp_hhi, guard_ok=guard_ok,
                crashes=dict(MB=mb_cr, HHi=hhi_cr, RES=res_cr, CDK=cdk_cr))


def sample(best):
    """30% local jitter around best (exploit), else uniform (explore)."""
    if best is not None and rng.random() < 0.30:
        out = {}
        for k, (lo, hi) in RANGES.items():
            span = (hi - lo)
            out[k] = float(np.clip(best[k] + rng.normal(0, 0.12 * span), lo, hi))
        return out
    return {k: float(rng.uniform(lo, hi)) for k, (lo, hi) in RANGES.items()}


if __name__ == "__main__":
    t0 = time.time()
    best = None; n = 0
    print(f"[seed {SEED}] rescue search -> {OUT}  (limit {TLIMIT}h)  target 100/25/75/0", flush=True)
    with open(OUT, "a") as fh:
        while time.time() - t0 < TLIMIT * 3600:
            p = sample(best)
            try:
                rec = evaluate(p['p16'], p['sigma'], p['kez'], p['ktl0'])
            except Exception as e:
                rec = dict(error=str(e)[:120], **p)
            rec['t_min'] = round((time.time() - t0) / 60, 1); rec['seed'] = SEED; rec['n'] = n
            fh.write(json.dumps(rec) + "\n"); fh.flush()
            if 'obj' in rec and (best is None or rec['obj'] < best['obj']):
                best = rec
                print(f"[seed {SEED}] n={n} obj={rec['obj']} "
                      f"MB={rec['MB']} HHi={rec['HHi']} RES={rec['RES']} CDK={rec['CDK']} "
                      f"p16={p['p16']:.2f} sig={p['sigma']:.2f} kez={p['kez']:.2f} ktl={p['ktl0']:.2f}",
                      flush=True)
            n += 1
    print(f"[seed {SEED}] DONE {n} evals in {(time.time()-t0)/3600:.2f}h. BEST: {json.dumps(best)}", flush=True)
