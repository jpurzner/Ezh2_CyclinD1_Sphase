"""Quick check: the model's TRUE quiescent (G0) fraction with the CALIBRATED absolute pRb gate (pRb<1.5,
same as validate_v44), MB vs GNP, and the EZH2 level of quiescent vs cycling cells. The 35th-percentile gate
used in the distribution sim is circular (forces ~35% G0); this uses the real threshold to measure how far
the model is from the measured IF MB transient-G0 (~34% Ki67-low). DMSO only, warmup, CDKI heterogeneity."""
import sys, os, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = 120; DT = 30.0; WARM = 5400.0
SIG_CDKI = float(os.environ.get('SIG_CDKI', '0.55')); TAIL_FRAC = float(os.environ.get('TAIL_FRAC', '0.12')); TAIL_X = float(os.environ.get('TAIL_X', '4.0'))
P18_MB = float(os.environ.get('P18_MB', '1.73'))   # MB CDKI level (push this to raise transient G0)
PRB_G0 = 1.5                                        # calibrated absolute G0 threshold (validate_v44)
rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('absolute_tolerance', 1e-8)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']
rng = np.random.default_rng(7)
mult = np.where(rng.random(N) < TAIL_FRAC, np.clip(rng.lognormal(np.log(TAIL_X), 0.25, N), 0.3, 8.0),
                np.clip(rng.lognormal(0, SIG_CDKI, N), 0.3, 8.0))
warm = rng.uniform(600.0, WARM, N)
SEL = ['EZH2', 'pRb', 'Dna']


def run(ctx, i):
    rr.reset(); rr['SHH'] = 0.5
    if ctx == 'GNP':
        rr['Ptch1_copy_number'] = 1.0; rr['p18'] = 0.464 * mult[i]; rr['p19'] = 0.36 * mult[i]
    else:
        rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = P18_MB * mult[i]; rr['p19'] = 0.58 * mult[i]
    rr['kSyP21'] = KS0 * mult[i]; rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    try:
        r = rr.simulate(0, warm[i], max(2, int(warm[i] / DT)), selections=SEL)
    except Exception:
        return None
    return {k: float(r[k][-1]) for k in SEL}


for ctx in ['MB', 'GNP']:
    cells = [c for c in (run(ctx, i) for i in range(N)) if c and np.isfinite(c['pRb'])]
    pr = np.array([c['pRb'] for c in cells]); ez = np.array([c['EZH2'] for c in cells])
    q = pr < PRB_G0                                  # quiescent = pRb below the committed threshold
    qf = q.mean()
    ezq = np.median(ez[q]) if q.sum() else np.nan; ezc = np.median(ez[~q]) if (~q).sum() else np.nan
    print(f"{ctx:4s} (P18_MB={P18_MB if ctx=='MB' else 0.464}): quiescent(pRb<1.5) {qf:.0%}  | EZH2 quiescent {ezq:.2f} vs cycling {ezc:.2f} (fold {ezq/ezc:.2f})  n={len(cells)}")
print("(IF: MB ~34% Ki67-low with EZH2 ~16x below cycling; GNP should be ~cycling)")
