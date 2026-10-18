"""DMSO-only calibration variant of sim_ezh2_bursting.py (JP 2026-10-09). Sweeps the bursting amplitude
to hit a target EZH2 abundance CV (~0.45, the IF bright-population value). Env-parameterized so a workflow
can fan out candidates in parallel without editing the file:
    CALIB_N (default 90), CALIB_TAU (200), CALIB_DUTY (0.5), CALIB_SIG_EXT (0.15)
Prints one parseable line:  CALIB sig_ext=.. duty=.. tau=.. n=.. cv=.. g1_cv=.. median=..
DMSO only (the CV target is defined on DMSO); the full 3-treatment run + figure is done separately once
the winning amplitude is chosen.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('CALIB_N', '90'))
DT = 30.0
WARM = 5400.0; TREAT = 1440.0
TAU_BURST = float(os.environ.get('CALIB_TAU', '200.0'))
DUTY = float(os.environ.get('CALIB_DUTY', '0.5'))
SIG_EXT = float(os.environ.get('CALIB_SIG_EXT', '0.15'))
KDEEZM = float(os.environ.get('CALIB_KDEEZM', '0.0098'))   # EZH2 mRNA decay; 0.0098=t1/2 1.2h, 0.00385=t1/2 3h (measured)
SIG_CDKI = 0.55; TAIL_FRAC, TAIL_X = 0.12, 4.0
SEL = ['EZH2', 'Dna', 'aRc', 'vfork', 'pRb']
K_ON = DUTY / TAU_BURST; K_OFF = (1 - DUTY) / TAU_BURST

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']; KPH0 = rr['kPhRbCd']
COMP = KDEEZM / rr['kDeEZm']                               # hold mRNA AND protein mean fixed: scale tx by kDeEZm ratio
B0 = {p: rr[p] * COMP for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}
rng = np.random.default_rng(2)
mult = np.where(rng.random(N) < TAIL_FRAC, np.clip(rng.lognormal(np.log(TAIL_X), 0.25, N), 0.3, 8.0),
                np.clip(rng.lognormal(0, SIG_CDKI, N), 0.3, 8.0))
fExt = np.exp(rng.normal(0, SIG_EXT, N))
warm = rng.uniform(600.0, WARM, N)


def run_cell(i):
    rr.reset()
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult[i]; rr['p19'] = 0.58 * mult[i]; rr['kSyP21'] = KS0 * mult[i]
    rr['EZH2'] = 6.0 * fExt[i]; rr['EZH2m'] = 0.3 * fExt[i]
    rr['kDeEZm'] = KDEEZM
    rr['kPhRbCd'] = KPH0; rr['HU'] = 0.0
    g = 1 if rng.random() < DUTY else 0
    t = 0.0; T_END = warm[i] + TREAT
    while t < T_END:
        g = 0 if (g and rng.random() < K_OFF * DT) else (1 if (not g and rng.random() < K_ON * DT) else g)
        eta = fExt[i] * (g / DUTY)
        for p in B0: rr[p] = B0[p] * eta
        try:
            r = rr.simulate(t, t + DT, 2, selections=SEL)
        except Exception:
            return None
        t += DT
    return {k: float(r[k][-1]) for k in SEL}


cells = []
for i in range(N):
    c = run_cell(i)
    if c and np.isfinite(c['EZH2']) and c['EZH2'] > 0:
        cells.append(c)

allpr = np.array([c['pRb'] for c in cells])
PRB_THR = float(np.quantile(allpr[allpr > 0], 0.35))
def phase(c):
    if c['Dna'] >= 0.95: return 'G2'
    if c['vfork'] * c['aRc'] > 0.03: return 'S'
    return 'G0' if c['pRb'] < PRB_THR else 'G1'
ez = np.array([c['EZH2'] for c in cells]); ph = np.array([phase(c) for c in cells])
g1 = ez[ph == 'G1']
cv = float(ez.std(ddof=1) / ez.mean())
g1cv = float(g1.std(ddof=1) / g1.mean()) if len(g1) > 5 else float('nan')
print(f"CALIB sig_ext={SIG_EXT} duty={DUTY} tau={TAU_BURST} n={len(cells)} "
      f"cv={cv:.4f} g1_cv={g1cv:.4f} median={np.median(ez):.3f}", flush=True)
