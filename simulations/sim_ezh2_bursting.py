"""Transcriptional-bursting EZH2 noise over time (JP 2026-10-09). Replaces the frozen per-cell factor
with a VARIABLE PRODUCTION RATE: EZH2 transcription fires in bursts (a telegraph gene on/off), so the
production rate itself is a stochastic, time-varying quantity. The stable EZH2 protein (t1/2 ~42 h)
integrates the bursts into SLOW drift (correlation ~ the protein lifetime, ~2 cell cycles).

Explicit, STATED assumption about the split the IF snapshot cannot resolve (temporal vs cell-to-cell):
  - INTRINSIC (bursting, within-cell, time-varying): telegraph on EZH2 transcription (kEZbas/kEZbas_Cd/
    kEZE2f scaled by g/duty; g in {0,1}). Keeps the E2f coupling, so CDK4/6i still depletes EZH2.
  - EXTRINSIC (frozen per-cell, cell-to-cell): a per-cell lognormal factor on the same transcription.
  total CV^2 = phase^2 + intrinsic^2 + extrinsic^2. Calibrated (DMSO CV ~0.45): abundance variance is
  BURSTING(intrinsic)-DOMINATED (~3/4 of non-phase variance); the frozen extrinsic is a small floor
  (sig_ext 0.10). The snapshot cannot resolve this split; it is a stated, defensible guess.

Hybrid integrator: chunk the deterministic ODE (dt=24 min); update the telegraph each chunk and set the
transcription params, apply treatment past onset, integrate. Async population via per-cell CDKI mult +
random warmup. Reads EZH2 by phase (model's clean pRb). Treatments: DMSO / CDK4/6i (kPhRbCd 0) / HU (1).

Run:  ./venv/bin/python simulations/sim_ezh2_bursting.py   (~18 min)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

N = 120
DT = 30.0                 # min per chunk
WARM = 5400.0; TREAT = 1440.0     # 90 h warmup (EZH2 t1/2 ~42 h, start at stationary), then 24 h treatment
TAU_BURST = 200.0         # min; telegraph correlation time (hours-scale bursts)
DUTY = 0.62               # gene ON fraction; calibrated so DMSO total CV ~0.45 (IF bright-pop)
SIG_EXT = 0.10            # frozen extrinsic ln-sd (cell-to-cell); calibrated to DMSO CV ~0.45 (IF bright-pop)
SIG_CDKI = 0.55; TAIL_FRAC, TAIL_X = 0.12, 4.0
SEL = ['EZH2', 'Dna', 'aRc', 'vfork', 'pRb']
K_ON = DUTY / TAU_BURST; K_OFF = (1 - DUTY) / TAU_BURST

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']; KPH0 = rr['kPhRbCd']; EZ0 = rr['EZH2']
B0 = {p: rr[p] for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}
rng = np.random.default_rng(2)
mult = np.where(rng.random(N) < TAIL_FRAC, np.clip(rng.lognormal(np.log(TAIL_X), 0.25, N), 0.3, 8.0),
                np.clip(rng.lognormal(0, SIG_CDKI, N), 0.3, 8.0))
fExt = np.exp(rng.normal(0, SIG_EXT, N))      # frozen extrinsic per cell
warm = rng.uniform(600.0, WARM, N)

TREAT_APPLY = {'DMSO': lambda: None,
               'CDK4/6i': lambda: rr.__setitem__('kPhRbCd', 0.0),
               'HU': lambda: rr.__setitem__('HU', 1.0)}


def run_cell(i, treat, record_trace=False):
    rr.reset()
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult[i]; rr['p19'] = 0.58 * mult[i]; rr['kSyP21'] = KS0 * mult[i]
    rr['EZH2'] = 6.0 * fExt[i]; rr['EZH2m'] = 0.3 * fExt[i]   # start near stationary (t1/2 ~42 h)
    rr['kPhRbCd'] = KPH0; rr['HU'] = 0.0
    g = 1 if rng.random() < DUTY else 0
    t = 0.0; T_END = warm[i] + TREAT; applied = False
    tr_t, tr_e = [], []
    while t < T_END:
        g = 0 if (g and rng.random() < K_OFF * DT) else (1 if (not g and rng.random() < K_ON * DT) else g)
        eta = fExt[i] * (g / DUTY)                            # extrinsic (frozen) x burst (telegraph)
        for p in B0: rr[p] = B0[p] * eta
        if (not applied) and t >= warm[i]:
            TREAT_APPLY[treat](); applied = True
        try:
            r = rr.simulate(t, t + DT, 2, selections=SEL)
        except Exception:
            return (None, None) if not record_trace else (None, None)
        if record_trace and t >= warm[i] - 2880:             # last ~48h+treat of trace
            tr_t.append(t); tr_e.append(float(r['EZH2'][-1]))
        t += DT
    end = {k: float(r[k][-1]) for k in SEL}
    return end, (np.array(tr_t), np.array(tr_e)) if record_trace else None


res = {t: [] for t in TREAT_APPLY}; traces = []
for treat in TREAT_APPLY:
    for i in range(N):
        end, tr = run_cell(i, treat, record_trace=(treat == 'DMSO' and i < 6))
        if end and np.isfinite(end['EZH2']) and end['EZH2'] > 0:
            res[treat].append(end)
        if tr is not None and len(tr[0]):
            traces.append(tr)
    print(f'{treat}: {len(res[treat])}/{N}', flush=True)

allpr = np.array([c['pRb'] for t in res for c in res[t]])
PRB_THR = float(np.quantile(allpr[allpr > 0], 0.35))
def phase(c):
    if c['Dna'] >= 0.95: return 'G2'
    if c['vfork'] * c['aRc'] > 0.03: return 'S'
    return 'G0' if c['pRb'] < PRB_THR else 'G1'

out = {'PRB_THR': PRB_THR, 'params': dict(TAU_BURST=TAU_BURST, DUTY=DUTY, SIG_EXT=SIG_EXT),
       'treatments': {}, 'traces': [{'t': (tr[0] / 60).round(2).tolist(), 'ez': tr[1].round(2).tolist()} for tr in traces]}
for t, cells in res.items():
    ez = np.array([c['EZH2'] for c in cells]); ph = np.array([phase(c) for c in cells])
    g1 = ez[ph == 'G1']
    out['treatments'][t] = {'log10_ezh2': np.log10(ez).round(4).tolist(), 'phase': ph.tolist(),
                            'median': float(np.median(ez)), 'cv': float(ez.std(ddof=1) / ez.mean()),
                            'g1_cv': float(g1.std(ddof=1) / g1.mean()) if len(g1) > 5 else None}
sc = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/model_ezh2_bursting.json"
json.dump(out, open(sc, 'w')); print('saved', sc)
for t in res:
    d = out['treatments'][t]
    fr = {p: round(float(np.mean(np.array(d['phase']) == p)), 2) for p in ['G0', 'G1', 'S', 'G2']}
    print(f"{t:8s}: median {d['median']:7.2f}  CV {d['cv']:.2f}  G1_CV {d['g1_cv']}  phase {fr}")
