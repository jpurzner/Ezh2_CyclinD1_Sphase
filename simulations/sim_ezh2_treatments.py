"""Replicate the single-cell EZH2 distributions under DMSO / CDK4/6i / HU with the model (JP 2026-10-09).

Async MB population: each cell gets (i) a per-cell CDK-inhibitor multiplier m ~ lognormal (desynchronizes
cells and makes a high-CDKI arrest/EZH2-low tail = the G0 low mode), and (ii) a per-cell lognormal EZH2
factor fE (the extrinsic abundance noise that sets the bright-peak width). Each cell is warmed up a random
time (async), then TREATED for 24 h, then EZH2 + phase markers are read — mirroring the 24 h IF.

Treatments (validated harness): DMSO (none); CDK4/6i = kPhRbCd=0 (block CyclinD-CDK4/6->Rb); HU = HU=1
(fork stall). Phases use the model's clean pRb (G0 vs G1), synthesis flux vfork*aRc (S), DNA (G2).

Compares to IF (AMOS AC_MB55_CC_Ezh2, if_EZH2 histograms). Units arbitrary -> align DMSO bright-peak
median only. Tests: bimodality, CDK4/6i EZH2 depletion, HU EZH2 preservation, bright-peak CV ~0.4.

Run:  ./venv/bin/python simulations/sim_ezh2_treatments.py   (~15 min)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

N = 400
SIG_CDKI = 0.55          # per-cell CDKI log-sd (bulk variance -> fate spread + arrest tail)
TAIL_FRAC, TAIL_X = 0.12, 4.0   # a high-CDKI tail (Atoh1+ G0 fraction) -> DMSO EZH2-low mode
SIG_EZ = 0.33            # per-cell EZH2 factor ln-sd (extrinsic abundance noise; bright-peak width knob)
WARM = 1800.0; TREAT = 1440.0    # warmup up to 30 h (async), then 24 h treatment
SEL = ['EZH2', 'Dna', 'aRc', 'vfork', 'pRb']

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']; KTLEZ0 = rr['kTlEZ']; KPH0 = rr['kPhRbCd']; EZ0 = rr['EZH2']
rng = np.random.default_rng(1)
mult = np.where(rng.random(N) < TAIL_FRAC, np.clip(rng.lognormal(np.log(TAIL_X), 0.25, N), 0.3, 8.0),
                np.clip(rng.lognormal(0, SIG_CDKI, N), 0.3, 8.0))
fE = np.exp(rng.normal(0, SIG_EZ, N))
warm = rng.uniform(200.0, WARM, N)

TREAT_SET = {'DMSO': lambda: None,
             'CDK4/6i': lambda: rr.__setitem__('kPhRbCd', 0.0),
             'HU': lambda: rr.__setitem__('HU', 1.0)}


def run_cell(i, treat):
    rr.reset()
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult[i]; rr['p19'] = 0.58 * mult[i]; rr['kSyP21'] = KS0 * mult[i]
    rr['kTlEZ'] = KTLEZ0 * fE[i]; rr['EZH2'] = EZ0 * fE[i]
    rr['kPhRbCd'] = KPH0; rr['HU'] = 0.0
    try:
        rr.simulate(0, warm[i], 2, selections=['time'])     # async warmup, untreated
        TREAT_SET[treat]()                                  # apply drug
        r = rr.simulate(warm[i], warm[i] + TREAT, 2, selections=SEL)
    except Exception:
        return None
    return {k: float(r[k][-1]) for k in SEL}


res = {t: [] for t in TREAT_SET}
for treat in TREAT_SET:
    for i in range(N):
        c = run_cell(i, treat)
        if c and np.isfinite(c['EZH2']) and c['EZH2'] > 0:
            res[treat].append(c)
    print(f'{treat}: {len(res[treat])}/{N} cells', flush=True)

# phase gate (model): G2=4N DNA; S=synthesis flux; else G0 (pRb-low) vs G1 (pRb-high)
allpr = np.array([c['pRb'] for t in res for c in res[t]])
PRB_THR = float(np.quantile(allpr[allpr > 0], 0.35))      # pRb G0/G1 split
def phase(c):
    if c['Dna'] >= 0.95: return 'G2'
    if c['vfork'] * c['aRc'] > 0.03: return 'S'
    return 'G0' if c['pRb'] < PRB_THR else 'G1'

out = {'PRB_THR': PRB_THR, 'SIG_EZ': SIG_EZ, 'treatments': {}}
for t, cells in res.items():
    ez = np.array([c['EZH2'] for c in cells]); ph = np.array([phase(c) for c in cells])
    out['treatments'][t] = {'log10_ezh2': np.log10(ez).round(4).tolist(),
                            'phase': ph.tolist(),
                            'median': float(np.median(ez)),
                            'cv': float(ez.std(ddof=1) / ez.mean())}
sc = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/model_ezh2_treatments.json"
json.dump(out, open(sc, 'w'))
print('saved', sc)
for t in res:
    ez = np.array([c['EZH2'] for c in res[t]]); ph = np.array([phase(c) for c in res[t]])
    fr = {p: round(float(np.mean(ph == p)), 2) for p in ['G0', 'G1', 'S', 'G2']}
    print(f'{t:8s}: median EZH2 {np.median(ez):8.1f}  CV {ez.std(ddof=1)/ez.mean():.2f}  phase-frac {fr}')
