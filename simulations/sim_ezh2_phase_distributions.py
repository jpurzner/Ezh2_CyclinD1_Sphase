"""Population sim for the EZH2 per-phase DISTRIBUTION test (JP 2026-10-17). Many cells over time (bursting),
recording EZH2 mRNA + PROTEIN + the phase markers at a snapshot, per treatment. Outputs per-phase (G0/G1/
S/G2) distributions for BOTH species to compare vs data/if_ezh2_phase_targets.json (protein: magnitude+
shape; transcript: fold-over-cycle vs scRNA). The S-phase boost (G0-low -> cycling-high, real in protein &
RNA, JP) should emerge from the model's G0 Cdh1-clearance + E2f-driven cycling expression.
Phase gate (model's own markers, parallels the IF gate): G2 Dna>=0.95; S vfork*aRc>0.03; else G0 pRb<thr / G1.
Treatments: DMSO / CDK4/6i (kPhRbCd 0) / HU (1). Run:  ./venv/bin/python simulations/sim_ezh2_phase_distributions.py
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

N = 160; DT = 30.0
WARM = 5400.0; TREAT = 1440.0
TAU = 200.0; DUTY = 0.62; SIG_EXT = 0.10
SIG_CDKI = 0.55; TAIL_FRAC, TAIL_X = 0.12, 4.0          # CDKI tail -> a real G0 (arrested) fraction
SEL = ['EZH2', 'EZH2m', 'Dna', 'aRc', 'vfork', 'pRb']
BINS = np.linspace(-1.0, 1.7, 40)                       # log10 of model EZH2 (median-free; aligned at scoring)
K_ON = DUTY / TAU; K_OFF = (1 - DUTY) / TAU

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']; KPH0 = rr['kPhRbCd']
B0 = {p: rr[p] for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}
rng = np.random.default_rng(5)
mult = np.where(rng.random(N) < TAIL_FRAC, np.clip(rng.lognormal(np.log(TAIL_X), 0.25, N), 0.3, 8.0),
                np.clip(rng.lognormal(0, SIG_CDKI, N), 0.3, 8.0))
fExt = np.exp(rng.normal(0, SIG_EXT, N))
warm = rng.uniform(600.0, WARM, N)
TREAT_APPLY = {'DMSO': lambda: None, 'CDK4/6i': lambda: rr.__setitem__('kPhRbCd', 0.0),
               'HU': lambda: rr.__setitem__('HU', 1.0)}


def run_cell(i, treat):
    rr.reset(); rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult[i]; rr['p19'] = 0.58 * mult[i]; rr['kSyP21'] = KS0 * mult[i]
    rr['EZH2'] = 6.0 * fExt[i]; rr['EZH2m'] = 0.3 * fExt[i]
    rr['kPhRbCd'] = KPH0; rr['HU'] = 0.0
    g = 1 if rng.random() < DUTY else 0
    t = 0.0; T_END = warm[i] + TREAT; applied = False
    while t < T_END:
        g = 0 if (g and rng.random() < K_OFF * DT) else (1 if (not g and rng.random() < K_ON * DT) else g)
        eta = fExt[i] * (g / DUTY)
        for p in B0: rr[p] = B0[p] * eta
        if (not applied) and t >= warm[i]:
            TREAT_APPLY[treat](); applied = True
        try:
            r = rr.simulate(t, t + DT, 2, selections=SEL)
        except Exception:
            return None
        t += DT
    return {k: float(r[k][-1]) for k in SEL}


res = {t: [] for t in TREAT_APPLY}
for treat in TREAT_APPLY:
    for i in range(N):
        e = run_cell(i, treat)
        if e and np.isfinite(e['EZH2']) and e['EZH2'] > 0 and np.isfinite(e['EZH2m']) and e['EZH2m'] > 0:
            res[treat].append(e)
    print(f'{treat}: {len(res[treat])}/{N}', flush=True)

allpr = np.array([c['pRb'] for t in res for c in res[t]])
PRB_THR = float(np.quantile(allpr[allpr > 0], 0.35))
def phase(c):
    if c['Dna'] >= 0.95: return 'G2'
    if c['vfork'] * c['aRc'] > 0.03: return 'S'
    return 'G0' if c['pRb'] < PRB_THR else 'G1'


def per_phase(cells, key):
    ph = np.array([phase(c) for c in cells]); v = np.log10(np.array([c[key] for c in cells]))
    out = {}
    for p in ['G0', 'G1', 'S', 'G2']:
        x = v[ph == p]
        if len(x) < 8: out[p] = {'n': int(len(x))}; continue
        h, _ = np.histogram(x, bins=BINS, density=True)
        out[p] = {'n': int(len(x)), 'median_log10': float(np.median(x)), 'log10_sd': float(x.std(ddof=1)),
                  'hist_density': h.round(5).tolist()}
    return out, {p: round(float(np.mean(ph == p)), 3) for p in ['G0', 'G1', 'S', 'G2']}


out = {'PRB_THR': PRB_THR, 'bins_log10': BINS.round(4).tolist(),
       'params': dict(N=N, tau=TAU, duty=DUTY, sig_ext=SIG_EXT), 'treatments': {}}
for t, cells in res.items():
    prot, frac = per_phase(cells, 'EZH2'); mrna, _ = per_phase(cells, 'EZH2m')
    out['treatments'][t] = {'phase_frac': frac, 'protein': prot, 'mRNA': mrna}
    def med(d, p): return 10 ** d[p]['median_log10'] if 'median_log10' in d.get(p, {}) else np.nan
    g0p, g0m = med(prot, 'G0'), med(mrna, 'G0')
    boost = f"prot S/G0 {med(prot,'S')/g0p:.1f} G2/G0 {med(prot,'G2')/g0p:.1f} | mRNA S/G0 {med(mrna,'S')/g0m:.1f} G2/G0 {med(mrna,'G2')/g0m:.1f}"
    print(f"{t:8s} frac {frac} | {boost}", flush=True)

SCO = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/model_ezh2_phase_dist.json"
os.makedirs(os.path.dirname(SCO), exist_ok=True)
json.dump(out, open(SCO, 'w')); print('saved', SCO)
