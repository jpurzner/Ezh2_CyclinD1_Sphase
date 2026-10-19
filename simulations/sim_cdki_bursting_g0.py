"""CDKI mRNA bursting -> transient-G0 variance (JP 2026-10-18). FINAL config.

JP's hypothesis: MB has transient-G0 excursions (GNPs do not) because MB has substantially higher CDKI,
and the SOURCE of the variance is bursting of the CDKI mRNAs. This sim implements per-cell transcriptional
bursting (bounded telegraph) on the synthesis of the CDKIs and measures the transient-G0 fraction.

ATTRIBUTION (sim_cdki_bursting_sweep.py): only p21 + p18 bursting cross the committed-Rb threshold; p27 is
buffered (and Skp2-degraded) and Skp2 bursting is anti-arrest -- so the crossing engine is p21 + p18.
Correlating the bursts breaks MB-vs-GNP specificity, so bursts are INDEPENDENT (rho=0); the separation comes
from MB's higher CDKI baseline + lower mitogen drive sitting near threshold while proliferating GNPs (high
SHH, low CDKI) stay clear. Config below gives MB ~28% transient-G0 (0% permanent) vs GNP ~7%; pushing past
~30% loses GNP specificity and starts producing PERMANENT arrest -> the last ~5-10% to the IF 34% Ki67-low
is the deeper (p57/OLIG2) dormancy population, not CDKI bursting (two-quiescence roadmap).

Model naming: P21 = biological p27; p21a = biological p21; p18_prot/p19_prot = INK4.
"""
import sys, os, json, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N        = int(os.environ.get('PS_N', '120'))
DT       = 30.0
WARM     = 7000.0
REC      = float(os.environ.get('REC', '6000'))      # ~100h window -> stable snapshot occupancy
CONTEXT  = os.environ.get('CONTEXT', 'MB')
SHH      = float(os.environ.get('SHH', '0.5' if CONTEXT == 'MB' else '1.0'))
P18_BASE = float(os.environ.get('P18_BASE', '1.73' if CONTEXT == 'MB' else '0.464'))
DUTY     = float(os.environ.get('DUTY', '0.45'))     # ~ fraction of time in a high-CDKI burst
TAU      = float(os.environ.get('TAU', '1500'))      # burst correlation time (min); ~1 cycle so protein tracks
A_ON     = float(os.environ.get('A_ON', '3.0'))      # bounded burst amplitude (ON)
A_OFF    = float(os.environ.get('A_OFF', '0.15'))    #                          (OFF); mRNA dyn range A_ON/A_OFF
P21_LVL  = float(os.environ.get('P21_LVL', '0.015')) # MB kSyp21a target (p21 is the strongest crossing driver)
SIG_FROZ = float(os.environ.get('SIG_FROZ', '0.35'))
PRB_G0   = 1.5                                       # calibrated committed-Rb threshold (validate_v44)
SEED     = int(os.environ.get('SEED', '11'))
GENES    = os.environ.get('BURST', 'p18,p21').split(',')
CAP      = {'p18': 6.0, 'p19': 3.0, 'kSyp21a': 0.12, 'kSyP21': 0.01, 'kSySkp2': 0.2}  # integrability ceilings
GENE_PARAMS = {'p18': ['p18', 'p19'], 'p21': ['kSyp21a'], 'p27': ['kSyP21'], 'skp2': ['kSySkp2']}
SEL = ['pRb', 'EZH2', 'p18_prot', 'p21a', 'P21', 'Dna']

rng = np.random.default_rng(SEED)
K_ON = DUTY / TAU; K_OFF = (1 - DUTY) / TAU
froz = np.clip(rng.lognormal(0, SIG_FROZ, N), 0.5, 2.0)
warm = rng.uniform(600.0, WARM, N)

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('absolute_tolerance', 1e-8)
rr.integrator.setValue('maximum_num_steps', 20000000)
KS0_P27 = rr['kSyP21']; KS0_SKP2 = rr['kSySkp2']


def base_params(i):
    rr.reset(); rr['SHH'] = SHH; rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if CONTEXT == 'GNP':
        rr['Ptch1_copy_number'] = 1.0
        base = {'p18': P18_BASE*froz[i], 'p19': 0.36*froz[i], 'kSyP21': KS0_P27*froz[i],
                'kSySkp2': KS0_SKP2, 'kSyp21a': P21_LVL*0.3*froz[i]}
    else:
        rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        base = {'p18': P18_BASE*froz[i], 'p19': 0.58*froz[i], 'kSyP21': KS0_P27*froz[i],
                'kSySkp2': KS0_SKP2, 'kSyp21a': P21_LVL*froz[i]}
    for k, v in base.items(): rr[k] = v
    return base


def _sim(t):
    for rtol, atol in [(1e-6, 1e-8), (1e-5, 1e-7), (1e-4, 1e-6)]:
        try:
            rr.integrator.setValue('relative_tolerance', rtol); rr.integrator.setValue('absolute_tolerance', atol)
            return rr.simulate(t, t+DT, 2, selections=SEL)
        except Exception:
            continue
    return None


def run_cell(i):
    base = base_params(i)
    g = {gene: (rng.random() < DUTY) for gene in GENES}
    def step(t0, t1, rec):
        traj = {k: [] for k in SEL}
        t = t0
        while t < t1:
            for gene in GENES:
                if g[gene] and rng.random() < K_OFF*DT: g[gene] = False
                elif (not g[gene]) and rng.random() < K_ON*DT: g[gene] = True
                eta = A_ON if g[gene] else A_OFF
                for p in GENE_PARAMS[gene]: rr[p] = min(base[p]*eta, CAP.get(p, 1e9))
            r = _sim(t)
            if r is None: return None
            if rec:
                for k in SEL: traj[k].append(float(r[k][-1]))
            t += DT
        return traj
    if step(0, warm[i], False) is None: return None
    tr = step(0, REC, True)
    if tr is None or len(tr['pRb']) < REC/DT*0.8: return None
    out = {k: np.array(v) for k, v in tr.items()}
    if not np.all(np.isfinite(out['pRb'])): return None
    return out


cells = [c for c in (run_cell(i) for i in range(N)) if c is not None]
n = len(cells)
allpr = np.concatenate([c['pRb'] for c in cells])
occ = (allpr < PRB_G0).mean()
transient = np.mean([(c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any() for c in cells])
perm = np.mean([(c['pRb'] < PRB_G0).all() for c in cells])
cyc = np.mean([(c['pRb'] >= PRB_G0).all() for c in cells])
def tcv(key): return float(np.median([c[key].std()/(c[key].mean()+1e-9) for c in cells]))
ez = np.concatenate([c['EZH2'] for c in cells]); q = allpr < PRB_G0
ezq = float(np.median(ez[q])) if q.any() else float('nan'); ezc = float(np.median(ez[~q])) if (~q).any() else float('nan')

print(f"=== CONTEXT={CONTEXT} SHH={SHH} P18={P18_BASE} n={n}/{N}  burst={GENES} duty={DUTY} tau={TAU:.0f} A_ON={A_ON} A_OFF={A_OFF} ===")
print(f"  snapshot G0 occupancy (pRb<{PRB_G0}) : {occ:.1%}")
print(f"  pure-cycling {cyc:.0%} | transient-G0 (dip&recover) {transient:.0%} | permanent-G0 {perm:.0%}")
print(f"  CDKI temporal CV (within-cell): p18_prot {tcv('p18_prot'):.2f}  p21a {tcv('p21a'):.2f}  p27(P21) {tcv('P21'):.2f}")
print(f"  EZH2 quiescent {ezq:.2f} vs cycling {ezc:.2f}  (fold {ezq/ezc:.2f})")

# save summary + a few example trajectories (subsampled) for the figure
SC = os.environ.get('PS_OUT') or "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
ex_idx = [i for i, c in enumerate(cells) if (c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any()][:4] or list(range(min(4, n)))
sub = slice(None, None, 2)
out = {'context': CONTEXT, 'shh': SHH, 'p18': P18_BASE, 'n': n, 'config': {'duty': DUTY, 'tau': TAU, 'A_ON': A_ON, 'A_OFF': A_OFF, 'genes': GENES},
       'occ': occ, 'transient': float(transient), 'perm': float(perm), 'cyc': float(cyc),
       'ezh2_q': ezq, 'ezh2_c': ezc, 'cv': {'p18': tcv('p18_prot'), 'p21': tcv('p21a'), 'p27': tcv('P21')},
       'dt': DT, 'prb_thr': PRB_G0,
       'examples': [{k: cells[i][k][sub].round(4).tolist() for k in ['pRb', 'p18_prot', 'p21a', 'EZH2']} for i in ex_idx],
       'prb_pool': allpr[::7].round(3).tolist(), 'ezh2_pool': ez[::7].round(3).tolist(), 'q_pool': q[::7].tolist()}
os.makedirs(SC, exist_ok=True)
json.dump(out, open(f"{SC}/cdki_burst_{CONTEXT}.json", 'w'))
print(f"  saved {SC}/cdki_burst_{CONTEXT}.json")
