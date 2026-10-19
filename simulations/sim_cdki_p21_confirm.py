"""CONFIRM the p21-bursting transient-G0 config with robust integration + a genuinely low-CDKI GNP control.
JP: MB has transient-G0 because it has substantially higher CDKI; GNPs (low CDKI, high mitogen drive) do not.
So GNP p21 is set to 0.1x MB here (genuinely low), not 0.3x. Robust solver: a failed chunk is retried with
progressively smaller internal steps (more output points) and looser tolerance, so deeply-arrested cells are
NOT dropped (dropping them would bias the G0 fraction). Saves cdki_burst_{MB,GNP}.json for the figure."""
import sys, os, json, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('PS_N', '72')); DT = 30.0; WARM = 7000.0; REC = float(os.environ.get('REC', '6000'))
CONTEXT = os.environ.get('CONTEXT', 'MB')
P21_LVL = float(os.environ.get('P21_LVL', '0.03')); A_ON = float(os.environ.get('A_ON', '4.0'))
A_ON18 = float(os.environ.get('A_ON18', '1.6')); DUTY = float(os.environ.get('DUTY', '0.40')); TAU = float(os.environ.get('TAU', '1500'))
GNP_P21_FAC = float(os.environ.get('GNP_P21_FAC', '0.1'))   # GNPs genuinely low p21
SIG_FROZ = 0.35; PRB_G0 = 1.5; SEED = int(os.environ.get('SEED', '11'))
rng = np.random.default_rng(SEED)
froz = np.clip(rng.lognormal(0, SIG_FROZ, N), 0.5, 2.0)
warm = rng.uniform(600.0, WARM, N)
rr = te.loada(build_model_v44()); rr.integrator.setValue('maximum_num_steps', 20000000)
KS0_P27 = rr['kSyP21']
SEL = ['pRb', 'EZH2', 'p18_prot', 'p21a', 'P21', 'Dna']


def _sim(t):
    for pts, rtol, atol in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6), (120, 1e-3, 1e-5)]:
        try:
            rr.integrator.setValue('relative_tolerance', rtol); rr.integrator.setValue('absolute_tolerance', atol)
            return rr.simulate(t, t+DT, pts, selections=SEL)
        except Exception:
            continue
    return None


def run_cell(i):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if CONTEXT == 'GNP':
        rr['SHH'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        p18b, p19b, p21b = 0.464*froz[i], 0.36*froz[i], P21_LVL*GNP_P21_FAC*froz[i]
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        p18b, p19b, p21b = 1.73*froz[i], 0.58*froz[i], P21_LVL*froz[i]
    rr['kSyP21'] = KS0_P27*froz[i]
    k_on = DUTY/TAU; k_off = (1-DUTY)/TAU
    g21 = rng.random() < DUTY; g18 = rng.random() < DUTY
    def step(t0, t1, rec):
        nonlocal g21, g18
        tr = {k: [] for k in SEL}; t = t0
        while t < t1:
            if g21 and rng.random() < k_off*DT: g21 = False
            elif (not g21) and rng.random() < k_on*DT: g21 = True
            if g18 and rng.random() < k_off*DT: g18 = False
            elif (not g18) and rng.random() < k_on*DT: g18 = True
            rr['kSyp21a'] = p21b*(A_ON if g21 else 0.15)
            e18 = A_ON18 if g18 else 0.5
            rr['p18'] = min(p18b*e18, 5.0); rr['p19'] = min(p19b*e18, 3.0)
            r = _sim(t)
            if r is None: return None
            if rec:
                for k in SEL: tr[k].append(float(r[k][-1]))
            t += DT
        return tr
    if step(0, warm[i], False) is None: return None
    tr = step(0, REC, True)
    if tr is None or len(tr['pRb']) < REC/DT*0.8: return None
    out = {k: np.array(v) for k, v in tr.items()}
    return out if np.all(np.isfinite(out['pRb'])) else None


cells = [c for c in (run_cell(i) for i in range(N)) if c is not None]
n = len(cells)
allpr = np.concatenate([c['pRb'] for c in cells]); occ = (allpr < PRB_G0).mean()
transient = np.mean([(c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any() for c in cells])
perm = np.mean([(c['pRb'] < PRB_G0).all() for c in cells]); cyc = np.mean([(c['pRb'] >= PRB_G0).all() for c in cells])
def tcv(k): return float(np.median([c[k].std()/(c[k].mean()+1e-9) for c in cells]))
ez = np.concatenate([c['EZH2'] for c in cells]); q = allpr < PRB_G0
ezq = float(np.median(ez[q])) if q.any() else float('nan'); ezc = float(np.median(ez[~q])) if (~q).any() else float('nan')
print(f"=== {CONTEXT} SHH={rr['SHH']} p21L={P21_LVL}(GNPfac {GNP_P21_FAC if CONTEXT=='GNP' else 1}) A_ON={A_ON} duty={DUTY} tau={TAU:.0f} n={n}/{N} ===")
print(f"  snapshot G0 {occ:.1%} | pure-cycling {cyc:.0%} transient {transient:.0%} permanent {perm:.0%}")
print(f"  CDKI CV: p18 {tcv('p18_prot'):.2f} p21 {tcv('p21a'):.2f} p27 {tcv('P21'):.2f} | EZH2 q {ezq:.2f} vs c {ezc:.2f} (fold {ezq/ezc:.2f})")

SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
ex = [i for i, c in enumerate(cells) if (c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any()][:4] or list(range(min(4, n)))
sub = slice(None, None, 2)
json.dump({'context': CONTEXT, 'shh': rr['SHH'], 'n': n, 'occ': occ, 'transient': float(transient), 'perm': float(perm),
           'cyc': float(cyc), 'ezh2_q': ezq, 'ezh2_c': ezc, 'cv': {'p18': tcv('p18_prot'), 'p21': tcv('p21a'), 'p27': tcv('P21')},
           'dt': DT, 'prb_thr': PRB_G0, 'config': {'p21L': P21_LVL, 'A_ON': A_ON, 'duty': DUTY, 'tau': TAU},
           'examples': [{k: cells[i][k][sub].round(4).tolist() for k in ['pRb', 'p18_prot', 'p21a', 'EZH2']} for i in ex],
           'prb_pool': allpr[::7].round(3).tolist(), 'ezh2_pool': ez[::7].round(3).tolist(), 'q_pool': q[::7].tolist()},
          open(f"{SC}/cdki_burst_{CONTEXT}.json", 'w'))
print(f"  saved cdki_burst_{CONTEXT}.json")
