"""Transient-G0 via CDKI mRNA bursting on the TWO-STEP model (JP 2026-10-18).

Uses build_model_v44(with_cdki_mrna=True): transcription -> mRNA species -> translation -> protein, with
independent measured mRNA + protein half-lives. The burst is now applied to the TRANSCRIPTION knobs
(p21a_tx, p18_tx; default 1), so it flows transcription -> mRNA (t1/2 ~5.2h p21 / 3.7h p18) -> protein and
NEVER leaks into mu_eff (which reads the static p18 level). This replaces the earlier param-bursting
(sim_cdki_p21_confirm.py), where the burst hit the synthesis parameter directly and leaked to growth.

p21 is the dominant crossing driver (short protein tracks its mRNA); p18 secondary. GNP is genuinely
low-CDKI + high-SHH (proliferating) -> stays cycling. Records mRNA AND protein so we can see the two-stage
filtering. Bonus: the smooth mRNA low-pass REDUCES the deep-arrest integrator mortality (GNP 0%, MB ~50%->~39%
with the baked config) -- it does NOT eliminate it; residual MB failures remain concentrated in deep arrest.
"""
import sys, os, json, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44
from src.bursting import CdkiBurstDriver        # canonical driver; burst SPEC comes from the model

CACHE = os.environ.get('PS_OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '_cache'))
os.makedirs(CACHE, exist_ok=True)
N = int(os.environ.get('PS_N', '72')); DT = 30.0; WARM = 7000.0; REC = float(os.environ.get('REC', '6000'))
CONTEXT = os.environ.get('CONTEXT', 'MB')
P21_LVL = float(os.environ.get('P21_LVL', '0.03'))     # baseline kSyp21a (transcription); burst multiplies on top
GNP_P21_FAC = float(os.environ.get('GNP_P21_FAC', '0.1'))
SIG_FROZ = 0.35; PRB_G0 = 1.5; SEED = int(os.environ.get('SEED', '11'))
rng = np.random.default_rng(SEED)
froz = np.clip(rng.lognormal(0, SIG_FROZ, N), 0.5, 2.0)
warm = rng.uniform(600.0, WARM, N)

rr = te.loada(build_model_v44())                # with_cdki_mrna + with_cdki_bursting are baked defaults
rr.integrator.setValue('maximum_num_steps', 20000000)
# Burst configuration is INHERITED from the model (burst_duty/burst_tau/burst_a_on/...); env vars override.
_ovr = {k: float(os.environ[e]) for k, e in
        [('burst_a_on', 'A_ON'), ('burst_a_off', 'A_OFF'), ('burst_a_on_p18', 'A_ON18'),
         ('burst_duty', 'DUTY'), ('burst_tau', 'TAU')] if e in os.environ}
burst = CdkiBurstDriver(rr, seed=SEED, p21_only=(os.environ.get('P21_ONLY') == '1'), **_ovr)
A_ON, A_OFF, DUTY, TAU = burst.amp['p21a_tx'][0], burst.amp['p21a_tx'][1], burst.duty, burst.tau
A_ON18 = burst.amp['p18_tx'][0]
KS0_P27 = rr['kSyP21']
SEL = ['pRb', 'EZH2', 'p18_prot', 'p21a', 'p21a_m', 'P21', 'Dna']


def _sim(t):
    for pts, rtol, atol in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6), (120, 1e-3, 1e-5)]:
        try:
            rr.integrator.setValue('relative_tolerance', rtol); rr.integrator.setValue('absolute_tolerance', atol)
            return rr.simulate(t, t + DT, pts, selections=SEL)
        except Exception:
            continue
    return None


def run_cell(i):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if CONTEXT == 'GNP':
        rr['SHH'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        rr['p18'] = 0.464 * froz[i]; rr['p19'] = 0.36 * froz[i]
        rr['kSyp21a'] = P21_LVL * GNP_P21_FAC * froz[i]
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = 1.73 * froz[i]; rr['p19'] = 0.58 * froz[i]
        rr['kSyp21a'] = P21_LVL * froz[i]
    rr['kSyP21'] = KS0_P27 * froz[i]
    burst.reset()                                   # independent telegraph realization for this cell

    def step(t0, t1, rec):
        tr = {k: [] for k in SEL}; t = t0
        while t < t1:
            burst.step(DT)                          # drives p21a_tx / p18_tx (-> mRNA -> protein)
            r = _sim(t)
            if r is None: return None
            if rec:
                for k in SEL: tr[k].append(float(r[k][-1]))
            t += DT
        return tr
    if step(0, warm[i], False) is None: return None
    tr = step(0, REC, True)
    if tr is None or len(tr['pRb']) < REC / DT * 0.8: return None
    out = {k: np.array(v) for k, v in tr.items()}
    return out if np.all(np.isfinite(out['pRb'])) else None


cells = [c for c in (run_cell(i) for i in range(N)) if c is not None]
n = len(cells)
allpr = np.concatenate([c['pRb'] for c in cells]); occ = (allpr < PRB_G0).mean()
transient = np.mean([(c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any() for c in cells])
perm = np.mean([(c['pRb'] < PRB_G0).all() for c in cells]); cyc = np.mean([(c['pRb'] >= PRB_G0).all() for c in cells])
def tcv(k): return float(np.median([c[k].std() / (c[k].mean() + 1e-9) for c in cells]))
ez = np.concatenate([c['EZH2'] for c in cells]); q = allpr < PRB_G0
ezq = float(np.median(ez[q])) if q.any() else float('nan'); ezc = float(np.median(ez[~q])) if (~q).any() else float('nan')
print(f"=== TWO-STEP {CONTEXT} SHH={rr['SHH']} p21L={P21_LVL} A_ON={A_ON} duty={DUTY} tau={TAU:.0f} n={n}/{N} ===")
print(f"  snapshot G0 {occ:.1%} | pure-cycling {cyc:.0%} transient {transient:.0%} permanent {perm:.0%}  [mortality {100*(1-n/N):.0f}%]")
print(f"  temporal CV: p21a_mRNA {tcv('p21a_m'):.2f} -> p21a_protein {tcv('p21a'):.2f} | p18_prot {tcv('p18_prot'):.2f} | p27 {tcv('P21'):.2f}")
print(f"  EZH2 quiescent {ezq:.2f} vs cycling {ezc:.2f} (fold {ezq/ezc:.2f})")

SC = CACHE
ex = [i for i, c in enumerate(cells) if (c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any()][:4] or list(range(min(4, n)))
sub = slice(None, None, 2)
json.dump({'context': CONTEXT, 'shh': rr['SHH'], 'n': n, 'N': N, 'occ': occ, 'transient': float(transient),
           'perm': float(perm), 'cyc': float(cyc), 'ezh2_q': ezq, 'ezh2_c': ezc,
           'cv': {'p21a_m': tcv('p21a_m'), 'p21a': tcv('p21a'), 'p18_prot': tcv('p18_prot'), 'p27': tcv('P21')},
           'dt': DT, 'prb_thr': PRB_G0, 'config': {'p21L': P21_LVL, 'A_ON': A_ON, 'duty': DUTY, 'tau': TAU, 'two_step': True},
           'examples': [{k: cells[i][k][sub].round(4).tolist() for k in ['pRb', 'p21a_m', 'p21a', 'p18_prot', 'EZH2']} for i in ex],
           'prb_pool': allpr[::7].round(3).tolist(), 'ezh2_pool': ez[::7].round(3).tolist(), 'q_pool': q[::7].tolist()},
          open(f"{SC}/cdki_mrna_burst_{CONTEXT}.json", 'w'))
print(f"  saved cdki_mrna_burst_{CONTEXT}.json")
