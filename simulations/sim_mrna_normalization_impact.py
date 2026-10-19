"""Does normalizing the CDKI mRNA half-lives (outlier-trimmed central values instead of the single Serdar
2025 chosen arm) change the model? (JP 2026-10-19). Two-step model, MB transient-G0 at the locked burst
config, baseline mRNA half-lives vs the log-IQR-trimmed-geomean 'normalized' set (p21 5.2->2.1h, p27
3.6->1.2h, p18 3.7->2.6h, etc.). Steady state is unchanged by construction (translation kTl=kDe_mRNA), so
only the burst-filtering DYNAMICS can move — p21 is the driver, and its mRNA t1/2 sets how its bursts
propagate, so this isolates whether the normalization matters for the transient-G0 result."""
import sys, os, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('PS_N', '100')); DT = 30.0; WARM = 7000.0; REC = 6000.0; PRB_G0 = 1.5; SIG_FROZ = 0.35
A_ON = 2.9; A_OFF = 0.15; A_ON18 = 1.6; DUTY = 0.28; TAU = 700.0; P21_LVL = 0.03   # locked reference config
BASE = {'kDe_p18m': None}  # placeholder; baseline = builder defaults
NORM = {'kDe_p18m': 0.004443, 'kDe_p19m': 0.004263, 'kDe_p21am': 0.005527, 'kDe_P21m': 0.009317, 'kDe_Skp2m': 0.003032}
rng = np.random.default_rng(21)
froz = np.clip(rng.lognormal(0, SIG_FROZ, N), 0.5, 2.0); warm = rng.uniform(600.0, WARM, N)
rr = te.loada(build_model_v44()); rr.integrator.setValue('maximum_num_steps', 20000000)
KS0 = rr['kSyP21']
# record the builder-default mRNA decay constants (the baseline / chosen Serdar values)
DEFAULTS = {k: rr[k] for k in NORM}


def _sim(t):
    for pts, rt, at in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6), (120, 1e-3, 1e-5)]:
        try:
            rr.integrator.setValue('relative_tolerance', rt); rr.integrator.setValue('absolute_tolerance', at)
            return rr.simulate(t, t + DT, pts, selections=['pRb'])
        except Exception:
            continue
    return None


def run_cell(i, mrna):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * froz[i]; rr['p19'] = 0.58 * froz[i]; rr['kSyp21a'] = P21_LVL * froz[i]; rr['kSyP21'] = KS0 * froz[i]
    for k, v in mrna.items(): rr[k] = v        # mRNA decay set (baseline defaults or normalized)
    k_on = DUTY / TAU; k_off = (1 - DUTY) / TAU
    g21 = rng.random() < DUTY; g18 = rng.random() < DUTY
    def step(t0, t1, rec):
        nonlocal g21, g18
        tr = []; t = t0
        while t < t1:
            if g21 and rng.random() < k_off * DT: g21 = False
            elif (not g21) and rng.random() < k_on * DT: g21 = True
            if g18 and rng.random() < k_off * DT: g18 = False
            elif (not g18) and rng.random() < k_on * DT: g18 = True
            rr['p21a_tx'] = A_ON if g21 else A_OFF; rr['p18_tx'] = A_ON18 if g18 else 0.5
            r = _sim(t)
            if r is None: return None
            if rec: tr.append(float(r['pRb'][-1]))
            t += DT
        return tr
    if step(0, warm[i], False) is None: return None
    pr = step(0, REC, True)
    if pr is None or len(pr) < REC / DT * 0.8: return None
    pr = np.array(pr)
    return pr if np.all(np.isfinite(pr)) else None


def summarize(label, mrna):
    prs = [run_cell(i, mrna) for i in range(N)]; prs = [p for p in prs if p is not None]; n = len(prs)
    occ = float((np.concatenate(prs) < PRB_G0).mean())
    trans = float(np.mean([(p < PRB_G0).any() and (p >= PRB_G0).any() for p in prs]))
    perm = float(np.mean([(p < PRB_G0).all() for p in prs]))
    print(f"  {label:34s}: snapshot-G0 {occ:5.1%} | transient {trans:4.0%} | permanent {perm:4.0%} | n={n}/{N}")
    return occ, trans, perm


print(f"MB transient-G0, locked config (A_ON={A_ON} duty={DUTY} tau={TAU:.0f}), N={N}")
print(f"  baseline mRNA t1/2 (Serdar): p21={0.6931/DEFAULTS['kDe_p21am']/60:.2f}h p27={0.6931/DEFAULTS['kDe_P21m']/60:.2f}h p18={0.6931/DEFAULTS['kDe_p18m']/60:.2f}h")
print(f"  normalized mRNA t1/2 (trimmed): p21={0.6931/NORM['kDe_p21am']/60:.2f}h p27={0.6931/NORM['kDe_P21m']/60:.2f}h p18={0.6931/NORM['kDe_p18m']/60:.2f}h\n")
summarize('BASELINE (Serdar chosen)', DEFAULTS)
summarize('NORMALIZED (outlier-trimmed geomean)', NORM)
