"""Quick sweep on the TWO-STEP model (with_cdki_mrna) to land MB ~34% transient-G0 with LOW permanent and
clean GNP specificity. One build; loops (A_ON, duty, tau). Integration is now stable (mRNA low-pass), so no
mortality concern. Bursts the transcription knobs p21a_tx (+ mild p18_tx)."""
import sys, os, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('PS_N', '40')); DT = 30.0; WARM = 7000.0; REC = float(os.environ.get('REC', '5000'))
P21_LVL = 0.03; A_OFF = 0.15; A_ON18 = 1.6; GNP_P21_FAC = 0.1; SIG_FROZ = 0.35; PRB_G0 = 1.5
rng0 = np.random.default_rng(11)
froz = np.clip(rng0.lognormal(0, SIG_FROZ, N), 0.5, 2.0); warm = rng0.uniform(600.0, WARM, N)
rr = te.loada(build_model_v44(with_cdki_mrna=True)); rr.integrator.setValue('maximum_num_steps', 20000000)
KS0 = rr['kSyP21']; SEL = ['pRb']


def _sim(t):
    for pts, rt, at in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6)]:
        try:
            rr.integrator.setValue('relative_tolerance', rt); rr.integrator.setValue('absolute_tolerance', at)
            return rr.simulate(t, t + DT, pts, selections=SEL)
        except Exception:
            continue
    return None


def run_cell(ctx, i, a_on, duty, tau):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if ctx == 'GNP':
        rr['SHH'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        rr['p18'] = 0.464 * froz[i]; rr['p19'] = 0.36 * froz[i]; rr['kSyp21a'] = P21_LVL * GNP_P21_FAC * froz[i]
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = 1.73 * froz[i]; rr['p19'] = 0.58 * froz[i]; rr['kSyp21a'] = P21_LVL * froz[i]
    rr['kSyP21'] = KS0 * froz[i]
    k_on = duty / tau; k_off = (1 - duty) / tau
    g21 = rng0.random() < duty; g18 = rng0.random() < duty

    def step(t0, t1, rec):
        nonlocal g21, g18
        tr = []; t = t0
        while t < t1:
            if g21 and rng0.random() < k_off * DT: g21 = False
            elif (not g21) and rng0.random() < k_on * DT: g21 = True
            if g18 and rng0.random() < k_off * DT: g18 = False
            elif (not g18) and rng0.random() < k_on * DT: g18 = True
            rr['p21a_tx'] = a_on if g21 else A_OFF
            rr['p18_tx'] = A_ON18 if g18 else 0.5
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


def summ(ctx, a_on, duty, tau):
    prs = [run_cell(ctx, i, a_on, duty, tau) for i in range(N)]
    prs = [p for p in prs if p is not None]; n = len(prs)
    if not n: return None
    occ = (np.concatenate(prs) < PRB_G0).mean()
    trans = np.mean([(p < PRB_G0).any() and (p >= PRB_G0).any() for p in prs])
    perm = np.mean([(p < PRB_G0).all() for p in prs])
    return occ, trans, perm, n


print(f"N={N} REC={REC:.0f} TWO-STEP (burst transcription -> mRNA -> protein)\n")
print(f"{'ctx':4s} {'A_ON':>5s} {'duty':>5s} {'tau':>5s} | {'occ%':>5s} {'trans%':>7s} {'perm%':>6s}  n")
for a_on in [2.0, 2.5, 3.0]:
    for duty in [0.30, 0.40]:
        for tau in [700, 1100]:
            r = summ('MB', a_on, duty, tau)
            if r: print(f"{'MB':4s} {a_on:5.1f} {duty:5.2f} {tau:5.0f} | {r[0]*100:5.1f} {r[1]*100:7.0f} {r[2]*100:6.0f}  {r[3]:3d}")
    print()
# GNP specificity at a mid config
for a_on, duty, tau in [(2.5, 0.40, 1100), (3.0, 0.40, 1100)]:
    r = summ('GNP', a_on, duty, tau)
    if r: print(f"{'GNP':4s} {a_on:5.1f} {duty:5.2f} {tau:5.0f} | {r[0]*100:5.1f} {r[1]*100:7.0f} {r[2]*100:6.0f}  {r[3]:3d}")
