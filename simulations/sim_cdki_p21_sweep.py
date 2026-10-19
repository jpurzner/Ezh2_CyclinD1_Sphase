"""p21-dominant bursting sweep (JP 2026-10-18): p21a is the strongest transient-G0 driver AND numerically
tiny, so bursting it (vs the stiff p18 level) should integrate cleanly (n~N) while still crossing the
committed-Rb threshold. Tests whether p21 bursting alone gives MB ~25-34% transient-G0 with GNP (high SHH,
low p21) staying clear. p18 kept at a modest fixed burst (A_ON18) to add some INK4 variance without stiffness."""
import sys, os, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('PS_N', '48')); DT = 30.0; WARM = 7000.0; REC = float(os.environ.get('REC', '6000'))
SIG_FROZ = 0.35; PRB_G0 = 1.5
rng0 = np.random.default_rng(11)
froz = np.clip(rng0.lognormal(0, SIG_FROZ, N), 0.5, 2.0)
warm = rng0.uniform(600.0, WARM, N)
rr = te.loada(build_model_v44())
rr.integrator.setValue('maximum_num_steps', 20000000)
KS0_P27 = rr['kSyP21']
SEL = ['pRb']


def _sim(t):
    for rtol, atol in [(1e-6, 1e-8), (1e-5, 1e-7), (1e-4, 1e-6)]:
        try:
            rr.integrator.setValue('relative_tolerance', rtol); rr.integrator.setValue('absolute_tolerance', atol)
            return rr.simulate(t, t+DT, 2, selections=SEL)
        except Exception:
            continue
    return None


def run_cell(ctx, i, p21_lvl, a_on, duty, tau, a_on18=1.6):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if ctx == 'GNP':
        rr['SHH'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        p18b, p19b, p21b = 0.464*froz[i], 0.36*froz[i], p21_lvl*0.3*froz[i]
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        p18b, p19b, p21b = 1.73*froz[i], 0.58*froz[i], p21_lvl*froz[i]
    rr['kSyP21'] = KS0_P27*froz[i]
    k_on = duty/tau; k_off = (1-duty)/tau
    g21 = rng0.random() < duty; g18 = rng0.random() < duty
    def step(t0, t1, rec):
        nonlocal g21, g18
        tr = []; t = t0
        while t < t1:
            if g21 and rng0.random() < k_off*DT: g21 = False
            elif (not g21) and rng0.random() < k_on*DT: g21 = True
            if g18 and rng0.random() < k_off*DT: g18 = False
            elif (not g18) and rng0.random() < k_on*DT: g18 = True
            rr['kSyp21a'] = p21b * (a_on if g21 else 0.15)
            e18 = a_on18 if g18 else 0.5
            rr['p18'] = min(p18b*e18, 5.0); rr['p19'] = min(p19b*e18, 3.0)
            r = _sim(t)
            if r is None: return None
            if rec: tr.append(float(r['pRb'][-1]))
            t += DT
        return tr
    if step(0, warm[i], False) is None: return None
    pr = step(0, REC, True)
    if pr is None or len(pr) < REC/DT*0.8: return None
    pr = np.array(pr)
    return pr if np.all(np.isfinite(pr)) else None


def summ(ctx, p21_lvl, a_on, duty, tau):
    prs = [run_cell(ctx, i, p21_lvl, a_on, duty, tau) for i in range(N)]
    prs = [p for p in prs if p is not None]; n = len(prs)
    if not n: return None
    occ = (np.concatenate(prs) < PRB_G0).mean()
    trans = np.mean([(p < PRB_G0).any() and (p >= PRB_G0).any() for p in prs])
    perm = np.mean([(p < PRB_G0).all() for p in prs])
    return occ, trans, perm, n


print(f"N={N} REC={REC:.0f}  p21-dominant bursting (p18 modest), genes=p21(+p18 mild)\n")
print(f"{'ctx':4s} {'p21L':>5s} {'A_ON':>5s} {'duty':>5s} {'tau':>5s} | {'occ%':>5s} {'trans%':>7s} {'perm%':>6s}  n")
for p21_lvl in [0.03, 0.06]:
    for a_on in [4.0, 6.0]:
        for duty in [0.40, 0.50]:
            r = summ('MB', p21_lvl, a_on, duty, 1500)
            if r: print(f"{'MB':4s} {p21_lvl:5.3f} {a_on:5.1f} {duty:5.2f} {1500:5.0f} | {r[0]*100:5.1f} {r[1]*100:7.0f} {r[2]*100:6.0f}  {r[3]:3d}")
    print()
# GNP control at the two strongest MB configs
for p21_lvl, a_on, duty in [(0.06, 6.0, 0.50), (0.06, 6.0, 0.40)]:
    r = summ('GNP', p21_lvl, a_on, duty, 1500)
    if r: print(f"{'GNP':4s} {p21_lvl:5.3f} {a_on:5.1f} {duty:5.2f} {1500:5.0f} | {r[0]*100:5.1f} {r[1]*100:7.0f} {r[2]*100:6.0f}  {r[3]:3d}")
