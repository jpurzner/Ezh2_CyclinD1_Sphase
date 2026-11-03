"""Sweep: how high must MB CDKI sit, and how deep must the bursts be, for MB to show ~34% transient-G0
while GNP stays cycling? One model build; loops (P18_MB x DUTY x TAU). Same bursting engine as
sim_cdki_bursting_g0.py (telegraph p18/p19 level, kSyp21a, kSyP21, kSySkp2 independently per cell)."""
import sys, os, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('PS_N', '40')); DT = 30.0; WARM = 7000.0; REC = float(os.environ.get('REC', '6000'))
SIG_FROZ = 0.35; P21_LVL = float(os.environ.get('P21_LVL', '0.015')); PRB_G0 = 1.5
A_ON = float(os.environ.get('A_ON', '2.5'))       # bounded burst: synth x A_ON when ON
A_OFF = float(os.environ.get('A_OFF', '0.2'))     #               x A_OFF when OFF  (dynamic range A_ON/A_OFF)
SHH_GNP = float(os.environ.get('SHH_GNP', '1.0')) # proliferating-GNP mitogen level (vs MB partial 0.5)
GENES = ['p18', 'p21', 'p27', 'skp2']
GENE_PARAMS = {'p18': ['p18', 'p19'], 'p21': ['kSyp21a'], 'p27': ['kSyP21'], 'skp2': ['kSySkp2']}
SEL = ['pRb']
rng0 = np.random.default_rng(11)
froz = np.clip(rng0.lognormal(0, SIG_FROZ, N), 0.5, 2.0)   # cap tail so base*A_ON stays integrable
warm = rng0.uniform(600.0, WARM, N)

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('absolute_tolerance', 1e-8)
rr.integrator.setValue('maximum_num_steps', 20000000)
KS0_P27 = rr['kSyP21']; KS0_SKP2 = rr['kSySkp2']


def base_params(ctx, i, p18_mb):
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    if ctx == 'GNP':
        rr['SHH'] = SHH_GNP; rr['Ptch1_copy_number'] = 1.0
        base = {'p18': 0.464*froz[i], 'p19': 0.36*froz[i], 'kSyP21': KS0_P27*froz[i],
                'kSySkp2': KS0_SKP2, 'kSyp21a': P21_LVL*0.3*froz[i]}
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        base = {'p18': p18_mb*froz[i], 'p19': 0.58*froz[i], 'kSyP21': KS0_P27*froz[i],
                'kSySkp2': KS0_SKP2, 'kSyp21a': P21_LVL*froz[i]}
    for k, v in base.items(): rr[k] = v
    return base


def _sim(t, dt):
    """one chunk with a tolerance-retry fallback for stiff bursts"""
    try:
        return rr.simulate(t, t+dt, 2, selections=SEL)
    except Exception:
        rr.integrator.setValue('relative_tolerance', 1e-5); rr.integrator.setValue('absolute_tolerance', 1e-7)
        r = rr.simulate(t, t+dt, 2, selections=SEL)
        rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('absolute_tolerance', 1e-8)
        return r


def run_cell(ctx, i, p18_mb, duty, tau, rng, rho=0.0, genes=None):
    genes = genes or GENES
    base = base_params(ctx, i, p18_mb)
    k_on = duty/tau; k_off = (1-duty)/tau
    g = {gene: (rng.random() < duty) for gene in genes}
    def step(t0, t1, rec=False):
        traj = []
        t = t0
        while t < t1:
            for gene in genes:
                if g[gene] and rng.random() < k_off*DT: g[gene] = False
                elif (not g[gene]) and rng.random() < k_on*DT: g[gene] = True
                eta = A_ON if g[gene] else A_OFF          # bounded two-state burst (dynamic range A_ON/A_OFF)
                for p in GENE_PARAMS[gene]: rr[p] = base[p]*eta
            try:
                r = _sim(t, DT)
            except Exception:
                return traj if rec else None
            if rec: traj.append(float(r['pRb'][-1]))
            t += DT
        return traj
    w = step(0, warm[i])
    if w is None: return None
    pr = np.array(step(0, REC, rec=True))
    if len(pr) < REC/DT*0.6 or not np.all(np.isfinite(pr)): return None
    return pr


def summarize(ctx, p18_mb, duty, tau, rho=0.0, genes=None):
    rng = np.random.default_rng(100 + int(p18_mb*10) + int(duty*100) + int(tau) + int(rho*100))
    prs = [run_cell(ctx, i, p18_mb, duty, tau, rng, rho, genes) for i in range(N)]
    prs = [p for p in prs if p is not None]; n = len(prs)
    if n == 0: return None
    occ = np.concatenate(prs); occ = (occ < PRB_G0).mean()
    trans = np.mean([(p < PRB_G0).any() and (p >= PRB_G0).any() for p in prs])
    perm = np.mean([(p < PRB_G0).all() for p in prs])
    return occ, trans, perm, n


def row(ctx, p18, duty, tau, rho, genes=None, tag=''):
    r = summarize(ctx, p18, duty, tau, rho, genes)
    if r: print(f"{ctx:4s} {p18:5.2f} {duty:5.2f} {tau:5.0f} {rho:4.1f} | {r[0]*100:5.1f} {r[1]*100:7.0f} {r[2]*100:6.0f}  {r[3]:3d}  {tag}")

print(f"N={N} REC={REC:.0f} P21_LVL={P21_LVL} A_ON={A_ON} A_OFF={A_OFF} SHH_GNP={SHH_GNP}")
print("bounded bursts, genes=p18+p21, rho=0  (occ=snapshot G0; trans=dip&recover; perm=stuck)\n")
print(f"{'ctx':4s} {'P18':>5s} {'duty':>5s} {'tau':>5s} | {'occ%':>5s} {'trans%':>7s} {'perm%':>6s}  n   note")
PG = ['p18', 'p21']
# duty ~= target quiescent fraction (for slow bursts where protein tracks the burst state); push toward 34%.
for duty in [0.45, 0.55, 0.65]:
    row('MB', 1.73, duty, 1500, 0.0, genes=PG, tag='MB (SHH0.5)')
print()
for duty in [0.45, 0.55, 0.65]:
    row('GNP', 0.464, duty, 1500, 0.0, genes=PG, tag=f'GNP (SHH{SHH_GNP})')
