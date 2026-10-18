"""Transcript-MATCHED withdrawal experiment (JP 2026-07-24): give ON / OFF / EZH2i enough Shh to reach the SAME
Ccnd1 transcript level, then withdraw Shh at MULTIPLE decline rates, over a POPULATION of cells. Question: if they
start at the same CyclinD1, does the withdrawal look the same — or does the EZH2 regulatory structure / Shh headroom
make them diverge?

Matched (target Cd_mRNA = ON@0.45 ≈ 1.61):  ON needs SHH 0.45, EZH2i 0.177, OFF 0.082  (same transcript/protein/rate).
Withdrawal: each cell's Shh mean ramps LINEARLY from its matched level to 0 over duration D (the decline rate),
with persistent cell-to-cell noise. Multiple D = multiple decline rates. Caches to _cache.npz.
Run: ./venv/bin/python simulations/sim_matched_withdrawal.py   (figure: fig_matched_withdrawal.py)"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MODEL = build_model_v44()
def divs(dna): dna = np.asarray(dna); return np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
def cond(name): return {} if name == 'ON' else ({'f0_prc2': 1.0} if name == 'OFF' else {'EZH2i': 1.0})
def sim_robust(rr, t0, t1, n, sel):
    for atol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        try: rr.integrator.setValue('absolute_tolerance', atol); return rr.simulate(t0, t1, n, selections=sel)
        except Exception: rr.integrator.setValue('maximum_time_step', 1.5)
    return None

MATCHED = {'ON': 0.450, 'OFF': 0.082, 'EZH2i': 0.177}   # Shh giving Cd_mRNA ≈ 1.61 in each condition
THRESH  = {'ON': 0.176, 'OFF': 0.012, 'EZH2i': 0.059}
RATES = [12.0, 48.0, 120.0]                              # withdrawal ramp durations (h): fast / medium / slow
NCELL = 6
TWARM = 150.0; TPOST = 260.0                             # warm (desynchronize) then withdraw + observe
GRID = np.arange(-30.0, TPOST + 0.1, 1.0)               # common time grid, 0 = withdrawal onset

def ou_noise(n, dt, tau=6.0, cv=0.10, seed=0):
    rng = np.random.default_rng(seed); s = cv * np.sqrt(2 / tau); x = np.zeros(n)
    for i in range(1, n): x[i] = x[i-1] - x[i-1] * dt / tau + s * np.sqrt(dt) * rng.normal()
    return x

def run_cell(c, D, seed):
    """warm at matched Shh (noisy) to desynchronize, then ramp mean -> 0 over D; return Cd_mRNA(GRID), Cd(GRID), div times (rel to onset)."""
    mean0 = MATCHED[c]
    tfull = np.arange(-TWARM, TPOST, 0.1)
    noise = ou_noise(len(tfull), 0.1, seed=seed)
    def shh_of(th):
        base = mean0 if th < 0 else (mean0 * max(0.0, 1 - th / D) if th < D else 0.0)
        i = int(np.clip((th + TWARM) / 0.1, 0, len(noise) - 1))
        return max(0.0, base * (1 + noise[i]))
    rr = te.loada(MODEL); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    for k, v in cond(c).items(): rr[k] = v
    rr.integrator.setValue('absolute_tolerance', 1e-8); rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('maximum_num_steps', 4000000)
    T_, MR, CD, DNA = [], [], [], []
    t = -TWARM * UPH; chunk = 8.0; Tend = TPOST * UPH
    while t < Tend - 1e-6:
        rr['SHH'] = shh_of(t / UPH)
        r = sim_robust(rr, t, t + chunk, 4, ['time', 'Cd_mRNA', 'Cd', 'Dna'])
        if r is not None:
            T_ += list(np.asarray(r['time'])[1:]); MR += list(np.asarray(r['Cd_mRNA'])[1:]); CD += list(np.asarray(r['Cd'])[1:]); DNA += list(np.asarray(r['Dna'])[1:])
        t += chunk
    T_ = np.array(T_) / UPH; MR = np.array(MR); CD = np.array(CD); DNA = np.array(DNA)
    dvt = T_[divs(DNA)]
    mr_g = np.interp(GRID, T_, MR); cd_g = np.interp(GRID, T_, CD)
    return mr_g, cd_g, dvt

if __name__ == '__main__':
    res = dict(grid=GRID, rates=np.array(RATES), matched=MATCHED, thresh=THRESH, conds=['ON', 'OFF', 'EZH2i'])
    for c in ('ON', 'OFF', 'EZH2i'):
        for D in RATES:
            MRs, CDs, DVs = [], [], []
            for cell in range(NCELL):
                mr, cd, dv = run_cell(c, D, seed=1000 + cell + int(D))
                MRs.append(mr); CDs.append(cd); DVs.append(dv)
            res[f'mr_{c}_{int(D)}'] = np.array(MRs)          # (NCELL, len(GRID))
            res[f'cd_{c}_{int(D)}'] = np.array(CDs)
            res[f'dv_{c}_{int(D)}'] = np.array(DVs, dtype=object)
            # divisions after onset + time of last division (arrest time), averaged over cells
            post = [len(dv[dv > 0]) for dv in DVs]; last = [dv[dv > 0][-1] if np.any(dv > 0) else 0.0 for dv in DVs]
            res[f'post_{c}_{int(D)}'] = np.array(post); res[f'last_{c}_{int(D)}'] = np.array(last)
            print(f"{c:5s} D={D:5.0f}h: post-withdrawal ÷ (mean)={np.mean(post):.1f}  last ÷ @ (mean)={np.mean(last):.0f}h  "
                  f"start Cd_mRNA={np.mean([m[np.argmin(np.abs(GRID-(-2)))] for m in MRs]):.2f}", flush=True)
    np.savez(os.path.join(HERE, 'sim_matched_withdrawal_cache.npz'), **res)
    print('-> cached; figure: fig_matched_withdrawal.py', flush=True)
