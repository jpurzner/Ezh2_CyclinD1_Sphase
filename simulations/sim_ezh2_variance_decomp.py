"""Many-cell, TIME-RESOLVED EZH2 simulation to decompose the variance into cell-to-cell vs within-cell-
over-time (JP 2026-10-16). The IF snapshot is blind to this split; we must guess it. This measures it for
a given parameterization via the law of total variance:
    V_total(snapshot)  =  V_between (Var of per-cell time-means = CELL-TO-CELL)
                       +  V_within  (mean of per-cell temporal Var = WITHIN-CELL OVER TIME)
A no-bursting control (VD_BURST=0: gene always on at the mean rate, only the frozen factor + cell cycle)
isolates the cell-cycle part of V_within, so V_within(burst) - V_within(noburst) = the BURSTING part.
Config via env: VD_N, VD_SIG_EXT (frozen cell-to-cell knob), VD_DUTY/VD_TAU (bursting temporal knob),
VD_BURST (1/0), VD_OUT (json path), VD_LABEL. No arrest tail (modest CDKI) so all cells cycle -> clean.
Run e.g.:  VD_OUT=/tmp/vd_baked.json VD_SIG_EXT=0.10 VD_DUTY=0.62 venv/bin/python simulations/sim_ezh2_variance_decomp.py
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

N = int(os.environ.get('VD_N', '120'))
DT = 30.0
WARM = 5400.0; REC = 12000.0          # 90 h warmup (protein near stationary), 200 h recorded
TAU = float(os.environ.get('VD_TAU', '200.0'))
DUTY = float(os.environ.get('VD_DUTY', '0.62'))
SIG_EXT = float(os.environ.get('VD_SIG_EXT', '0.10'))
BURST = os.environ.get('VD_BURST', '1') == '1'
LABEL = os.environ.get('VD_LABEL', 'baked')
OUT = os.environ.get('VD_OUT', '/tmp/vd.json')
SIG_CDKI = 0.30                        # modest CDKI spread, NO arrest tail -> all cells cycle
REC_EVERY = 4                          # record every 4 chunks = 2 h
SEL = ['EZH2']
K_ON = DUTY / TAU; K_OFF = (1 - DUTY) / TAU

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']
B0 = {p: rr[p] for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}
rng = np.random.default_rng(3)
mult = np.clip(rng.lognormal(0, SIG_CDKI, N), 0.4, 3.0)      # frozen CDKI (cell-to-cell), no tail
fExt = np.exp(rng.normal(0, SIG_EXT, N))                     # frozen extrinsic (cell-to-cell)


def run_cell(i):
    rr.reset(); rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult[i]; rr['p19'] = 0.58 * mult[i]; rr['kSyP21'] = KS0 * mult[i]
    rr['EZH2'] = 6.0 * fExt[i]; rr['EZH2m'] = 0.3 * fExt[i]
    g = 1 if rng.random() < DUTY else 0
    t = 0.0; T_END = WARM + REC; tr = []; k = 0
    while t < T_END:
        if BURST:
            g = 0 if (g and rng.random() < K_OFF * DT) else (1 if (not g and rng.random() < K_ON * DT) else g)
            eta = fExt[i] * (g / DUTY)
        else:
            eta = fExt[i]                                    # no telegraph: gene on at the mean rate
        for p in B0: rr[p] = B0[p] * eta
        try:
            r = rr.simulate(t, t + DT, 2, selections=SEL)
        except Exception:
            break
        if t >= WARM and k % REC_EVERY == 0:
            tr.append(float(r['EZH2'][-1]))
        t += DT; k += 1
    return np.array(tr)


traces = [run_cell(i) for i in range(N)]
minL = min(len(t) for t in traces if len(t) > 10)
X = np.array([t[:minL] for t in traces if len(t) >= minL])   # cells x time
tvec = (np.arange(minL) * REC_EVERY * DT / 60.0)             # hours

mu_i = X.mean(axis=1); mu = float(mu_i.mean())
V_between = float(mu_i.var(ddof=1))                          # cell-to-cell (variance of per-cell time-means)
V_within = float(X.var(axis=1, ddof=1).mean())              # within-cell over time (avg temporal variance)
V_total = float(X.reshape(-1).var(ddof=1))                  # pooled = snapshot variance
cv_total = float(np.sqrt(V_total) / mu)
tcv_percell = float(np.mean(X.std(axis=1, ddof=1) / mu_i))  # mean within-cell temporal CV

out = {'label': LABEL, 'burst': BURST, 'sig_ext': SIG_EXT, 'duty': DUTY, 'tau': TAU, 'n': int(X.shape[0]),
       'mu': mu, 'V_between': V_between, 'V_within': V_within, 'V_total': V_total,
       'cv_total': cv_total, 'frac_cell2cell': V_between / V_total, 'frac_temporal': V_within / V_total,
       'mean_within_cell_temporal_cv': tcv_percell,
       'tvec': tvec.round(1).tolist(),
       'traces': X[:24].round(3).tolist()}                  # 24 traces for plotting
json.dump(out, open(OUT, 'w'))
print(f"[{LABEL}] burst={BURST} n={X.shape[0]} mu={mu:.2f} cv_total={cv_total:.3f} | "
      f"cell2cell={V_between/V_total:.0%} temporal={V_within/V_total:.0%} | within-cell CV {tcv_percell:.3f}")
print('saved', OUT)
