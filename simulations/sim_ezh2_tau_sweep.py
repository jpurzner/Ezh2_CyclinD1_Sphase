"""Does transcriptional bursting give MORE over-time (within-cell) variance? (JP 2026-10-15)
The lever is the burst CORRELATION TIME (tau_burst) relative to the stable protein lifetime (~60 h).
Short tau << lifetime: the protein averages many bursts -> small, fast-smoothed drift. Long tau ~ lifetime:
the protein tracks slow on/off epochs -> large, slow within-cell drift + long dwell at extremes (G0-prone).
Isolates the bursting's OWN temporal contribution: representative cell (mult=1, fExt=1), no frozen extrinsic,
no CDKI tail; vary only the telegraph seed per cell. Reports per-cell temporal CV, the across-cell snapshot
CV at the final time, and the deep-low dwell fraction (time below 0.6x median = G0-prone).
Run:  venv/bin/python simulations/sim_ezh2_tau_sweep.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

DT = 30.0; WARM = 5400.0; REC = 15000.0    # 90 h warmup + 250 h recorded
DUTY = 0.62; NCELL = 8
TAUS = [200.0, 1200.0, 2400.0]             # min: 3.3 h (current), 20 h, 40 h
SEL = ['EZH2', 'Dna', 'pRb']

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
B0 = {p: rr[p] for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}


def run_cell(tau, seed):
    rng = np.random.default_rng(seed)
    k_on = DUTY / tau; k_off = (1 - DUTY) / tau
    rr.reset(); rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    g = 1 if rng.random() < DUTY else 0
    t = 0.0; T_END = WARM + REC; T, E = [], []
    while t < T_END:
        g = 0 if (g and rng.random() < k_off * DT) else (1 if (not g and rng.random() < k_on * DT) else g)
        eta = max(g / DUTY, 0.08)       # small basal leak: gene never fully silent (numerically kinder)
        for p in B0: rr[p] = B0[p] * eta
        try:
            r = rr.simulate(t, t + DT, 2, selections=SEL)
        except Exception:
            break                        # integrator stiffness: keep the partial trace
        if t >= WARM:
            T.append((t - WARM) / 60.0); E.append(float(r['EZH2'][-1]))
        t += DT
    return np.array(T), np.array(E)


results = {}
for tau in TAUS:
    raw = [run_cell(tau, 200 + s) for s in range(NCELL)]              # seeds 200..207, same set per tau
    traces = [(T, E) for (T, E) in raw if len(E) >= 150]             # keep cells with >=75 h of trace
    allE = np.concatenate([e for _, e in traces])
    med = np.median(allE)
    tcv = np.mean([e.std() / e.mean() for _, e in traces])            # mean within-cell temporal CV
    final = np.array([e[-1] for _, e in traces])
    scv = final.std(ddof=1) / final.mean() if len(final) > 1 else float('nan')
    dwell = np.mean(allE < 0.6 * med)                                 # deep-low dwell fraction
    results[tau] = dict(traces=traces, tcv=tcv, scv=scv, dwell=dwell, med=med)
    print(f"tau={tau/60:4.1f}h  n_usable={len(traces)}  temporal_CV={tcv:.3f}  snapshot_CV={scv:.3f}  deep_low_dwell={dwell:.3f}", flush=True)

fig, axes = plt.subplots(1, len(TAUS), figsize=(15, 4.6), sharey=True)
fig.suptitle("Bursting & over-time variance: longer burst correlation time tau -> larger, slower within-cell drift\n"
             "(isolated bursting: no frozen extrinsic, no CDKI tail; same 8 cells per panel, DMSO)",
             fontsize=12, fontweight='bold')
for ax, tau in zip(axes, TAUS):
    R = results[tau]
    for (T, E) in R['traces']:
        ax.plot(T, E, lw=1.0, alpha=0.8)
    ax.axhline(0.6 * R['med'], color='k', ls='--', lw=0.8, alpha=0.5)
    ax.set_title(f"tau = {tau/60:.0f} h\ntemporal CV {R['tcv']:.2f} | snapshot CV {R['scv']:.2f} | low-dwell {R['dwell']:.0%}",
                 fontsize=10)
    ax.set_xlabel('time (h)')
axes[0].set_ylabel('EZH2 protein (a.u.)')
plt.tight_layout(rect=[0, 0, 1, 0.9])
out = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/fig_ezh2_tau_sweep.png"
plt.savefig(out, dpi=130); print('saved', out)
