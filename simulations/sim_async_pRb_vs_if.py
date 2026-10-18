"""Asynchronous-population snapshot of model pRb vs the measured IF pRb distribution (JP 2026-10-06).

Goal (JP): take a marker we have lots of IF data for (pRb, 36 experiments, bimodal G0/committed),
simulate many cells dividing out of sequence with each other, and look at the distribution OVER TIME
of the simulation. The deterministic model already produces a pRb spread purely from cell-cycle-phase
asynchrony — this asks how much of the measured spread that accounts for, and how much is MISSING
(the G0 mode + extra width) that intrinsic bursting / noise would have to supply.

Method: one long deterministic trajectory of the baked MB model = cells sampled at all cycle phases
(an asynchronous population). Classify each time sample by DNA content (2N / S / G2) like the IF gate.
Compare to IF (AMOS imaging_distribution, experiment ABV, MB55/56 DMSO). Units are arbitrary on both
sides, so we align ONLY the median (one multiplicative scale) and do NOT rescale the width — the width
and shape gap is the result.

IF target JSON produced by amos.imaging_distribution (see scratchpad if_pRb_target.json).
Run:  ./venv/bin/python simulations/sim_async_pRb_vs_if.py  /path/to/if_pRb_target.json   (~2 min)
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

IF_JSON = sys.argv[1] if len(sys.argv) > 1 else "/tmp/if_pRb_target.json"
C_IF, C_MOD = '#b0402f', '#1f6f4a'
plt.rcParams.update({'font.size': 10.5, 'axes.titlesize': 11, 'axes.labelsize': 10.5,
                     'axes.titleweight': 'bold', 'font.family': 'DejaVu Sans'})


def shape_stats(x):
    x = np.asarray(x, float); x = x[np.isfinite(x) & (x > 0)]
    lx = np.log10(x)
    from scipy.stats import skew
    return dict(n=len(x), cv=float(x.std(ddof=1) / x.mean()), log10_sd=float(lx.std(ddof=1)),
                skew=float(skew(x)), p05=float(np.quantile(x, .05)), p50=float(np.quantile(x, .50)),
                p95=float(np.quantile(x, .95)))


# ---------- model: one long asynchronous trajectory ----------
rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
except Exception: pass
rr['p18'] = 1.73; rr['p19'] = 0.58
rr.reset()
T = 90000.0; NP = 30000                      # ~60 MB cycles, sampled every ~3 min
r = rr.simulate(0, T, NP, selections=['time', 'pRb', 'Dna'])
warm = r['time'] >= 0.15 * T                  # drop initial transient
t = r['time'][warm]; prb = r['pRb'][warm]; dna = r['Dna'][warm]
keep = np.isfinite(prb) & (prb > 0)
t, prb, dna = t[keep], prb[keep], dna[keep]
# DNA-content state (model Dna: 0=G1/2N, 0->1 across S, 1=G2/4N)
state = np.where(dna <= 0.05, '2N', np.where(dna >= 0.95, 'G2', 'S'))
mod = {s: prb[state == s] for s in ('2N', 'S', 'G2')}
mod_pool = prb
print('MODEL async pRb snapshot (deterministic, phase-mixing only):')
for s in ('2N', 'S', 'G2'):
    st = shape_stats(mod[s]); print(f'  {s}: n={st["n"]:6d} cv={st["cv"]:.3f} log10_sd={st["log10_sd"]:.3f} skew={st["skew"]:+.2f}')
sp = shape_stats(mod_pool); print(f'  pooled: cv={sp["cv"]:.3f} log10_sd={sp["log10_sd"]:.3f} skew={sp["skew"]:+.2f}')

# ---------- IF target ----------
IF = json.load(open(IF_JSON))
if_pool = IF['pooled']
print(f'\nIF pRb (ABV DMSO): pooled cv={if_pool["cv"]:.3f} log10_sd={if_pool["log10_sd"]:.3f} skew={if_pool["skew"]:.2f}')

# ---------- align medians only (one arbitrary-unit scale); DO NOT rescale width ----------
if_log_med = float(if_pool['p50']); if_log_med = np.log10(if_log_med)
mod_log_med = float(np.median(np.log10(mod_pool)))
shift = if_log_med - mod_log_med               # additive in log10 = multiplicative in raw
def mod_log(x): return np.log10(np.asarray(x, float)) + shift

# ---------- figure ----------
fig = plt.figure(figsize=(14, 8))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.28)

# (A) pooled overlay
axA = fig.add_subplot(gs[0, :2])
if_raw = np.asarray(IF['raw_log10_sample'], float)
bins = np.linspace(min(if_raw.min(), mod_log(mod_pool).min()), max(if_raw.max(), mod_log(mod_pool).max()), 60)
axA.hist(if_raw, bins=bins, density=True, color=C_IF, alpha=0.55, label=f'IF pRb (measured, CV {if_pool["cv"]:.2f})')
axA.hist(mod_log(mod_pool), bins=bins, density=True, histtype='step', color=C_MOD, lw=2.2,
         label=f'model async (deterministic, CV {sp["cv"]:.2f})')
axA.set_xlabel('log10 pRb (median-aligned, arb. units)'); axA.set_ylabel('density')
axA.set_title('(A) Pooled: measured IF vs deterministic async model (median aligned, width NOT rescaled)')
axA.legend(fontsize=9)

# (B) stats comparison
axB = fig.add_subplot(gs[0, 2])
labels = ['CV', 'log10_sd', 'skew']
ifv = [if_pool['cv'], if_pool['log10_sd'], if_pool['skew']]
mov = [sp['cv'], sp['log10_sd'], sp['skew']]
x = np.arange(3); w = 0.36
axB.bar(x - w/2, ifv, w, color=C_IF, label='IF')
axB.bar(x + w/2, mov, w, color=C_MOD, label='model')
for xi, (a, b) in enumerate(zip(ifv, mov)):
    axB.text(xi - w/2, a, f'{a:.2f}', ha='center', va='bottom', fontsize=8)
    axB.text(xi + w/2, b, f'{b:.2f}', ha='center', va='bottom', fontsize=8)
axB.set_xticks(x); axB.set_xticklabels(labels); axB.set_title('(B) Shape: IF vs model (pooled)')
axB.legend(fontsize=9)

# (C-E) by-state overlays
if_state = {d['cell_state']: d for d in IF['by_state']}
for j, s in enumerate(('2N', 'S', 'G2')):
    ax = fig.add_subplot(gs[1, j])
    if s in if_state and 'hist_counts' in if_state[s]:
        edges = np.asarray(if_state[s]['hist_edges_log10'], float)
        cnt = np.asarray(if_state[s]['hist_counts'], float); cnt = cnt / cnt.sum() / np.diff(edges).mean()
        ax.bar(edges[:-1], cnt, width=np.diff(edges), align='edge', color=C_IF, alpha=0.55,
               label=f'IF {s} (CV {if_state[s]["cv"]:.2f})')
    ml = mod_log(mod[s])
    ax.hist(ml, bins=30, density=True, histtype='step', color=C_MOD, lw=2.0,
            label=f'model {s} (CV {shape_stats(mod[s])["cv"]:.2f})')
    ax.set_title(f'({chr(67+j)}) {s}'); ax.set_xlabel('log10 pRb (aligned)'); ax.legend(fontsize=8)
    if j == 0: ax.set_ylabel('density')

fig.suptitle('Asynchronous model pRb (phase-mixing only) vs measured IF pRb — how much spread is '
             'cell-cycle asynchrony, and how much is missing (→ bursting/noise + G0)', fontsize=12, fontweight='bold', y=0.98)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_async_pRb_vs_if')
fig.savefig(out + '.png', dpi=150, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('\nwrote', out + '.png/.pdf')
print(f'VERDICT: IF pooled CV {if_pool["cv"]:.2f} vs model async CV {sp["cv"]:.2f} '
      f'-> asynchrony explains ~{100*sp["cv"]/if_pool["cv"]:.0f}% of the CV; '
      f'IF log10_sd {if_pool["log10_sd"]:.2f} vs model {sp["log10_sd"]:.2f}.')
