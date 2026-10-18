"""Asynchronous-population snapshot of a model species vs its measured IF distribution, for a
TOTAL-PROTEIN marker (JP 2026-10-08). Generalizes sim_async_pRb_vs_if.py.

JP's point: pRb is a PHOSPHOSPECIFIC antibody — its distribution reflects the commitment switch, not
abundance. EZH2 is a NORMAL (total-protein) antibody, so its distribution IS the protein-abundance
distribution (lognormal), the right target for calibrating abundance noise. Unlike pRb, the model's
EZH2 genuinely varies over the cycle (builds up, halves at division), so phase-mixing contributes.

One long deterministic MB trajectory = cells at all phases (async, dividing out of sequence). Sample
the species over time, classify by DNA content like the IF gate, compare shape to IF (median aligned,
width NOT rescaled, since units are arbitrary). The width/shape gap = what abundance noise must supply.

Run:  ./venv/bin/python simulations/sim_async_marker_vs_if.py <if_json> <model_species> [out_tag]
e.g.  ... sim_async_marker_vs_if.py .../if_EZH2_target.json EZH2 ezh2
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import skew as _skew
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

IF_JSON = sys.argv[1]
SPECIES = sys.argv[2] if len(sys.argv) > 2 else 'EZH2'
TAG = sys.argv[3] if len(sys.argv) > 3 else SPECIES.lower()
C_IF, C_MOD = '#2472a4', '#1f6f4a'
plt.rcParams.update({'font.size': 10.5, 'axes.titlesize': 11, 'axes.labelsize': 10.5,
                     'axes.titleweight': 'bold', 'font.family': 'DejaVu Sans'})


def shape(x):
    x = np.asarray(x, float); x = x[np.isfinite(x) & (x > 0)]
    lx = np.log10(x)
    return dict(n=len(x), cv=float(x.std(ddof=1) / x.mean()), log10_sd=float(lx.std(ddof=1)),
                skew=float(_skew(x)), log10_skew=float(_skew(lx)))


# ---- model: one long asynchronous MB trajectory ----
rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
except Exception: pass
rr['p18'] = 1.73; rr['p19'] = 0.58
rr.reset()
T, NP = 90000.0, 30000
r = rr.simulate(0, T, NP, selections=['time', SPECIES, 'Dna'])
warm = r['time'] >= 0.15 * T
x = r[SPECIES][warm]; dna = r['Dna'][warm]
keep = np.isfinite(x) & (x > 0); x, dna = x[keep], dna[keep]
state = np.where(dna <= 0.05, '2N', np.where(dna >= 0.95, 'G2', 'S'))
mod = {s: x[state == s] for s in ('2N', 'S', 'G2')}
sm = shape(x)
print(f'MODEL async {SPECIES} (deterministic, phase-mixing only): pooled cv={sm["cv"]:.3f} '
      f'log10_sd={sm["log10_sd"]:.3f} log10_skew={sm["log10_skew"]:+.2f}')
for s in ('2N', 'S', 'G2'):
    ss = shape(mod[s]); print(f'  {s}: n={ss["n"]:6d} cv={ss["cv"]:.3f} log10_sd={ss["log10_sd"]:.3f}')

IF = json.load(open(IF_JSON)); ifp = IF['pooled']
print(f'IF {IF["marker"]} ({IF["eid"]} {IF["cond"]}): cv={ifp["cv"]:.3f} log10_sd={ifp["log10_sd"]:.3f} skew={ifp["skew"]:.2f}')

# ---- align medians only; do NOT rescale width ----
shift = np.log10(float(ifp['p50'])) - float(np.median(np.log10(x)))
def mlog(v): return np.log10(np.asarray(v, float)) + shift

frac_cv = 100 * sm['cv'] / ifp['cv']
frac_sd = 100 * sm['log10_sd'] / ifp['log10_sd']

# ---- figure ----
fig = plt.figure(figsize=(14, 7.5))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.28)

axA = fig.add_subplot(gs[0, :2])
ifraw = np.asarray(IF['raw_log10_sample'], float)
bins = np.linspace(min(ifraw.min(), mlog(x).min()), max(ifraw.max(), mlog(x).max()), 70)
axA.hist(ifraw, bins=bins, density=True, color=C_IF, alpha=0.5, label=f'IF {IF["marker"]} (measured, CV {ifp["cv"]:.2f}, log-sd {ifp["log10_sd"]:.2f})')
axA.hist(mlog(x), bins=bins, density=True, histtype='step', color=C_MOD, lw=2.2,
         label=f'model async (deterministic, CV {sm["cv"]:.2f}, log-sd {sm["log10_sd"]:.2f})')
axA.set_xlabel(f'log10 {IF["marker"]} (median-aligned, arb. units)'); axA.set_ylabel('density')
axA.set_title(f'(A) Pooled: measured IF {IF["marker"]} (total-protein Ab) vs deterministic async model')
axA.legend(fontsize=8.5)

axB = fig.add_subplot(gs[0, 2])
labels = ['CV', 'log10_sd']
ifv = [ifp['cv'], ifp['log10_sd']]; mov = [sm['cv'], sm['log10_sd']]
xb = np.arange(2); w = 0.36
axB.bar(xb - w/2, ifv, w, color=C_IF, label='IF'); axB.bar(xb + w/2, mov, w, color=C_MOD, label='model')
for xi, (a, b) in enumerate(zip(ifv, mov)):
    axB.text(xi - w/2, a, f'{a:.2f}', ha='center', va='bottom', fontsize=8)
    axB.text(xi + w/2, b, f'{b:.2f}', ha='center', va='bottom', fontsize=8)
axB.set_xticks(xb); axB.set_xticklabels(labels); axB.legend(fontsize=9)
axB.set_title(f'(B) Spread: async explains\n~{frac_cv:.0f}% of CV, ~{frac_sd:.0f}% of log-width')

ifst = {d['cell_state']: d for d in IF['by_state']}
for j, s in enumerate(('2N', 'S', 'G2')):
    ax = fig.add_subplot(gs[1, j])
    if s in ifst and 'hist_counts' in ifst[s]:
        e = np.asarray(ifst[s]['hist_edges_log10'], float)
        c = np.asarray(ifst[s]['hist_counts'], float); c = c / c.sum() / np.diff(e).mean()
        ax.bar(e[:-1], c, width=np.diff(e), align='edge', color=C_IF, alpha=0.5, label=f'IF {s} (CV {ifst[s]["cv"]:.2f})')
    ax.hist(mlog(mod[s]), bins=30, density=True, histtype='step', color=C_MOD, lw=2.0, label=f'model {s} (CV {shape(mod[s])["cv"]:.2f})')
    ax.set_title(f'({chr(67+j)}) {s}'); ax.set_xlabel(f'log10 {IF["marker"]} (aligned)'); ax.legend(fontsize=8)
    if j == 0: ax.set_ylabel('density')

fig.suptitle(f'Asynchronous model {IF["marker"]} (phase-mixing only) vs measured IF {IF["marker"]} '
             f'(total-protein antibody, lognormal) — the width gap is the abundance-noise budget',
             fontsize=11.5, fontweight='bold', y=0.98)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'sim_async_{TAG}_vs_if')
fig.savefig(out + '.png', dpi=150, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
print(f'VERDICT: async phase-mixing explains ~{frac_cv:.0f}% of the measured {IF["marker"]} CV '
      f'({sm["cv"]:.2f} vs {ifp["cv"]:.2f}); log-width {sm["log10_sd"]:.2f} vs {ifp["log10_sd"]:.2f}. '
      f'IF log-skew ~{ifp["skew"]:.2f} (lognormal).')
