"""G0 arrest is MB-specific — the high-CDK-inhibitor tail of the population transiently arrests
(JP 2026-09-22, revised). Current baked model (integrated governor + CDK4/6-escape + proximal/distal).

Three panels (per JP):
  (A) G0 fraction: GNP vs MB (each a population with inherited CDK-inhibitor variance + an MB Atoh1+ tail).
  (B) DISTRIBUTION of the CDK-inhibitor proteins (p18/p19/p27) across the population — the cells that fall
      into transient G0 are the high-CDKI tail (coloured by fate); GNP sits well below the arrest zone.
  (C) The simulated cells that FALL INTO G0 — CyclinB traces of tail cells, dwelling in G0 between divisions.

Efficient: fate depends (to good approximation) only on a cell's CDK-inhibitor level, so we classify vs
CDKI fold ONCE per cell type, then map sampled populations through it.

Run:  ./venv/bin/python simulations/fig_v44_g0_mb_specific.py    (~10 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

C_CYC = '#5aa9e6'; C_G0 = '#8e2f2f'; C_GNP = '#1b9e77'; C_MB = '#762A83'
PER_TG0 = 40.0     # cycle period (h) above this = transient-G0 dwell
GNP_TOT, MB_TOT = 0.464 + 0.36, 1.73 + 0.58    # baseline total INK4 (p18+p19), a.u.
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
                     'axes.titleweight': 'bold', 'font.family': 'DejaVu Sans'})

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
KS0 = rr['kSyP21']


def set_cell(kind, fold):
    rr.reset()
    if kind == 'GNP':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0
        rr['p18'] = 0.464 * fold; rr['p19'] = 0.36 * fold; rr['kSyP21'] = KS0 * fold
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = 1.73 * fold; rr['p19'] = 0.58 * fold; rr['kSyP21'] = KS0 * fold


def classify(kind, fold, T=24000):
    set_cell(kind, fold)
    r = rr.simulate(0, T, T // 5, selections=['time', 'Dna'])
    d = r['time'][np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]]
    if len(d) < 2:
        return 'G0'                                   # permanent arrest
    return 'G0' if float(np.mean(np.diff(d)) / 60) > PER_TG0 else 'cycling'   # transient-G0 vs cycling


FOLDS = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0])
grids = {k: [classify(k, f) for f in FOLDS] for k in ('GNP', 'MB')}
print('MB  grid:', list(zip(FOLDS, grids['MB'])), flush=True)
print('GNP grid:', list(zip(FOLDS, grids['GNP'])), flush=True)
fate = lambda kind, mult: grids[kind][int(np.argmin(np.abs(FOLDS - mult)))]

# ---- populations: inherited lognormal CDKI variance (+ MB Atoh1+ high-CDKI tail) ----
rng = np.random.default_rng(3)
N = 5000
CV = 0.5; sig = np.sqrt(np.log(1 + CV**2))
bulk = lambda n: np.clip(rng.lognormal(-0.5 * sig**2, sig, n), 0.2, 8.0)
TAIL_FRAC, TAIL_X = 0.18, 5.0
gnp_m = bulk(N)
tail_mask = rng.random(N) < TAIL_FRAC
mb_m = np.where(tail_mask, np.clip(rng.lognormal(np.log(TAIL_X), 0.2, N), 0.2, 8.0), bulk(N))

gnp_fate = np.array([fate('GNP', m) for m in gnp_m])
mb_fate = np.array([fate('MB', m) for m in mb_m])
gnp_g0 = 100 * np.mean(gnp_fate == 'G0'); mb_g0 = 100 * np.mean(mb_fate == 'G0')
print(f'G0 fraction: GNP {gnp_g0:.1f}%  MB {mb_g0:.1f}%', flush=True)

# absolute CDK-inhibitor protein (p18+p19, a.u.) per cell
gnp_cdkn = gnp_m * GNP_TOT
mb_cdkn = mb_m * MB_TOT
G0_THR_ABS = 3.0 * MB_TOT     # MB transient-G0 onset (~3x MB mean) in absolute units

# ---- example cells that fall into G0 (MB tail cells) ----
def trace(kind, fold, T=21000):
    set_cell(kind, fold)
    r = rr.simulate(0, T, 14000, selections=['time', 'Cb', 'Dna'])
    dd = r['time'][np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]] / 60
    return r['time'] / 60, r['Cb'], dd

tails = [(3.0, trace('MB', 3.0)), (5.0, trace('MB', 5.0)), (7.0, trace('MB', 7.0))]

# ============================ figure ============================
fig = plt.figure(figsize=(15.5, 5.4))
gs = GridSpec(1, 3, figure=fig, wspace=0.34, left=0.05, right=0.985, top=0.82, bottom=0.15, width_ratios=[0.75, 1.2, 1.25])

# (A) G0 fraction: GNP vs MB (no center bar)
axA = fig.add_subplot(gs[0, 0])
axA.bar([0, 1], [gnp_g0, mb_g0], color=[C_GNP, C_MB], edgecolor='k', width=0.62)
for x, v in zip([0, 1], [gnp_g0, mb_g0]):
    axA.text(x, v + 0.5, f'{v:.0f}%', ha='center', fontweight='bold', fontsize=13)
axA.set_xticks([0, 1]); axA.set_xticklabels(['GNP', 'MB']); axA.set_xlim(-0.6, 1.6)
axA.set_ylabel('% of cells in transient G0'); axA.set_ylim(0, max(gnp_g0, mb_g0) * 1.35 + 3)
axA.set_title('(A) G0 arrest is\nMB-specific', fontsize=11.5)

# (B) distribution of the CDK-inhibitor proteins, coloured by fate
axB = fig.add_subplot(gs[0, 1])
bins = np.linspace(0, mb_cdkn.max() * 1.02, 45)
axB.hist(mb_cdkn[mb_fate == 'cycling'], bins=bins, color=C_CYC, alpha=0.85, label='MB cycling', edgecolor='white', lw=0.3)
axB.hist(mb_cdkn[mb_fate == 'G0'], bins=bins, color=C_G0, alpha=0.9, label='MB transient G0', edgecolor='white', lw=0.3)
axB.hist(gnp_cdkn, bins=bins, histtype='step', color=C_GNP, lw=1.8, label='GNP')
axB.axvspan(G0_THR_ABS, bins[-1], color=C_G0, alpha=0.06)
axB.axvline(G0_THR_ABS, color=C_G0, ls='--', lw=1.2)
axB.text(G0_THR_ABS + 0.15, axB.get_ylim()[1] * 0.9, 'transient-G0\nthreshold', fontsize=8.5, color=C_G0, va='top')
axB.set_xlabel('CDK-inhibitor protein  (p18+p19+p27, a.u.)'); axB.set_ylabel('number of cells')
axB.set_title('(B) The cells that fall into G0 are the\nhigh-CDK-inhibitor tail', fontsize=11.5)
axB.legend(fontsize=9.5, loc='upper right')

# (C) simulated cells that fall into G0 (MB tail cells) — CyclinB, G0 dwells shaded
axC = fig.add_subplot(gs[0, 2])
for k, (fold, (t, cb, dd)) in enumerate(tails):
    off = k * 3.1
    axC.plot(t, cb + off, color=C_MB, lw=1.1)
    edges = np.concatenate([[0.0], dd, [t[-1]]])
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a > PER_TG0 * 0.5:                     # shade the long inter-division G0 dwells
            axC.axvspan(a + 3, b - 3, ymin=(off) / 12.5, ymax=(off + 2.7) / 12.5, color=C_G0, alpha=0.06)
    per = np.mean(np.diff(dd)) if len(dd) >= 2 else np.nan
    axC.text(3, off + 2.75, f'CDK-inhibitor {fold:.0f}× mean  —  G0 dwell, ~{per:.0f} h cycle',
             fontsize=9, color=C_MB, fontweight='bold')
axC.set_xlim(0, 350); axC.set_ylim(0, 12.5); axC.set_yticks([])
axC.set_xlabel('time (h)'); axC.set_ylabel('CyclinB–CDK1  (each spike = a division)')
axC.set_title('(C) The simulated MB cells that fall into G0\ndwell between divisions', fontsize=11.5)

fig.suptitle('Transient-G0 arrest is MB-specific — the high-CDK-inhibitor tail of the population dwells in '
             'reversible G0 between divisions', fontsize=13, fontweight='bold', y=0.965)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_g0_mb_specific')
fig.savefig(out + '.png', dpi=190, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
