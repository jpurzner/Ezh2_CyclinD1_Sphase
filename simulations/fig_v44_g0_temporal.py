"""Transient G0 as a TEMPORAL, reversible state: each cell's CDK-inhibitor level drifts up and down across
divisions (a mean-reverting random walk, inherited but fluctuating), so a lineage divides for a while, drifts
into G0, then drifts back down and re-enters (JP 2026-09-22). Current baked model.

Mechanism: the CDKN expression targets (p18/p19 -> p18_prot/p19_prot which relax to them; kSyP21 & birth-p27
P21_div) are RE-DRAWN at every division from an AR(1) process in log-space:
    x_{n+1} = rho*x_n + sqrt(1-rho^2)*sigma*z ,   CDKI level = mean * exp(x)   (clipped to a reversible range)
rho = division-to-division inheritance (sisters/mother-daughter correlation); sigma = temporal spread. The
population cell-to-cell spread is the stationary distribution of this same process (one mechanism, both scales).
When a cell's level drifts above ~3x the MB mean it enters a reversible transient-G0 dwell (CDK4/6-escape);
GNP's mean is far below its (permanent) arrest threshold, so GNP just cycles.

  (A) G0 fraction (time-averaged), GNP vs MB.
  (B) DISTRIBUTION of the CDK-inhibitor protein (p18_prot+p19_prot+p27), pooled over time x cells, coloured by
      instantaneous fate -- the G0 time-points are the high-CDKI tail.
  (C) Example MB lineages: CyclinB divisions + the fluctuating CDKN protein -- divide, drift into G0, re-enter.

Run:  ./venv/bin/python simulations/fig_v44_g0_temporal.py    (~6 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, multiprocessing as mp
import tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

C_CYC = '#5aa9e6'; C_G0 = '#8e2f2f'; C_GNP = '#1b9e77'; C_MB = '#762A83'; C_CDKN = '#e0902f'
SEL = ['time', 'Cb', 'Dna', 'pRb', 'aRc', 'p18_prot', 'p19_prot', 'P21']
T_TOTAL, CHUNK, NP = 1500 * 60.0, 180.0, 7          # 1500 h, 3 h chunks
RHO, SIGMA, XMAX = 0.5, 0.8, np.log(6.5)            # AR(1) inheritance / temporal spread / reversible ceiling
G0_AGE = 28 * 60.0                                  # min; pre-S longer than a normal cycle = transient-G0 dwell
DS = 6                                               # downsample (record every 6th point ~ 30 min)
BASE = {'GNP': dict(Ptch1=1.0, p18=0.464, p19=0.36), 'MB': dict(Ptch1=0.3, p18=1.73, p19=0.58)}
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
                     'axes.titleweight': 'bold', 'font.family': 'DejaVu Sans'})

_RR = None; _KS0 = None
def _init(_):
    global _RR, _KS0
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 400000)
    except Exception: pass
    _KS0 = _RR['kSyP21']


def _apply(kind, mult):
    b = BASE[kind]
    _RR['p18'] = b['p18'] * mult; _RR['p19'] = b['p19'] * mult
    _RR['kSyP21'] = _KS0 * mult; _RR['P21_div'] = 0.6 * mult


def run_lineage(args):
    kind, seed = args
    rng = np.random.default_rng(seed)
    _RR.reset()
    _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = BASE[kind]['Ptch1']
    if kind == 'MB':
        try: _RR['Cd2_expr'] = _RR['CD2_EXPR_MB']
        except Exception: pass
    x = float(rng.normal(0, SIGMA))
    _apply(kind, np.exp(np.clip(x, -XMAX, XMAX)))
    T, CB, PRB, ARC, DNA, CDKN = [], [], [], [], [], []
    t = 0.0
    for atol in (1e-8,):
        _RR.integrator.setValue('absolute_tolerance', atol)
    while t < T_TOTAL:
        try:
            r = _RR.simulate(t, t + CHUNK, NP, selections=SEL)
        except Exception:
            break
        T.append(r['time'][1:]); CB.append(r['Cb'][1:]); PRB.append(r['pRb'][1:])
        ARC.append(r['aRc'][1:]); DNA.append(r['Dna'][1:])
        CDKN.append(r['p18_prot'][1:] + r['p19_prot'][1:] + r['P21'][1:])
        dna = r['Dna']
        if np.any((dna[:-1] > 0.9) & (dna[1:] < 0.1)):          # a division happened -> re-draw the CDKN level
            x = RHO * x + np.sqrt(1 - RHO**2) * SIGMA * float(rng.normal())
            _apply(kind, np.exp(np.clip(x, -XMAX, XMAX)))
        t += CHUNK
    if not T:
        return kind, np.nan, None
    T = np.concatenate(T); CB = np.concatenate(CB); PRB = np.concatenate(PRB)
    ARC = np.concatenate(ARC); DNA = np.concatenate(DNA); CDKN = np.concatenate(CDKN)
    warm = T >= 200 * 60.0                                       # drop the initial transient
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    dvt = T[np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]]     # division times
    if len(dvt) > 0:
        j = np.searchsorted(dvt, T, side='right') - 1
        age = np.where(j >= 0, T - dvt[np.clip(j, 0, len(dvt) - 1)], T)   # time since last division
    else:
        age = T.copy()
    last_div = dvt[-1] if len(dvt) else -1.0
    # TRANSIENT G0 only: a prolonged pre-commitment dwell that is later RESOLVED by a division (reversible).
    # Permanent exits (high-CDKI cells that never re-divide, e.g. GNP differentiation) are excluded.
    g0 = preS & (age > G0_AGE) & (T < last_div)
    g0_frac = 100.0 * float(np.mean(g0[warm]))
    d = slice(None, None, DS)
    return kind, g0_frac, dict(t=T[d] / 60.0, cb=CB[d], cdkn=CDKN[d], g0=g0[d])


def main():
    N = 24
    tasks = [(k, 7000 + i * 131 + (0 if k == 'GNP' else 999)) for k in ('GNP', 'MB') for i in range(N)]
    print(f'temporal G0: N={N}/type, T={T_TOTAL/60:.0f} h, chunk={CHUNK/60:.0f} h ...', flush=True)
    with mp.get_context('spawn').Pool(max(1, mp.cpu_count() - 1), initializer=_init, initargs=(None,)) as pool:
        res = pool.map(run_lineage, tasks, chunksize=1)

    gnp = [r for r in res if r[0] == 'GNP' and r[2] is not None]
    mb = [r for r in res if r[0] == 'MB' and r[2] is not None]
    gnp_g0 = float(np.mean([r[1] for r in gnp])); mb_g0 = float(np.mean([r[1] for r in mb]))
    print(f'G0 time-fraction: GNP {gnp_g0:.1f}%  MB {mb_g0:.1f}%', flush=True)
    print('MB per-lineage G0%:', sorted(round(r[1], 1) for r in mb), flush=True)

    # pooled CDKN protein distribution (MB), split by instantaneous fate
    mb_cdkn = np.concatenate([r[2]['cdkn'] for r in mb]); mb_g0mask = np.concatenate([r[2]['g0'] for r in mb])
    gnp_cdkn = np.concatenate([r[2]['cdkn'] for r in gnp])
    thr = np.percentile(mb_cdkn[mb_g0mask], 5) if mb_g0mask.any() else np.nan

    # example MB lineages with clear intermittent G0 (g0% in a mid range)
    ex = sorted([r for r in mb if 8 < r[1] < 55], key=lambda r: abs(r[1] - 25))[:2] or mb[:2]

    # ============================ figure ============================
    fig = plt.figure(figsize=(15.5, 5.6))
    gs = GridSpec(2, 3, figure=fig, wspace=0.34, hspace=0.55, left=0.05, right=0.985, top=0.80, bottom=0.12,
                  width_ratios=[0.72, 1.15, 1.35], height_ratios=[1, 1])

    axA = fig.add_subplot(gs[:, 0])
    axA.bar([0, 1], [gnp_g0, mb_g0], color=[C_GNP, C_MB], edgecolor='k', width=0.62)
    for xx, v in zip([0, 1], [gnp_g0, mb_g0]):
        axA.text(xx, v + 0.4, f'{v:.0f}%', ha='center', fontweight='bold', fontsize=13)
    axA.set_xticks([0, 1]); axA.set_xticklabels(['GNP', 'MB']); axA.set_xlim(-0.6, 1.6)
    axA.set_ylabel('% of time in transient G0'); axA.set_ylim(0, max(gnp_g0, mb_g0) * 1.35 + 3)
    axA.set_title('(A) Time in G0 is\nMB-specific', fontsize=11.5)

    axB = fig.add_subplot(gs[:, 1])
    bins = np.linspace(0, np.percentile(mb_cdkn, 99.5), 45)
    axB.hist(mb_cdkn[~mb_g0mask], bins=bins, color=C_CYC, alpha=0.85, label='MB cycling', edgecolor='white', lw=0.3)
    axB.hist(mb_cdkn[mb_g0mask], bins=bins, color=C_G0, alpha=0.9, label='MB in G0', edgecolor='white', lw=0.3)
    axB.hist(gnp_cdkn, bins=bins, histtype='step', color=C_GNP, lw=1.8, label='GNP')
    if not np.isnan(thr):
        axB.axvline(thr, color=C_G0, ls='--', lw=1.2)
    axB.set_xlabel('CDK-inhibitor protein  (p18+p19+p27, a.u.)'); axB.set_ylabel('cell-time (samples)')
    axB.set_title('(B) The high-CDK-inhibitor time-points\nare the cells in G0', fontsize=11.5)
    axB.legend(fontsize=9.5, loc='upper right')

    for k, r in enumerate(ex):
        axc = fig.add_subplot(gs[k, 2])
        d = r[2]
        axc.plot(d['t'], d['cb'], color=C_MB, lw=0.9, zorder=3)
        axr = axc.twinx()
        axr.plot(d['t'], d['cdkn'], color=C_CDKN, lw=1.5, zorder=2)
        axr.fill_between(d['t'], 0, d['cdkn'], where=d['g0'], color=C_G0, alpha=0.12, zorder=0, step='mid')
        axr.set_ylim(0, np.percentile(mb_cdkn, 99.5)); axr.tick_params(axis='y', labelcolor=C_CDKN, labelsize=8)
        axc.set_ylim(0, 2.7); axc.set_xlim(0, d['t'][-1])
        axc.set_zorder(axr.get_zorder() + 1); axc.patch.set_visible(False)
        axc.set_yticks([]); axc.tick_params(axis='x', labelsize=8)
        axc.text(0.01, 0.86, f'MB lineage ({r[1]:.0f}% of time in G0)', transform=axc.transAxes,
                 fontsize=8.5, color=C_MB, fontweight='bold')
        if k == 0:
            axc.set_title('(C) A cell divides, drifts into G0, then re-enters\n(green=CDK-inhibitor, spikes=divisions, shaded=G0)',
                          fontsize=10.5)
        if k == 1:
            axc.set_xlabel('time (h)')
        axr.set_ylabel('CDKN', color=C_CDKN, fontsize=9)

    fig.suptitle('Transient G0 is a reversible, time-varying state — the CDK-inhibitor level fluctuates across '
                 'divisions, so MB cells divide, pause in G0, and re-enter', fontsize=12.5, fontweight='bold', y=0.965)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_g0_temporal')
    fig.savefig(out + '.png', dpi=190, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
    print('wrote', out + '.png/.pdf')


if __name__ == '__main__':
    main()
