"""Pick the simulated MB cells that fall into a G0 arrest and show each one DIVIDING and then ARRESTING
(JP 2026-09-22). Current baked model, temporal CDK-inhibitor model: each cell's CDKN expression (p18/p19 +
kSyP21 + birth-p27) is re-drawn at every division from a mean-reverting random walk, so a cell divides for a
while and then, when it draws a high CDKN level, stops dividing (falls into G0). We simulate many MB lineages,
keep the ones with a clear arrest (a long non-dividing gap after >=4 divisions), and plot them.

Run:  ./venv/bin/python simulations/fig_v44_g0_arrest_examples.py    (~5 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, multiprocessing as mp
import tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

C_DIV = '#762A83'; C_CDKN = '#e0902f'; C_ARR = '#8e2f2f'
SEL = ['time', 'Cb', 'Dna', 'p18_prot', 'p19_prot', 'P21']
T_TOTAL, CHUNK, NP = 900 * 60.0, 180.0, 7
RHO, SIGMA, XMAX = 0.7, 0.95, np.log(12.0)      # inheritance / temporal spread / ceiling (allows sustained arrest)
MIN_STOP = 45.0                                 # h; a >45 h no-division stretch = the cell stopped (G0)
DS = 5
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 11.5, 'axes.labelsize': 10.5,
                     'axes.titleweight': 'bold', 'font.family': 'DejaVu Sans'})

_RR = None; _KS0 = None
def _init(_):
    global _RR, _KS0
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 400000)
    except Exception: pass
    _KS0 = _RR['kSyP21']


def _apply(mult):
    _RR['p18'] = 1.73 * mult; _RR['p19'] = 0.58 * mult
    _RR['kSyP21'] = _KS0 * mult; _RR['P21_div'] = 0.6 * mult


def run_lineage(seed):
    rng = np.random.default_rng(seed)
    _RR.reset()
    _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 0.3
    try: _RR['Cd2_expr'] = _RR['CD2_EXPR_MB']
    except Exception: pass
    _RR.integrator.setValue('absolute_tolerance', 1e-8)
    x = float(rng.normal(0, SIGMA)); _apply(np.exp(np.clip(x, -XMAX, XMAX)))
    T, CB, DNA, CDKN = [], [], [], []
    t = 0.0
    while t < T_TOTAL:
        try:
            r = _RR.simulate(t, t + CHUNK, NP, selections=SEL)
        except Exception:
            break
        T.append(r['time'][1:]); CB.append(r['Cb'][1:]); DNA.append(r['Dna'][1:])
        CDKN.append(r['p18_prot'][1:] + r['p19_prot'][1:] + r['P21'][1:])
        dna = r['Dna']
        if np.any((dna[:-1] > 0.9) & (dna[1:] < 0.1)):
            x = RHO * x + np.sqrt(1 - RHO**2) * SIGMA * float(rng.normal())
            _apply(np.exp(np.clip(x, -XMAX, XMAX)))
        t += CHUNK
    if not T:
        return None
    T = np.concatenate(T) / 60.0; CB = np.concatenate(CB); DNA = np.concatenate(DNA); CDKN = np.concatenate(CDKN)
    dvt = T[np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]]
    d = slice(None, None, DS)
    return dict(t=T[d], cb=CB[d], cdkn=CDKN[d], dvt=dvt)


def main():
    NLIN = 60
    with mp.get_context('spawn').Pool(max(1, mp.cpu_count() - 1), initializer=_init, initargs=(None,)) as pool:
        res = [r for r in pool.map(run_lineage, [4000 + i * 97 for i in range(NLIN)], chunksize=1) if r]

    # keep cells that DIVIDE then STOP: >=4 divisions, then the single longest no-division stretch (a G0 arrest).
    cand = []
    for r in res:
        dvt = r['dvt']
        if len(dvt) < 4:
            continue
        gaps = np.diff(dvt); internal = gaps.max() if len(gaps) else 0.0
        term = r['t'][-1] - dvt[-1]
        if term >= internal:
            a0, a1, reenters = float(dvt[-1]), float(r['t'][-1]), False   # stops and stays (within window)
        else:
            gi = int(np.argmax(gaps)); a0, a1, reenters = float(dvt[gi]), float(dvt[gi + 1]), True
        stop = a1 - a0
        if stop >= MIN_STOP and np.sum(dvt <= a0 + 1) >= 3:              # a real dividing phase before the stop
            cand.append((stop, a0, a1, reenters, r))
    cand.sort(key=lambda c: -c[0])
    pick = cand[:6]
    print(f'{len(res)} lineages, {len(cand)} that divide then stop; showing {len(pick)}', flush=True)
    for stop, a0, a1, re, r in pick:
        print(f'  divisions={len(r["dvt"])}, stop={stop:.0f} h, {"re-enters" if re else "stays arrested"}', flush=True)

    ncol = 2; nrow = int(np.ceil(len(pick) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 2.35 * nrow), squeeze=False)
    cdkn_max = max(r['cdkn'].max() for _, _, _, _, r in pick) * 1.05
    for k, (stop, a0, a1, reenters, r) in enumerate(pick):
        ax = axes[k // ncol][k % ncol]; d = r
        ax.plot(d['t'], d['cb'], color=C_DIV, lw=1.0, zorder=3)
        axr = ax.twinx()
        axr.plot(d['t'], d['cdkn'], color=C_CDKN, lw=1.6, zorder=2)
        axr.set_ylim(0, cdkn_max); axr.set_ylabel('CDK-inhibitor', color=C_CDKN, fontsize=8.5)
        axr.tick_params(axis='y', labelcolor=C_CDKN, labelsize=8)
        axr.axvspan(a0, a1, color=C_ARR, alpha=0.10, zorder=0)
        ax.annotate('G0 arrest', xy=((a0 + a1) / 2, 2.55), ha='center', va='top', fontsize=8.5,
                    color=C_ARR, fontweight='bold')
        ax.set_ylim(0, 2.75); ax.set_xlim(0, d['t'][-1])
        ax.set_zorder(axr.get_zorder() + 1); ax.patch.set_visible(False)
        ax.set_yticks([]); ax.tick_params(axis='x', labelsize=8.5)
        tail = 're-enters' if reenters else 'stays arrested'
        ax.text(0.01, 0.9, f'{len(d["dvt"])} divisions, then stops {stop:.0f} h ({tail})', transform=ax.transAxes,
                fontsize=8.5, color=C_DIV, fontweight='bold', va='top')
        if k % ncol == 0:
            ax.set_ylabel('CyclinB\n(divisions)', fontsize=9)
        if k // ncol == nrow - 1:
            ax.set_xlabel('time (h)')
    for j in range(len(pick), nrow * ncol):
        axes[j // ncol][j % ncol].axis('off')

    fig.suptitle('MB cells that divide, then arrest — when a cell’s CDK-inhibitor level (orange) drifts up, '
                 'it stops dividing (G0)', fontsize=12.5, fontweight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_g0_arrest_examples')
    fig.savefig(out + '.png', dpi=190, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
    print('wrote', out + '.png/.pdf')


if __name__ == '__main__':
    main()
