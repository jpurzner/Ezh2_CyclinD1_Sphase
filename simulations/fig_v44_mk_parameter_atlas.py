"""Figure: parameter ATLAS of H3K27me3 (Mk) plateau at Ccnd1 over divisions (v44.1).

28 plateau plots, each a REAL sweep of one model parameter (proliferating GNP baseline; the mark is
followed across the replicative sawtooth and summarized by its settled mean + min/max band). Markers
are COLORED BY CYCLE PERIOD so that division-timing-mediated accumulation is visible: any parameter
that lengthens the cycle (or arrests it) cuts the dilution frequency and lets Mk accumulate,
independent of the writer/eraser balance.

Grouped: I read-write / PRC2 kinetics · II EZH2 abundance & protein stability · III turnover & eraser
· IV cell-cycle timing (→ dilution frequency) · V mitogen load · VI CDK inhibitors.

Run:  ./venv/bin/python simulations/fig_v44_mk_parameter_atlas.py [--workers N]
"""
import os, sys, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

T_END, N_PTS = 18000, 6000
MK_SETTLE = 9000     # window for the settled-mean/plateau Mk
DIV_SETTLE = 2000    # divisions counted after this (period estimate)
SEL = ['time', 'Mk', 'Dna', 'MPF', 'EZH2', 'Gli1']

# swept parameters that must be RESET to baseline every pooled call (reset() restores species only)
SWEPT = ['a_rw_prc2', 'a0_prc2', 'g_prc2', 'K_prc2', 'n_prc2', 'f0_prc2', 'K_tx_mk', 'p_tx_mk',
         'kTlEZ', 'kDeEZ', 'kEZbas', 'kEZbas_Cd', 'kEZE2f', 'Kez_cd', 'K_E2f_EZ',
         'del_mk', 'k_jmjd3_gli', 'mu', 'M_commit', 'M_size', 'kSyDna', 'MPF_div',
         'SHH', 'k_Cd_translation', 'MYCN_amplification', 'Ptch1_copy_number', 'p18', 'kSyP21']

_RR = None
_BASE = None
def _init(_):
    global _RR, _BASE
    import tellurium as te
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-8)
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000)
        _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass
    _BASE = {n: _RR[n] for n in SWEPT}

def _apply(pname, val):
    _RR.reset()
    for k, v in _BASE.items():
        _RR[k] = v
    _RR[pname] = val

def _simulate(pname, val):
    for maxstep in (20.0, 5.0, 2.0):
        try:
            _apply(pname, val)
            _RR.integrator.setValue('maximum_time_step', maxstep)
            return _RR.simulate(0, T_END, N_PTS, selections=SEL)
        except Exception:
            continue
    return None

def _divtimes(d, settle):
    t = d['time']; dna = d['Dna']
    idx = np.where((dna[:-1] > 0.5) & (dna[1:] < 0.5) & (t[1:] >= settle))[0]
    return t[idx]

def work(task):
    pidx, jx, pname, val = task
    d = _simulate(pname, val)
    if d is None:
        return (pidx, jx, None)
    t = d['time']; s = t >= MK_SETTLE; mk = d['Mk']
    dts = _divtimes(d, DIV_SETTLE)
    ndiv = int(len(dts))
    period = float(np.mean(np.diff(dts)) / 60.0) if ndiv >= 2 else np.nan
    return (pidx, jx, dict(mean=float(mk[s].mean()), lo=float(mk[s].min()), hi=float(mk[s].max()),
                           ndiv=ndiv, period=period))


# ---------------- parameter sweep specification ----------------
# cat, name, human label, kind ('mul' factor of baseline | 'abs'), lo, hi, n, note
CAT = {
    'RW':  ('#6c3483', 'I. read-write / PRC2 kinetics'),
    'EZ':  ('#1f618d', 'II. EZH2 abundance & protein stability'),
    'ER':  ('#148f77', 'III. turnover & eraser'),
    'CY':  ('#b9770e', 'IV. cell-cycle timing (→ dilution frequency)'),
    'MI':  ('#1e8449', 'V. mitogen load'),
    'CK':  ('#a93226', 'VI. CDK inhibitors'),
}
SPEC = [
    ('RW', 'a_rw_prc2', 'read-write recruit  a_rw', 'mul', 0.0, 2.5),
    ('RW', 'a0_prc2', 'accessory recruit  a0', 'mul', 0.2, 4.0),
    ('RW', 'g_prc2', 'nascent-tx eviction  g', 'mul', 0.0, 6.0),
    ('RW', 'K_prc2', 'repression K  (Ccnd1)', 'mul', 0.3, 3.0),
    ('RW', 'n_prc2', 'repression Hill n', 'abs', 1.0, 6.0),
    ('RW', 'f0_prc2', 'leaky floor  f0', 'abs', 0.02, 0.5),
    ('RW', 'K_tx_mk', 'eviction threshold  K_tx', 'mul', 0.4, 3.0),
    ('RW', 'p_tx_mk', 'eviction Hill  p_tx', 'abs', 1.0, 6.0),
    ('EZ', 'kTlEZ', 'EZH2 translation (abundance)', 'mul', 0.2, 2.2),
    ('EZ', 'kDeEZ', 'EZH2 degradation (stability⁻¹)', 'mul', 0.3, 4.0),
    ('EZ', 'kEZbas', 'EZH2 basal transcription', 'mul', 0.0, 8.0),
    ('EZ', 'kEZbas_Cd', 'EZH2 mitogen-scaled basal', 'mul', 0.0, 5.0),
    ('EZ', 'kEZE2f', 'EZH2 E2F-driven tx', 'mul', 0.2, 3.0),
    ('EZ', 'Kez_cd', 'EZH2 mitogen Km', 'mul', 0.3, 3.0),
    ('EZ', 'K_E2f_EZ', 'EZH2 E2F Km', 'mul', 0.3, 3.0),
    ('ER', 'del_mk', 'mark turnover  del_mk', 'mul', 0.2, 12.0),
    ('ER', 'k_jmjd3_gli', 'Gli→Kdm6b eraser', 'mul', 0.0, 8.0),
    ('CY', 'mu', 'growth rate  mu  (period)', 'mul', 0.3, 2.5),
    ('CY', 'M_commit', 'commitment size  M_commit', 'mul', 0.6, 1.8),
    ('CY', 'M_size', 'S-entry size  M_size', 'mul', 0.6, 1.6),
    ('CY', 'kSyDna', 'fork speed  kSyDna (S-dur⁻¹)', 'mul', 0.3, 2.5),
    ('CY', 'MPF_div', 'mitosis threshold  MPF_div', 'mul', 0.5, 1.8),
    ('MI', 'SHH', 'SHH (mitogen dose)', 'abs', 0.15, 1.2),
    ('MI', 'k_Cd_translation', 'CyclinD1 translation', 'mul', 0.3, 2.0),
    ('MI', 'MYCN_amplification', 'MYCN amplification', 'abs', 1.0, 8.0),
    ('MI', 'Ptch1_copy_number', 'Ptch1 copy (1→MB)', 'abs', 1.0, 0.06),
    ('CK', 'p18', 'p18 INK4 (CDK4/6 brake)', 'mul', 0.3, 4.0),
    ('CK', 'kSyP21', 'p21/p27 synthesis (CIP/KIP)', 'mul', 0.5, 5.0),
]
NP = 10  # points per sweep


import pickle
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mk_parameter_atlas_cache.pkl')

def main():
    workers = 8
    fresh = '--fresh' in sys.argv
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    if not fresh and os.path.exists(CACHE):
        print('using cached results (pass --fresh to re-simulate)', flush=True)
        with open(CACHE, 'rb') as f:
            plot(pickle.load(f))
        return
    import tellurium as te
    r0 = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    base = {n: r0[n] for n in SWEPT}
    # build absolute sweep values + plot-x per parameter
    sweeps = []
    for (cat, name, label, kind, lo, hi) in SPEC:
        if kind == 'mul':
            facs = np.linspace(lo, hi, NP)
            vals = base[name] * facs
            xs = facs; xlab = '× baseline'
        else:
            vals = np.linspace(lo, hi, NP)
            xs = vals; xlab = 'value'
        sweeps.append(dict(cat=cat, name=name, label=label, kind=kind, vals=vals, xs=xs,
                           xlab=xlab, base=base[name],
                           base_x=(1.0 if kind == 'mul' else base[name])))
    tasks = []
    for pidx, sw in enumerate(sweeps):
        for jx, v in enumerate(sw['vals']):
            tasks.append((pidx, jx, sw['name'], float(v)))
    print(f'running {len(tasks)} real simulations on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=4)
    for pidx, jx, out in res:
        sweeps[pidx].setdefault('out', [None] * NP)[jx] = out
    with open(CACHE, 'wb') as f:
        pickle.dump(sweeps, f)
    plot(sweeps)


def plot(sweeps):
    PMIN, PMAX = 16.0, 42.0
    norm = Normalize(vmin=PMIN, vmax=PMAX)
    cmap = plt.get_cmap('viridis')
    ncol, nrow = 7, 4
    fig = plt.figure(figsize=(26, 15.5))
    gs = GridSpec(nrow, ncol, figure=fig, hspace=0.5, wspace=0.32,
                  left=0.045, right=0.93, top=0.865, bottom=0.055)
    fig.suptitle("H3K27me3 (Mk) plateau atlas at $Ccnd1$ over divisions — v44.1  (every point = a real simulation; markers colored by cycle period)",
                 fontsize=17, fontweight="bold", y=0.965)
    fig.text(0.5, 0.937, "settled mean Mk (line) + sawtooth min–max (band) vs each parameter, GNP baseline.  "
             "KEY: a longer/arrested cycle (yellow markers, ✕) cuts dilution frequency → Mk accumulates even without changing the writer/eraser balance.",
             ha="center", fontsize=11, color="#555", style="italic")
    # category color key (one colored label per group)
    keys = list(CAT.keys())
    x0, x1 = 0.10, 0.90
    for i, k in enumerate(keys):
        xk = x0 + (x1 - x0) * i / (len(keys) - 1)
        fig.text(xk, 0.910, CAT[k][1], fontsize=9.5, color=CAT[k][0], fontweight='bold', ha='center')

    for pidx, sw in enumerate(sweeps):
        r, c = divmod(pidx, ncol)
        ax = fig.add_subplot(gs[r, c])
        col, catname = CAT[sw['cat']]
        xs = sw['xs']; out = sw.get('out', [None] * NP)
        mean = np.array([o['mean'] if o else np.nan for o in out])
        lo = np.array([o['lo'] if o else np.nan for o in out])
        hi = np.array([o['hi'] if o else np.nan for o in out])
        per = np.array([o['period'] if o else np.nan for o in out])
        ndv = np.array([o['ndiv'] if o else 0 for o in out])
        ax.fill_between(xs, lo, hi, color=col, alpha=0.13, lw=0)
        ax.plot(xs, mean, '-', color=col, lw=1.4, alpha=0.7, zorder=2)
        # cycling points: color by period
        cyc = ndv >= 2
        if cyc.any():
            ax.scatter(xs[cyc], mean[cyc], c=per[cyc], cmap=cmap, norm=norm, s=42,
                       edgecolor='k', linewidth=0.5, zorder=4)
        # arrested / non-cycling: red X (Mk saturates -> no dilution)
        arr = ndv < 2
        if arr.any():
            ax.scatter(xs[arr], mean[arr], marker='x', c='#c0392b', s=46, linewidth=1.8, zorder=5)
        ax.axvline(sw['base_x'], color='k', ls=':', lw=0.9, alpha=0.5)
        ax.set_ylim(-0.03, 1.02)
        ax.set_title(sw['label'], fontsize=9.0, color=col, fontweight='bold')
        ax.set_xlabel(sw['xlab'], fontsize=7.6); ax.tick_params(labelsize=7.2)
        ax.grid(alpha=0.2)
        if c == 0:
            ax.set_ylabel("mean Mk", fontsize=8.5)

    # shared colorbar (period)
    cax = fig.add_axes([0.945, 0.12, 0.013, 0.66])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend='both')
    cb.set_label("cycle period (h)  —  longer = fewer dilutions", fontsize=10)
    fig.text(0.945, 0.085, "✕ = arrested\n(no dilution →\nMk saturates)", fontsize=8, color='#c0392b', ha='left')

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mk_parameter_atlas")
    fig.savefig(out + ".png", dpi=135, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")


if __name__ == '__main__':
    main()
