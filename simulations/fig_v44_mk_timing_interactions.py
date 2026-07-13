"""Figure: does division timing gate the H3K27me3 plateau? + top-lever 2-D interaction maps (v44.1).

Two questions, all REAL v44.1 sims (proliferating GNP baseline, mark followed across the sawtooth):

(1) TIMING TEST.  A runtime k_meth multiplier is injected on the methylation reaction so DEPOSITION
    speed is tunable independently of repression.  We compare two ways to change the mark's kinetics:
      - slow DEPOSITION (k_meth down): lowers the plateau but leaves timing-sensitivity ~0.
      - make the mark PERSISTENT (loss = del_mk & k_jmjd3_gli down): lengthens the relaxation
        timescale -> division timing (cycle period, swept via mu) becomes a real lever.
    Conclusion: the plateau is buffered against cycle period UNLESS the mark is persistent (slow loss)
    -- which is the strong-chromatin-memory regime the MB/GNP data exclude.

(2) INTERACTION MAPS.  Pairwise phase diagrams (mean Mk) for the top continuous levers:
    EZH2 stability (kDeEZ), eraser (k_jmjd3_gli), read-write (a_rw), turnover (del_mk), abundance (kTlEZ).

Run:  ./venv/bin/python simulations/fig_v44_mk_timing_interactions.py [--workers N] [--fresh]
"""
import os, sys, pickle, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

T_END, N_PTS, MK_SETTLE, DIV_SETTLE = 20000, 6000, 10000, 2000
SEL = ['time', 'Mk', 'Dna']
BASE = dict(del_mk=0.0015, k_jmjd3_gli=0.047832059241740936, a_rw_prc2=0.00429701,
            kDeEZ=0.00116852, kTlEZ=0.0150989, mu=0.0005, k_meth=1.0, SHH=0.5)
SWEPT = list(BASE.keys())
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mk_timing_interactions_cache.pkl')

_RR = None
def _init(_):
    global _RR
    import tellurium as te
    m = build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_chain=False, with_two_step_rb=False)
    m = m.replace("Mk_methylation: => Mk; Cell*PRC2*(1 - Mk);",
                  "k_meth = 1.0;\n  Mk_methylation: => Mk; Cell*k_meth*PRC2*(1 - Mk);")
    _RR = te.loada(m)
    _RR.integrator.setValue('absolute_tolerance', 1e-8)
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000)
        _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass

def work(task):
    group, key, ov = task
    for maxstep in (20.0, 5.0, 2.0):
        try:
            _RR.reset()
            for k, v in BASE.items():
                _RR[k] = v
            for k, v in ov.items():
                _RR[k] = v
            _RR.integrator.setValue('maximum_time_step', maxstep)
            d = _RR.simulate(0, T_END, N_PTS, selections=SEL)
            break
        except Exception:
            d = None
    if d is None:
        return (group, key, None)
    t = d['time']; s = t >= MK_SETTLE; mk = d['Mk']
    dna = d['Dna']; idx = np.where((dna[:-1] > 0.5) & (dna[1:] < 0.5) & (t[1:] >= DIV_SETTLE))[0]
    per = float(np.mean(np.diff(t[idx])) / 60.0) if len(idx) >= 2 else np.nan
    return (group, key, dict(mean=float(mk[s].mean()), ndiv=int(len(idx)), period=per))


# sweep specs
PERSIST_LEVELS = [2.0, 1.0, 0.5, 0.25, 0.1, 0.04]
KMETH_LEVELS = [2.0, 1.3, 1.0, 0.6, 0.3, 0.15, 0.08]
MUFAC = np.linspace(0.35, 2.3, 8)
MAPN = 11
MAPS = [
    ('map_stab_eraser', 'kDeEZ', (0.3, 3.5), 'k_jmjd3_gli', (0.0, 6.0),
     "EZH2 stability (1/kDeEZ) × eraser", "EZH2 degradation kDeEZ (× base; →right = less stable)", "eraser k_jmjd3_gli (× base)", False, False),
    ('map_rw_turn', 'a_rw_prc2', (0.0, 2.5), 'del_mk', (0.2, 10.0),
     "read-write × turnover", "read-write a_rw (× base)", "turnover del_mk (× base)", False, False),
    ('map_abund_eraser', 'kTlEZ', (0.3, 2.2), 'k_jmjd3_gli', (0.0, 6.0),
     "EZH2 abundance × eraser", "EZH2 abundance kTlEZ (× base)", "eraser k_jmjd3_gli (× base)", False, False),
    ('map_stab_rw', 'kDeEZ', (0.3, 3.5), 'a_rw_prc2', (0.0, 2.5),
     "EZH2 stability × read-write", "EZH2 degradation kDeEZ (× base; →right = less stable)", "read-write a_rw (× base)", False, False),
]


def build_tasks():
    tasks = []
    # (1) persistence x period
    for pl in PERSIST_LEVELS:
        for j, mf in enumerate(MUFAC):
            tasks.append(('persist', (pl, j), dict(del_mk=BASE['del_mk'] * pl,
                          k_jmjd3_gli=BASE['k_jmjd3_gli'] * pl, mu=BASE['mu'] * mf)))
    # (1b) deposition (k_meth) x period-extremes  (contrast curve)
    for kl in KMETH_LEVELS:
        for j, mf in enumerate([MUFAC[0], MUFAC[-1]]):
            tasks.append(('kmeth', (kl, j), dict(k_meth=kl, mu=BASE['mu'] * mf)))
    # (2) interaction maps
    for (name, px, prx, py, pry, *_rest) in MAPS:
        XF = np.linspace(prx[0], prx[1], MAPN); YF = np.linspace(pry[0], pry[1], MAPN)
        for ix, xf in enumerate(XF):
            for iy, yf in enumerate(YF):
                tasks.append((name, (ix, iy), {px: BASE[px] * xf, py: BASE[py] * yf}))
    return tasks


def main():
    workers = 8
    fresh = '--fresh' in sys.argv
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    if not fresh and os.path.exists(CACHE):
        print('using cached results (--fresh to re-simulate)', flush=True)
        with open(CACHE, 'rb') as f:
            plot(pickle.load(f))
        return
    tasks = build_tasks()
    print(f'running {len(tasks)} real simulations on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=4)
    R = {}
    for g, k, o in res:
        R.setdefault(g, {})[k] = o
    with open(CACHE, 'wb') as f:
        pickle.dump(R, f)
    plot(R)


def plot(R):
    fig = plt.figure(figsize=(19, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.34, wspace=0.3, left=0.06, right=0.955, top=0.9, bottom=0.07)
    fig.suptitle("Does division timing set the H3K27me3 plateau?  Timing test + top-lever interaction maps — v44.1 (real sims)",
                 fontsize=15, fontweight="bold", y=0.965)
    fig.text(0.5, 0.925, "TIMING is buffered unless the mark is made PERSISTENT (slow loss) — the strong-memory regime the MB/GNP data exclude.  "
             "Continuous plateau is set by writer/eraser/turnover KINETICS (maps).", ha="center", fontsize=10, color="#555", style="italic")

    # ---- A1: Mk vs cycle period, one line per persistence level ----
    axA = fig.add_subplot(gs[0, 0])
    cmapP = plt.get_cmap('plasma')
    for i, pl in enumerate(PERSIST_LEVELS):
        pts = [(R['persist'][(pl, j)]['period'], R['persist'][(pl, j)]['mean'])
               for j in range(len(MUFAC)) if R['persist'].get((pl, j))]
        pts = [p for p in pts if not np.isnan(p[0])]
        if not pts:
            continue
        pts.sort()
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        lab = f"loss ×{pl:g}" + (" (baseline)" if pl == 1.0 else (" (persistent)" if pl < 1 else ""))
        axA.plot(xs, ys, "-o", ms=3.5, color=cmapP(i / (len(PERSIST_LEVELS) - 1)), lw=1.8, label=lab)
    axA.set_xlabel("cycle period (h)   ← faster    slower →"); axA.set_ylabel("settled mean Mk")
    axA.set_ylim(-0.03, 1.02); axA.grid(alpha=0.25)
    axA.set_title("A.  Mk vs cycle period, by mark persistence", fontsize=10.5, fontweight="bold")
    axA.legend(fontsize=7.0, loc="center right", framealpha=0.9)
    axA.text(0.03, 0.04, "baseline: nearly flat (buffered);\npersistent mark: steep (timing matters)",
             transform=axA.transAxes, fontsize=7.2, va="bottom", color="#555", style="italic")

    # ---- A2: timing-sensitivity vs factor -- persistence vs deposition ----
    axB = fig.add_subplot(gs[0, 1])
    pfac, psens = [], []
    for pl in PERSIST_LEVELS:
        a = R['persist'].get((pl, 0)); b = R['persist'].get((pl, len(MUFAC) - 1))
        if a and b:
            pfac.append(pl); psens.append(a['mean'] - b['mean'])   # slow-cycle minus fast-cycle
    kfac, ksens = [], []
    for kl in KMETH_LEVELS:
        a = R['kmeth'].get((kl, 0)); b = R['kmeth'].get((kl, 1))
        if a and b:
            kfac.append(kl); ksens.append(a['mean'] - b['mean'])
    op = np.argsort(pfac); ok = np.argsort(kfac)
    axB.plot(np.array(pfac)[op], np.array(psens)[op], "-o", color="#c0392b", lw=2, ms=5,
             label="make mark PERSISTENT (↓ loss: del_mk & eraser)")
    axB.plot(np.array(kfac)[ok], np.array(ksens)[ok], "-s", color="#2471a3", lw=2, ms=5,
             label="slow DEPOSITION (↓ k_meth)")
    axB.axhline(0, color="k", lw=0.8, ls=":")
    axB.axvline(1.0, color="k", lw=0.8, ls=":", alpha=0.5)
    axB.set_xscale("log"); axB.set_xlabel("kinetic factor (× baseline, log)")
    axB.set_ylabel("timing sensitivity  ΔMk (slow-cycle − fast-cycle)")
    axB.set_title("B.  What unlocks division-timing sensitivity", fontsize=10.5, fontweight="bold")
    axB.grid(alpha=0.25, which="both"); axB.legend(fontsize=7.6, loc="upper right", framealpha=0.9)
    axB.annotate("persistent →\ntiming matters", xy=(0.05, 0.15), fontsize=7.4, color="#c0392b", fontweight="bold")
    axB.annotate("slower deposition\njust lowers level\n(no timing effect)", xy=(0.1, 0.005), fontsize=7.0, color="#2471a3")

    # ---- interaction maps (B1 top-right, then bottom row) ----
    map_axes = [gs[0, 2], gs[1, 0], gs[1, 1], gs[1, 2]]
    for (spec, gax) in zip(MAPS, map_axes):
        name, px, prx, py, pry, title, xlab, ylab, xflip, _ = spec
        ax = fig.add_subplot(gax)
        XF = np.linspace(prx[0], prx[1], MAPN); YF = np.linspace(pry[0], pry[1], MAPN)
        Z = np.full((MAPN, MAPN), np.nan)
        for (ix, iy), o in R[name].items():
            if o:
                Z[iy, ix] = o['mean']
        im = ax.pcolormesh(XF, YF, Z, cmap="RdYlGn", vmin=0, vmax=1, shading="auto")
        cs = ax.contour(XF, YF, Z, levels=[0.25, 0.5], colors="k", linewidths=0.9, alpha=0.6)
        ax.clabel(cs, fmt="%.2f", fontsize=6.5)
        ax.plot(1.0, 1.0, "*", color="k", ms=13, mec="w", mew=1.0)
        ax.text(1.05, 1.05, "GNP", fontsize=7, fontweight="bold")
        if xflip:
            ax.invert_xaxis()
        ax.set_xlabel(xlab, fontsize=8); ax.set_ylabel(ylab, fontsize=8)
        ax.set_title(title, fontsize=9.6, fontweight="bold"); ax.tick_params(labelsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).set_label("mean Mk", fontsize=7.5)

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mk_timing_interactions")
    fig.savefig(out + ".png", dpi=145, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")


if __name__ == '__main__':
    main()
