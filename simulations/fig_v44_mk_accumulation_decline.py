"""Figure: what drives ACCUMULATION vs DECLINE of H3K27me3 (Mk) at Ccnd1 over divisions (v44.1).

Every value is a REAL simulation of the baked v44.1 model (build_model_v44, with_prc2 default) run in
a proliferating GNP-baseline condition (SHH 0.5, Ptch1 1.0); the mark is followed across the
replicative sawtooth (Mk halved each S) and summarized by its settled mean + sawtooth min/max.

Panels:
  A  representative Mk(t) trajectories over divisions (accumulate vs decline)
  B  settled mean Mk by biological condition (GNP / MB / perturbations)
  C  writer sweep      : EZH2 abundance (kTlEZ)          -> Mk up
  D  eraser sweep      : Gli->Kdm6b rate (k_jmjd3_gli)   -> Mk down
  E  turnover sweep    : del_mk                          -> Mk down
  F  read-write sweep  : a_rw_prc2                       -> Mk up (amplifier)
  G  mitogen sweep     : SHH (raises Gli1 eraser AND cycling/dilution) + Gli1 + period
  H  writer x eraser PHASE DIAGRAM: mean Mk heatmap, accumulate vs decline basin

Run:  ./venv/bin/python simulations/fig_v44_mk_accumulation_decline.py [--workers N]
"""
import os, sys, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

# ---- baked v44.1 baselines (reset every call: reset() restores species only, not parameters) ----
KTLEZ0 = 0.0150989
KJMJD30 = 0.047832059241740936
DELMK0 = 0.0015
ARW0 = 0.004297010112265872
BASE_PARAMS = dict(kTlEZ=KTLEZ0, k_jmjd3_gli=KJMJD30, del_mk=DELMK0, a_rw_prc2=ARW0,
                   SHH=0.5, Ptch1_copy_number=1.0, EZH2i=0.0, MYCN_amplification=1.0,
                   p16=0.0, p18=0.4, kSyP21=0.002)
T_END, N_PTS, SETTLE = 16000, 6400, 8000
SEL = ['time', 'Mk', 'Dna', 'MPF', 'EZH2', 'Gli1']

_RR = None
def _init(_):
    global _RR
    _RR = build_model_v44(with_ezh2=True, with_hh=True)
    import tellurium as te
    _RR = te.loada(_RR)
    _RR.integrator.setValue('absolute_tolerance', 1e-8)
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000)
        _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass

def _simulate(ov):
    """Full-baseline reset + task overrides -> simulate. Returns the SEL result dict."""
    _RR.reset()
    for k, v in BASE_PARAMS.items():
        _RR[k] = v
    for k, v in (ov or {}).items():
        _RR[k] = v
    for maxstep in (20.0, 5.0, 2.0):
        try:
            _RR.integrator.setValue('maximum_time_step', maxstep)
            return _RR.simulate(0, T_END, N_PTS, selections=SEL)
        except Exception:
            _RR.reset()
            for k, v in BASE_PARAMS.items():
                _RR[k] = v
            for k, v in (ov or {}).items():
                _RR[k] = v
    return None

def _ndiv(d, settle=SETTLE):
    t = d['time']; dna = d['Dna']
    return int(np.sum((dna[:-1] > 0.5) & (dna[1:] < 0.5) & (t[1:] >= settle)))

def _period_h(d, settle=SETTLE):
    t = d['time']; dna = d['Dna']
    idx = np.where((dna[:-1] > 0.5) & (dna[1:] < 0.5) & (t[1:] >= settle))[0]
    if len(idx) < 2:
        return np.nan
    return float(np.mean(np.diff(t[idx])) / 60.0)

def work(task):
    kind, key, ov = task
    d = _simulate(ov)
    if d is None:
        return (kind, key, None)
    t = d['time']; s = t >= SETTLE; mk = d['Mk']
    out = dict(mean=float(mk[s].mean()), lo=float(mk[s].min()), hi=float(mk[s].max()),
               ndiv=_ndiv(d), period=_period_h(d),
               ezh2=float(d['EZH2'][s].mean()), gli1=float(d['Gli1'][s].mean()))
    if kind == 'traj':
        # downsample the trajectory for plotting
        step = max(1, len(t) // 1600)
        out['t'] = (t[::step] / 60.0).tolist(); out['mk'] = mk[::step].tolist()
        out['dna'] = d['Dna'][::step].tolist()
    return (kind, key, out)


def build_tasks():
    tasks = []
    # A: representative trajectories (accumulate -> decline)
    TRAJ = {
        'GNP (accumulate)':        {},
        'high EZH2 (2x writer)':   dict(kTlEZ=KTLEZ0 * 2.0),
        'MB (high Gli1 eraser)':   dict(Ptch1_copy_number=0.09, MYCN_amplification=6.0, p16=0.306, p18=1.553, kSyP21=0.004),
        'high eraser (4x)':        dict(k_jmjd3_gli=KJMJD30 * 4.0),
        'EZH2i (writer off)':      dict(EZH2i=1.0),
    }
    for k, ov in TRAJ.items():
        tasks.append(('traj', k, ov))
    # B: condition bars
    BAR = {
        'GNP': {}, 'MB': dict(Ptch1_copy_number=0.09, MYCN_amplification=6.0, p16=0.306, p18=1.553, kSyP21=0.004),
        '+EZH2i': dict(EZH2i=1.0), 'EZH2 x2': dict(kTlEZ=KTLEZ0 * 2.0), 'EZH2 x0.5': dict(kTlEZ=KTLEZ0 * 0.5),
        'eraser x4': dict(k_jmjd3_gli=KJMJD30 * 4.0), 'turnover x6': dict(del_mk=DELMK0 * 6.0),
        'no read-write': dict(a_rw_prc2=0.0),
    }
    for k, ov in BAR.items():
        tasks.append(('bar', k, ov))
    # C: writer sweep
    for i, f in enumerate(np.linspace(0.2, 2.2, 13)):
        tasks.append(('writer', f, dict(kTlEZ=KTLEZ0 * f)))
    # D: eraser sweep
    for f in np.linspace(0.0, 8.0, 13):
        tasks.append(('eraser', f, dict(k_jmjd3_gli=KJMJD30 * f)))
    # E: turnover sweep
    for f in np.linspace(0.2, 12.0, 13):
        tasks.append(('turnover', f, dict(del_mk=DELMK0 * f)))
    # F: read-write sweep
    for f in np.linspace(0.0, 2.5, 13):
        tasks.append(('readwrite', f, dict(a_rw_prc2=ARW0 * f)))
    # G: broken-Hh-feedback sweep (Ptch1 copy number 1.0 -> 0.05 = GNP -> MB): raises Gli1 eraser
    for pt in np.linspace(1.0, 0.06, 12):
        tasks.append(('ptch', float(pt), dict(Ptch1_copy_number=float(pt))))
    # H: writer x eraser phase diagram
    WX = np.linspace(0.3, 2.0, 14); EY = np.linspace(0.0, 6.0, 14)
    for ix, wf in enumerate(WX):
        for iy, ef in enumerate(EY):
            tasks.append(('phase', (ix, iy), dict(kTlEZ=KTLEZ0 * wf, k_jmjd3_gli=KJMJD30 * ef)))
    return tasks, WX, EY


def main():
    workers = 8
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    tasks, WX, EY = build_tasks()
    print(f'running {len(tasks)} real simulations on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=4)
    R = {}
    for kind, key, out in res:
        R.setdefault(kind, {})[key] = out
    plot(R, WX, EY)


# ----------------------------------------------------------------------------------------
GREEN, RED, PURP, TEAL, GREY = "#1e8449", "#c0392b", "#8e44ad", "#148f77", "#7f8c8d"

def _sweep_xy(dct):
    xs = sorted(dct.keys())
    mean = np.array([dct[x]['mean'] if dct[x] else np.nan for x in xs])
    lo = np.array([dct[x]['lo'] if dct[x] else np.nan for x in xs])
    hi = np.array([dct[x]['hi'] if dct[x] else np.nan for x in xs])
    ndv = np.array([dct[x]['ndiv'] if dct[x] else 0 for x in xs])
    return np.array(xs), mean, lo, hi, ndv

def _panel_sweep(ax, dct, xlabel, title, color, updown, xscale=1.0, xunit=""):
    xs, mean, lo, hi, ndv = _sweep_xy(dct)
    x = xs * xscale
    ax.fill_between(x, lo, hi, color=color, alpha=0.16, lw=0, label="sawtooth min–max")
    ax.plot(x, mean, "-o", color=color, ms=4, lw=2, label="settled mean Mk")
    arr = ndv == 0  # arrested (no divisions -> no dilution)
    if arr.any():
        ax.plot(x[arr], mean[arr], "s", color="#e67e22", ms=7, mfc="none", mew=1.6,
                label="arrested (no dilution)")
    ax.set_xlabel(xlabel + (f"  {xunit}" if xunit else "")); ax.set_ylabel("Mk (H3K27me3 @ $Ccnd1$)")
    ax.set_ylim(-0.03, 1.02); ax.set_title(title, fontsize=10.5, fontweight="bold")
    ax.grid(alpha=0.25); ax.axhline(0.0, color=GREY, lw=0.8, ls=":")
    ax.text(0.02, 0.95, updown, transform=ax.transAxes, fontsize=8.5, va="top", color=color,
            fontweight="bold")

def plot(R, WX, EY):
    fig = plt.figure(figsize=(16.5, 13.2))
    gs = GridSpec(3, 4, figure=fig, height_ratios=[1.05, 1.0, 1.05], hspace=0.42, wspace=0.34,
                  left=0.06, right=0.975, top=0.9, bottom=0.06)
    fig.suptitle("What drives accumulation vs decline of H3K27me3 at $Ccnd1$ over divisions  —  v44.1 (all points = real simulations)",
                 fontsize=15, fontweight="bold", y=0.965)
    fig.text(0.5, 0.925, "proliferating GNP baseline (SHH 0.5, Ptch1 1.0); mark followed across the replicative sawtooth (Mk ÷2 each S).  "
             "ACCUMULATE when writer flux > (turnover + erase + dilution); DECLINE when losses win.",
             ha="center", fontsize=9.5, color="#555", style="italic")

    # ---- A: trajectories ----
    axA = fig.add_subplot(gs[0, 0:2])
    order = ['high EZH2 (2x writer)', 'GNP (accumulate)', 'MB (high Gli1 eraser)', 'high eraser (4x)', 'EZH2i (writer off)']
    cols = ['#0e6655', GREEN, '#b9770e', RED, '#7b241c']
    for k, c in zip(order, cols):
        o = R['traj'].get(k)
        if not o:
            continue
        t = np.array(o['t']); mk = np.array(o['mk'])
        m = t >= 60  # skip first hour transient
        axA.plot(t[m], mk[m], color=c, lw=1.6, label=f"{k}  (⟨Mk⟩={o['mean']:.2f}, {o['ndiv']} div)")
    axA.set_xlabel("time (h)"); axA.set_ylabel("Mk (H3K27me3 @ $Ccnd1$)")
    axA.set_ylim(-0.03, 1.02); axA.set_title("A.  Mk trajectories over divisions", fontsize=10.5, fontweight="bold")
    axA.grid(alpha=0.25); axA.legend(fontsize=7.4, loc="center right", framealpha=0.9)
    axA.text(0.015, 0.03, "each drop = replicative dilution (÷2 at S);\nrise between = PRC2 re-methylation",
             transform=axA.transAxes, fontsize=7.6, va="bottom", color=GREY, style="italic")

    # ---- B: condition bars ----
    axB = fig.add_subplot(gs[0, 2:4])
    border = ['GNP', 'EZH2 x2', 'EZH2 x0.5', 'no read-write', 'turnover x6', 'eraser x4', 'MB', '+EZH2i']
    vals = [R['bar'][k]['mean'] if R['bar'].get(k) else 0 for k in border]
    ndv = [R['bar'][k]['ndiv'] if R['bar'].get(k) else 0 for k in border]
    base = R['bar']['GNP']['mean']
    bcol = [GREEN if v >= base - 0.02 else (RED if v < base * 0.7 else "#e67e22") for v in vals]
    bars = axB.bar(range(len(border)), vals, color=bcol, alpha=0.85, edgecolor="#2c3e50", lw=0.7)
    axB.axhline(base, color=GREEN, ls="--", lw=1.2, alpha=0.7)
    axB.text(len(border) - 0.4, base + 0.01, "GNP baseline", fontsize=7.6, color=GREEN, ha="right")
    for i, (v, nd) in enumerate(zip(vals, ndv)):
        axB.text(i, v + 0.015, f"{v:.2f}", ha="center", fontsize=7.6)
        if nd == 0:
            axB.text(i, 0.03, "arrest", ha="center", fontsize=6.8, color="#e67e22", rotation=90, va="bottom")
    axB.set_xticks(range(len(border))); axB.set_xticklabels(border, rotation=30, ha="right", fontsize=8)
    axB.set_ylabel("settled mean Mk"); axB.set_ylim(0, 1.02)
    axB.set_title("B.  Mk by condition / perturbation", fontsize=10.5, fontweight="bold")
    axB.grid(alpha=0.2, axis="y")

    # ---- C-F: 1-D sweeps ----
    axC = fig.add_subplot(gs[1, 0]); _panel_sweep(axC, R['writer'], "EZH2 abundance  (× baseline kTlEZ)",
        "C.  WRITER  ↑ → accumulate", PURP, "writer ↑  →  Mk ↑")
    axD = fig.add_subplot(gs[1, 1]); _panel_sweep(axD, R['eraser'], "Gli→Kdm6b rate  (× baseline k_jmjd3_gli)",
        "D.  ERASER  ↑ → decline", TEAL, "eraser ↑  →  Mk ↓")
    axE = fig.add_subplot(gs[1, 2]); _panel_sweep(axE, R['turnover'], "turnover  (× baseline del_mk)",
        "E.  TURNOVER  ↑ → decline", "#2471a3", "turnover ↑  →  Mk ↓")
    axF = fig.add_subplot(gs[1, 3]); _panel_sweep(axF, R['readwrite'], "read-write  (× baseline a_rw_prc2)",
        "F.  READ-WRITE  ↑ → accumulate", "#7d3c98", "a_rw ↑  →  Mk ↑")
    for ax in (axC, axD, axE, axF):
        ax.legend(fontsize=6.6, loc="upper right", framealpha=0.9)

    # ---- G: broken-Hh-feedback (Ptch1) sweep = GNP -> MB, raises Gli1 eraser ----
    axG = fig.add_subplot(gs[2, 0:2])
    xs, mean, lo, hi, ndv = _sweep_xy(R['ptch'])
    gli = np.array([R['ptch'][x]['gli1'] for x in xs])
    ez = np.array([R['ptch'][x]['ezh2'] for x in xs])
    axG.fill_between(xs, lo, hi, color=GREEN, alpha=0.15, lw=0)
    axG.plot(xs, mean, "-o", color=GREEN, ms=4, lw=2, label="mean Mk")
    arr = ndv == 0
    if arr.any():
        axG.plot(xs[arr], mean[arr], "s", color="#e67e22", ms=7, mfc="none", mew=1.6, label="arrested (no dilution)")
    axG.set_xlabel("Ptch1 copy number   (1.0 = GNP, intact feedback  →  ~0.1 = MB, broken)")
    axG.set_ylabel("Mk", color=GREEN); axG.set_ylim(-0.03, 1.02)
    axG.tick_params(axis="y", labelcolor=GREEN); axG.invert_xaxis()
    axG.set_title("G.  Broken Hh feedback (Ptch1↓, GNP→MB): Gli1 eraser ↑ beats writer ↑ → Mk ↓", fontsize=9.8, fontweight="bold")
    axG.grid(alpha=0.22)
    axG2 = axG.twinx()
    axG2.plot(xs, gli, "-^", color="#b9770e", ms=3.5, lw=1.6, label="Gli1 (eraser drive)")
    axG2.plot(xs, ez, ":d", color=PURP, ms=3.5, lw=1.5, label="EZH2 (writer)")
    axG2.set_ylabel("Gli1  /  EZH2", color="#b9770e"); axG2.tick_params(axis="y", labelcolor="#b9770e")
    l1, la1 = axG.get_legend_handles_labels(); l2, la2 = axG2.get_legend_handles_labels()
    axG.legend(l1 + l2, la1 + la2, fontsize=7.2, loc="upper right", framealpha=0.9)

    # ---- H: writer x eraser phase diagram ----
    axH = fig.add_subplot(gs[2, 2:4])
    Z = np.full((len(EY), len(WX)), np.nan)
    for (ix, iy), o in R['phase'].items():
        if o:
            Z[iy, ix] = o['mean']
    im = axH.pcolormesh(WX, EY, Z, cmap="RdYlGn", vmin=0, vmax=1, shading="auto")
    cs = axH.contour(WX, EY, Z, levels=[0.25, 0.5], colors="k", linewidths=1.0, alpha=0.6)
    axH.clabel(cs, fmt="Mk=%.2f", fontsize=7)
    axH.plot(1.0, 1.0, "*", color="k", ms=15, mec="w", mew=1.0)
    axH.text(1.05, 1.05, "GNP\nbaseline", fontsize=7.6, fontweight="bold")
    axH.set_xlabel("WRITER  →  EZH2 abundance (× kTlEZ)")
    axH.set_ylabel("ERASER  →  Gli→Kdm6b rate (× k_jmjd3_gli)")
    axH.set_title("H.  Writer × eraser phase diagram (mean Mk)", fontsize=10.5, fontweight="bold")
    cb = fig.colorbar(im, ax=axH, fraction=0.046, pad=0.02); cb.set_label("settled mean Mk")
    axH.text(0.03, 0.95, "ACCUMULATE", transform=axH.transAxes, fontsize=9, fontweight="bold",
             color="#145a32", va="top")
    axH.text(0.62, 0.06, "DECLINE", transform=axH.transAxes, fontsize=9, fontweight="bold",
             color="#7b241c", va="bottom")

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mk_accumulation_decline")
    fig.savefig(out + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")


if __name__ == '__main__':
    main()
