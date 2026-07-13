"""Figure: the POSSIBILITY SPACE of the Ccnd1 H3K27me3 mark — extended bounds + structural knobs (v44.1).

NOT a calibration exercise. The point is to map what behaviors are REACHABLE when the mark's
kinetic AND structural variables are pushed well beyond the data-constrained regime. All points are
real dynamic sims. Two runtime structural knobs are injected on the read-write recruitment:
    PRC2 read-write term  a_rw_prc2 * Mk        ->   a_rw_prc2 * Mk^n_rw / (Krw_prc2^n_rw + Mk^n_rw)
so cooperativity (n_rw) and half-saturation (Krw) are tunable, and a0_prc2 (mark-independent
nucleation) can be driven to 0. n_rw=1 with a0 at baseline ~ the calibrated (monostable) model.

Findings mapped here:
  * calibrated architecture (linear read-write, a0>0): MONOSTABLE at all parameter values -- the
    low-init and high-init trajectories converge; no latching, no hysteresis.
  * cooperative read-write (n_rw>=2) + suppressed nucleation (a0->0): BISTABLE -- a latching
    epigenetic switch. The high-mark state is self-sustaining and ARRESTS the cycle; the low-mark
    state PROLIFERATES. Same parameters, two history-dependent fates. Driver ramps trace hysteresis.

Run:  ./venv/bin/python simulations/fig_v44_mk_possibility_space.py [--workers N] [--fresh]
"""
import os, sys, pickle, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, BoundaryNorm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

T_END, N_PTS, SETTLE, DIVS = 20000, 5000, 12000, 2000
SEL = ['time', 'Mk', 'Dna']
A0_0, ARW0, DEL0, KJ0, KTL0 = 0.000377635, 0.00429701, 0.0015, 0.047832059241740936, 0.0150989
BASE = dict(a0_prc2=A0_0, a_rw_prc2=ARW0, n_rw=1.0, Krw_prc2=0.4,
            del_mk=DEL0, k_jmjd3_gli=KJ0, kTlEZ=KTL0, SHH=0.5)
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mk_possibility_space_cache.pkl')

_RR = None
def _init(_):
    global _RR
    import tellurium as te
    m = build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_chain=False, with_two_step_rb=False)
    m = m.replace("a0_prc2 + a_rw_prc2*Mk)", "a0_prc2 + a_rw_prc2*Mk^n_rw/(Krw_prc2^n_rw + Mk^n_rw))")
    m = m.replace("Mk_methylation:", "n_rw = 1.0;\n  Krw_prc2 = 0.4;\n  Mk_methylation:", 1)
    _RR = te.loada(m)
    _RR.integrator.setValue('absolute_tolerance', 1e-8)
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000)
        _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass

def _set(ov):
    _RR.reset()
    for k, v in BASE.items():
        _RR[k] = v
    mk0 = ov.get('Mk0', 0.02)
    for k, v in ov.items():
        if k == 'Mk0':
            continue
        _RR[k] = v
    _RR['Mk'] = mk0

def _sim(ov, traj=False):
    for maxstep in (20.0, 5.0, 2.0):
        try:
            _set(ov); _RR.integrator.setValue('maximum_time_step', maxstep)
            d = _RR.simulate(0, T_END, N_PTS, selections=SEL); break
        except Exception:
            d = None
    if d is None:
        return None
    t = d['time']; s = t >= SETTLE; dna = d['Dna']
    ndiv = int(np.sum((dna[:-1] > 0.5) & (dna[1:] < 0.5) & (t[1:] >= DIVS)))
    out = dict(mean=float(d['Mk'][s].mean()), ndiv=ndiv)
    if traj:
        step = max(1, len(t) // 1400)
        out['t'] = (t[::step] / 60.0).tolist(); out['mk'] = d['Mk'][::step].tolist()
    return out

def _hyst(regime):
    """Quasi-static kTlEZ (EZH2) ramp up then down, state carried over -> hysteresis loop."""
    _set(regime); _RR.reset()
    for k, v in BASE.items():
        _RR[k] = v
    for k, v in regime.items():
        if k != 'Mk0':
            _RR[k] = v
    _RR['Mk'] = 0.02
    facs = np.concatenate([np.linspace(0.08, 3.2, 28), np.linspace(3.2, 0.08, 28)])
    up = np.ones(56, dtype=bool); up[28:] = False
    drv, mkv = [], []
    t0 = 0.0; last = 0.02
    for f in facs:
        _RR['kTlEZ'] = KTL0 * f
        val = last
        for ms in (20.0, 5.0, 2.0):
            try:
                _RR.integrator.setValue('maximum_time_step', ms)
                d = _RR.simulate(t0, t0 + 3500, 200, selections=['time', 'Mk'])
                val = float(np.mean(d['Mk'][-60:])); break
            except Exception:
                continue
        t0 += 3500; last = val
        drv.append(float(f)); mkv.append(val)
    return dict(drv=drv, mk=mkv, up=up.tolist())

def work(task):
    group, key, payload = task
    try:
        if group == 'hyst':
            return (group, key, _hyst(payload))
        if group == 'traj':
            return (group, key, _sim(payload, traj=True))
        return (group, key, _sim(payload))
    except Exception:
        return (group, key, None)


NMAP = 11
def build_tasks():
    tasks = []
    # A: extended a_rw two-IC sweep, monostable (n_rw=1) vs cooperative (n_rw=4,a0=0)
    for f in np.linspace(0.05, 12.0, 12):
        for mk0 in (0.02, 0.98):
            tasks.append(('A_mono', (round(f, 3), mk0), dict(a_rw_prc2=ARW0 * f, Mk0=mk0)))
            tasks.append(('A_coop', (round(f, 3), mk0), dict(a_rw_prc2=ARW0 * f, n_rw=4.0, a0_prc2=0.0, Mk0=mk0)))
    # B: bistability map  n_rw x a_rw  (a0=0)
    NR = np.linspace(1.0, 6.0, NMAP); AR = np.linspace(0.2, 20.0, NMAP)
    for ix, nr in enumerate(NR):
        for iy, ar in enumerate(AR):
            for mk0 in (0.02, 0.98):
                tasks.append(('B', (ix, iy, mk0), dict(n_rw=nr, a_rw_prc2=ARW0 * ar, a0_prc2=0.0, Mk0=mk0)))
    # C: nucleation map  a0 x a_rw  (n_rw=4)
    A0 = np.linspace(0.0, 3.0, NMAP); AR2 = np.linspace(0.2, 20.0, NMAP)
    for ix, a0f in enumerate(A0):
        for iy, ar in enumerate(AR2):
            for mk0 in (0.02, 0.98):
                tasks.append(('C', (ix, iy, mk0), dict(a0_prc2=A0_0 * a0f, n_rw=4.0, a_rw_prc2=ARW0 * ar, Mk0=mk0)))
    # D: two-fate trajectories at a bistable point
    biz = dict(n_rw=4.0, a0_prc2=0.0, a_rw_prc2=ARW0 * 4.0)
    tasks.append(('traj', 'low', {**biz, 'Mk0': 0.02}))
    tasks.append(('traj', 'high', {**biz, 'Mk0': 0.98}))
    # E: hysteresis loops (bistable needs a tiny seed a0>0 so EZH2 can nucleate the switch)
    tasks.append(('hyst', 'bistable', dict(n_rw=4.0, a0_prc2=A0_0 * 0.05, a_rw_prc2=ARW0 * 6.0)))
    tasks.append(('hyst', 'monostable', dict(n_rw=1.0)))
    return tasks, NR, AR, A0, AR2


def main():
    workers = 8
    fresh = '--fresh' in sys.argv
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    tasks, NR, AR, A0, AR2 = build_tasks()
    axes = dict(NR=NR, AR=AR, A0=A0, AR2=AR2)
    if not fresh and os.path.exists(CACHE):
        print('cached (--fresh to re-simulate)', flush=True)
        with open(CACHE, 'rb') as f:
            R, axes = pickle.load(f)
        plot(R, axes); return
    print(f'running {len(tasks)} real simulations on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=4)
    R = {}
    for g, k, o in res:
        R.setdefault(g, {})[k] = o
    with open(CACHE, 'wb') as f:
        pickle.dump((R, axes), f)
    plot(R, axes)


def _grid_gap(R, group, NX, NY):
    lo = np.full((NY, NX), np.nan); hi = np.full((NY, NX), np.nan); dvh = np.full((NY, NX), np.nan)
    for (ix, iy, mk0), o in R[group].items():
        if not o:
            continue
        if mk0 < 0.5:
            lo[iy, ix] = o['mean']
        else:
            hi[iy, ix] = o['mean']; dvh[iy, ix] = o['ndiv']
    return lo, hi, dvh


def plot(R, axes):
    NR, AR, A0, AR2 = axes['NR'], axes['AR'], axes['A0'], axes['AR2']
    fig = plt.figure(figsize=(19, 12))
    gs = GridSpec(2, 3, figure=fig, hspace=0.36, wspace=0.32, left=0.055, right=0.965, top=0.875, bottom=0.07)
    fig.suptitle("Possibility space of the $Ccnd1$ H3K27me3 mark — extended bounds + structural knobs (v44.1 dynamics; not data-constrained)",
                 fontsize=14.5, fontweight="bold", y=0.965)
    fig.text(0.5, 0.928, "Calibrated architecture is monostable (buffered). Cooperative read-write (n_rw≥2) + suppressed nucleation (a0→0) unlocks a BISTABLE latching switch "
             "coupled to cell-cycle arrest — a reachable behavior the data-fit model does not use.", ha="center", fontsize=9.5, color="#555", style="italic")

    # ---- A: extended two-IC sweep, monostable vs cooperative ----
    axA = fig.add_subplot(gs[0, 0])
    for grp, col, lab in [('A_mono', '#2471a3', 'linear read-write (calibrated form)'),
                          ('A_coop', '#c0392b', 'cooperative n_rw=4, a0→0')]:
        fs = sorted({k[0] for k in R[grp]})
        loy = [R[grp].get((f, 0.02), {}).get('mean', np.nan) for f in fs]
        hiy = [R[grp].get((f, 0.98), {}).get('mean', np.nan) for f in fs]
        axA.plot(fs, loy, '-o', color=col, ms=3.5, lw=1.6, label=f"{lab} — from Mk₀=0")
        axA.plot(fs, hiy, '--s', color=col, ms=3.5, lw=1.6, mfc='none', label=f"{lab} — from Mk₀=1")
    axA.set_xlabel("read-write strength a_rw  (× baseline)"); axA.set_ylabel("settled mean Mk")
    axA.set_ylim(-0.03, 1.03); axA.grid(alpha=0.25)
    axA.set_title("A.  Extended read-write sweep: two initial conditions", fontsize=10.2, fontweight="bold")
    axA.legend(fontsize=6.6, loc="center right", framealpha=0.9)
    axA.text(0.03, 0.55, "linear: Mk₀=0 and Mk₀=1\ncollapse together (MONOSTABLE)", transform=axA.transAxes,
             fontsize=6.8, color='#2471a3', style='italic')
    axA.text(0.30, 0.14, "cooperative: split =\nBISTABLE / latching", transform=axA.transAxes,
             fontsize=6.8, color='#c0392b', style='italic', fontweight='bold')

    # ---- B: bistability map n_rw x a_rw ----
    axB = fig.add_subplot(gs[0, 1])
    lo, hi, _ = _grid_gap(R, 'B', len(NR), len(AR))
    gap = np.abs(hi - lo)
    im = axB.pcolormesh(NR, AR, gap, cmap='magma', vmin=0, vmax=1, shading='auto')
    cs = axB.contour(NR, AR, gap, levels=[0.15], colors='cyan', linewidths=1.6)
    axB.set_xlabel("read-write cooperativity  n_rw"); axB.set_ylabel("read-write strength a_rw (× base)")
    axB.set_title("B.  Latching (hysteresis gap) vs cooperativity × strength\n(a0=0)", fontsize=9.8, fontweight="bold")
    fig.colorbar(im, ax=axB, fraction=0.046, pad=0.02).set_label("|Mk(hi-init) − Mk(lo-init)|", fontsize=7.5)
    axB.text(0.04, 0.94, "cyan = bistable boundary", transform=axB.transAxes, fontsize=7, color='cyan', va='top')
    axB.text(0.5, 0.06, "n_rw=1 (linear): never latches", transform=axB.transAxes, fontsize=7, color='w', ha='center')

    # ---- C: nucleation map a0 x a_rw ----
    axC = fig.add_subplot(gs[0, 2])
    lo, hi, _ = _grid_gap(R, 'C', len(A0), len(AR2))
    gap = np.abs(hi - lo)
    im = axC.pcolormesh(A0, AR2, gap, cmap='magma', vmin=0, vmax=1, shading='auto')
    axC.contour(A0, AR2, gap, levels=[0.15], colors='cyan', linewidths=1.6)
    axC.axvline(1.0, color='w', ls=':', lw=1); axC.text(1.05, AR2[-2], 'baseline a0', color='w', fontsize=7)
    axC.set_xlabel("mark-independent nucleation  a0 (× base)"); axC.set_ylabel("read-write strength a_rw (× base)")
    axC.set_title("C.  Latching needs low nucleation (n_rw=4)", fontsize=9.8, fontweight="bold")
    fig.colorbar(im, ax=axC, fraction=0.046, pad=0.02).set_label("hysteresis gap", fontsize=7.5)

    # ---- D: two-fate trajectories ----
    axD = fig.add_subplot(gs[1, 0])
    for key, col, lab in [('low', '#1e8449', 'Mk₀=0 → PROLIFERATES (Mk stays low)'),
                          ('high', '#8e44ad', 'Mk₀=1 → LATCHES + ARRESTS (Mk high)')]:
        o = R['traj'].get(key)
        if o:
            t = np.array(o['t']); mk = np.array(o['mk'])
            axD.plot(t, mk, color=col, lw=1.8, label=f"{lab}  ({o['ndiv']} div)")
    axD.set_xlabel("time (h)"); axD.set_ylabel("Mk"); axD.set_ylim(-0.03, 1.03); axD.grid(alpha=0.25)
    axD.set_title("D.  Same parameters, two fates (bistable point)", fontsize=10.2, fontweight="bold")
    axD.legend(fontsize=7.2, loc="center right", framealpha=0.9)

    # ---- E: hysteresis loops ----
    axE = fig.add_subplot(gs[1, 1])
    for key, col, lab in [('monostable', '#2471a3', 'monostable mark (n_rw=1): small loop = cell-cycle memory'),
                          ('bistable', '#c0392b', 'bistable mark (n_rw=4): large latch = epigenetic memory')]:
        o = R['hyst'].get(key)
        if not o:
            continue
        drv = np.array(o['drv']); mk = np.array(o['mk']); up = np.array(o['up'])
        axE.plot(drv[up], mk[up], '-', color=col, lw=1.8, label=f"{lab}")
        axE.plot(drv[~up], mk[~up], '--', color=col, lw=1.8, alpha=0.85)
    axE.annotate("", xy=(0.7, 0.9), xytext=(1.6, 0.9), arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
    axE.text(1.15, 0.93, "latched ON\n(EZH2 removed,\nmark persists)", fontsize=6.6, color="#c0392b", ha="center", va="bottom")
    axE.set_xlabel("EZH2 abundance kTlEZ (× base)  — ramped up (—) then down (– –)"); axE.set_ylabel("Mk")
    axE.set_ylim(-0.03, 1.06); axE.grid(alpha=0.25)
    axE.set_title("E.  Driver-ramp hysteresis (memory)", fontsize=10.2, fontweight="bold")
    axE.legend(fontsize=6.6, loc="center left", framealpha=0.9)

    # ---- F: regime classification (n_rw x a_rw grid from B) ----
    axF = fig.add_subplot(gs[1, 2])
    lo, hi, dvh = _grid_gap(R, 'B', len(NR), len(AR))
    Z = np.full_like(lo, np.nan)
    for iy in range(len(AR)):
        for ix in range(len(NR)):
            l, h, dv = lo[iy, ix], hi[iy, ix], dvh[iy, ix]
            if np.isnan(l):
                continue
            if abs(h - l) > 0.15:
                Z[iy, ix] = 2                      # bistable
            elif h < 0.12:
                Z[iy, ix] = 0                      # collapsed
            elif h > 0.85 and (dv == 0):
                Z[iy, ix] = 3                      # saturated/arrested
            else:
                Z[iy, ix] = 1                      # monostable graded
    cmap = ListedColormap(['#d5dbdb', '#7dcea0', '#e74c3c', '#5b2c6f'])
    im = axF.pcolormesh(NR, AR, Z, cmap=cmap, norm=BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5], cmap.N), shading='auto')
    axF.set_xlabel("read-write cooperativity  n_rw"); axF.set_ylabel("read-write strength a_rw (× base)")
    axF.set_title("F.  Regime map (a0=0)", fontsize=10.2, fontweight="bold")
    cb = fig.colorbar(im, ax=axF, fraction=0.046, pad=0.02, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(['collapsed', 'monostable\ngraded', 'BISTABLE\nlatching', 'saturated\n/arrested'], fontsize=7)

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mk_possibility_space")
    fig.savefig(out + ".png", dpi=145, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")


if __name__ == '__main__':
    main()
