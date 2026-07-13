"""Figure: DEFINING transient G0 vs immediate G1/early-S entry, and its CyclinD1/p27/CDKi control (v44).

Operational definition (per daughter, at each division E_div):
  * born-state = (CyclinD1 setpoint, birth p27 = P21_div) -> the CyclinD1/p27 ratio (Fan-Meyer).
  * G0 dwell  = time from birth until COMMITMENT = first time p27 (P21) drops below P27_THR (the
    Skp2-p27 toggle flips: E2f/Skp2 up, CDK2 freed). Marked by p27, NOT pRb (pRb saturates ~20min
    post-division, long before functional commitment).
  * IMMEDIATE G1/early-S : dwell <= TAU_MIN (born committed / CDK2-inc; Spencer carryover keeps pRb/CycE/low-p27).
  * TRANSIENT G0         : TAU_MIN < dwell < inf  (dwells in the p27-high/low-CDK2 basin, then commits).
  * SUSTAINED G0 / ARREST: never commits (no S-entry / <2 divisions) -- CyclinD1 below the commitment threshold.

Panels: (A) three single-cell exemplars showing p27, CDK2-activity proxy (free Ce+Ca), pRb (bad marker),
Dna; G0 window shaded. (B) CyclinD1 x birth-p27 phase diagram of outcome (the Fan-Meyer decision boundary).
Run: ./venv/bin/python simulations/fig_v44_transient_g0_definition.py [--workers N]
"""
import os, sys, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, BoundaryNorm
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

KTL0 = 0.801                      # baseline k_Cd_translation (CyclinD1 setpoint)
P27_THR = 0.1                     # p27 (P21) commitment threshold (same as validate classify)
TAU_MIN = 2.0                     # h: dwell below this = "immediate" (no dwell beyond minimal G1)
T_END, N_PTS, SETTLE = 12000, 24000, 4000
SEL = ['time', 'P21', 'Ce', 'Ca', 'CeP21', 'CaP21', 'pRb', 'Dna', 'aRc', 'MPF', 'mass', 'Cd']

GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002)
MB = dict(SHH=0.5, Ptch1_copy_number=0.09, MYCN_amplification=6.0, p16=0.306, p18=1.553, kSyP21=0.002)

_RR = None
def _init(_):
    global _RR
    import tellurium as te
    _RR = te.loada(build_model_v44())     # default = chain; cell-cycle machinery unaffected
    _RR.integrator.setValue('absolute_tolerance', 1e-9)
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 2000000)
        _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass

def run_cell(cond, traj=False):
    _RR.reset()
    for k, v in cond.items():
        _RR[k] = v
    for maxstep in (20.0, 5.0, 2.0):
        try:
            _RR.integrator.setValue('maximum_time_step', maxstep)
            d = _RR.simulate(0, T_END, N_PTS, selections=SEL); break
        except Exception:
            d = None
    if d is None:
        return None
    t = d['time']; dt = t[1] - t[0]
    Dna = d['Dna']; P21 = d['P21']; aRc = d['aRc']
    # division events = Dna reset 1->0 (E_div)
    div_idx = np.where((Dna[:-1] > 0.5) & (Dna[1:] < 0.1))[0] + 1
    div_idx = div_idx[t[div_idx] >= SETTLE]
    # per-daughter G0 dwell: from birth to first commitment (p27 < thr) before next division
    dwells = []
    for i, b in enumerate(div_idx[:-1]):
        nxt = div_idx[i + 1]
        seg = slice(b, nxt)
        p = P21[seg]
        below = np.where(p < P27_THR)[0]
        dwell_h = (below[0] * dt) / 60.0 if len(below) else (nxt - b) * dt / 60.0
        dwells.append(dwell_h)
    ndiv = len(div_idx)
    if ndiv < 2:
        cls, dwell = 'arrest', np.nan
    else:
        dwell = float(np.mean(dwells)) if dwells else np.nan
        cls = 'immediate' if dwell <= TAU_MIN else 'transient_G0'
    out = dict(cls=cls, dwell=dwell, ndiv=int(ndiv))
    if traj:
        m = t >= SETTLE
        step = max(1, int(np.sum(m)) // 2400)
        idx = np.where(m)[0][::step]
        out['t'] = (t[idx] / 60.0).tolist()
        for s in ['P21', 'pRb', 'Dna', 'aRc']:
            out[s] = d[s][idx].tolist()
        out['CDK2'] = (d['Ce'][idx] + d['Ca'][idx]).tolist()   # free (active) CDK2 proxy
        out['div_t'] = (t[div_idx] / 60.0).tolist()
    return out

def work(task):
    kind, key, cond = task
    return (kind, key, run_cell(cond, traj=(kind == 'exemplar')))


def build_tasks():
    tasks = []
    # A: exemplars -- pick conditions that give the 3 classes
    # transient G0 is a near-threshold / high-birth-p27 phenomenon: a deterministic cell at the
    # setpoint commits fast (carryover keeps p27 low), so we force each class via birth p27 / mitogen.
    EX = {
        'immediate':    {**GNP, 'EZH2i': 1.0, 'P21_div': 0.2},           # low birth p27 + de-repressed CyclinD1
        'transient_G0': {**GNP, 'P21_div': 1.6},                         # high birth p27 -> dwells, then commits
        'arrest':       {**GNP, 'SHH': 0.0, 'P21_div': 1.6},             # no mitogen -> CyclinD1 collapses -> never commits
    }
    for k, c in EX.items():
        tasks.append(('exemplar', k, c))
    # B: CyclinD1 x birth-p27 phase diagram (GNP background)
    CDF = np.linspace(0.25, 1.6, 16)      # x k_Cd_translation baseline
    P21D = np.linspace(0.1, 2.6, 16)      # birth p27
    for ix, cf in enumerate(CDF):
        for iy, p in enumerate(P21D):
            tasks.append(('phase', (ix, iy), {**GNP, 'k_Cd_translation': KTL0 * cf, 'P21_div': float(p)}))
    return tasks, CDF, P21D


def main():
    workers = 8
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    tasks, CDF, P21D = build_tasks()
    print(f'running {len(tasks)} single-cell sims on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=4)
    R = {}
    for kind, key, o in res:
        R.setdefault(kind, {})[key] = o
    plot(R, CDF, P21D)


def plot(R, CDF, P21D):
    fig = plt.figure(figsize=(16, 9))
    gs = GridSpec(3, 3, figure=fig, hspace=0.5, wspace=0.34, left=0.06, right=0.965, top=0.88, bottom=0.08)
    fig.suptitle("Transient G0 vs immediate G1/early-S entry — definition, mechanism, and CyclinD1/p27 control (v44)",
                 fontsize=14, fontweight="bold", y=0.965)
    fig.text(0.5, 0.925, "G0 = p27-high / pre-S dwell after birth (NOT pRb, which saturates early).  "
             "immediate: dwell ≤ 2 h · transient G0: commits after a dwell · arrest: never commits.",
             ha="center", fontsize=9.5, color="#555", style="italic")

    order = [('immediate', '#1e8449', 'A1. IMMEDIATE (born committed)'),
             ('transient_G0', '#b9770e', 'A2. TRANSIENT G0 (dwell, then commit)'),
             ('arrest', '#c0392b', 'A3. ARREST (never commits)')]
    for col, (key, c, title) in enumerate(order):
        ax = fig.add_subplot(gs[0, col])
        o = R['exemplar'].get(key)
        if not o or 't' not in o:
            ax.set_title(title + ' (n/a)', fontsize=9); continue
        t = np.array(o['t'])
        ax.plot(t, o['P21'], color='#8e44ad', lw=1.6, label='p27 (=P21)')
        ax.plot(t, o['CDK2'], color='#2471a3', lw=1.5, label='CDK2 (free Ce+Ca)')
        ax.plot(t, np.array(o['pRb']) / max(np.max(o['pRb']), 1e-6), color='#95a5a6', lw=1.2, ls=':', label='pRb (norm, bad marker)')
        ax.plot(t, o['Dna'], color='#17a2b8', lw=1.1, alpha=0.8, label='Dna (S/G2)')
        ax.axhline(P27_THR, color='#8e44ad', lw=0.8, ls='--', alpha=0.5)
        # shade p27-high pre-S (G0) windows
        p = np.array(o['P21']); dna = np.array(o['Dna']); arc = np.array(o['aRc'])
        g0 = (p > P27_THR) & (dna < 0.05) & (arc < 0.05)
        ax.fill_between(t, 0, 1.05, where=g0, color='#f9e79f', alpha=0.5, lw=0, transform=ax.get_xaxis_transform())
        ax.set_ylim(-0.03, 1.08); ax.set_xlabel('time (h)'); ax.set_title(title, fontsize=9.2, fontweight='bold', color=c)
        dw = 'no division' if o['cls'] == 'arrest' else f"dwell ≈ {o['dwell']:.1f} h, {o['ndiv']} div"
        ax.text(0.02, 0.02, dw, transform=ax.transAxes, fontsize=7.6, va='bottom', color=c, fontweight='bold')
        if col == 0:
            ax.legend(fontsize=6.4, loc='upper right', framealpha=0.9)
        ax.grid(alpha=0.2)
    fig.text(0.06, 0.60, "shaded = G0 (p27-high, pre-S). Note pRb (grey dotted) jumps to max right after "
             "division in ALL three — it can't distinguish the classes; p27 and CDK2-activity can.",
             fontsize=8, color="#555", style="italic")

    # B: phase diagram
    axB = fig.add_subplot(gs[1:, 0:2])
    Z = np.full((len(P21D), len(CDF)), np.nan)
    code = {'immediate': 0, 'transient_G0': 1, 'arrest': 2}
    for (ix, iy), o in R['phase'].items():
        if o:
            Z[iy, ix] = code[o['cls']]
    cmap = ListedColormap(['#7dcea0', '#f5b041', '#e74c3c'])
    im = axB.pcolormesh(CDF, P21D, Z, cmap=cmap, norm=BoundaryNorm([-.5, .5, 1.5, 2.5], cmap.N), shading='auto')
    axB.plot(1.0, 0.6, '*', color='k', ms=16, mec='w', mew=1.2); axB.text(1.03, 0.6, ' GNP\n setpoint', fontsize=8, fontweight='bold')
    axB.set_xlabel("CyclinD1 setpoint  (× baseline k_Cd_translation)  →", fontsize=9.5)
    axB.set_ylabel("birth p27  (P21_div)  →", fontsize=9.5)
    axB.set_title("B. Commitment phase diagram: the CyclinD1/p27 ratio sets the fate (Fan–Meyer boundary)", fontsize=10, fontweight='bold')
    cb = fig.colorbar(im, ax=axB, fraction=0.045, pad=0.02, ticks=[0, 1, 2])
    cb.ax.set_yticklabels(['immediate', 'transient G0', 'arrest'], fontsize=8)
    axB.text(0.03, 0.95, "low CyclinD1 / high p27 → arrest", transform=axB.transAxes, fontsize=8, color='#7b241c', va='top', fontweight='bold')
    axB.text(0.62, 0.06, "high CyclinD1 / low p27\n→ immediate", transform=axB.transAxes, fontsize=8, color='#145a32', va='bottom', ha='left', fontweight='bold')

    # legend/definition panel (bottom-right)
    axL = fig.add_subplot(gs[1:, 2]); axL.axis('off')
    axL.text(0, 1.0, "DEFINITION (per daughter, at birth):", fontsize=9.5, fontweight='bold', va='top')
    axL.text(0, 0.86, "born-state = CyclinD1/p27 ratio\n  (birth p27 = P21_div; CyclinD1 setpoint)\n\n"
             "G0 dwell = time from birth until p27\n  drops below threshold (Skp2–p27 toggle\n  flips; CDK2 freed). Marked by p27, not pRb.\n\n"
             "• IMMEDIATE  dwell ≤ 2 h\n   (Spencer carryover keeps pRb/CycE/low-p27)\n\n"
             "• TRANSIENT G0  dwells, then commits\n   (grows to M_commit / CyclinD1 rises)\n\n"
             "• ARREST  never commits\n   (CyclinD1 below the threshold)",
             fontsize=8.2, va='top', family='monospace')

    out = os.path.join(os.path.dirname(__file__), "fig_v44_transient_g0_definition")
    fig.savefig(out + ".png", dpi=150, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")
    for k in ['immediate', 'transient_G0', 'arrest']:
        o = R['exemplar'].get(k)
        if o:
            print(f"  exemplar {k}: class={o['cls']} dwell={o['dwell']} ndiv={o['ndiv']}")


if __name__ == '__main__':
    main()
