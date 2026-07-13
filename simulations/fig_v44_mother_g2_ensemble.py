"""Figure: mother-G2 p27 integrator -> the Overton graded->binary fate bifurcation (ensemble, v44).

Exploratory. Builds on the TWO-STEP Rb base (CyclinD1/p27-ratio commitment), where the inherited birth
p27 is the fate-setter -- the mother-G2 integrator is directionally correct there (low mother-G2 mitogen
-> arrest/quiescent; high -> immediate), unlike on the single-step model where the daughter's own CyclinD1
dominates. An ENSEMBLE of cells heterogeneous in mitogen response (k_Cd_translation, lognormal) is swept
across mitogen (SHH); at each level the population splits into cycling (CDK2inc/immediate) vs quiescent
(CDK2low/transient-G0 + arrest). Spencer/Overton: physiological mitogen sits INSIDE the bistable window.

Run: ./venv/bin/python simulations/fig_v44_mother_g2_ensemble.py [--workers N]
"""
import os, sys, multiprocessing as mp
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

KTL0 = 0.801
MG2 = dict(k_sg2=0.008, g_moth=3.0, K_mit_g2=1.5)   # mother-G2 calibration (right direction on two-step)
N_CELLS = 60
SHH_SWEEP = np.array([0.25, 0.35, 0.45, 0.6, 0.8, 1.0, 1.3])
T_END, N_PTS, SETTLE, DWELL_CUT = 14000, 28000, 5000, 2.0

_RR_ON = _RR_OFF = None
def _init(_):
    global _RR_ON, _RR_OFF
    import tellurium as te
    _RR_ON = te.loada(build_model_v44(with_two_step_rb=True, with_mother_g2=True))
    _RR_OFF = te.loada(build_model_v44(with_two_step_rb=True, with_mother_g2=False))
    for r in (_RR_ON, _RR_OFF):
        r.integrator.setValue('absolute_tolerance', 1e-8); r.integrator.setValue('relative_tolerance', 1e-6)
        try:
            r.integrator.setValue('maximum_num_steps', 2000000); r.integrator.setValue('maximum_time_step', 20.0)
        except Exception:
            pass

def classify_cell(r, shh, ktl, mother):
    r.reset()
    r['SHH'] = shh; r['k_Cd_translation'] = ktl
    if mother:
        for k, v in MG2.items():
            r[k] = v
    for maxstep in (20.0, 5.0, 2.0):
        try:
            r.integrator.setValue('maximum_time_step', maxstep)
            d = r.simulate(0, T_END, N_PTS, selections=['time', 'P21', 'Dna']); break
        except Exception:
            d = None
    if d is None:
        return 'fail', np.nan
    t = d['time']; dt = t[1] - t[0]; Dna = d['Dna']
    div = np.where((Dna[:-1] > 0.5) & (Dna[1:] < 0.1))[0] + 1
    div = div[t[div] >= SETTLE]
    if len(div) < 2:
        return 'arrest', np.nan
    dws = []
    for i, b in enumerate(div[:-1]):
        below = np.where(d['P21'][b:div[i + 1]] < 0.1)[0]
        dws.append((below[0] * dt / 60.0) if len(below) else (div[i + 1] - b) * dt / 60.0)
    dw = float(np.mean(dws))
    return ('immediate' if dw < DWELL_CUT else 'transient_G0'), dw

def work(task):
    ish, icell, shh, ktl, mother = task
    r = _RR_ON if mother else _RR_OFF
    cls, dw = classify_cell(r, shh, ktl, mother)
    return (ish, icell, mother, cls, dw)


def main():
    workers = 8
    for a in sys.argv:
        if a.startswith('--workers='):
            workers = int(a.split('=')[1])
    rng = np.random.default_rng(20260713)
    KTL = KTL0 * np.exp(rng.normal(0, 0.4, N_CELLS))   # heterogeneous mitogen response
    tasks = []
    for ish, shh in enumerate(SHH_SWEEP):
        for icell, ktl in enumerate(KTL):
            for mother in (True, False):
                tasks.append((ish, icell, float(shh), float(ktl), mother))
    print(f'running {len(tasks)} single-cell sims on {workers} workers ...', flush=True)
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(work, tasks, chunksize=6)
    # tally
    fates = {True: {i: [] for i in range(len(SHH_SWEEP))}, False: {i: [] for i in range(len(SHH_SWEEP))}}
    dwell = {True: {i: [] for i in range(len(SHH_SWEEP))}, False: {i: [] for i in range(len(SHH_SWEEP))}}
    for ish, icell, mother, cls, dw in res:
        fates[mother][ish].append(cls)
        if not np.isnan(dw):
            dwell[mother][ish].append(dw)
    plot(fates, dwell, KTL)


def _frac(cls_list, kinds):
    n = len(cls_list) or 1
    return sum(1 for c in cls_list if c in kinds) / n

def plot(fates, dwell, KTL):
    fig = plt.figure(figsize=(16, 5.5))
    gs = GridSpec(1, 3, figure=fig, wspace=0.3, left=0.06, right=0.975, top=0.82, bottom=0.14)
    fig.suptitle("Mother-G2 p27 integrator → graded fate bifurcation across mitogen (ensemble, two-step base)",
                 fontsize=13.5, fontweight="bold", y=0.97)
    fig.text(0.5, 0.885, f"N={N_CELLS} cells heterogeneous in mitogen response (k_Cd_translation, lognormal). "
             "Quiescent = arrest + transient-G0 (CDK2low); cycling = immediate (CDK2inc). Exploratory / first calibration.",
             ha="center", fontsize=9, color="#555", style="italic")

    # A: quiescent fraction vs mitogen, mother-G2 ON vs OFF
    axA = fig.add_subplot(gs[0, 0])
    for mother, c, lab in [(True, '#c0392b', 'mother-G2 ON'), (False, '#7f8c8d', 'OFF (baseline)')]:
        q = [_frac(fates[mother][i], ('arrest', 'transient_G0')) for i in range(len(SHH_SWEEP))]
        axA.plot(SHH_SWEEP, q, '-o', color=c, lw=2, ms=5, label=lab)
    axA.set_xlabel("SHH (mitogen dose)"); axA.set_ylabel("quiescent fraction (CDK2low)")
    axA.set_ylim(-0.03, 1.03); axA.grid(alpha=0.25); axA.legend(fontsize=8.5)
    axA.set_title("A. Graded bifurcation: quiescent fraction ↓ with mitogen", fontsize=10, fontweight='bold')
    axA.axvspan(0.4, 0.7, color='#f9e79f', alpha=0.4); axA.text(0.55, 0.05, 'physiological\n(bistable) window', fontsize=7, ha='center', color='#7d6608')

    # B: stacked fate fractions vs mitogen (mother-G2 ON)
    axB = fig.add_subplot(gs[0, 1])
    imm = [_frac(fates[True][i], ('immediate',)) for i in range(len(SHH_SWEEP))]
    tg0 = [_frac(fates[True][i], ('transient_G0',)) for i in range(len(SHH_SWEEP))]
    arr = [_frac(fates[True][i], ('arrest',)) for i in range(len(SHH_SWEEP))]
    axB.stackplot(SHH_SWEEP, imm, tg0, arr, labels=['immediate (CDK2inc)', 'transient G0', 'arrest'],
                  colors=['#7dcea0', '#f5b041', '#e74c3c'], alpha=0.9)
    axB.set_xlabel("SHH (mitogen dose)"); axB.set_ylabel("population fraction"); axB.set_ylim(0, 1)
    axB.set_xlim(SHH_SWEEP[0], SHH_SWEEP[-1]); axB.legend(fontsize=7.5, loc='center right')
    axB.set_title("B. Fate composition (mother-G2 ON)", fontsize=10, fontweight='bold')

    # C: dwell distribution at an intermediate mitogen (the split)
    axC = fig.add_subplot(gs[0, 2])
    imid = int(np.argmin(np.abs(SHH_SWEEP - 0.5)))
    dw = dwell[True][imid]
    if dw:
        axC.hist(dw, bins=np.linspace(0, max(20, max(dw) + 1), 22), color='#8e44ad', alpha=0.8, edgecolor='w')
    axC.axvline(DWELL_CUT, color='k', ls='--', lw=1); axC.text(DWELL_CUT + 0.3, axC.get_ylim()[1] * 0.9, 'immediate | transient-G0', fontsize=7)
    axC.set_xlabel("per-daughter G0 dwell (h)"); axC.set_ylabel("cells")
    axC.set_title(f"C. Dwell distribution at SHH={SHH_SWEEP[imid]:.2f} (the split)", fontsize=10, fontweight='bold')

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mother_g2_ensemble")
    fig.savefig(out + ".png", dpi=150, bbox_inches="tight"); fig.savefig(out + ".pdf", bbox_inches="tight")
    print("wrote", out + ".png /.pdf")
    for i, shh in enumerate(SHH_SWEEP):
        q = _frac(fates[True][i], ('arrest', 'transient_G0'))
        print(f"  SHH {shh:.2f}: quiescent {q:.2f} (ON) vs {_frac(fates[False][i],('arrest','transient_G0')):.2f} (OFF)")


if __name__ == '__main__':
    main()
