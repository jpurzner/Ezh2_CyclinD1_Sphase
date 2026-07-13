"""Population layer: GNP vs MB transient-G0 fractions + the Overton graded->binary split + Spencer sister
concordance, from the DATA-GROUNDED distribution machinery (not a placeholder, not the shelved mitogen tracker).

Reuses the sim_inheritance_partition engine (per-lineage sim, Dna-reset division detection, per-cycle G0 dwell)
but runs it for BOTH cell types with the calibrated per-lineage draws (scale-free, medians only):
    k_Cd_translation = KTL0 * LogN(0, 0.633)   (CyclinD1 abundance CV 0.70, data)
    kTlEZ            = KTLEZ0 * LogN(0, 0.51)   (EZH2 abundance CV ~0.55, data)
    P21_div (birth p27) = MED[cond] * LogN(0, 0.32)  (CV 0.33, data), redrawn each division = partition noise
GNP birth-p27 median 0.6, MB 1.8 (P21_DIV_MB, the fold-safe G0 restore). MB adds its INK4/mitogen identity.

Three readouts:
  (A) population transient-G0 fraction, GNP vs MB (frac of lineage-cycles with G0 dwell >= DWELL_CUT).
  (B) Overton: G0 fraction vs mitogen (SHH) -- graded across the population, binarized per-cell by the Skp2 toggle.
  (C) Spencer: sister concordance -- two daughters of one mother share the lineage draws + inherited state, differ
      only by an independent partition redraw; fraction that land in the SAME fate.

Run:  ./venv/bin/python simulations/sim_population_g0.py [--n=120] [--workers=9] [--pilot] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import multiprocessing as mp
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']
KTL0, KTLEZ0 = 0.801, 0.004
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.633, 0.32, 0.51     # CyclinD1 CV0.70, p27 CV0.33, EZH2 CV0.55 (data)
PART = 0.30                                           # partition-noise sigma on birth-p27
SETTLE, MEAS, CHUNK, NP = 5000.0, 12000.0, 60.0, 9
DWELL_CUT, P27_HI = 2.0, 0.1                          # transient if mean G0 dwell >= 2h; G0 = pre-S & p27>0.1
# cell-type identities (per-condition params) + birth-p27 median
GNP = dict(params=dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002), med=0.6)
MB  = dict(params=dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002), med=1.8)

_RR = None
def _init(_):
    global _RR
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000); _RR.integrator.setValue('maximum_time_step', 5.0)
    except Exception: pass


def _run_lineage(shh, params, ktl, ktlez, med, part, seed):
    """One lineage: return per-cycle G0 dwells (hours). Empty -> arrest."""
    rng = np.random.default_rng(seed)
    _RR.reset()
    _RR['SHH'] = shh
    for k, v in params.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez)
    _RR['P21_div'] = float(med)
    draw = lambda: med * np.exp(rng.normal(-part * part / 2.0, part)) if part > 0 else med
    T, P21, ARC, DNA = [], [], [], []
    t = 0.0; TEND = SETTLE + MEAS
    while t < TEND:
        if part > 0: _RR['P21_div'] = float(draw())         # partition redraw for the NEXT division
        r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(t, t + CHUNK, NP, selections=SEL); break
            except Exception:
                r = None
        if r is None: break
        if t + CHUNK > SETTLE:
            for A, s in [(T, 'time'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna')]:
                A.append(r[s][1:])
        t += CHUNK
    if not T: return np.array([])
    T = np.concatenate(T); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA)
    dt = np.median(np.diff(T))
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]   # canonical division = Dna reset
    if len(div) < 2: return np.array([])
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = preS & (P21 > P27_HI)
    return np.array([float(g0[a:b].sum() * dt / 60.0) for a, b in zip(div[:-1], div[1:])])


def _work(task):
    i, shh, tag, params, ktl, ktlez, med, rep = task
    dwells = _run_lineage(shh, params, ktl, ktlez, med, PART, (hash((i, round(shh, 3), tag, rep)) & 0xFFFFFFFF))
    if len(dwells) == 0:
        return (i, shh, tag, rep, 'arrest', 0.0, 0)
    cls = 'transient' if np.mean(dwells) >= DWELL_CUT else 'immediate'
    frac_tr = float(np.mean(dwells >= DWELL_CUT))          # this lineage's cycle-level transient fraction
    return (i, shh, tag, rep, cls, frac_tr, len(dwells))


def main():
    def arg(f, d):
        for a in sys.argv:
            if a.startswith(f + '='): return a.split('=', 1)[1]
        return d
    N = int(arg('--n', 120)); WORKERS = int(arg('--workers', max(1, mp.cpu_count() - 1)))
    if '--pilot' in sys.argv: N = 24
    SHH_SWEEP = [0.2, 0.3, 0.4, 0.5, 0.7]
    rng = np.random.default_rng(11)
    ktl_i = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)     # shared per-lineage draws (paired
    ktlez_i = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)  # GNP/MB use the SAME abundance cell)
    tasks = []
    for tag, cond in (('GNP', GNP), ('MB', MB)):
        for shh in SHH_SWEEP:                              # Overton: mitogen sweep
            for i in range(N):
                tasks.append((i, shh, tag, cond['params'], float(ktl_i[i]), float(ktlez_i[i]), cond['med'], 0))
        for i in range(N):                                 # Spencer: rep=1 at SHH=0.5 (paired sister -- same draws, diff noise)
            tasks.append((i, 0.5, tag, cond['params'], float(ktl_i[i]), float(ktlez_i[i]), cond['med'], 1))
    print(f"population G0: N={N}/type x {len(SHH_SWEEP)} SHH + sisters on {WORKERS} workers ({len(tasks)} runs) ...", flush=True)
    with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(_work, tasks, chunksize=4)

    def cell(tag, shh, rep):   # -> list of rows for a condition
        return [r for r in res if r[2] == tag and abs(r[1] - shh) < 1e-6 and r[3] == rep]
    frac = lambda rows: np.mean([r[5] for r in rows]) if rows else np.nan          # cycle-level transient fraction
    cdk2low = lambda rows: np.mean([r[4] in ('transient', 'arrest') for r in rows]) if rows else np.nan  # per-lineage CDK2low

    print("\n(A) POPULATION TRANSIENT-G0 (SHH=0.5):")
    for tag in ('GNP', 'MB'):
        rows = cell(tag, 0.5, 0)
        comp = {c: np.mean([r[4] == c for r in rows]) for c in ('immediate', 'transient', 'arrest')}
        print(f"  {tag}: transient-G0 cycle-level {frac(rows)*100:.0f}% | CDK2low(lineage) {cdk2low(rows)*100:.0f}% | "
              f"immediate {comp['immediate']*100:.0f}% / transient {comp['transient']*100:.0f}% / arrest {comp['arrest']*100:.0f}%")

    print("\n(B) OVERTON graded split -- CDK2low(lineage) fraction vs mitogen (SHH):")
    print(f"  {'SHH':>5} | " + ' '.join(f'{s:>5}' for s in SHH_SWEEP))
    for tag in ('GNP', 'MB'):
        print(f"  {tag:>5} | " + ' '.join(f'{cdk2low(cell(tag,s,0))*100:4.0f}%' for s in SHH_SWEEP))

    print("\n(C) SPENCER sister concordance (same lineage draws + inherited identity, independent partition noise):")
    for tag in ('GNP', 'MB'):
        a = {r[0]: r[4] for r in cell(tag, 0.5, 0)}; b = {r[0]: r[4] for r in cell(tag, 0.5, 1)}
        both = [i for i in a if i in b]
        conc = np.mean([a[i] == b[i] for i in both]) if both else np.nan
        print(f"  {tag}: {conc*100:.0f}% same-fate sisters (n={len(both)} pairs; target ~98%)")


if __name__ == '__main__':
    main()
