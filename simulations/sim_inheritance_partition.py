"""PROTOTYPE: calibrated distributions + real division-to-division inheritance with PARTITIONING NOISE.

Two upgrades over sim_steadystate_g0_mark (which used placeholder widths and a FIXED per-lineage birth-p27):
  1. CALIBRATED per-lineage abundance distributions (shape/CV only, never absolute means):
       k_Cd_translation = 0.801 · LogNormal(0, 0.633)   (CyclinD1 CV 0.70, data)
       P21_div (birth p27 median) = 0.6 · LogNormal(0, 0.32)   (p27 CV 0.33, data; unimodal input)
       kTlEZ = 0.004 · LogNormal(0, 0.51)               (EZH2 abundance CV ~0.55, data)
     (cell-cycle period mu fixed at default to ISOLATE the calibrated abundance axes.)
  2. PARTITIONING NOISE at division: instead of resetting p27 to a fixed lineage value every cycle, each
     division redraws birth-p27 = P21_div_lineage · exp(N(-s²/2, s))  (mean-preserving, s = PART_SIG).
     This is injected in the Python loop (per-chunk P21_div update) — the calibrated model file is untouched.
     Effect: within a lineage the birth-p27 (the Spencer/Fan-Meyer bifurcation variable) now varies cycle
     to cycle, so a near-threshold cell can FLIP between transient-G0 and immediate re-entry — the genuine
     stochastic-inheritance behaviour the deterministic model cannot show (each lineage otherwise has ONE fate).

     NB CyclinD1 is NOT given partition noise: in v44 Cd is mitogen-clamped (continuous, not halved at
     division), so a division-time perturbation relaxes back within the fate-decision window — partition
     noise on p27 (which persists until the Skp2-p27 feedforward clears it) is what actually moves fate.

Classification (probe_transient_g0 conventions): divisions=MPF peaks; per-cycle dwell = time in (pre-S ∧
p27>0.1); lineage class from mean dwell (<2h immediate / ≥2h transient / <2 divisions arrest).

Run:  ./venv/bin/python simulations/sim_inheritance_partition.py [--n=200] [--part=0.30] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.signal import find_peaks
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SHH_LEVELS = [0.40, 0.55, 0.70, 0.85, 1.00]
SHH_TRACE = 0.55                                   # SHH used for the within-lineage panels
SETTLE, MEAS, CHUNK, NP = 5000.0, 11000.0, 60.0, 9
DWELL_CUT, P27_HI = 2.0, 0.1
KTL0, KTLEZ0, P21_MED = 0.801, 0.004, 0.42   # birth-p27 median = OPERATING POINT (a.u., not data-constrained;
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.633, 0.32, 0.51

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, shh, f0, part, ktl, p21d, ktlez, want_trace = task
    rng = np.random.default_rng((hash((i, round(shh, 3), f0, round(part, 3))) & 0xFFFFFFFF))
    _RR.reset()
    _RR['SHH'] = shh
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['f0_mk'] = f0
    _RR['P21_div'] = float(p21d)
    lognorm = lambda: p21d * np.exp(rng.normal(-part * part / 2.0, part)) if part > 0 else p21d  # mean-preserving
    T, MPF, P21, ARC, DNA = [], [], [], [], []
    t = 0.0; TEND = SETTLE + MEAS
    while t < TEND:
        if part > 0: _RR['P21_div'] = float(lognorm())        # fresh birth-p27 for the NEXT division
        r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(t, t + CHUNK, NP, selections=SEL); break
            except Exception:
                r = None
        if r is None: break
        if t + CHUNK > SETTLE:
            for A, s in [(T, 'time'), (MPF, 'MPF'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna')]:
                A.append(r[s][1:])
        t += CHUNK
    if not T:
        return (i, shh, f0, part, 'arrest', np.nan, np.nan, np.nan, 0, None)
    T = np.concatenate(T); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA)
    dt = np.median(np.diff(T))
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]      # canonical division = Dna reset (robust at coarse dt; MPF peaks are too sharp to sample)
    if len(div) < 2:
        return (i, shh, f0, part, 'arrest', np.nan, np.nan, np.nan, 0, None)
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = preS & (P21 > P27_HI)
    dwells = []                                                # per-cycle G0 dwell (h) between consecutive divisions
    for a, b in zip(div[:-1], div[1:]):
        dwells.append(float(g0[a:b].sum() * dt / 60.0))
    dwells = np.array(dwells)
    dwell_mean = float(np.mean(dwells))
    cls = 'transient' if dwell_mean >= DWELL_CUT else 'immediate'
    dcv = float(np.std(dwells) / (np.mean(dwells) + 1e-9)) if len(dwells) >= 2 else 0.0
    frac_tr = float(np.mean(dwells >= DWELL_CUT))              # mixedness: fraction of this lineage's cycles that are transient-G0
    seq = dwells.astype(np.float32) if want_trace else None
    return (i, shh, f0, part, cls, dwell_mean, dcv, frac_tr, len(dwells), seq)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 200)); PART = float(_arg('--part', 0.30))
    WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_inheritance_partition_cache.npz'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl_i = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21d_i = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktlez_i = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        tasks = []
        for part in (0.0, PART):
            for shh in SHH_LEVELS:
                for f0 in (0.233, 1.0):
                    for i in range(N):
                        wt = (shh == SHH_TRACE)
                        tasks.append((i, shh, f0, part, ktl_i[i], p21d_i[i], ktlez_i[i], wt))
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            results = []
            for k, res in enumerate(pool.imap_unordered(_work, tasks, chunksize=8)):
                results.append(res)
                if (k + 1) % 500 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        cc = {'immediate': 0, 'transient': 1, 'arrest': 2}
        store = {}
        for part in (0.0, PART):
            ptag = 'nn' if part == 0 else 'pn'
            for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
                frac = np.zeros((len(SHH_LEVELS), 3))
                for si, shh in enumerate(SHH_LEVELS):
                    sub = [r for r in results if abs(r[1] - shh) < 1e-6 and r[2] == f0 and r[3] == part]
                    codes = np.array([cc[r[4]] for r in sub])
                    for c in range(3): frac[si, c] = np.mean(codes == c) if len(codes) else np.nan
                store[f'{ptag}_{tag}_frac'] = frac
            # within-lineage stats at SHH_TRACE, WITH mark
            sub = [r for r in results if abs(r[1] - SHH_TRACE) < 1e-6 and r[2] == 0.233 and r[3] == part and r[4] != 'arrest']
            store[f'{ptag}_dcv'] = np.array([r[6] for r in sub])
            store[f'{ptag}_fractr'] = np.array([r[7] for r in sub])
            seqs = [r[9] for r in sub if r[9] is not None and len(r[9]) >= 6]
            # example lineages: pick ones that STRADDLE the 2h line (recur through transient G0), for both nn & pn
            if seqs:
                flippy = [s for s in seqs if s.min() < DWELL_CUT < s.max()]
                pool = flippy if len(flippy) >= 3 else seqs
                means = np.array([np.mean(s) for s in pool])
                idx = np.argsort(np.abs(means - DWELL_CUT))[:3]
                for j, ii in enumerate(idx): store[f'ex_{ptag}_{j}'] = pool[ii]
        np.savez(CACHE, shh=np.array(SHH_LEVELS), N=N, PART=PART, SHH_TRACE=SHH_TRACE, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    shh = z['shh']; Nc = int(z['N']); PART = float(z['PART']); ST = float(z['SHH_TRACE'])
    labels = ['immediate re-entry', 'transient G0', 'arrest']; cols = ['#2e8b57', '#e08a1e', '#8b1a1a']

    fig, ax = plt.subplots(2, 2, figsize=(15, 10.5))

    def stacked(a, ptag, ttl):
        w = 0.020
        for grp, tag, off, hatch in [('WITH', 'W', -w, None), ('WITHOUT', 'N', +w, '//')]:
            frac = z[f'{ptag}_{tag}_frac']; bottom = np.zeros(len(shh))
            for c in range(3):
                a.bar(shh + off, frac[:, c], width=w * 1.8, bottom=bottom, color=cols[c], hatch=hatch,
                      edgecolor='white', lw=0.4, label=(labels[c] if grp == 'WITH' else None))
                bottom += frac[:, c]
        a.set_xlabel('steady-state Hh (SHH)  (left=WITH mark, right hatched=WITHOUT)'); a.set_ylabel('fraction of cells')
        a.set_title(ttl, fontweight='bold', fontsize=11.5); a.set_ylim(0, 1); a.legend(fontsize=8.5, loc='lower right')

    stacked(ax[0, 0], 'nn', '(A) Calibrated distributions, NO partition noise\n(fixed birth-p27 per lineage — each cell has ONE fate)')
    stacked(ax[0, 1], 'pn', f'(B) Calibrated + PARTITION NOISE (σ={PART}) at each division\n(birth-p27 redrawn per cycle)')

    def seqpanel(a, ptag, ttl):
        ylo, yhi = 0.02, 60.0
        a.axhspan(DWELL_CUT, yhi, color='#e08a1e', alpha=0.10); a.axhspan(ylo, DWELL_CUT, color='#2e8b57', alpha=0.10)
        a.axhline(DWELL_CUT, color='k', ls='--', lw=1.2)
        any_ = False
        for j, c in enumerate(['#1f77b4', '#9467bd', '#2ca02c']):
            key = f'ex_{ptag}_{j}'
            if key in z.files:
                s = np.clip(z[key], ylo, yhi); a.plot(np.arange(1, len(s) + 1), s, '-o', color=c, lw=1.8, ms=5, label=f'lineage {j+1}'); any_ = True
        a.set_yscale('log'); a.set_ylim(ylo, yhi)
        a.text(0.98, 0.965, 'transient G0  (dwell ≥ 2 h)', transform=a.transAxes, ha='right', va='top', fontsize=8.5, color='#7a5210')
        a.text(0.98, 0.03, 'immediate re-entry (< 2 h)', transform=a.transAxes, ha='right', va='bottom', fontsize=8.5, color='#1f5c3a')
        a.set_xlabel('division number (within one lineage)'); a.set_ylabel("that cycle's G0 dwell (h, log)")
        a.set_title(ttl, fontweight='bold', fontsize=11)
        if any_: a.legend(fontsize=8.5, loc='center right')

    # (C) deterministic near-threshold lineages — REGULAR multi-periodic recurrence through transient G0
    seqpanel(ax[1, 0], 'nn', f'(C) NO noise: near-threshold lineages already recur through G0\non a REGULAR (deterministic multi-periodic) schedule (SHH={ST}, WITH mark)')
    # (D) same cells + partition noise — the recurrence becomes IRREGULAR / stochastic
    seqpanel(ax[1, 1], 'pn', f'(D) + PARTITION NOISE (σ={PART}): the recurrence becomes\nIRREGULAR / stochastic — sister divisions diverge (SHH={ST}, WITH mark)')

    fig.suptitle(f'Calibrated distributions + division-inheritance partition noise (N={Nc}/level/condition) — prototype',
                 fontsize=13.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_inheritance_partition.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_inheritance_partition.pdf', bbox_inches='tight')
    plt.close()
    for ptag, nm in [('nn', 'NO noise'), ('pn', f'partition σ={PART}')]:
        print(f'{nm}:  transient-G0 frac (WITH mark) by SHH = {z[f"{ptag}_W_frac"][:,1].round(2)}   arrest = {z[f"{ptag}_W_frac"][:,2].round(2)}')
    print(f'within-lineage dwell CV (SHH={ST}, WITH):  no-noise median {np.median(z["nn_dcv"]):.2f}  vs  noise median {np.median(z["pn_dcv"]):.2f}')
    print('Saved sim_inheritance_partition.png')
