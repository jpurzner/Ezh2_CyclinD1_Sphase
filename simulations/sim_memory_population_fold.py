"""Does the memory-mark CyclinD1 MB/GNP conflict SURVIVE in a realistic sustained-division POPULATION?

The 5.07 conflict was derived from the DETERMINISTIC settled MB trajectory, where a strong self-reinforcing
mark rebuilds after each division and SATURATES -> loses mitogen-dose sensitivity -> MB/GNP snaps to ~6.75.
But the measured 5.07 is a POPULATION average over asynchronously dividing cells, where replicative dilution
(Mk halved every S) fights the rebuild and the emergent transient-G0 fraction accumulates the mark separately.
JP: test it in the real regime -- the model's single trajectory is not necessarily accurate.

Ensembles of GNP and MB cells (population-of-models heterogeneity: CV in k_Cd_translation, P21_div, kTlEZ ->
emergent transient G0), SUSTAINED proliferation, BASELINE mark vs MEMORY mark. Measures population-average
CyclinD1 (time-avg over asynchronous cyclers = snapshot population mean), the mark split dividing-vs-G0, and
the G0 fraction. Verdict = does population memory MB/GNP stay near 5.07 (conflict = artifact) or ~6.75 (robust)?

Run:  ./venv/bin/python simulations/sim_memory_population_fold.py [--n=40] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'Cd_mRNA', 'Mk', 'aRc', 'Dna', 'P21']
# genotypes (from validate_v44: GNP + SHH vs MB)
GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
MB  = dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002, HHi=0, EZH2i=0)
GENO = {'GNP': GNP, 'MB': MB}
# mark regimes: baseline (memoryless) + a k_w_mk sweep at durable del_mk (EZH2i-fold-matched K_mk).
# Tests whether ANY k_w keeps the DIVIDING population graded (low mark -> 5.07) while G0 accumulates.
REG = {
    'base': dict(k_w_mk=0.00449887, del_mk=0.00328353, K_mk=0.303, f0_mk=0.144970),
    'k008': dict(k_w_mk=0.008,      del_mk=0.001,      K_mk=0.40,  f0_mk=0.144970),
    'k014': dict(k_w_mk=0.014,      del_mk=0.001,      K_mk=0.55,  f0_mk=0.144970),
    'k025': dict(k_w_mk=0.025,      del_mk=0.001,      K_mk=0.72,  f0_mk=0.144970),
}
KW_OF = {'base': 0.0045, 'k008': 0.008, 'k014': 0.014, 'k025': 0.025}
SETTLE, RECORD, CHUNK, NP = 5000.0, 15000.0, 60.0, 9
KTL0, P21_MED, KTLEZ0, MU0 = 0.801, 0.42, 0.004, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, PART, P27_HI = 0.633, 0.32, 0.51, 0.30, 0.1
CACHE = 'simulations/sim_memory_population_fold_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000); _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass


def _work(task):
    i, geno, reg, ktl, p21d, ktlez = task
    rng = np.random.default_rng((hash((i, geno, reg)) & 0xFFFFFFFF))
    _RR.reset()
    for k, v in GENO[geno].items():
        _RR[k] = v
    for k, v in REG[reg].items():
        _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = MU0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try:
            _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception:
            return None
        tm += CHUNK
    CD, MK, ARC, DNA, P21 = [], [], [], [], []
    t = 0.0
    while t < RECORD:
        _RR['P21_div'] = draw(); r = None
        for atol, ms in [(1e-9, 20.0), (1e-8, 20.0), (1e-7, 8.0), (1e-6, 3.0), (1e-5, 1.5)]:
            try:
                _RR.integrator.setValue('absolute_tolerance', atol); _RR.integrator.setValue('maximum_time_step', ms)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception:
                r = None
        if r is None:
            break
        for A, s in [(CD, 'Cd_mRNA'), (MK, 'Mk'), (ARC, 'aRc'), (DNA, 'Dna'), (P21, 'P21')]:
            A.append(r[s][1:])
        t += CHUNK; tm += CHUNK
    if not CD:
        return None
    CD = np.concatenate(CD); MK = np.concatenate(MK); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA); P21 = np.concatenate(P21)
    div = int(np.sum((DNA[:-1] > 0.9) & (DNA[1:] < 0.1)))
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = preS & (P21 > P27_HI)
    return dict(geno=geno, reg=reg, cd_mean=float(np.mean(CD)), mk_mean=float(np.mean(MK)),
                mk_cyc=float(np.mean(MK[~g0])) if np.any(~g0) else np.nan,
                mk_g0=float(np.mean(MK[g0])) if np.any(g0) else np.nan,
                g0_frac=float(np.mean(g0)), ndiv=div,
                mk_hist=np.histogram(MK, bins=np.linspace(0, 1, 26))[0].astype(float))


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='):
            return type(d)(a.split('=', 1)[1])
    return d


if __name__ == '__main__':
    N = _arg('--n', 40); WORKERS = _arg('--workers', max(1, mp.cpu_count() - 1))
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        tasks = [(i, g, rg, ktl[i], p21[i], ktz[i]) for g in GENO for rg in REG for i in range(N)]
        print(f'{len(tasks)} cells ({len(GENO)} geno x {len(REG)} regime x N={N}), sustained {RECORD/60:.0f}h', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
            res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]

        def grp(g, rg):
            return [r for r in res if r['geno'] == g and r['reg'] == rg]
        agg = {}
        for g in GENO:
            for rg in REG:
                sub = grp(g, rg)
                agg[(g, rg)] = dict(
                    cd=float(np.mean([r['cd_mean'] for r in sub])),
                    mk=float(np.nanmean([r['mk_mean'] for r in sub])),
                    mk_cyc=float(np.nanmean([r['mk_cyc'] for r in sub])),
                    mk_g0=float(np.nanmean([r['mk_g0'] for r in sub])),
                    g0=float(np.mean([r['g0_frac'] for r in sub])),
                    ndiv=float(np.mean([r['ndiv'] for r in sub])),
                    hist=np.sum([r['mk_hist'] for r in sub], 0), n=len(sub))
        np.savez(CACHE, **{f'{g}_{rg}_{k}': v for (g, rg), d in agg.items() for k, v in d.items()})
        print('cached ->', CACHE, flush=True)

    z = dict(np.load(CACHE))
    def G(g, rg, k):
        return float(z[f'{g}_{rg}_{k}']) if z[f'{g}_{rg}_{k}'].ndim == 0 else z[f'{g}_{rg}_{k}']
    order = list(REG)
    fold = {rg: G('MB', rg, 'cd') / G('GNP', rg, 'cd') for rg in order}
    print('\n===== POPULATION CyclinD1 MB/GNP + mark saturation (sustained division) =====')
    print(f'  measured target = 5.07 ; raw cascade (feedback off) = 6.86')
    print(f'  {"regime":8s} {"k_w_mk":>7s} {"MB/GNP":>7s} | {"MB mk_cyc":>9s} {"MB mk_G0":>8s} {"GNP mk_cyc":>10s} {"MB G0frac":>9s} {"MBdiv":>6s}')
    for rg in order:
        print(f'  {rg:8s} {KW_OF[rg]:7.4f} {fold[rg]:7.2f} | {G("MB",rg,"mk_cyc"):9.2f} {G("MB",rg,"mk_g0"):8.2f} {G("GNP",rg,"mk_cyc"):10.2f} {G("MB",rg,"g0"):9.2f} {G("MB",rg,"ndiv"):6.0f}')

    kws = [KW_OF[rg] for rg in order]
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5))
    # (A) population fold vs k_w_mk
    a = ax[0]
    a.plot(kws, [fold[rg] for rg in order], '-o', color='#e45756', lw=2, ms=7)
    a.axhline(5.07, color='k', ls='--', lw=1.5); a.text(kws[-1], 5.13, 'measured 5.07', fontsize=9, ha='right')
    a.axhline(6.86, color='gray', ls=':', lw=1.2); a.text(kws[-1], 6.92, 'raw cascade 6.86', fontsize=8, ha='right', color='gray')
    for rg in order:
        a.annotate(f'{fold[rg]:.2f}', (KW_OF[rg], fold[rg]), textcoords='offset points', xytext=(0, 8), ha='center', fontsize=8)
    a.set_xlabel('read-write k_w_mk (mark strength)'); a.set_ylabel('population CyclinD1 MB/GNP'); a.set_ylim(4.5, 7.2)
    a.set_title('(A) Population fold climbs to raw ~7 as the\nmark strengthens -- conflict SURVIVES the population', fontweight='bold', fontsize=10.5); a.grid(alpha=0.15)
    # (B) mark in DIVIDING cells vs k_w_mk -- does dilution keep it graded?
    a = ax[1]
    a.plot(kws, [G('MB', rg, 'mk_cyc') for rg in order], '-o', color='#54a24b', lw=2, ms=6, label='MB dividing/cyc')
    a.plot(kws, [G('MB', rg, 'mk_g0') for rg in order], '-s', color='#b279a2', lw=2, ms=6, label='MB transient-G0')
    a.plot(kws, [G('GNP', rg, 'mk_cyc') for rg in order], '--o', color='#7aa457', lw=1.5, ms=5, label='GNP dividing/cyc')
    a.axhspan(0.6, 1.0, color='red', alpha=0.05); a.text(kws[1], 0.94, 'saturated -> dose-sensitivity lost', fontsize=8, color='#b03030')
    a.set_xlabel('read-write k_w_mk'); a.set_ylabel('mean mark Mk'); a.set_ylim(0, 1)
    a.set_title('(B) The mark SATURATES the dividing population too\n(rebuild beats division-dilution) -- no graded window', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=8.5); a.grid(alpha=0.15)
    # (C) mark distribution (k025 memory) MB vs GNP
    a = ax[2]
    ctr = 0.5 * (np.linspace(0, 1, 26)[:-1] + np.linspace(0, 1, 26)[1:])
    for g, c in [('MB', '#e45756'), ('GNP', '#4c78a8')]:
        for rg, ls, lab in [('base', ':', 'baseline'), ('k025', '-', 'memory')]:
            h = G(g, rg, 'hist'); h = h / h.sum()
            a.plot(ctr, h, color=c, lw=2, ls=ls, label=f'{g} {lab}')
    a.set_xlabel('H3K27me3 mark Mk'); a.set_ylabel('population density')
    a.set_title('(C) Mark distribution: memory mark piles up\nat Mk~1 (saturated), baseline stays low/graded', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=8); a.grid(alpha=0.15)
    fig.suptitle('Memory-mark CyclinD1 MB/GNP conflict in a REALISTIC sustained-division population (GNP vs MB ensembles, N=%d)' % int(_arg('--n', 40)), fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_memory_population_fold.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_memory_population_fold.pdf', bbox_inches='tight')
    print('\nSaved sim_memory_population_fold.png')
