"""Whole GNP developmental division period — a ramp-UP → sustain → ramp-DOWN Hh trajectory giving a
target MEAN of ~8 divisions, comparing WITH vs WITHOUT the H3K27me3 mark at MATCHED PROLIFERATIVE OUTPUT.

Rationale (JP): comparing WITH/WITHOUT at the same Hh is confounded — the mark represses CyclinD1, so a
marked cell needs a HIGHER Hh to divide at all. The fair comparison matches OUTPUT (division count): give
each condition the Hh plateau that yields the target mean over the developmental window, then compare how
they get there (entry/exit timing, transient G0, CyclinD1).

Two levers, both matter (2-D Hh × window sweep, JP): the Hh plateau LEVEL sets the mean/entry-floor, and
the plateau DURATION (window) caps the fast-cycler TAIL — the 25-division tail of the earlier 140h window
was the high-mu/low-p27 cells cycling ~10h for the full window. A 100h window caps the tail at p95~14 /
max~18 while the doses (WITH 0.75, WITHOUT 0.22) hold both means at ~8. Short 24h tail so the count is the
Hh-ON developmental period, not post-withdrawal coasting (unmarked cells coast/divide at Hh=0.05; marked
cells arrest — the panel-C exit asymmetry). NB the mark leaves an irreducible ~13% never-divide fraction
at ANY Hh (the repressed low-CyclinD1 tail); the unmarked condition can be pushed to frac-zero 0.

Trajectory: settle arrested at SHH=0.05 → ramp up to plateau over T_UP → plateau T_PLAT → ramp down to
0.05 over T_DOWN → short tail. Ensemble (calibrated abundance + mu + partition). Division = Dna 1→0 reset;
per-time state = cycling / transient-G0 / arrest.

--calib: sweep the plateau level for each condition and report MEAN (median)[IQR] divisions (pick PLAT_W/N).
main:    run at PLAT_W (mark on) and PLAT_N (mark off) and make the developmental figure.

Run:  ./venv/bin/python simulations/sim_gnp_developmental.py --calib   (find the two Hh levels)
      ./venv/bin/python simulations/sim_gnp_developmental.py [--n=300] [--fresh]   (main figure)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'SHH', 'Dna', 'P21', 'aRc', 'Cd', 'EZH2']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH_FLOOR = 0.05
T_UP, T_PLAT, T_DOWN, TAIL = 48 * 60.0, 100 * 60.0, 48 * 60.0, 24 * 60.0  # short tail: count the Hh-ON developmental period, not post-withdrawal coasting
SETTLE, CHUNK, NP = 3000.0, 60.0, 9
DT = CHUNK / (NP - 1)
DWELL_CUT, P27_HI, PART = 2.0, 0.1, 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, MU_SDLOG = 0.633, 0.32, 0.51, 0.22
PLAT_W, PLAT_N = 0.75, 0.22                              # matched mean~8 Hh plateaus on the 100h window (2-D Hh×window sweep; tail capped p95~14/max~18)
T_TRACE = np.arange(0.0, T_UP + T_PLAT + T_DOWN + TAIL + 1e-6, 12.0)

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 1000000)
    except Exception: pass
    try: _RR.integrator.setValue('maximum_time_step', 20.0)   # faster kDeEZ (2026-07 recal) -> stiff ramp-down arrests
    except Exception: pass


def shh_of(t, plat):
    if t < T_UP: return HH_FLOOR + (plat - HH_FLOOR) * t / T_UP
    if t < T_UP + T_PLAT: return plat
    td = t - T_UP - T_PLAT
    if td < T_DOWN: return plat + (HH_FLOOR - plat) * td / T_DOWN
    return HH_FLOOR


def _work(task):
    i, plat, f0, ktl, p21d, ktlez, mu, wt = task
    rng = np.random.default_rng((hash((i, round(plat, 3), f0)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = HH_FLOOR
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu); _RR['f0_prc2'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    NULL = (plat, f0, 0, None, None, None, None)
    tm = 0.0
    while tm < SETTLE:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return NULL
        tm += CHUNK
    T, SH, DNA, P21, ARC, CD, EZ = [], [], [], [], [], [], []
    t = 0.0; TEND = T_UP + T_PLAT + T_DOWN + TAIL
    while t < TEND:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t, plat)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        T.append(r['time'][1:] - tm + t); SH.append(r['SHH'][1:]); DNA.append(r['Dna'][1:]); P21.append(r['P21'][1:]); ARC.append(r['aRc'][1:]); CD.append(r['Cd'][1:]); EZ.append(r['EZH2'][1:])
        t += CHUNK; tm += CHUNK
    if not T:
        return NULL
    T = np.concatenate(T); DNA = np.concatenate(DNA); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); CD = np.concatenate(CD); EZ = np.concatenate(EZ)
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    ndiv = len(div)
    if not wt:
        return (plat, f0, ndiv, None, None, None, None)
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; g0inst = (~inS & ~inG2) & (P21 > 0.1)
    state = np.zeros(len(T), np.int8)
    dd = np.diff(np.concatenate([[0], g0inst.astype(np.int8), [0]]))
    for s, e in zip(np.where(dd == 1)[0], np.where(dd == -1)[0]):
        ee = min(e, len(T) - 1)
        if (ee - s) * DT < DWELL_CUT * 60: continue
        if np.any(div > ee): state[s:e] = 1
        else: state[s:] = 2; break
    st = np.clip(np.round(np.interp(T_TRACE, T, state)), 0, 2).astype(np.int8)
    cd = np.interp(T_TRACE, T, CD).astype(np.float32)
    ez = np.interp(T_TRACE, T, EZ).astype(np.float32)
    return (plat, f0, ndiv, st, cd, T[div].astype(np.float32), ez)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv; CALIB = '--calib' in sys.argv
    N = int(_arg('--n', 40 if CALIB else 300)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    rng = np.random.default_rng(11)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)

    if CALIB:
        PLATS = [0.18, 0.22, 0.26, 0.30, 0.36, 0.45, 0.55, 0.70, 0.85]
        tasks = [(i, pl, f0, ktl[i], p21[i], ktz[i], mu[i], False) for f0 in (0.05054636367014487, 1.0) for pl in PLATS for i in range(N)]
        print(f'CALIB: {len(tasks)} cells on {WORKERS} workers ... (target MEAN divisions = 8)', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=6))
        print('plateau   WITHOUT mean (median)[IQR]        WITH mean (median)[IQR]')
        for pl in PLATS:
            row = []
            for f0 in (1.0, 0.05054636367014487):
                nd = np.array([r[2] for r in res if abs(r[0] - pl) < 1e-6 and r[1] == f0])
                row.append(f'{nd.mean():.1f} ({np.median(nd):.0f})[{np.percentile(nd,25):.0f}-{np.percentile(nd,75):.0f}]')
            print(f'  {pl:.2f}      {row[0]:26s}   {row[1]}')
        sys.exit(0)

    CACHE = 'simulations/sim_gnp_developmental_cache.npz'
    if FRESH or not os.path.exists(CACHE):
        tasks = ([(i, PLAT_W, 0.05054636367014487, ktl[i], p21[i], ktz[i], mu[i], True) for i in range(N)] +
                 [(i, PLAT_N, 1.0, ktl[i], p21[i], ktz[i], mu[i], True) for i in range(N)])
        print(f'MAIN: {len(tasks)} cells on {WORKERS} workers (PLAT_W={PLAT_W}, PLAT_N={PLAT_N}) ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=6))
        store = {}
        for f0, tag in [(0.05054636367014487, 'W'), (1.0, 'N')]:
            sub = [r for r in res if r[1] == f0 and r[3] is not None]
            nd = np.array([r[2] for r in sub])
            sts = np.array([r[3] for r in sub]); cds = np.array([r[4] for r in sub]); ezs = np.array([r[6] for r in sub])
            store[f'{tag}_ndiv'] = nd
            store[f'{tag}_cyc'] = (sts == 0).mean(0); store[f'{tag}_tg0'] = (sts == 1).mean(0); store[f'{tag}_arr'] = (sts == 2).mean(0)
            store[f'{tag}_cd'] = np.nanmean(cds, 0); store[f'{tag}_ez'] = np.nanmean(ezs, 0)
            ex = np.argsort(nd)[len(nd) // 2]                # median-count example cell trace
            store[f'{tag}_ex_cd'] = cds[ex]; store[f'{tag}_ex_dt'] = np.array([r[5] for r in sub], dtype=object)[ex]
        shh_W = np.array([shh_of(tt, PLAT_W) for tt in T_TRACE]); shh_N = np.array([shh_of(tt, PLAT_N) for tt in T_TRACE])  # T_TRACE in minutes
        np.savez(CACHE, Ttr=T_TRACE / 60.0, N=N, PLAT_W=PLAT_W, PLAT_N=PLAT_N, shh_W=shh_W, shh_N=shh_N,
                 W_ex_dt=store.pop('W_ex_dt'), N_ex_dt=store.pop('N_ex_dt'), **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE, allow_pickle=True)
    T = z['Ttr']; Nc = int(z['N']); PW = float(z['PLAT_W']); PN = float(z['PLAT_N'])

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    # (A) the Hh trajectories (matched output requires WITH mark higher) + regions
    a = ax[0, 0]
    a.plot(T, [shh_of(t * 60, PW) for t in T], color='#8b1a1a', lw=2.6, label=f'WITH mark — Hh plateau {PW:.2f}')
    a.plot(T, [shh_of(t * 60, PN) for t in T], color='#3b7bbf', lw=2.6, label=f'WITHOUT mark — Hh plateau {PN:.2f}')
    a.set_xlabel('developmental time (h)'); a.set_ylabel('Hh (SHH)')
    a.set_title('(A) Developmental Hh trajectory — ramp up → sustain → ramp down\nMATCHED OUTPUT needs a HIGHER plateau WITH the mark', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15)
    # (B) division-count distributions (should both be ~8-12)
    a = ax[0, 1]
    bins = np.arange(-0.5, max(z['W_ndiv'].max(), z['N_ndiv'].max()) + 1.5)
    a.hist(z['N_ndiv'], bins=bins, alpha=0.6, color='#3b7bbf', label=f'WITHOUT (mean {z["N_ndiv"].mean():.1f}, p95 {np.percentile(z["N_ndiv"],95):.0f})')
    a.hist(z['W_ndiv'], bins=bins, alpha=0.6, color='#8b1a1a', label=f'WITH (mean {z["W_ndiv"].mean():.1f}, p95 {np.percentile(z["W_ndiv"],95):.0f})')
    a.set_xlabel('divisions over the developmental period'); a.set_ylabel('cells')
    a.set_title('(B) Division-count distribution — MATCHED at mean~8, tail CAPPED by the 100h window', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15)
    # (C) population dynamics: fraction cycling / transient-G0 / arrest over the period
    a = ax[1, 0]
    for tag, col, nm in [('W', '#8b1a1a', 'WITH'), ('N', '#3b7bbf', 'WITHOUT')]:
        a.plot(T, z[f'{tag}_cyc'], '-', color=col, lw=2.5, label=f'{nm} — cycling')
        a.plot(T, z[f'{tag}_tg0'], '--', color=col, lw=1.7, alpha=0.8, label=f'{nm} — transient-G0')
    a.axvspan(0, T_UP / 60, color='#2e8b57', alpha=0.05); a.axvspan((T_UP + T_PLAT) / 60, (T_UP + T_PLAT + T_DOWN) / 60, color='#c0392b', alpha=0.05)
    a.set_xlabel('developmental time (h)'); a.set_ylabel('fraction of cells')
    a.set_title('(C) Population dynamics: cycling (solid) & transient-G0 (dashed)\nentry on ramp-up, exit on ramp-down', fontweight='bold', fontsize=11); a.legend(fontsize=8, loc='center right'); a.grid(alpha=0.15)
    # (D) mean CyclinD1 over the period + a median example cell
    a = ax[1, 1]
    for tag, col, nm in [('W', '#8b1a1a', 'WITH'), ('N', '#3b7bbf', 'WITHOUT')]:
        a.plot(T, z[f'{tag}_cd'], '-', color=col, lw=2.4, label=f'{nm} — mean CyclinD1')
    a.set_xlabel('developmental time (h)'); a.set_ylabel('CyclinD1 (population mean)')
    a.set_title(f'(D) CyclinD1 over the period — WITH mark actually HIGHER\n(~{PW/PN:.1f}× Hh over-compensates the repression; same divisions at higher CyclinD1)', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15)

    fig.suptitle(f'Whole GNP developmental period — WITH vs WITHOUT H3K27me3 at MATCHED output '
                 f'(mean {z["W_ndiv"].mean():.1f} vs {z["N_ndiv"].mean():.1f} div; Hh {PW:.2f} vs {PN:.2f}; 100h window caps the tail), N={Nc}/cond',
                 fontsize=12, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_gnp_developmental.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_gnp_developmental.pdf', bbox_inches='tight')
    plt.close()

    # ---- companion figure: EZH2 abundance over the developmental period ----
    fig2, ax2 = plt.subplots(1, 2, figsize=(15, 5.5))
    a = ax2[0]
    for tag, col, nm, pl in [('W', '#8b1a1a', 'WITH', PW), ('N', '#3b7bbf', 'WITHOUT', PN)]:
        a.plot(T, z[f'{tag}_ez'], '-', color=col, lw=2.6, label=f'{nm} — mean EZH2 (Hh {pl:.2f})')
    a.set_xlabel('developmental time (h)'); a.set_ylabel('EZH2 (population mean)')
    a.set_title('(A) EZH2 over the developmental period — tracks the cycle/mitogen\nAt MATCHED output the arms CONVERGE (steady EZH2 ≈1.0 both, despite 3.4× Hh)', fontweight='bold', fontsize=11)
    a.legend(fontsize=9, loc='upper right'); a.grid(alpha=0.15)
    a2 = a.twinx()
    a2.plot(T, [shh_of(t * 60, PW) for t in T], ':', color='#8b1a1a', lw=1.3, alpha=0.5)
    a2.plot(T, [shh_of(t * 60, PN) for t in T], ':', color='#3b7bbf', lw=1.3, alpha=0.5)
    a2.set_ylabel('Hh (SHH) — dotted', color='grey'); a2.set_ylim(0, 1.0)
    # (B) EZH2 (solid) vs CyclinD1 (dashed) — the negative feedback in the same cells
    a = ax2[1]
    for tag, col, nm in [('W', '#8b1a1a', 'WITH'), ('N', '#3b7bbf', 'WITHOUT')]:
        a.plot(T, z[f'{tag}_ez'], '-', color=col, lw=2.4, label=f'{nm} — EZH2')
    a.set_xlabel('developmental time (h)'); a.set_ylabel('EZH2 (solid)')
    ab = a.twinx()
    for tag, col in [('W', '#8b1a1a'), ('N', '#3b7bbf')]:
        ab.plot(T, z[f'{tag}_cd'], '--', color=col, lw=1.8, alpha=0.8)
    ab.set_ylabel('CyclinD1 (dashed)')
    a.set_title('(B) EZH2 (solid) vs CyclinD1 (dashed) in the same cells\nthe EZH2 ⊣ CyclinD1 feedback; within a ramp both rise (cycle-driven)', fontweight='bold', fontsize=11)
    a.legend(fontsize=9, loc='upper right'); a.grid(alpha=0.15)
    fig2.suptitle(f'EZH2 abundance over the GNP developmental period — WITH vs WITHOUT H3K27me3 at matched output '
                  f'(Hh {PW:.2f} vs {PN:.2f}), N={Nc}/cond', fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_gnp_developmental_ezh2.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_gnp_developmental_ezh2.pdf', bbox_inches='tight')
    plt.close()

    for tag, nm in [('W', 'WITH   '), ('N', 'WITHOUT')]:
        nd = z[f'{tag}_ndiv']
        print(f'{nm} mean {nd.mean():.1f}  median {np.median(nd):.0f}  IQR {np.percentile(nd,25):.0f}-{np.percentile(nd,75):.0f}  p95 {np.percentile(nd,95):.0f}  max {nd.max():.0f}  frac0 {np.mean(nd==0):.2f}')
    print('Saved sim_gnp_developmental.png')
