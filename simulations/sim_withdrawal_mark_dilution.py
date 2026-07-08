"""Hh WITHDRAWAL: the H3K27me3 dilution-vs-restoration race + is transient G0 more/less WITH the mark?
Adds a CELL-CYCLE-DURATION (mu) distribution — the pivot of the dilution race — on top of the calibrated
abundance distributions + division-inheritance partition noise + Dna-reset division detection.

Per-lineage draws (shape/CV only): k_Cd_translation=0.801·LogN(0,0.633) (CyclinD1 CV 0.70),
P21_div=0.42·LogN(0,0.32) (birth-p27 CV 0.33), kTlEZ=0.004·LogN(0,0.51) (EZH2 CV 0.55),
mu=0.0005·LogN(0,0.22) (cell-cycle-duration axis; CV not GNP/MB-calibrated — plausible ~0.2). Birth-p27
redrawn per division (partition noise σ=0.30). Paired across mark conditions.

Protocol: settle cycling at SHH=1.0 (random phase) → measure period & mark at high Hh → withdraw to 0.15
over 100 h → hold. Divisions = Dna 1→0 reset; per-cycle G0 = pre-S∧p27>0.1 (≥2 h run = transient-G0 if a
later division follows, else arrest).

Panels: (A) transient-G0 occupancy vs Hh, WITH vs WITHOUT mark — the answer. (B) arrest occupancy vs Hh.
(C) mean H3K27me3 mark vs Hh (it rebuilds as the cell slows on withdrawal). (D) the RACE: per-cell mark
vs cell-cycle period at high Hh — fast/short period → diluted-low mark (regime 1), slow/long → plateau (2/3).

Run:  ./venv/bin/python simulations/sim_withdrawal_mark_dilution.py [--n=500] [--workers=9] [--fresh]
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

SEL = ['time', 'SHH', 'P21', 'aRc', 'Dna', 'Mk']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SHH_HI, SHH_LO = 1.00, 0.15
SETTLE, MEAS_HI, D_RAMP, HOLD, CHUNK, NP = 5000.0, 4000.0, 6000.0, 6000.0, 60.0, 9
DT = CHUNK / (NP - 1)
DWELL_CUT, P27_HI, PART = 2.0, 0.1, 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, MU_SDLOG = 0.633, 0.32, 0.51, 0.22
SHH_BINS = np.linspace(0.12, 1.03, 17)

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _state(T, DNA, P21, ARC):
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = (preS & (P21 > P27_HI)).astype(np.int8)
    st = np.zeros(len(T), np.int8)
    d = np.diff(np.concatenate([[0], g0, [0]]))
    for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
        ee = min(e, len(T) - 1)
        if T[ee] - T[s] < DWELL_CUT * 60: continue
        if np.any(div > ee): st[s:e] = 1
        else: st[s:] = 2; break
    return div, st


def _work(task):
    i, f0, ktl, p21d, ktlez, mu = task
    rng = np.random.default_rng((hash((i, f0)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = SHH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu)
    _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    nb = len(SHH_BINS) - 1; z0 = np.zeros(nb)
    # settle (desync)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return (f0, np.nan, np.nan, z0, z0.copy(), z0.copy(), z0.copy(), z0.copy())
        tm += CHUNK
    # measure period & mark at high Hh
    Th, Dh, Mh = [], [], []
    t = 0.0
    while t < MEAS_HI:
        _RR['P21_div'] = draw()
        try: r = _RR.simulate(tm, tm + CHUNK, NP, selections=['time', 'Dna', 'Mk'])
        except Exception: break
        Th.append(r['time'][1:]); Dh.append(r['Dna'][1:]); Mh.append(r['Mk'][1:]); tm += CHUNK; t += CHUNK
    period_hi, mk_hi = np.nan, np.nan
    if Th:
        Th = np.concatenate(Th); Dh = np.concatenate(Dh); Mh = np.concatenate(Mh)
        dv = np.where((Dh[:-1] > 0.9) & (Dh[1:] < 0.1))[0]
        if len(dv) >= 2: period_hi = float(np.mean(np.diff(Th[dv])) / 60.0)
        mk_hi = float(np.mean(Mh))
    # withdrawal + hold
    shh_of = lambda tt: SHH_HI + (SHH_LO - SHH_HI) * min(tt / D_RAMP, 1.0)
    T, SH, P21, ARC, DNA, MK = [], [], [], [], [], []
    t = 0.0
    while t < D_RAMP + HOLD:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        for A, s in [(T, 'time'), (SH, 'SHH'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna'), (MK, 'Mk')]:
            A.append(r[s][1:])
        t += CHUNK; tm += CHUNK
    if not T:
        return (f0, period_hi, mk_hi, z0, z0.copy(), z0.copy(), z0.copy(), z0.copy())
    T = np.concatenate(T); SH = np.concatenate(SH); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA); MK = np.concatenate(MK)
    div, st = _state(T, DNA, P21, ARC)
    bi = np.digitize(SH, SHH_BINS)
    cnt = np.array([np.sum(bi == b) for b in range(1, len(SHH_BINS))], float)
    tg0 = np.array([np.sum((bi == b) & (st == 1)) for b in range(1, len(SHH_BINS))], float)
    arr = np.array([np.sum((bi == b) & (st == 2)) for b in range(1, len(SHH_BINS))], float)
    dvh = np.histogram(SH[div], SHH_BINS)[0].astype(float) if len(div) else z0.copy()
    mks = np.array([np.sum(MK[bi == b]) for b in range(1, len(SHH_BINS))], float)
    return (f0, period_hi, mk_hi, dvh, cnt, tg0, arr, mks)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 500)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_withdrawal_mark_dilution_cache.npz'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        tasks = [(i, f0, ktl[i], p21[i], ktz[i], mu[i]) for f0 in (0.233, 1.0) for i in range(N)]
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = []
            for k, r in enumerate(pool.imap_unordered(_work, tasks, chunksize=6)):
                res.append(r)
                if (k + 1) % 200 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        store = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            sub = [r for r in res if r[0] == f0]
            CNT = np.sum([r[4] for r in sub], 0); TG0 = np.sum([r[5] for r in sub], 0)
            ARR = np.sum([r[6] for r in sub], 0); DVH = np.sum([r[3] for r in sub], 0); MKS = np.sum([r[7] for r in sub], 0)
            store[f'{tag}_tg0'] = np.where(CNT > 0, TG0 / CNT, np.nan)
            store[f'{tag}_arr'] = np.where(CNT > 0, ARR / CNT, np.nan)
            store[f'{tag}_mk'] = np.where(CNT > 0, MKS / CNT, np.nan)
            store[f'{tag}_rate'] = np.where(CNT > 0, DVH / (CNT * DT) * 1440.0, np.nan)
            store[f'{tag}_period'] = np.array([r[1] for r in sub], float)
            store[f'{tag}_mkhi'] = np.array([r[2] for r in sub], float)
        np.savez(CACHE, bins=SHH_BINS, N=N, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    b = z['bins']; xc = 0.5 * (b[:-1] + b[1:]); Nc = int(z['N'])
    COND = [('W', '#8b1a1a', 'WITH mark (repressed)'), ('N', '#3b7bbf', 'WITHOUT mark')]

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    a = ax[0, 0]
    for tag, col, nm in COND: a.plot(xc, z[f'{tag}_tg0'], '-o', color=col, lw=2.6, ms=5, label=nm)
    a.set_xlabel('Hh (SHH) during withdrawal'); a.set_ylabel('transient-G0 time-occupancy')
    a.set_title('(A) Transient-G0 occupancy vs Hh — WITH vs WITHOUT mark', fontweight='bold', fontsize=12); a.grid(alpha=0.15); a.legend(fontsize=9)
    a = ax[0, 1]
    for tag, col, nm in COND: a.plot(xc, z[f'{tag}_arr'], '-o', color=col, lw=2.6, ms=5, label=nm)
    a.set_xlabel('Hh (SHH) during withdrawal'); a.set_ylabel('deep-arrest time-occupancy')
    a.set_title('(B) Arrest occupancy vs Hh (context: where the transient-G0 goes)', fontweight='bold', fontsize=12); a.grid(alpha=0.15); a.legend(fontsize=9)
    a = ax[1, 0]
    for tag, col, nm in COND: a.plot(xc, z[f'{tag}_mk'], '-o', color=col, lw=2.6, ms=5, label=nm)
    a.set_xlabel('Hh (SHH) during withdrawal'); a.set_ylabel('mean H3K27me3 mark (Mk)')
    a.set_title('(C) Mean mark FALLS as Hh withdraws — restoration is CYCLE-GATED\n(EZH2 writer is E2f-gated → collapses as cells slow/arrest; WITH mark lower = more arrest → more EZH2 loss)', fontweight='bold', fontsize=10.5); a.grid(alpha=0.15); a.legend(fontsize=9)
    a = ax[1, 1]
    for tag, col, nm in COND:
        p = z[f'{tag}_period']; m = z[f'{tag}_mkhi']; ok = np.isfinite(p) & np.isfinite(m)
        a.plot(p[ok], m[ok], 'o', color=col, ms=4, alpha=0.4, label=nm)
    a.set_xlabel('cell-cycle period at high Hh (h)'); a.set_ylabel('H3K27me3 mark (Mk) at high Hh')
    a.set_title('(D) The dilution-vs-restoration RACE (mu-distributed cells)\nfast/short period → diluted-LOW mark (regime 1); slow → PLATEAU (regime 2/3)', fontweight='bold', fontsize=11); a.grid(alpha=0.15); a.legend(fontsize=9)

    fig.suptitle(f'Hh withdrawal: H3K27me3 dilution race + transient-G0 WITH vs WITHOUT mark  (calibrated + mu-distribution + partition noise, N={Nc}/cond)',
                 fontsize=12.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_withdrawal_mark_dilution.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_withdrawal_mark_dilution.pdf', bbox_inches='tight')
    plt.close()
    print('transient-G0 occ WITH  =', np.round(z['W_tg0'], 3))
    print('transient-G0 occ WITHOUT=', np.round(z['N_tg0'], 3))
    print('arrest occ WITH  =', np.round(z['W_arr'], 3), ' WITHOUT =', np.round(z['N_arr'], 3))
    print('Saved sim_withdrawal_mark_dilution.png')
