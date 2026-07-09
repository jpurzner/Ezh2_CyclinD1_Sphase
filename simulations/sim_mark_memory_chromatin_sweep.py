"""Which CHROMATIN parameters create the withdrawal memory, at the DATA-SUPPORTED fast EZH2 (kDeEZ t1/2
10h, the 25/27 value)? Map the two chromatin axes: read-write re-accumulation (k_w_mk) x mark durability/
turnover (del_mk). For each (k_w_mk, del_mk) runs ramp-up (entry) + ramp-down (withdrawal) ensembles WITH
(f0_mk=0.145) vs WITHOUT (f0_mk=1.0; mark can't repress -> independent of both -> single baseline).
Heatmaps of mark level, withdrawal-arrest memory, and entry suppression -> the memory-competent region.

Run:  ./venv/bin/python simulations/sim_mark_memory_chromatin_sweep.py [--n=35] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'SHH', 'P21', 'aRc', 'Dna', 'Mk']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SHH_LO, SHH_HI = 0.15, 1.00
SETTLE, D_RAMP, HOLD, CHUNK, NP = 5000.0, 6000.0, 5000.0, 60.0, 9
RECORD = D_RAMP + HOLD; DT = CHUNK / (NP - 1)
DWELL_CUT, P27_HI, PART = 2.0, 0.1, 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.633, 0.32, 0.51
SHH_BINS = np.linspace(0.12, 1.03, 17)
KDEEZ_FIXED = 0.0012                              # DATA-SUPPORTED fast EZH2, t1/2 ~10h (25/27)
KWMK = [0.0016, 0.004, 0.010, 0.025, 0.060]      # read-write re-accumulation (1x .. 40x)
DELMK = [0.0002, 0.0005, 0.001, 0.002, 0.004]    # mark turnover / durability (mark t1/2 58 .. 3 h)
LOW = 0.5 * (SHH_BINS[:-1] + SHH_BINS[1:]) < 0.35
CACHE = 'simulations/sim_mark_memory_chromatin_sweep_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 1000000)
    except Exception: pass
    try: _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception: pass


def shh_of(t, proto):
    if proto == 'up': return SHH_LO + (SHH_HI - SHH_LO) * min(t / D_RAMP, 1.0)
    return SHH_HI + (SHH_LO - SHH_HI) * min(t / D_RAMP, 1.0)


def _work(task):
    i, proto, f0, kwmk, delmk, ktl, p21d, ktlez, mu = task
    rng = np.random.default_rng((hash((i, proto, f0, round(kwmk, 6), round(delmk, 7))) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = SHH_LO if proto == 'up' else SHH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu)
    _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d); _RR['kDeEZ'] = KDEEZ_FIXED
    _RR['k_w_mk'] = float(kwmk); _RR['del_mk'] = float(delmk)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    nb = len(SHH_BINS) - 1; z = np.zeros(nb)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return (kwmk, delmk, proto, f0, z, z.copy(), z.copy(), z.copy())
        tm += CHUNK
    T, SH, P21, ARC, DNA, MK = [], [], [], [], [], []
    t = 0.0
    while t < RECORD:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t, proto)); r = None
        for atol, ms in [(1e-9, 20.0), (1e-8, 20.0), (1e-7, 8.0), (1e-6, 3.0), (1e-5, 1.5)]:
            try:
                _RR.integrator.setValue('absolute_tolerance', atol); _RR.integrator.setValue('maximum_time_step', ms)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        for A, s in [(T, 'time'), (SH, 'SHH'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna'), (MK, 'Mk')]:
            A.append(r[s][1:])
        t += CHUNK; tm += CHUNK
    if not T:
        return (kwmk, delmk, proto, f0, z, z.copy(), z.copy(), z.copy())
    T = np.concatenate(T); SH = np.concatenate(SH); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA); MK = np.concatenate(MK)
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = (preS & (P21 > P27_HI)).astype(np.int8)
    state = np.zeros(len(T), np.int8)
    d = np.diff(np.concatenate([[0], g0, [0]]))
    for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
        ee = min(e, len(T) - 1)
        if T[ee] - T[s] < DWELL_CUT * 60: continue
        if np.any(div > ee): state[s:e] = 1
        else: state[s:] = 2; break
    bi = np.digitize(SH, SHH_BINS)
    nsamp = np.array([np.sum(bi == b) for b in range(1, len(SHH_BINS))], float)
    ndiv = np.histogram(SH[div], SHH_BINS)[0].astype(float) if len(div) else z.copy()
    narr = np.array([np.sum((bi == b) & (state == 2)) for b in range(1, len(SHH_BINS))], float)
    mksum = np.array([np.sum(MK[bi == b]) for b in range(1, len(SHH_BINS))], float)
    return (kwmk, delmk, proto, f0, ndiv, nsamp, narr, mksum)


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    N = int(_arg('--n', 35)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        rng = np.random.default_rng(7)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.full(N, MU0)
        tasks = []
        for proto in ('up', 'down'):
            for kw in KWMK:
                for dm in DELMK:
                    for i in range(N): tasks.append((i, proto, 0.145, kw, dm, ktl[i], p21[i], ktz[i], mu[i]))
            for i in range(N): tasks.append((i, proto, 1.0, KWMK[0], DELMK[0], ktl[i], p21[i], ktz[i], mu[i]))  # WITHOUT baseline
        print(f'{len(tasks)} cell-runs ({len(KWMK)}x{len(DELMK)} WITH + 1 WITHOUT, 2 protos, N={N}) at kDeEZ t1/2 10h', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=6))

        def agg(sub, kind):
            NS = np.sum([r[5] for r in sub], 0)
            if kind == 'rate': v = np.where(NS > 0, np.sum([r[4] for r in sub], 0) / (NS * DT) * 1440.0, np.nan)
            elif kind == 'arr': v = np.where(NS > 0, np.sum([r[6] for r in sub], 0) / NS, np.nan)
            else: v = np.where(NS > 0, np.sum([r[7] for r in sub], 0) / NS, np.nan)
            return float(np.nanmean(v[LOW]))
        wo_arr = agg([r for r in res if r[3] == 1.0 and r[2] == 'down'], 'arr')
        wo_rate = agg([r for r in res if r[3] == 1.0 and r[2] == 'up'], 'rate')
        nk, nd = len(KWMK), len(DELMK)
        MKm = np.full((nd, nk), np.nan); WDmem = np.full((nd, nk), np.nan); ENTsup = np.full((nd, nk), np.nan)
        for b, kw in enumerate(KWMK):
            for a, dm in enumerate(DELMK):
                dn = [r for r in res if abs(r[0]-kw) < 1e-9 and abs(r[1]-dm) < 1e-12 and r[3] == 0.145 and r[2] == 'down']
                up = [r for r in res if abs(r[0]-kw) < 1e-9 and abs(r[1]-dm) < 1e-12 and r[3] == 0.145 and r[2] == 'up']
                MKm[a, b] = agg(dn, 'mk'); WDmem[a, b] = agg(dn, 'arr') - wo_arr; ENTsup[a, b] = wo_rate - agg(up, 'rate')
        np.savez(CACHE, KWMK=KWMK, DELMK=DELMK, MKm=MKm, WDmem=WDmem, ENTsup=ENTsup, wo_arr=wo_arr, wo_rate=wo_rate)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    KW = list(z['KWMK']); DM = list(z['DELMK']); MKm, WDmem, ENTsup = z['MKm'], z['WDmem'], z['ENTsup']
    mkt = [np.log(2)/d/60 for d in DM]
    def _fpos(val, grid):                            # fractional (log-interp) grid index for the marker
        lg = np.log(np.array(grid, float)); lv = np.log(val)
        if lv <= lg[0]: return 0.0
        if lv >= lg[-1]: return float(len(grid) - 1)
        j = int(np.searchsorted(lg, lv)) - 1
        return j + (lv - lg[j]) / (lg[j + 1] - lg[j])
    cur_b = _fpos(0.00449887, KW)                    # current baked k_w_mk
    cur_a = _fpos(0.00328353, DM)                    # current baked del_mk (mark t½ 3.5h)

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    def heat(a, M, title, cmap, vmin=None, vmax=None):
        im = a.imshow(M, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        a.set_xticks(range(len(KW))); a.set_xticklabels([f'{k:g}\n({k/KW[0]:.0f}x)' for k in KW], fontsize=8)
        a.set_yticks(range(len(DM))); a.set_yticklabels([f'{d:g}\n(t½ {t:.0f}h)' for d, t in zip(DM, mkt)], fontsize=8)
        a.set_xlabel('read-write re-accumulation  k_w_mk'); a.set_ylabel('mark durability  del_mk (turnover)')
        a.plot(cur_b, cur_a, 'o', ms=13, mfc='none', mec='#00e5ff', mew=2.2)
        for ii in range(len(DM)):
            for jj in range(len(KW)):
                if np.isfinite(M[ii, jj]): a.text(jj, ii, f'{M[ii,jj]:.2f}', ha='center', va='center', fontsize=7.5,
                                                  color='white' if (vmax and M[ii, jj] > 0.55*vmax) else 'black')
        fig.colorbar(im, ax=a); a.set_title(title, fontweight='bold', fontsize=10.5)
    heat(ax[0], MKm, '(A) mark level Mk @low Hh (withdrawal)', 'viridis', 0)
    heat(ax[1], WDmem, '(B) WITHDRAWAL memory (WITH − WITHOUT arrest)\nthe memory-competent region', 'magma', 0)
    heat(ax[2], ENTsup, '(C) ENTRY suppression (WITHOUT − WITH div rate)', 'cividis')
    fig.suptitle('Which CHROMATIN parameters create the memory — at the DATA-SUPPORTED fast EZH2 (kDeEZ t½ 10h). cyan ring = current calibration (GNP, N=%d)' % int(_arg('--n', 35)),
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_memory_chromatin_sweep.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_memory_chromatin_sweep.pdf', bbox_inches='tight')
    print('Saved sim_mark_memory_chromatin_sweep.png')
    print('rows=del_mk', DM, ' cols=k_w_mk', KW)
    print('withdrawal memory:\n', np.round(WDmem, 3)); print('mark level:\n', np.round(MKm, 3))
