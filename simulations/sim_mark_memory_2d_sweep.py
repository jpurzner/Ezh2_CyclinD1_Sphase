"""2-D memory sweep: EZH2 half-life (kDeEZ, the WRITER's persistence) x H3K27me3 read-write re-accumulation
rate (k_w_mk, the CHROMATIN layer's self-reinforcement / rebuild-after-S-dilution). Tests JP's point:
the mark's memory need not ride on EZH2's lifetime -- a fast-re-accumulating, self-sustaining mark could
hold memory even with a short-lived EZH2. For each (kDeEZ, k_w_mk) runs ramp-UP (entry) + ramp-DOWN
(withdrawal) WITH mark (f0_mk=0.145); WITHOUT (f0_mk=1.0, mark can't repress -> independent of k_w_mk) is
the per-kDeEZ baseline. Heatmaps of mark level, withdrawal-arrest memory, and entry suppression.

Run:  ./venv/bin/python simulations/sim_mark_memory_2d_sweep.py [--n=35] [--workers=9] [--fresh]
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
KDEEZ = [1.5e-4, 3e-4, 6e-4, 1.2e-3, 2.4e-3]     # t1/2 77, 39, 19, 10(current), 5 h
KWMK = [0.0016, 0.004, 0.010, 0.025, 0.060]      # read-write re-accumulation (default 0.0016 -> 40x)
LOW = 0.5 * (SHH_BINS[:-1] + SHH_BINS[1:]) < 0.35
CACHE = 'simulations/sim_mark_memory_2d_sweep_cache.npz'
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
    i, proto, f0, kdeez, kwmk, ktl, p21d, ktlez, mu = task
    rng = np.random.default_rng((hash((i, proto, f0, round(kdeez, 7), round(kwmk, 6))) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = SHH_LO if proto == 'up' else SHH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu)
    _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d); _RR['kDeEZ'] = float(kdeez); _RR['k_w_mk'] = float(kwmk)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    nb = len(SHH_BINS) - 1; z = np.zeros(nb)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return (kdeez, kwmk, proto, f0, z, z.copy(), z.copy(), z.copy())
        tm += CHUNK
    T, SH, P21, ARC, DNA, MK = [], [], [], [], [], []
    t = 0.0
    while t < RECORD:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t, proto)); r = None
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
        return (kdeez, kwmk, proto, f0, z, z.copy(), z.copy(), z.copy())
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
    return (kdeez, kwmk, proto, f0, ndiv, nsamp, narr, mksum)


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
        for kd in KDEEZ:
            for proto in ('up', 'down'):
                for kw in KWMK:                              # WITH mark: full 2-D grid
                    for i in range(N): tasks.append((i, proto, 0.145, kd, kw, ktl[i], p21[i], ktz[i], mu[i]))
                for i in range(N): tasks.append((i, proto, 1.0, kd, KWMK[0], ktl[i], p21[i], ktz[i], mu[i]))  # WITHOUT: 1 per kDeEZ
        print(f'{len(tasks)} cell-runs ({len(KDEEZ)} kDeEZ x [{len(KWMK)} k_w_mk WITH + 1 WITHOUT] x 2 protos x N={N})', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=6))

        def agg(sub, kind):
            ND = np.sum([r[4] for r in sub], 0); NS = np.sum([r[5] for r in sub], 0)
            NA = np.sum([r[6] for r in sub], 0); MK = np.sum([r[7] for r in sub], 0)
            if kind == 'rate': v = np.where(NS > 0, ND / (NS * DT) * 1440.0, np.nan)
            elif kind == 'arr': v = np.where(NS > 0, NA / NS, np.nan)
            else: v = np.where(NS > 0, MK / NS, np.nan)
            return float(np.nanmean(v[LOW]))
        nk, nw = len(KDEEZ), len(KWMK)
        MKm = np.full((nk, nw), np.nan); WDmem = np.full((nk, nw), np.nan); ENTsup = np.full((nk, nw), np.nan)
        for a, kd in enumerate(KDEEZ):
            wo_dn = [r for r in res if abs(r[0]-kd) < 1e-12 and r[3] == 1.0 and r[2] == 'down']
            wo_up = [r for r in res if abs(r[0]-kd) < 1e-12 and r[3] == 1.0 and r[2] == 'up']
            wo_arr = agg(wo_dn, 'arr'); wo_rate = agg(wo_up, 'rate')
            for b, kw in enumerate(KWMK):
                w_dn = [r for r in res if abs(r[0]-kd) < 1e-12 and abs(r[1]-kw) < 1e-9 and r[3] == 0.145 and r[2] == 'down']
                w_up = [r for r in res if abs(r[0]-kd) < 1e-12 and abs(r[1]-kw) < 1e-9 and r[3] == 0.145 and r[2] == 'up']
                MKm[a, b] = agg(w_dn, 'mk')
                WDmem[a, b] = agg(w_dn, 'arr') - wo_arr
                ENTsup[a, b] = wo_rate - agg(w_up, 'rate')
        np.savez(CACHE, KDEEZ=KDEEZ, KWMK=KWMK, MKm=MKm, WDmem=WDmem, ENTsup=ENTsup)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    KD = list(z['KDEEZ']); KW = list(z['KWMK']); th = [np.log(2)/k/60 for k in KD]
    MKm, WDmem, ENTsup = z['MKm'], z['WDmem'], z['ENTsup']
    cur_a = int(np.argmin([abs(t-10) for t in th])); cur_b = 0     # current model: t1/2~10h, k_w_mk default

    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    def heat(a, M, title, cmap, vmin=None, vmax=None):
        im = a.imshow(M, origin='lower', aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
        a.set_xticks(range(len(KW))); a.set_xticklabels([f'{k:g}\n({k/KW[0]:.0f}x)' for k in KW], fontsize=8)
        a.set_yticks(range(len(th))); a.set_yticklabels([f'{t:.0f}h' for t in th], fontsize=9)
        a.set_xlabel('H3K27me3 read-write re-accumulation  k_w_mk'); a.set_ylabel('EZH2 half-life (kDeEZ)')
        a.plot(cur_b, cur_a, 'o', ms=13, mfc='none', mec='#00e5ff', mew=2.2)                # current 25/27 model
        for ii in range(len(th)):
            for jj in range(len(KW)):
                if np.isfinite(M[ii, jj]): a.text(jj, ii, f'{M[ii,jj]:.2f}', ha='center', va='center', fontsize=7.5,
                                                  color='white' if (vmax and M[ii,jj] > 0.55*vmax) else 'black')
        fig.colorbar(im, ax=a); a.set_title(title, fontweight='bold', fontsize=10.5)
    heat(ax[0], MKm, '(A) mark level Mk @low Hh (withdrawal)\ndoes faster re-accumulation restore the mark?', 'viridis', 0)
    heat(ax[1], WDmem, '(B) WITHDRAWAL memory  (WITH − WITHOUT arrest)\ncan the mark hold memory at fast kDeEZ?', 'magma', 0)
    heat(ax[2], ENTsup, '(C) ENTRY suppression  (WITHOUT − WITH div rate)', 'cividis')
    fig.suptitle('2-D memory sweep: EZH2 half-life x H3K27me3 re-accumulation (k_w_mk) — cyan ring = current 25/27 model (GNP, N=%d)' % N,
                 fontsize=12.5, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_memory_2d_sweep.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_memory_2d_sweep.pdf', bbox_inches='tight')
    print('Saved sim_mark_memory_2d_sweep.png')
    print('rows = kDeEZ t1/2:', [f'{t:.0f}h' for t in th], ' cols = k_w_mk:', KW)
    print('mark Mk:\n', np.round(MKm, 3)); print('withdrawal memory:\n', np.round(WDmem, 3))
