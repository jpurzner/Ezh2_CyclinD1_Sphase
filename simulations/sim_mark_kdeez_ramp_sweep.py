"""Sweep EZH2 protein half-life (kDeEZ) and observe the MARK's ramp-UP (entry) and ramp-DOWN
(withdrawal) effects. Tests the tension the 25/27 recalibration exposed: a FASTER EZH2 (short t1/2)
makes the H3K27me3 mark a poor slow-integrator -> weak entry/withdrawal/hysteresis effects; a SLOWER
EZH2 (long t1/2) lets the mark integrate -> strong effects. For each kDeEZ, runs ramp-up + ramp-down
ensembles WITH (f0_mk=0.145) vs WITHOUT (f0_mk=1.0) and reports low-Hh division rate, arrest occupancy,
mark level, and the entry/exit hysteresis gap vs t1/2.

Run:  ./venv/bin/python simulations/sim_mark_kdeez_ramp_sweep.py [--n=50] [--workers=9] [--fresh]
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
KDEEZ = [9e-5, 1.5e-4, 3e-4, 6e-4, 1.2e-3, 2.4e-3, 5e-3]     # t1/2 128,77,39,19,10,5,2.3 h
LOW = 0.5 * (SHH_BINS[:-1] + SHH_BINS[1:]) < 0.35            # low-Hh bins (entry/exit regime)
CACHE = 'simulations/sim_mark_kdeez_ramp_sweep_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 1000000)
    except Exception: pass
    try: _RR.integrator.setValue('maximum_time_step', 20.0)      # faster kDeEZ -> stiff withdrawal arrests
    except Exception: pass


def shh_of(t, proto):
    if proto == 'up': return SHH_LO + (SHH_HI - SHH_LO) * min(t / D_RAMP, 1.0)
    return SHH_HI + (SHH_LO - SHH_HI) * min(t / D_RAMP, 1.0)


def _work(task):
    i, proto, f0, kdeez, ktl, p21d, ktlez, mu = task
    rng = np.random.default_rng((hash((i, proto, f0, round(kdeez, 7))) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = SHH_LO if proto == 'up' else SHH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu)
    _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d); _RR['kDeEZ'] = float(kdeez)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    nb = len(SHH_BINS) - 1; z = np.zeros(nb)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return (kdeez, proto, f0, z, z.copy(), z.copy(), z.copy())
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
        return (kdeez, proto, f0, z, z.copy(), z.copy(), z.copy())
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
    return (kdeez, proto, f0, ndiv, nsamp, narr, mksum)


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    N = int(_arg('--n', 50)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        rng = np.random.default_rng(7)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.full(N, MU0)
        tasks = [(i, proto, f0, kd, ktl[i], p21[i], ktz[i], mu[i])
                 for kd in KDEEZ for f0 in (0.145, 1.0) for proto in ('up', 'down') for i in range(N)]
        print(f'{len(tasks)} cell-runs ({len(KDEEZ)} kDeEZ x 2 f0 x 2 protos x N={N})', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=6))
        store = {}
        for kd in KDEEZ:
            for f0, tag in [(0.145, 'W'), (1.0, 'N')]:
                for proto in ('up', 'down'):
                    sub = [r for r in res if abs(r[0]-kd) < 1e-9 and r[2] == f0 and r[1] == proto]
                    ND = np.sum([r[3] for r in sub], 0); NS = np.sum([r[4] for r in sub], 0)
                    NA = np.sum([r[5] for r in sub], 0); MK = np.sum([r[6] for r in sub], 0)
                    key = f'{kd:.1e}_{tag}_{proto}'
                    store[key + '_rate'] = np.where(NS > 0, ND / (NS * DT) * 1440.0, np.nan)
                    store[key + '_arr'] = np.where(NS > 0, NA / NS, np.nan)
                    store[key + '_mk'] = np.where(NS > 0, MK / NS, np.nan)
        np.savez(CACHE, KDEEZ=KDEEZ, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    KD = list(z['KDEEZ']); th = [np.log(2) / k / 60 for k in KD]
    def lowmean(kd, tag, proto, kind):
        v = z[f'{kd:.1e}_{tag}_{proto}_{kind}']; return float(np.nanmean(v[LOW]))
    ent = {t: [lowmean(kd, t, 'up', 'rate') for kd in KD] for t in 'WN'}      # entry proliferation (low Hh)
    arr = {t: [lowmean(kd, t, 'down', 'arr') for kd in KD] for t in 'WN'}     # withdrawal arrest (low Hh)
    mk = {t: [lowmean(kd, t, 'down', 'mk') for kd in KD] for t in 'WN'}       # mark level (low Hh, withdrawal)
    hys = {t: [lowmean(kd, t, 'down', 'rate') - lowmean(kd, t, 'up', 'rate') for kd in KD] for t in 'WN'}

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))
    def mark_models(a):
        a.axvline(10, color='#8b1a1a', ls=':', lw=1.2, alpha=0.6); a.axvline(77, color='grey', ls=':', lw=1.2, alpha=0.6)
        a.text(10, a.get_ylim()[1], ' 25/27 model\n (t½10h)', color='#8b1a1a', fontsize=7, va='top')
        a.text(77, a.get_ylim()[1], ' old\n (77h)', color='grey', fontsize=7, va='top')
    a = ax[0, 0]
    a.plot(th, ent['W'], '-o', color='#8b1a1a', label='WITH mark'); a.plot(th, ent['N'], '-o', color='#3b7bbf', label='WITHOUT')
    a.set_xscale('log'); a.set_xlabel('EZH2 half-life (h)'); a.set_ylabel('entry division rate @low Hh (/day)')
    a.set_title('(A) RAMP-UP entry — the gap = mark ENTRY suppression', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both'); mark_models(a)
    a = ax[0, 1]
    a.plot(th, arr['W'], '-o', color='#8b1a1a', label='WITH mark'); a.plot(th, arr['N'], '-o', color='#3b7bbf', label='WITHOUT')
    a.set_xscale('log'); a.set_xlabel('EZH2 half-life (h)'); a.set_ylabel('withdrawal arrest occupancy @low Hh')
    a.set_title('(B) RAMP-DOWN withdrawal — the gap = mark WITHDRAWAL memory', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both'); mark_models(a)
    a = ax[1, 0]
    a.plot(th, mk['W'], '-o', color='#2e8b57')
    a.set_xscale('log'); a.set_xlabel('EZH2 half-life (h)'); a.set_ylabel('mean mark Mk @low Hh (withdrawal)')
    a.set_title('(C) MECHANISM — mark accumulates as EZH2 lives longer', fontweight='bold', fontsize=11); a.grid(alpha=0.15, which='both'); mark_models(a)
    a = ax[1, 1]
    a.plot(th, hys['W'], '-o', color='#8b1a1a', label='WITH mark'); a.plot(th, hys['N'], '-o', color='#3b7bbf', label='WITHOUT')
    a.axhline(0, color='k', lw=0.7, alpha=0.4)
    a.set_xscale('log'); a.set_xlabel('EZH2 half-life (h)'); a.set_ylabel('hysteresis gap (down − up rate) @low Hh')
    a.set_title('(D) HYSTERESIS loop width vs EZH2 half-life', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both'); mark_models(a)
    fig.suptitle('EZH2 half-life (kDeEZ) sweep — the H3K27me3 mark\'s ramp-up/down effects require a SLOW-turnover EZH2 (GNP, N=%d/cond)' % int(_arg('--n', 50)),
                 fontsize=12.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_kdeez_ramp_sweep.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_kdeez_ramp_sweep.pdf', bbox_inches='tight')
    print('Saved sim_mark_kdeez_ramp_sweep.png')
    for kd, t in zip(KD, th):
        print(f't½ {t:5.0f}h: entry W/N {lowmean(kd,"W","up","rate"):.2f}/{lowmean(kd,"N","up","rate"):.2f}  arrest W/N {lowmean(kd,"W","down","arr"):.2f}/{lowmean(kd,"N","down","arr"):.2f}  Mk {lowmean(kd,"W","down","mk"):.2f}')
