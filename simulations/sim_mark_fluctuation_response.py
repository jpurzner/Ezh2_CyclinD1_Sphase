"""How do WITH vs WITHOUT the H3K27me3 mark respond to FLUCTUATIONS in Hh signalling?
Frequency response of proliferation. Drive Hh(t) = Hh0 + A·sin(2π t/T) at a range of drive periods T
(fast→slow), with Hh0 near the commitment threshold. Measure the division-rate MODULATION (how much
proliferation tracks the fluctuation) as a function of drive frequency, WITH vs WITHOUT the mark.

Hypothesis: the mark is a slow integrator → an extra LOW-PASS stage. Fast Hh fluctuations should be
BUFFERED (small proliferation modulation) more WITH the mark; slow fluctuations tracked by both. The cell
cycle (~1 cycle response time) is itself a filter; comparing WITH vs WITHOUT isolates the mark's added
filtering and its cutoff.

Ensemble: calibrated abundance dists (CyclinD1 CV 0.70 / birth-p27 CV 0.33 / EZH2 CV 0.55) + growth-time
(mu) dist + partition noise; divisions = Dna 1→0 reset. Phase-binning over the drive period = the
proliferation waveform over one fluctuation cycle.

Panels: (A) rate vs drive phase at a FAST drive (buffering), (B) at a SLOW drive (tracking), (C) modulation
index vs drive period — the low-pass (WITH vs WITHOUT), (D) the mark's own amplitude vs drive period.

Run:  ./venv/bin/python simulations/sim_mark_fluctuation_response.py [--n=80] [--workers=9] [--fresh]
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

SEL = ['time', 'Dna', 'Mk', 'Cd', 'P21', 'aRc']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH0, AMP = 0.45, 0.30                                    # mean Hh (near threshold) + fluctuation amplitude
T_DRIVE_H = [12.0, 20.0, 33.0, 55.0, 90.0, 150.0]        # drive periods (hours), fast -> slow
NB = 12                                                  # phase bins over one drive period
SETTLE, N_PER, CHUNK, NP = 5000.0, 4, 60.0, 9            # settle + 4 drive periods recorded
DT = CHUNK / (NP - 1)
PART = 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, MU_SDLOG = 0.633, 0.32, 0.51, 0.22
TWO_PI = 2.0 * np.pi

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, T_h, f0, ktl, p21d, ktlez, mu = task
    T = T_h * 60.0
    rng = np.random.default_rng((hash((i, round(T_h, 1), f0)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = HH0
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu); _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    z0 = np.zeros(NB); NULL = (T_h, f0, z0, z0.copy(), z0.copy(), z0.copy(), z0.copy())
    # settle at the mean (random phase to desync lineages)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return NULL
        tm += CHUNK
    hh = lambda tt: float(np.clip(HH0 + AMP * np.sin(TWO_PI * tt / T), 0.02, 1.10))
    DNA, MK, CD, P21, ARC, PH = [], [], [], [], [], []
    t = 0.0; TREC = N_PER * T
    while t < TREC:
        _RR['P21_div'] = draw(); _RR['SHH'] = hh(t); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        tt = r['time'][1:]
        DNA.append(r['Dna'][1:]); MK.append(r['Mk'][1:]); CD.append(r['Cd'][1:]); P21.append(r['P21'][1:]); ARC.append(r['aRc'][1:])
        PH.append((TWO_PI * (t + (tt - tt[0])) / T) % TWO_PI)
        t += CHUNK; tm += CHUNK
    if not DNA:
        return NULL
    DNA = np.concatenate(DNA); MK = np.concatenate(MK); CD = np.concatenate(CD); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); PH = np.concatenate(PH)
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; g0inst = (~inS & ~inG2) & (P21 > 0.1)   # instantaneous p27-high G0 (commitment readout)
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    edges = np.linspace(0, TWO_PI, NB + 1)
    bi = np.clip(np.digitize(PH, edges) - 1, 0, NB - 1)
    cnt = np.array([np.sum(bi == b) for b in range(NB)], float)
    mks = np.array([np.sum(MK[bi == b]) for b in range(NB)], float)
    cds = np.array([np.sum(CD[bi == b]) for b in range(NB)], float)
    g0s = np.array([np.sum(g0inst[bi == b]) for b in range(NB)], float)
    divh = np.zeros(NB)
    if len(div):
        dbi = np.clip(np.digitize(PH[div], edges) - 1, 0, NB - 1)
        for b in dbi: divh[b] += 1
    return (T_h, f0, divh, cnt, mks, cds, g0s)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 80)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_mark_fluctuation_response_cache.npz'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        tasks = [(i, T_h, f0, ktl[i], p21[i], ktz[i], mu[i]) for f0 in (0.233, 1.0) for T_h in T_DRIVE_H for i in range(N)]
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = []
            for k, r in enumerate(pool.imap_unordered(_work, tasks, chunksize=6)):
                res.append(r)
                if (k + 1) % 200 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        store = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            rate = np.full((len(T_DRIVE_H), NB), np.nan); mk = np.full((len(T_DRIVE_H), NB), np.nan)
            cd = np.full((len(T_DRIVE_H), NB), np.nan); g0 = np.full((len(T_DRIVE_H), NB), np.nan)
            for ti, T_h in enumerate(T_DRIVE_H):
                sub = [r for r in res if r[0] == T_h and r[1] == f0]
                D = np.sum([r[2] for r in sub], 0); C = np.sum([r[3] for r in sub], 0); M = np.sum([r[4] for r in sub], 0)
                CDs = np.sum([r[5] for r in sub], 0); G0s = np.sum([r[6] for r in sub], 0)
                rate[ti] = np.where(C > 0, D / (C * DT) * 1440.0, np.nan)
                mk[ti] = np.where(C > 0, M / C, np.nan); cd[ti] = np.where(C > 0, CDs / C, np.nan); g0[ti] = np.where(C > 0, G0s / C, np.nan)
            store[f'{tag}_rate'] = rate; store[f'{tag}_mk'] = mk; store[f'{tag}_cd'] = cd; store[f'{tag}_g0'] = g0
        np.savez(CACHE, T=np.array(T_DRIVE_H), NB=NB, HH0=HH0, AMP=AMP, N=N, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    T = z['T']; NB = int(z['NB']); Nc = int(z['N'])
    ph = (np.arange(NB) + 0.5) / NB * TWO_PI
    hh_ph = HH0 + AMP * np.sin(ph)

    def fund(y):                                          # fundamental amplitude / mean = modulation index
        y = np.asarray(y, float)
        a = 2 * np.mean(y * np.cos(ph)); b = 2 * np.mean(y * np.sin(ph))
        return np.hypot(a, b) / (np.mean(y) + 1e-9)

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    ifast, islow = 0, len(T) - 1

    def phase_panel(a, ti, ttl):
        a2 = a.twinx(); a2.plot(np.degrees(ph), hh_ph, color='#aaa', lw=1.8, alpha=0.7); a2.plot(np.degrees(ph), hh_ph, 'o', color='#aaa', ms=3, alpha=0.7)
        a2.set_ylabel('Hh drive (grey)', color='#888'); a2.set_ylim(0, 1.1)
        for tag, col, nm in [('W', '#8b1a1a', 'WITH mark'), ('N', '#3b7bbf', 'WITHOUT mark')]:
            a.plot(np.degrees(ph), z[f'{tag}_cd'][ti], '-o', color=col, lw=2.6, ms=5, label=nm)
        a.set_xlabel('drive phase (deg)'); a.set_ylabel('cycle-averaged CyclinD1')
        a.set_title(ttl, fontweight='bold', fontsize=11.5); a.legend(fontsize=9, loc='upper left'); a.grid(alpha=0.15)

    phase_panel(ax[0, 0], ifast, f'(A) FAST fluctuation (T={T[ifast]:.0f} h): CyclinD1 tracks Hh in BOTH\nWITH mark = lower mean (repressed) but same relative swing — not buffered')
    phase_panel(ax[0, 1], islow, f'(B) SLOW fluctuation (T={T[islow]:.0f} h): CyclinD1 tracks the drive (both)')

    # (C) modulation index vs drive period: CyclinD1 (mechanism) + G0/commitment fraction (proliferation)
    a = ax[1, 0]
    for tag, col, nm in [('W', '#8b1a1a', 'WITH mark'), ('N', '#3b7bbf', 'WITHOUT mark')]:
        modcd = np.array([fund(z[f'{tag}_cd'][ti]) for ti in range(len(T))])
        modg0 = np.array([fund(z[f'{tag}_g0'][ti]) for ti in range(len(T))])
        a.plot(T, modcd, '-o', color=col, lw=2.7, ms=6, label=f'{nm} — CyclinD1 (relative)')
        a.plot(T, modg0, ':s', color=col, lw=2.0, ms=5, alpha=0.8, label=f'{nm} — G0/commitment')
    a.set_xscale('log'); a.set_xlabel('drive period (h)  — faster ← → slower'); a.set_ylabel('modulation index (0 = buffered, high = tracks Hh)')
    a.set_title('(C) CyclinD1 relative swing ≈ equal & flat (mark SCALES, no low-pass);\nG0/commitment (dotted) IS low-pass filtered (buffered vs fast) but ~EQUAL WITH/WITHOUT →\nthe buffering is the cell cycle + commitment switch, NOT the mark', fontweight='bold', fontsize=9.5); a.legend(fontsize=7.5, loc='upper left'); a.grid(alpha=0.15, which='both')

    # (D) the mark's own amplitude vs drive period (the filter itself)
    a = ax[1, 1]
    mka = np.array([(np.nanmax(z['W_mk'][ti]) - np.nanmin(z['W_mk'][ti])) for ti in range(len(T))])
    a.plot(T, mka, '-o', color='#8e44ad', lw=2.7, ms=6, label='H3K27me3 mark amplitude (WITH mark)')
    a.set_xscale('log'); a.set_xlabel('drive period (h)  — faster ← → slower'); a.set_ylabel('mark (Mk) peak-to-trough over drive')
    a.set_title('(D) The mark is the slow filter: its own oscillation is small for\nFAST drives (can\'t follow) and grows for SLOW drives', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both')

    fig.suptitle(f'Response to Hh FLUCTUATIONS — WITH vs WITHOUT H3K27me3 repression (Hh0={HH0}±{AMP} sine, calibrated+mu+partition, N={Nc}/cond)',
                 fontsize=12.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_fluctuation_response.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_fluctuation_response.pdf', bbox_inches='tight')
    plt.close()
    for tag, nm in [('W', 'WITH   '), ('N', 'WITHOUT')]:
        print(f'{nm} modulation index by drive period {list(T)} =', np.round([fund(z[f"{tag}_rate"][ti]) for ti in range(len(T))], 3))
    print('Saved sim_mark_fluctuation_response.png')
