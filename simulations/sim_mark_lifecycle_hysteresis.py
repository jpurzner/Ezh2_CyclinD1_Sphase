"""Ramp-UP (entry) + STEADY-STATE + ramp-DOWN (withdrawal), WITH vs WITHOUT H3K27me3 repression of CyclinD1,
using the UPGRADED single-cell machinery: data-calibrated abundance distributions, division-inheritance
PARTITION NOISE on birth-p27, and Dna-reset division detection (robust; MPF peaks under-sample at coarse dt).

Per-lineage calibrated draws (shape/CV only): k_Cd_translation=0.801·LogN(0,0.633) (CyclinD1 CV 0.70),
P21_div median 0.42·LogN(0,0.32) (birth-p27 CV 0.33), kTlEZ=0.004·LogN(0,0.51) (EZH2 abundance CV ~0.55).
At every division birth-p27 is redrawn mean-preserving (σ=0.30) = partitioning noise. Paired across all
conditions (same cells; only f0_mk and the Hh protocol differ).

Protocols (SHH 0.15..1.0):
  UP    : settle ARRESTED at 0.15 -> ramp to 1.0 over 100 h -> hold      (ENTRY)
  STEADY: constant Hh at each of several levels                          (EQUILIBRIUM)
  DOWN  : settle CYCLING at 1.0 -> ramp to 0.15 over 100 h -> hold       (WITHDRAWAL)

Per-time state (probe_transient_g0 conventions): divisions = Dna 1->0 reset; a pre-S∧p27>0.1 run ≥2 h that
is later followed by a division = TRANSIENT G0, one that never resolves = ARREST, else actively cycling.
Readouts vs the INSTANTANEOUS Hh: fraction cycling (=1-arrest) and transient-G0 occupancy -> the entry/
steady/exit curves form a hysteresis loop; comparing WITH vs WITHOUT the mark shows the mark's effect.

Run:  ./venv/bin/python simulations/sim_mark_lifecycle_hysteresis.py [--n=200] [--workers=9] [--fresh]
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

SEL = ['time', 'SHH', 'P21', 'aRc', 'Dna']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SHH_LO, SHH_HI = 0.15, 1.00
STEADY_LEVELS = [0.20, 0.35, 0.50, 0.65, 0.80, 1.00]
SETTLE, D_RAMP, HOLD, CHUNK, NP = 5000.0, 6000.0, 5000.0, 60.0, 9
RECORD = D_RAMP + HOLD
DWELL_CUT, P27_HI, PART = 2.0, 0.1, 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.633, 0.32, 0.51   # MU_SDLOG (growth-time / cell-cycle-duration) via --sdmu
NDS = 150                                                          # downsample points per cell
SHH_BINS = np.linspace(0.12, 1.03, 17)          # 16 Hh bins (coarse enough for stable ramp division-rate)

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def shh_of(t, proto, level):
    if proto == 'steady': return level
    if proto == 'up': return SHH_LO + (SHH_HI - SHH_LO) * min(t / D_RAMP, 1.0)
    return SHH_HI + (SHH_LO - SHH_HI) * min(t / D_RAMP, 1.0)       # down


def _work(task):
    i, proto, level, f0, ktl, p21d, ktlez, mu = task
    rng = np.random.default_rng((hash((i, proto, level, f0)) & 0xFFFFFFFF))
    settle_shh = {'up': SHH_LO, 'down': SHH_HI, 'steady': level}[proto]
    _RR.reset(); _RR['SHH'] = settle_shh
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu); _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    # settle with partition noise + RANDOM extra phase (desynchronize lineages so ramp divisions don't
    # alias into Hh bins as synchronized waves); track model time tm across settle+record.
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: break
        tm += CHUNK
    # record window
    T, SH, P21, ARC, DNA = [], [], [], [], []
    t = 0.0
    while t < RECORD:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t, proto, level)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        for A, s in [(T, 'time'), (SH, 'SHH'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna')]:
            A.append(r[s][1:])
        t += CHUNK; tm += CHUNK
    if not T:
        z0 = np.zeros(len(SHH_BINS) - 1)
        return (proto, f0, z0, z0.copy(), z0.copy(), z0.copy())
    T = np.concatenate(T); SH = np.concatenate(SH); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA)
    dt = float(np.median(np.diff(T)))
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]         # Dna reset = division
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = (preS & (P21 > P27_HI)).astype(np.int8)
    state = np.zeros(len(T), np.int8)                            # 0 cycling / 1 transient-G0 / 2 arrest
    d = np.diff(np.concatenate([[0], g0, [0]]))
    for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
        ee = min(e, len(T) - 1)
        if T[ee] - T[s] < DWELL_CUT * 60: continue
        if np.any(div > ee):
            state[s:e] = 1
        else:
            state[s:] = 2; break
    # per-Hh-bin histograms (proliferation as DIVISION RATE = divisions per cell-time in that Hh bin)
    bi = np.digitize(SH, SHH_BINS); nb = len(SHH_BINS) - 1
    time_h = np.array([np.sum(bi == b) for b in range(1, len(SHH_BINS))], float) * dt      # cell-minutes per bin
    tg0_h = np.array([np.sum((bi == b) & (state == 1)) for b in range(1, len(SHH_BINS))], float) * dt
    arr_h = np.array([np.sum((bi == b) & (state == 2)) for b in range(1, len(SHH_BINS))], float) * dt
    div_h = np.histogram(SH[div], SHH_BINS)[0].astype(float) if len(div) else np.zeros(nb)
    return (proto, f0, div_h, time_h, tg0_h, arr_h)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 200)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    SDMU = float(_arg('--sdmu', 0.0))                             # growth-time (cell-cycle-duration) spread; 0 = fixed
    SUF = '' if SDMU == 0 else f'_mu{SDMU:.2f}'
    CACHE = f'simulations/sim_mark_lifecycle_hysteresis_cache{SUF}.npz'
    OUT = f'simulations/sim_mark_lifecycle_hysteresis{SUF}'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl_i = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21d_i = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktlez_i = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu_i = np.clip(MU0 * np.exp(rng.normal(0, SDMU, N)), 0.00028, 0.00085) if SDMU > 0 else np.full(N, MU0)
        tasks = []
        for f0 in (0.233, 1.0):
            for i in range(N):
                tasks.append((i, 'up', 0.0, f0, ktl_i[i], p21d_i[i], ktlez_i[i], mu_i[i]))
                tasks.append((i, 'down', 0.0, f0, ktl_i[i], p21d_i[i], ktlez_i[i], mu_i[i]))
                for L in STEADY_LEVELS:
                    tasks.append((i, 'steady', L, f0, ktl_i[i], p21d_i[i], ktlez_i[i], mu_i[i]))
        print(f'running {len(tasks)} cell-runs on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            results = []
            for k, res in enumerate(pool.imap_unordered(_work, tasks, chunksize=8)):
                results.append(res)
                if (k + 1) % 1000 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)

        store = {}
        nb = len(SHH_BINS) - 1
        for proto in ('up', 'steady', 'down'):
            for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
                subs = [r for r in results if r[0] == proto and r[1] == f0]
                if not subs:
                    for k in ('rate', 'tg0', 'arr'): store[f'{proto}_{tag}_{k}'] = np.full(nb, np.nan)
                    continue
                D = np.sum([r[2] for r in subs], 0); TT = np.sum([r[3] for r in subs], 0)
                G = np.sum([r[4] for r in subs], 0); A = np.sum([r[5] for r in subs], 0)
                store[f'{proto}_{tag}_rate'] = np.where(TT > 0, D / TT * 1440.0, np.nan)   # divisions per cell per day
                store[f'{proto}_{tag}_tg0'] = np.where(TT > 0, G / TT, np.nan)             # transient-G0 time-occupancy
                store[f'{proto}_{tag}_arr'] = np.where(TT > 0, A / TT, np.nan)             # arrest time-occupancy
        np.savez(CACHE, bins=SHH_BINS, N=N, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    b = z['bins']; xc = 0.5 * (b[:-1] + b[1:]); Nc = int(z['N'])
    PROT = [('up', 'RAMP-UP (entry)', '#1f8a4c'), ('steady', 'STEADY STATE', '#333'), ('down', 'RAMP-DOWN (withdrawal)', '#c0392b')]
    rmax = np.nanmax([np.nanmax(z[f'{p}_{t}_rate']) for p, _, _ in PROT for t in 'WN']) * 1.08

    def sm(y):                                                   # 3-pt smooth (nan-safe) for the ramp curves
        y = np.asarray(y, float); ys = y.copy()
        for k in range(1, len(y) - 1):
            w = y[k - 1:k + 2]
            if not np.any(np.isnan(w)): ys[k] = 0.25 * w[0] + 0.5 * w[1] + 0.25 * w[2]
        return ys

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    # (A)-(C): division rate vs Hh (proliferation), WITH vs WITHOUT mark, per protocol
    for a, (proto, ttl, _c), lab in zip([ax[0, 0], ax[0, 1], ax[1, 0]], PROT, 'ABC'):
        for tag, col, nm in [('W', '#8b1a1a', 'WITH mark (repressed)'), ('N', '#3b7bbf', 'WITHOUT mark')]:
            if proto == 'steady':
                a.plot(xc, z[f'{proto}_{tag}_rate'], 'o-', color=col, lw=1.6, ms=8, label=nm)
            else:
                a.plot(xc, sm(z[f'{proto}_{tag}_rate']), '-', color=col, lw=2.8, label=nm)
                a.plot(xc, z[f'{proto}_{tag}_rate'], 'o', color=col, ms=3, alpha=0.35)
        a.set_xlabel('instantaneous Hh (SHH)'); a.set_ylabel('division rate (per cell per day)')
        a.set_title(f'({lab}) {ttl} — proliferation vs Hh', fontweight='bold', fontsize=12)
        a.set_ylim(-0.03, rmax); a.grid(alpha=0.15); a.legend(fontsize=9, loc='upper left')

    # overlay ARREST occupancy (dotted, right axis) on the ramps: with a growth-time spread the mark's
    # withdrawal effect is invisible in the RATE (fast low-mark cyclers, below the ~16h dilution crossover,
    # dominate divisions & coast) but clear in arrest occupancy (the slow high-mark cells).
    for a, proto in [(ax[0, 0], 'up'), (ax[1, 0], 'down')]:
        a2 = a.twinx()
        for tag, col in [('W', '#8b1a1a'), ('N', '#3b7bbf')]:
            a2.plot(xc, z[f'{proto}_{tag}_arr'], ':', color=col, lw=2.2, alpha=0.8)
        a2.set_ylim(0, 0.6); a2.set_ylabel('arrest occupancy (dotted)', color='#555', fontsize=9); a2.tick_params(labelsize=8, colors='#555')
    ax[1, 0].set_title('(C) RAMP-DOWN — division rate is FLAT (fast low-mark cyclers coast & dominate divisions),\nbut ARREST occupancy (dotted) shows the mark still arrests the slow cells', fontweight='bold', fontsize=10.5)

    # (D) hysteresis synthesis: division rate vs Hh — ramps smoothed lines, steady as points; WITH solid / WITHOUT faded
    a = ax[1, 1]
    for proto, ttl, col in PROT:
        yW, yN = z[f'{proto}_W_rate'], z[f'{proto}_N_rate']
        if proto == 'steady':
            a.plot(xc, yW, 'o', color=col, ms=9, label=f'{ttl} — WITH'); a.plot(xc, yN, 'o', mfc='none', mec=col, ms=9, label=f'{ttl} — WITHOUT')
        else:
            a.plot(xc, sm(yW), '-', color=col, lw=3.0, label=f'{ttl} — WITH mark')
            a.plot(xc, sm(yN), '--', color=col, lw=1.9, alpha=0.5, label=f'{ttl} — WITHOUT')
    a.set_xlabel('instantaneous Hh (SHH)'); a.set_ylabel('division rate (per cell per day)')
    a.set_title('(D) Hysteresis: at low Hh, entry (green, lags LOW) < steady (pts) < withdrawal (red, coasts HIGH)\nWITH mark (solid) opens the wider loop', fontweight='bold', fontsize=11)
    a.set_ylim(-0.03, rmax); a.grid(alpha=0.15); a.legend(fontsize=7.5, loc='upper left')

    mutxt = f'+ growth-time dist (mu CV {SDMU:.2f})' if SDMU > 0 else 'fixed growth time'
    fig.suptitle(f'Ramp-up / steady-state / ramp-down — WITH vs WITHOUT H3K27me3 repression of CyclinD1  '
                 f'(calibrated distributions + partition noise, {mutxt}, N={Nc}/condition)', fontsize=12.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig(f'{OUT}.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'{OUT}.pdf', bbox_inches='tight')
    plt.close()
    for proto, ttl, _c in PROT:
        print(f'{ttl:26s} div-rate/day @Hh (WITH) = {np.round(z[f"{proto}_W_rate"],2)}')
        print(f'{"":26s} div-rate/day @Hh (NO)   = {np.round(z[f"{proto}_N_rate"],2)}')
    print(f'Saved {OUT}.png')
