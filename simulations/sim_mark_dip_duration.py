"""How do WITH vs WITHOUT the H3K27me3 mark respond to a transient Hh DIP of varying DURATION?
Instead of a sine, Hh sits at a cycling baseline, DIPS to a low level for a duration D, then RETURNS —
and we measure the quiescence response and the RECOVERY after the dip. This probes the mark's MEMORY: a
slow integrator should let cells ride through a brief dip but convert a long dip into lasting arrest.

Protocol: settle cycling at SHH=0.8 (random phase) → dip to SHH=0.15 for duration D → return to 0.8 →
250 h recovery. Ensemble: calibrated abundance dists + mu (growth) dist + partition noise. Divisions =
Dna 1→0 reset; per-time state (run-based) = cycling / transient-G0 (p27-high dwell that re-divides) /
arrest (never re-divides in the window). Dip durations swept 4–250 h.

Readouts vs dip duration, WITH vs WITHOUT: (A) time-domain at a representative dip — fraction cycling /
transient-G0 / arrest vs time (the response & recovery). (B) fraction ARRESTED (never recovers) — the
critical dip duration for irreversible exit. (C) RECOVERY time (time to re-divide after the dip) — the
mark's memory. (D) G0 occupancy DURING the dip — how fast/deep the quiescence response is.

Run:  ./venv/bin/python simulations/sim_mark_dip_duration.py [--n=400] [--workers=9] [--fresh]
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

SEL = ['time', 'Dna', 'P21', 'aRc']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH_HI, HH_LO = 0.80, 0.15
D_DIP_H = [4.0, 10.0, 24.0, 48.0, 96.0, 160.0, 250.0]    # dip durations (hours)
D_TRACE = 48.0                                           # representative dip for the time-domain panel
SETTLE, RECOVERY, CHUNK, NP = 5000.0, 15000.0, 60.0, 9   # 250 h recovery
DT = CHUNK / (NP - 1)
DWELL_CUT, P27_HI, PART = 2.0, 0.1, 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, MU_SDLOG = 0.633, 0.32, 0.51, 0.22
T_TRACE = np.arange(0.0, D_TRACE * 60 + RECOVERY + 1e-6, 12.0)   # fixed grid for the representative-dip trace (min)

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, Ddip_h, f0, ktl, p21d, ktlez, mu, wt = task
    Ddip = Ddip_h * 60.0
    rng = np.random.default_rng((hash((i, round(Ddip_h, 1), f0)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = HH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu); _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    NULLTRACE = None
    NULL = (Ddip_h, f0, 0, np.nan, np.nan, NULLTRACE)
    # settle cycling at baseline (random phase)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return NULL
        tm += CHUNK
    hh = lambda tt: HH_LO if tt < Ddip else HH_HI
    T, DNA, P21, ARC = [], [], [], []
    t = 0.0; TREC = Ddip + RECOVERY
    while t < TREC:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(hh(t)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        T.append(r['time'][1:] - tm + t); DNA.append(r['Dna'][1:]); P21.append(r['P21'][1:]); ARC.append(r['aRc'][1:])
        t += CHUNK; tm += CHUNK
    if not T:
        return NULL
    T = np.concatenate(T); DNA = np.concatenate(DNA); P21 = np.concatenate(P21); ARC = np.concatenate(ARC)  # T in post-settle min
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; g0inst = (~inS & ~inG2) & (P21 > 0.1)
    state = np.zeros(len(T), np.int8)
    dd = np.diff(np.concatenate([[0], g0inst.astype(np.int8), [0]]))
    for s, e in zip(np.where(dd == 1)[0], np.where(dd == -1)[0]):
        ee = min(e, len(T) - 1)
        if (ee - s) * DT < DWELL_CUT * 60: continue
        if np.any(div > ee): state[s:e] = 1
        else: state[s:] = 2; break
    dtimes = T[div]
    after = dtimes[dtimes > Ddip]
    recovered = int(len(after) > 0)
    rec_time = float((after.min() - Ddip) / 60.0) if recovered else np.nan
    dipmask = T < Ddip
    g0occ = float(np.mean(g0inst[dipmask])) if dipmask.any() else np.nan
    trace = np.clip(np.round(np.interp(T_TRACE, T, state)), 0, 2).astype(np.int8) if wt else None
    return (Ddip_h, f0, recovered, rec_time, g0occ, trace)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 400)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_mark_dip_duration_cache.npz'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        tasks = [(i, Dh, f0, ktl[i], p21[i], ktz[i], mu[i], Dh == D_TRACE) for f0 in (0.233, 1.0) for Dh in D_DIP_H for i in range(N)]
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = []
            for k, r in enumerate(pool.imap_unordered(_work, tasks, chunksize=6)):
                res.append(r)
                if (k + 1) % 400 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        store = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            fa = []; rt = []; g0 = []
            for Dh in D_DIP_H:
                sub = [r for r in res if r[0] == Dh and r[1] == f0]
                rec = np.array([r[2] for r in sub]); rts = np.array([r[3] for r in sub]); g0o = np.array([r[4] for r in sub])
                fa.append(1.0 - np.mean(rec) if len(rec) else np.nan)
                rt.append(np.nanmedian(rts[rec == 1]) if np.any(rec == 1) else np.nan)
                g0.append(np.nanmean(g0o))
            store[f'{tag}_arr'] = np.array(fa); store[f'{tag}_rec'] = np.array(rt); store[f'{tag}_g0'] = np.array(g0)
            traces = np.array([r[5] for r in res if r[1] == f0 and r[5] is not None])
            store[f'{tag}_cyc_t'] = (traces == 0).mean(0); store[f'{tag}_tg0_t'] = (traces == 1).mean(0); store[f'{tag}_arr_t'] = (traces == 2).mean(0)
        np.savez(CACHE, D=np.array(D_DIP_H), Ttr=T_TRACE / 60.0, D_TRACE=D_TRACE, N=N, HH_HI=HH_HI, HH_LO=HH_LO, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    D = z['D']; Ttr = z['Ttr']; Nc = int(z['N']); DTR = float(z['D_TRACE'])
    COND = [('W', '#8b1a1a', 'WITH mark (repressed)'), ('N', '#3b7bbf', 'WITHOUT mark')]

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    # (A) time-domain at the representative dip
    a = ax[0, 0]
    a.axvspan(0, DTR, color='#cccccc', alpha=0.35, label=f'Hh dip ({DTR:.0f} h)')
    for tag, col, nm in COND:
        a.plot(Ttr, z[f'{tag}_cyc_t'], '-', color=col, lw=2.6, label=f'{nm} — cycling')
        a.plot(Ttr, z[f'{tag}_arr_t'], ':', color=col, lw=2.0, alpha=0.8, label=f'{nm} — arrested')
    a.set_xlabel('time (h)  (dip starts at 0)'); a.set_ylabel('fraction of cells'); a.set_xlim(0, min(Ttr.max(), DTR + 200))
    a.set_title(f'(A) Time-domain: a {DTR:.0f} h Hh dip — cycling (solid) & arrest (dotted)\nWITH mark recovers slower / arrests more after the dip', fontweight='bold', fontsize=10.5); a.legend(fontsize=7.5, loc='center right'); a.grid(alpha=0.15)
    # (B) fraction arrested (never recovers) vs dip duration
    a = ax[0, 1]
    for tag, col, nm in COND: a.plot(D, z[f'{tag}_arr'], '-o', color=col, lw=2.7, ms=6, label=nm)
    a.set_xscale('log'); a.set_xlabel('dip duration (h)'); a.set_ylabel('fraction ARRESTED (never re-divides)')
    a.set_title('(B) How long a dip causes IRREVERSIBLE exit?\nWITH mark = arrest at SHORTER dips (lower threshold)', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both'); a.set_ylim(-0.02, 1.02)
    # (C) recovery time after the dip vs dip duration
    a = ax[1, 0]
    for tag, col, nm in COND: a.plot(D, z[f'{tag}_rec'], '-o', color=col, lw=2.7, ms=6, label=nm)
    a.set_xscale('log'); a.set_xlabel('dip duration (h)'); a.set_ylabel('recovery time — dip end → first division (h)')
    a.set_title('(C) Among cells that RECOVER: recovery is slower WITH the mark at short/mid dips (memory);\nboth fall at long dips = survivorship (only the fast/robust cells still recover)', fontweight='bold', fontsize=9.5); a.legend(fontsize=9); a.grid(alpha=0.15, which='both')
    # (D) G0 occupancy during the dip vs dip duration
    a = ax[1, 1]
    for tag, col, nm in COND: a.plot(D, z[f'{tag}_g0'], '-o', color=col, lw=2.7, ms=6, label=nm)
    a.set_xscale('log'); a.set_xlabel('dip duration (h)'); a.set_ylabel('p27-high G0 occupancy DURING the dip')
    a.set_title('(D) How deep/fast is the quiescence response to the dip itself?', fontweight='bold', fontsize=11); a.legend(fontsize=9); a.grid(alpha=0.15, which='both')

    fig.suptitle(f'Response to an Hh DIP of varying DURATION — WITH vs WITHOUT H3K27me3 repression (baseline {HH_HI}, dip {HH_LO}, calibrated+mu+partition, N={Nc}/cond)',
                 fontsize=12, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_dip_duration.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_dip_duration.pdf', bbox_inches='tight')
    plt.close()
    print('dip durations (h)      =', list(D))
    for tag, nm in [('W', 'WITH   '), ('N', 'WITHOUT')]:
        print(f'{nm} frac arrested   =', np.round(z[f'{tag}_arr'], 3))
        print(f'{nm} recovery time h =', np.round(z[f'{tag}_rec'], 1))
    print('Saved sim_mark_dip_duration.png')
