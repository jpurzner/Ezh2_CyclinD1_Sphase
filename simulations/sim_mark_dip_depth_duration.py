"""Hh DIP: 2-D sweep of DURATION x DEPTH — how much Hedgehog loss, for how long, commits a cell to
permanent exit, WITH vs WITHOUT the H3K27me3 mark. Longer baseline & recovery on both ends.

Protocol: settle cycling at SHH=0.8 for ~150 h (random phase) → dip to SHH=Hh_lo for duration D → return
to 0.8 → ~360 h recovery. Ensemble (calibrated abundance + mu + partition noise). Division = Dna 1→0 reset;
a cell "ARRESTS" if it never re-divides in the (long) recovery window. Sweep dip duration D and dip DEPTH
(the level Hh drops to) as a 2-D grid.

Readout: fraction ARRESTED over (duration x depth), WITH vs WITHOUT the mark, and the mark-attributable
difference. The mark accumulates during the low-Hh dip (no divisions dilute it), so a deeper/longer dip is
remembered and a transient loss becomes lasting quiescence.

Panels: (A) arrest heatmap WITH mark, (B) WITHOUT, (C) WITH−WITHOUT (mark-attributable arrest),
(D) arrest vs duration at a few depths (WITH solid / WITHOUT dashed).

Run:  ./venv/bin/python simulations/sim_mark_dip_depth_duration.py [--n=200] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'Dna']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH_HI = 0.80
D_DIP_H = [12.0, 25.0, 50.0, 100.0, 200.0]              # dip durations (hours)
HH_LO_LIST = [0.05, 0.15, 0.25, 0.35, 0.45]             # dip DEPTHS (level Hh drops to; lower = deeper)
SETTLE, RECOVERY, CHUNK, NP = 9000.0, 21600.0, 60.0, 9  # ~150 h baseline, ~360 h recovery
DWELL_CUT = 2.0
PART = 0.30
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, EZ_SDLOG, MU_SDLOG = 0.633, 0.32, 0.51, 0.22

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, Ddip_h, hh_lo, f0, ktl, p21d, ktlez, mu = task
    Ddip = Ddip_h * 60.0
    rng = np.random.default_rng((hash((i, round(Ddip_h, 1), round(hh_lo, 3), f0)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = HH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez); _RR['mu'] = float(mu); _RR['f0_mk'] = f0; _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    NULL = (Ddip_h, hh_lo, f0, 0, np.nan)
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return NULL
        tm += CHUNK
    hh = lambda tt: hh_lo if tt < Ddip else HH_HI
    T, DNA = [], []
    t = 0.0; TREC = Ddip + RECOVERY
    while t < TREC:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(hh(t)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        T.append(r['time'][1:] - tm + t); DNA.append(r['Dna'][1:]); t += CHUNK; tm += CHUNK
    if not T:
        return NULL
    T = np.concatenate(T); DNA = np.concatenate(DNA)
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]
    after = T[div][T[div] > Ddip]
    recovered = int(len(after) > 0)
    rec_time = float((after.min() - Ddip) / 60.0) if recovered else np.nan
    return (Ddip_h, hh_lo, f0, recovered, rec_time)


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 200)); WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_mark_dip_depth_duration_cache.npz'
    nD, nH = len(D_DIP_H), len(HH_LO_LIST)

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        tasks = [(i, Dh, hl, f0, ktl[i], p21[i], ktz[i], mu[i]) for f0 in (0.233, 1.0) for Dh in D_DIP_H for hl in HH_LO_LIST for i in range(N)]
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            res = []
            for k, r in enumerate(pool.imap_unordered(_work, tasks, chunksize=4)):
                res.append(r)
                if (k + 1) % 500 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        store = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            arr = np.full((nH, nD), np.nan); rec = np.full((nH, nD), np.nan)
            for di, Dh in enumerate(D_DIP_H):
                for hi, hl in enumerate(HH_LO_LIST):
                    sub = [r for r in res if r[0] == Dh and abs(r[1] - hl) < 1e-6 and r[2] == f0]
                    rc = np.array([r[3] for r in sub]); rt = np.array([r[4] for r in sub])
                    arr[hi, di] = 1.0 - np.mean(rc) if len(rc) else np.nan
                    rec[hi, di] = np.nanmedian(rt[rc == 1]) if np.any(rc == 1) else np.nan
            store[f'{tag}_arr'] = arr; store[f'{tag}_rec'] = rec
        np.savez(CACHE, D=np.array(D_DIP_H), HL=np.array(HH_LO_LIST), N=N, HH_HI=HH_HI, **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    D = z['D']; HL = z['HL']; Nc = int(z['N'])

    def _edges(c, log=False):
        c = np.asarray(c, float); lc = np.log10(c) if log else c
        e = np.empty(len(c) + 1); e[1:-1] = 0.5 * (lc[:-1] + lc[1:]); e[0] = lc[0] - 0.5 * (lc[1] - lc[0]); e[-1] = lc[-1] + 0.5 * (lc[-1] - lc[-2])
        return 10 ** e if log else e
    XE = _edges(D, log=True); YE = _edges(HL)

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))

    def heat(a, M, ttl, cmap='magma', vmin=0, vmax=1, diff=False):
        pm = a.pcolormesh(XE, YE, M, cmap=cmap, shading='flat', vmin=vmin, vmax=vmax)
        a.set_xscale('log'); a.set_xlabel('dip duration (h)'); a.set_ylabel('dip depth — Hh floor (lower = deeper)')
        a.set_title(ttl, fontweight='bold', fontsize=11)
        for hi in range(len(HL)):
            for di in range(len(D)):
                v = M[hi, di]
                if np.isfinite(v): a.text(D[di], HL[hi], f'{v:+.2f}' if diff else f'{v:.2f}', ha='center', va='center', fontsize=7.5, color='white' if (diff and abs(v) > 0.25) or (not diff and v > 0.45) else 'black')
        fig.colorbar(pm, ax=a, label='fraction arrested' if not diff else 'Δ arrest (WITH − WITHOUT)')

    heat(ax[0, 0], z['W_arr'], '(A) Fraction ARRESTED — WITH mark (repressed)')
    heat(ax[0, 1], z['N_arr'], '(B) Fraction ARRESTED — WITHOUT mark')
    dcmap = LinearSegmentedColormap.from_list('d', ['#3b7bbf', '#f7f7f7', '#8b1a1a'])
    heat(ax[1, 0], z['W_arr'] - z['N_arr'], '(C) Mark-attributable arrest (WITH − WITHOUT)\nred = the mark makes this dip commit cells to exit', cmap=dcmap, vmin=-0.6, vmax=0.6, diff=True)

    a = ax[1, 1]
    depth_idx = [0, len(HL) // 2, len(HL) - 1]
    cols = ['#4a1486', '#c47f17', '#1f8a4c']
    for j, hi in enumerate(depth_idx):
        a.plot(D, z['W_arr'][hi], '-o', color=cols[j], lw=2.6, ms=6, label=f'Hh floor {HL[hi]:.2f} — WITH')
        a.plot(D, z['N_arr'][hi], '--s', color=cols[j], lw=1.8, ms=4, alpha=0.7, label=f'Hh floor {HL[hi]:.2f} — WITHOUT')
    a.set_xscale('log'); a.set_xlabel('dip duration (h)'); a.set_ylabel('fraction arrested')
    a.set_title('(D) Arrest vs dip duration at 3 depths — WITH (solid) vs WITHOUT (dashed)', fontweight='bold', fontsize=11); a.legend(fontsize=7.5, loc='upper left'); a.grid(alpha=0.15, which='both'); a.set_ylim(-0.02, 1.02)

    fig.suptitle(f'Hh DIP duration × depth — arrest WITH vs WITHOUT H3K27me3 (baseline {HH_HI}, ~150 h before / ~360 h after, calibrated+mu+partition, N={Nc}/cell)',
                 fontsize=12, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_mark_dip_depth_duration.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_mark_dip_depth_duration.pdf', bbox_inches='tight')
    plt.close()
    print('dip durations (h):', list(D), ' depths (Hh floor):', list(HL))
    for tag, nm in [('W', 'WITH   '), ('N', 'WITHOUT')]:
        print(f'{nm} arrest grid (rows=depth {list(HL)}, cols=duration):')
        print(np.round(z[f'{tag}_arr'], 2))
    print('Saved sim_mark_dip_depth_duration.png')
