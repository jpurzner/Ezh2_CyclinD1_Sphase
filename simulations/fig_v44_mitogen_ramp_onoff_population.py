"""GNP mitogen ramp-ON -> sustain -> ramp-OFF as a POPULATION ensemble (the specified sampling
distributions), WITH vs WITHOUT the H3K27me3 mark. Population-median EZH2 & CyclinD1 traces with IQR
bands + population division rate, and the ENTRY / EXIT mitogen-threshold DISTRIBUTIONS across cells.

Sampling (scale-free = default x lognormal, CV only): CyclinD1 abundance (k_Cd_translation) CV 0.70
[data-grounded]; birth p27 (P21_div) CV 0.33; EZH2 abundance (kTlEZ) CV 0.55 [p27/EZH2 CVs placeholder
pending data]. All cells share the SHH staircase; heterogeneity -> a DISTRIBUTION of entry/exit thresholds.

Run:  ./venv/bin/python simulations/fig_v44_mitogen_ramp_onoff_population.py [--n=200] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

GNP = dict(SHH=0.05, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SHH_LO, SHH_HI = 0.05, 0.60
SETTLE, SEG, NPSEG = 4000.0, 300.0, 12
NUP, NSUS, NDN = 30, 20, 40
# defaults (current model) x lognormal(sigma) -> CVs 0.70 / 0.33 / 0.55
KTL0, KTL_SD = 0.801, 0.633
P21_0, P21_SD = 0.6, 0.32
KTLEZ0, KTLEZ_SD = 0.0150989, 0.51
SEL = ['time', 'SHH', 'Cd_mRNA', 'EZH2', 'Dna']
SCHED = ([('up', s) for s in np.linspace(SHH_LO, SHH_HI, NUP)]
         + [('hold', SHH_HI)] * NSUS
         + [('dn', s) for s in np.linspace(SHH_HI, SHH_LO, NDN)])
NSEG = len(SCHED)
# common grid (segment-wise, drop duplicate endpoints): time (h from ramp onset) and SHH per point
_t = []; _shh = []
t0 = 0.0
for _, s in SCHED:
    pts = np.linspace(t0, t0 + SEG, NPSEG)[1:]
    _t.extend(pts); _shh.extend([s] * (NPSEG - 1)); t0 += SEG
TGRID = np.array(_t) / 60.0            # hours
SHHGRID = np.array(_shh)
_RR = None


def _init(mark):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True) if mark
                   else build_model_v44(with_ezh2=True, with_hh=True, params={'f0_mk': 1.0}))
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000); _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass


def _work(task):
    i, ktl, p21, ktz = task
    _RR.reset()
    for k, v in GNP.items():
        _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['P21_div'] = float(p21); _RR['kTlEZ'] = float(ktz)

    def seg(t0, t1, npts):
        for atol, ms in [(1e-9, 20.0), (1e-8, 20.0), (1e-7, 8.0), (1e-6, 3.0), (1e-5, 1.0), (1e-4, 0.5)]:
            try:
                _RR.integrator.setValue('absolute_tolerance', atol); _RR.integrator.setValue('maximum_time_step', ms)
                return _RR.simulate(t0, t1, npts, selections=SEL)
            except Exception:
                continue
        return None
    if seg(0, SETTLE, 4) is None:
        return None
    CD, EZ, DN = [], [], []
    holding = False; lcd = lez = ldn = 0.0
    t = SETTLE
    for _, s in SCHED:
        c = None
        if not holding:
            _RR['SHH'] = float(s); c = seg(t, t + SEG, NPSEG)
            if c is None:
                holding = True                       # stiff deep-arrest -> hold arrested values for the rest
        if c is None:
            CD.append(np.full(NPSEG - 1, lcd)); EZ.append(np.full(NPSEG - 1, lez)); DN.append(np.full(NPSEG - 1, ldn))
        else:
            CD.append(c['Cd_mRNA'][1:]); EZ.append(c['EZH2'][1:]); DN.append(c['Dna'][1:])
            lcd, lez, ldn = c['Cd_mRNA'][-1], c['EZH2'][-1], c['Dna'][-1]
        t += SEG
    CD = np.concatenate(CD); EZ = np.concatenate(EZ); DN = np.concatenate(DN)
    dv = np.where((DN[:-1] > 0.9) & (DN[1:] < 0.1))[0]
    div_t = TGRID[dv]; div_shh = SHHGRID[dv]
    ent = float(div_shh[np.argmin(div_t)]) if len(dv) else np.nan
    ex = float(div_shh[np.argmax(div_t)]) if len(dv) else np.nan
    return dict(cd=CD, ez=EZ, div_t=div_t, ent=ent, ex=ex, ndiv=int(len(dv)))


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='):
            return type(d)(a.split('=', 1)[1])
    return d


CACHE = 'simulations/fig_v44_mitogen_ramp_onoff_population_cache.npz'

if __name__ == '__main__':
    N = _arg('--n', 200); WORKERS = _arg('--workers', max(1, mp.cpu_count() - 1))
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        rng = np.random.default_rng(2026)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, KTL_SD, N)), 0.1, 6.0)
        p21 = np.clip(P21_0 * np.exp(rng.normal(0, P21_SD, N)), 0.1, 3.0)
        ktz = np.clip(KTLEZ0 * np.exp(rng.normal(0, KTLEZ_SD, N)), 0.002, 0.07)
        tasks = [(i, ktl[i], p21[i], ktz[i]) for i in range(N)]
        out = {}
        for tag, mark in [('with', True), ('without', False)]:
            print(f'running {tag} mark: {N} cells...', flush=True)
            with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(mark,)) as pool:
                res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]
            CD = np.vstack([r['cd'] for r in res]); EZ = np.vstack([r['ez'] for r in res])
            out[f'{tag}_cd'] = CD; out[f'{tag}_ez'] = EZ
            out[f'{tag}_ent'] = np.array([r['ent'] for r in res])
            out[f'{tag}_ex'] = np.array([r['ex'] for r in res])
            out[f'{tag}_ndiv'] = np.array([r['ndiv'] for r in res])
            out[f'{tag}_divt'] = np.concatenate([r['div_t'] for r in res]) if res else np.array([])
            out[f'{tag}_ncell'] = len(res)
        np.savez(CACHE, TGRID=TGRID, SHHGRID=SHHGRID, **out)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    TGRID = z['TGRID']; SHHGRID = z['SHHGRID']
    bins = np.linspace(0, TGRID[-1], 60); bw = (bins[1] - bins[0])

    fig = plt.figure(figsize=(14, 12))
    gs = GridSpec(3, 2, height_ratios=[1, 1, 0.85], hspace=0.32, wspace=0.22)

    def trace_panel(axrow, tag, title):
        ax = fig.add_subplot(gs[axrow, :])
        cd = z[f'{tag}_cd']; ez = z[f'{tag}_ez']; nc = int(z[f'{tag}_ncell'])
        for arr, col, lab in [(cd, '#117a65', 'CyclinD1'), (ez, '#8e44ad', 'EZH2')]:
            m = np.median(arr, 0); lo = np.percentile(arr, 25, 0); hi = np.percentile(arr, 75, 0)
            ax.fill_between(TGRID, lo, hi, color=col, alpha=0.18)
            ax.plot(TGRID, m, color=col, lw=2, label=f'{lab} (median ± IQR)')
        # population division rate (per cell per hour), scaled onto left axis
        dt = z[f'{tag}_divt']; h, _ = np.histogram(dt, bins=bins)
        rate = h / nc / bw
        ymax = max(np.percentile(cd, 97), 3) * 1.1
        sc = ymax * 0.42 / max(rate.max(), 1e-6)
        ax.fill_between(0.5 * (bins[:-1] + bins[1:]), 0, rate * sc, step='mid', color='#e67e22', alpha=0.28)
        ax.plot([], [], color='#e67e22', lw=5, alpha=0.4, label=f'pop. division rate (peak {rate.max():.03f}/cell/h)')
        a2 = ax.twinx(); a2.plot(TGRID, SHHGRID, color='#888', lw=1.5, ls=':', label='SHH')
        a2.set_ylim(0, 0.72); a2.set_ylabel('SHH (mitogen)', color='#666')
        ax.set_ylim(0, ymax); ax.set_ylabel('level (a.u.)')
        ax.set_title(f'{title}   (n={nc} cells; CyclinD1 CV0.70 / p27 CV0.33 / EZH2 CV0.55)', fontweight='bold', fontsize=10.5)
        ax.legend(fontsize=8, ncol=3, loc='upper left'); ax.grid(alpha=0.15)
    trace_panel(0, 'with', 'WITH H3K27me3 repression (real model, 25/27)')
    trace_panel(1, 'without', 'WITHOUT H3K27me3 repression (f0_mk=1.0)')

    def dist_panel(col, key, title, xlab):
        ax = fig.add_subplot(gs[2, col])
        for tag, c, lab in [('with', '#117a65', 'WITH mark'), ('without', '#c0392b', 'WITHOUT mark')]:
            v = z[f'{tag}_{key}']; v = v[np.isfinite(v)]
            ax.hist(v, bins=np.linspace(0, 0.62, 26), color=c, alpha=0.5, label=f'{lab} (median {np.median(v):.2f})', density=True)
            ax.axvline(np.median(v), color=c, ls='--', lw=1.5)
        ax.set_xlabel(xlab); ax.set_ylabel('density'); ax.set_title(title, fontweight='bold', fontsize=10.5)
        ax.legend(fontsize=8); ax.grid(alpha=0.15)
    dist_panel(0, 'ent', '(C) ENTRY threshold distribution', 'SHH at first division (entry)')
    dist_panel(1, 'ex', '(D) EXIT threshold distribution', 'SHH at last division (exit)')

    fig.suptitle('GNP mitogen ramp-ON → sustain → ramp-OFF — POPULATION ensemble (specified sampling), WITH vs WITHOUT H3K27me3 repression',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.savefig('simulations/fig_v44_mitogen_ramp_onoff_population.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/fig_v44_mitogen_ramp_onoff_population.pdf', bbox_inches='tight')
    print('Saved fig_v44_mitogen_ramp_onoff_population.png')
    for tag in ('with', 'without'):
        e = z[f'{tag}_ent']; x = z[f'{tag}_ex']; nd = z[f'{tag}_ndiv']
        print(f'  {tag:8s}: entry SHH {np.nanmedian(e):.3f} [{np.nanpercentile(e,25):.2f}-{np.nanpercentile(e,75):.2f}]  '
              f'exit {np.nanmedian(x):.3f}  divisions {np.median(nd):.0f} [{np.percentile(nd,10):.0f}-{np.percentile(nd,90):.0f}]  '
              f'never-divide {100*np.mean(nd==0):.0f}%')
