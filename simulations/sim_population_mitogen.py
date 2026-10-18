"""MANY-CELL population simulation: each cell samples its variable parameters from the calibrated distributions,
run over multiple divisions across a mitogen dose-response, GNP & MB, EZH2->CyclinD1 repression ON vs OFF. (JP 2026-07-20.)

Per-cell distributions (calibrated, scale-free -- same as sim_population_g0 / sim_g0_bifurcation):
  k_Cd_translation = 0.801  * LogN(0, 0.70)   -- CyclinD1 abundance, CV ~0.70 (data)   [the CyclinD1-setpoint axis]
  kTlEZ            = 0.0151 * LogN(0, 0.51)   -- EZH2 abundance,     CV ~0.55 (data)
  P21_div          = LogN(median 0.45, sig 0.55) -- INHERITED birth-p27 (shared by sisters), the G0-duration axis
  mu               = 0.0005 * LogN(0, 0.15)   -- modest growth-rate spread (cell-size heterogeneity)

Population readouts vs mitogen: FRACTION proliferating (graded threshold from the heterogeneity), and the median +
spread of division rate / CyclinD1 / EZH2 / H3K27me3. The single-cell sharp threshold (sim_ezh2_repression_mitogen)
becomes a GRADED population dose-response; the repression ON vs OFF shifts the whole distribution of thresholds.

Parallel over cells (mp.Pool, model built once per worker; per-cell params set at runtime). Run:
  ./venv/bin/python simulations/sim_population_mitogen.py [--n=120] [--quick]
"""
import os, sys, pickle, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, multiprocessing as mp, tellurium as te
from src.build_model_v44_heldt import build_model_v44

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
KTL0, KTLEZ0, MU0 = 0.801, 0.0150988626303029, 0.0005
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.006)
SEL = ['time', 'Dna', 'Cd', 'EZH2', 'Mk', 'Skp2']
_RR = None


def _arg(flag, d):
    for a in sys.argv:
        if a.startswith(flag + '='): return type(d)(a.split('=', 1)[1])
    return d


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6); _RR.integrator.setValue('maximum_num_steps', 2000000)


def _one(cond, f0, cellp, shh, settle=6000, run=14000):
    _RR.reset()
    for k, v in cond.items(): _RR[k] = v
    for k, v in cellp.items(): _RR[k] = v
    _RR['f0_prc2'] = f0; _RR['SHH'] = shh
    r = None
    for atol in (1e-8, 1e-7, 1e-6, 1e-5):
        try:
            _RR.integrator.setValue('absolute_tolerance', atol)
            r = _RR.simulate(0, run, int(run / 6), selections=SEL); break
        except Exception:
            try: _RR.integrator.setValue('maximum_time_step', 5.0)
            except Exception: pass
    if r is None: return (0.0, np.nan, np.nan, np.nan)
    t = r['time']; m = t >= settle; tt = t[m]; dna = np.asarray(r['Dna'])[m]
    dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    span = (tt[-1] - tt[0]) / UPH
    return (len(dv) / span if span > 0 else 0.0, np.asarray(r['Cd'])[m].mean(),
            np.asarray(r['EZH2'])[m].mean(), np.asarray(r['Mk'])[m].mean())


def _task(a):
    ct_name, cond, f0, ci, cellp, shh_list = a
    out = [_one(cond, f0, cellp, s) for s in shh_list]
    return (ct_name, f0, ci, np.array(out))    # out[:, (div_h, Cd, EZH2, Mk)]


if __name__ == '__main__':
    N = _arg('--n', 40 if '--quick' in sys.argv else 120)
    NSHH = 8 if '--quick' in sys.argv else 14
    SHH = np.round(np.geomspace(0.03, 1.3, NSHH), 4)
    rng = np.random.default_rng(2026)
    # per-cell draws (paired: a cell keeps its identity across conditions)
    cd_scale = np.clip(np.exp(rng.normal(0, 0.70, N)), 0.2, 5.0)
    ez_scale = np.clip(np.exp(rng.normal(0, 0.51, N)), 0.3, 3.0)
    p21 = np.clip(np.exp(rng.normal(np.log(0.45), 0.55, N)), 0.1, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, 0.15, N)), 3e-4, 8e-4)
    cellps = [dict(k_Cd_translation=KTL0 * cd_scale[i], kTlEZ=KTLEZ0 * ez_scale[i], P21_div=float(p21[i]), mu=float(mu[i])) for i in range(N)]
    print(f"Population: {N} cells x {NSHH} Shh x 2 celltypes x 2 repression (ON f0=0.175 / OFF f0=1.0). Multi-division.", flush=True)
    print(f"  draws: CyclinD1 CV {cd_scale.std()/cd_scale.mean():.2f}, EZH2 CV {ez_scale.std()/ez_scale.mean():.2f}, "
          f"birth-p27 median {np.median(p21):.2f} CV {p21.std()/p21.mean():.2f}", flush=True)

    tasks = []
    for ct_name, cond in (('GNP', GNP), ('MB', MB)):
        for f0 in (0.174825546661583, 1.0):
            for ci in range(N):
                tasks.append((ct_name, cond, f0, ci, cellps[ci], SHH))
    t0 = time.time()
    res = {}
    with mp.Pool(6, initializer=_init, initargs=(None,)) as pool:
        done = 0
        for ct_name, f0, ci, arr in pool.imap_unordered(_task, tasks, chunksize=2):
            res.setdefault((ct_name, f0), {})[ci] = arr
            done += 1
            if done % 60 == 0: print(f"  {done}/{len(tasks)} cell-sweeps ({(time.time()-t0)/60:.1f} min)", flush=True)

    # assemble arrays: for each (ct, f0) -> array [N cells, NSHH, 4 (div_h, Cd, EZH2, Mk)]
    out = dict(shh=SHH, cd_scale=cd_scale, ez_scale=ez_scale, p21=p21)
    for (ct, f0), d in res.items():
        A = np.stack([d[i] for i in range(N)])          # [N, NSHH, 4]
        tag = f"{ct}_{'ON' if f0 < 0.99 else 'OFF'}"
        out[tag] = A
    with open(os.path.join(HERE, 'sim_population_mitogen_cache.pkl'), 'wb') as fh:
        pickle.dump(out, fh)

    # ---- summary: population proliferating fraction + median readouts; population threshold (50% proliferating) ----
    def pop_thr(A, frac=0.5):
        prolif = (A[:, :, 0] * 24 > 0.3)                 # cell "proliferating" if >0.3 div/day
        f = prolif.mean(axis=0)                          # fraction proliferating vs Shh
        if f.max() < frac: return np.nan
        i = np.argmax(f >= frac)
        if i == 0: return float(SHH[0])
        x0, x1, y0, y1 = SHH[i-1], SHH[i], f[i-1], f[i]
        return float(x0 + (frac - y0) * (x1 - x0) / (y1 - y0 + 1e-9))
    print(f"\n{'condition':>10} | population mitogen threshold (50% proliferating)", flush=True)
    for ct in ('GNP', 'MB'):
        on, off = pop_thr(out[f'{ct}_ON']), pop_thr(out[f'{ct}_OFF'])
        print(f"  {ct:>8} | ON {on:.3f}   OFF {off:.3f}   -> repression shifts pop. threshold {on/off if off else float('nan'):.2f}x", flush=True)
    print(f"\n==== DONE in {(time.time()-t0)/60:.1f} min -> sim_population_mitogen_cache.pkl ; figure next ====", flush=True)
