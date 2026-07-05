"""High-resolution parameter sweep of the EZH2->CyclinD1 circuit against two readouts:
  (1) mitogen (SHH) THRESHOLD for the first division (cell-cycle entry), and
  (2) Hh WITHDRAWAL response — divisions, exit-Hh, and cumulative proliferation.

Three knobs:
  - EZH2 mitogen-responsiveness  = kEZE2f          (default 0.022)
  - CyclinD1 inhibition from EZH2 = f0_mk           (default 0.233; inhibition depth = 1 - f0_mk)
  - rate of Hh decline           = withdrawal ramp duration D (fast->slow)

Robust for a long unattended run: every sim wrapped (nan on failure), the cache is CHECKPOINTED after
every row so a crash never loses progress, and progress is logged (flush). Re-run without --fresh to
re-plot from the cache; --mini runs a tiny grid to smoke-test.

Run:  ./venv/bin/python simulations/sweep_param_exploration.py [--fresh] [--mini]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MINI = '--mini' in sys.argv
FRESH = '--fresh' in sys.argv or MINI
CACHE = 'simulations/sweep_param_exploration_cache.npz' if not MINI else 'simulations/sweep_param_exploration_mini.npz'
THR = 2.0                       # nominal CyclinD1 commitment threshold (a.u.)
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)

# ---- grids ----
if MINI:
    KEZ = np.array([0.011, 0.022, 0.044]); F0 = np.array([0.1, 0.35]); DRATES = np.array([9000.0])
else:
    KEZ = np.round(np.geomspace(0.005, 0.060, 50), 5)      # EZH2 responsiveness (default 0.022)
    F0 = np.round(np.linspace(0.02, 0.70, 50), 4)          # inhibition (depth = 1-f0; default f0 0.233)
    DRATES = np.round(np.geomspace(3000.0, 30000.0, 10), 0)  # Hh decline: 50 h .. 500 h


def _rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6); return r


def _set(r, kez, f0):
    r.reset()
    for k, v in GNP.items():
        r[k] = v
    r['kEZE2f'] = float(kez); r['f0_mk'] = float(f0)


def divrate(r, shh, T=8000, settle=3500):
    r['SHH'] = float(shh)
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            res = r.simulate(0, T, int(T / 0.5), selections=['time', 'MPF']); m = res['time'] >= settle; t = res['time'][m]
            pk, _ = find_peaks(res['MPF'][m], prominence=0.15, distance=int(200 / (t[1] - t[0])))
            return len(pk) / ((t[-1] - t[0]) / 60) * 168.0
        except Exception:
            continue
    return np.nan


def entry_threshold(r, kez, f0, lo=0.02, hi=2.2, it=7):
    _set(r, kez, f0)
    if not (divrate(r, hi) > 0.5):
        return np.nan
    _set(r, kez, f0)
    if divrate(r, lo) > 0.5:
        return lo
    for _ in range(it):
        mid = 0.5 * (lo + hi); _set(r, kez, f0)
        if divrate(r, mid) > 0.5:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def movavg(t, x, win=1320):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


def withdrawal(r, kez, f0, D, s0=1.0, s1=0.06, t0=1500.0):
    """ramp SHH s0->s1 over D (exit always occurs by SHH~0.5, so s1=0.06 avoids the deep-arrest stiff
    zone); return (ndiv after t0, exit-SHH via cycle-avg Cd<THR, cumulative divisions)."""
    _set(r, kez, f0); r['SHH'] = s0
    try:
        r.integrator.setValue('maximum_num_steps', 60000)
    except Exception:
        pass
    shh = lambda tt: s0 if tt < t0 else (s1 if tt > t0 + D else s0 + (s1 - s0) * (tt - t0) / D)
    chunk = max(150.0, D / 90.0)           # ~90 chunks regardless of decline rate (bounded overhead)
    T_end = t0 + D + 800
    T, MP, CD, SH = [], [], [], []; t = 0.0; fails = 0
    while t < T_end:
        r['SHH'] = float(shh(t)); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                r.integrator.setValue('absolute_tolerance', atol)
                res = r.simulate(t, t + chunk, max(40, int(chunk / 4)), selections=['time', 'MPF', 'Cd', 'SHH']); ok = True; break
            except Exception:
                pass
        if not ok:
            fails += 1
            if fails > 3:
                break
            t += chunk; continue
        T.append(res['time'][1:]); MP.append(res['MPF'][1:]); CD.append(res['Cd'][1:]); SH.append(res['SHH'][1:])
        t += chunk
    if not T:
        return np.nan, np.nan, np.nan
    T = np.concatenate(T); MP = np.concatenate(MP); CD = np.concatenate(CD); SH = np.concatenate(SH)
    pk, _ = find_peaks(MP, prominence=0.15, distance=int(200 / (T[1] - T[0])))
    ndiv = int(np.sum(T[pk] >= t0))
    cav = movavg(T, CD); above = np.where((cav >= THR) & (T >= t0))[0]
    exit_shh = float(SH[above[-1]]) if len(above) else float(SH[-1])
    return ndiv, exit_shh, ndiv


def save(**kw):
    np.savez(CACHE, KEZ=KEZ, F0=F0, DEPTH=1 - F0, DRATES=DRATES, THR=THR, **kw)


if FRESH or not os.path.exists(CACHE):
    r = _rr()
    nF, nK, nD = len(F0), len(KEZ), len(DRATES)
    ZA = np.full((nF, nK), np.nan)                  # entry threshold
    ZB_nd = np.full((nD, nF, nK), np.nan)           # withdrawal divisions
    ZB_ex = np.full((nD, nF, nK), np.nan)           # withdrawal exit-SHH
    import time as _time
    print(f'grids: KEZ={nK}, F0={nF}, D={nD}  -> entry {nF*nK} thresholds, withdrawal {nD*nF*nK} ramps', flush=True)
    print('=== A: entry threshold ===', flush=True)
    for i, f0 in enumerate(F0):
        for j, kz in enumerate(KEZ):
            ZA[i, j] = entry_threshold(r, kz, f0)
        save(ZA=ZA, ZB_nd=ZB_nd, ZB_ex=ZB_ex)      # checkpoint
        print(f'  [A] f0={f0:.3f}  ({i+1}/{nF})  thresholds={np.round(ZA[i],2)}', flush=True)
    print('=== B: withdrawal (per decline rate) ===', flush=True)
    for di, D in enumerate(DRATES):
        for i, f0 in enumerate(F0):
            for j, kz in enumerate(KEZ):
                nd, ex, _ = withdrawal(r, kz, f0, D)
                ZB_nd[di, i, j] = nd; ZB_ex[di, i, j] = ex
            save(ZA=ZA, ZB_nd=ZB_nd, ZB_ex=ZB_ex)  # checkpoint every row
            print(f'  [B] D={D/60:.0f}h  f0={f0:.3f}  ({di+1}/{nD}, {i+1}/{nF})', flush=True)
    print('cached ->', CACHE, flush=True)
print('DONE. Plot with plot_param_exploration.py', flush=True)
