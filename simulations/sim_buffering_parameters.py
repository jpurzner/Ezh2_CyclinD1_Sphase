"""Which parameters support BUFFERING? Systematic one-at-a-time sensitivity of the EZH2/H3K27me3 buffering of
CyclinD1/proliferation, on the current v44 default (GNP). (JP 2026-07-20.)

Three buffering metrics (GNP):
  shift      = mitogen threshold ON / OFF  (STATIC buffering: how much more mitogen the repression demands)
  gain       = d(log CyclinD1)/d(log Shh) over the proliferative band, ON  (dose-response flattening; LOWER = buffered)
  fluc_atten = mark CV / Shh CV under an Shh oscillation (period 33 h)  (DYNAMIC low-pass; LOWER = buffered)
plus  mark_slope = Mk(hi Shh) - Mk(lo Shh)  (>0 = mark peaks at proliferation, ChIP-seq-correct)

Sweeps the buffering-relevant parameters one at a time from the baked default:
  f0_prc2 (repression floor), K_prc2 (repression sensitivity), n_prc2 (Hill), a_rw_prc2 (mark read-write/level),
  del_mk (mark turnover/timescale), kEZbas_Cd (mitogen-dose EZH2), kEZE2f (cycle EZH2), kDeEZ (EZH2 turnover).
thr_OFF is mark-independent (f0=1 -> no repression) -> computed once.

Run: ./venv/bin/python simulations/sim_buffering_parameters.py [--quick]
"""
import os, sys, pickle, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, multiprocessing as mp, tellurium as te
from src.build_model_v44_heldt import build_model_v44

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
SEL = ['time', 'Dna', 'Cd', 'EZH2', 'Mk']
BASE = dict(f0_prc2=0.174826, K_prc2=0.00685695, n_prc2=4.01574, a_rw_prc2=0.00226571,
            del_mk=0.000440007, kEZbas_Cd=0.000318437, kEZE2f=0.0122639, kDeEZ=0.000276622)
SHH = np.round(np.geomspace(0.03, 1.3, 14), 4)
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6); _RR.integrator.setValue('maximum_num_steps', 2000000)


def _steady(over, shh, settle=6000, run=13000):
    _RR.reset()
    for k, v in GNP.items(): _RR[k] = v
    for k, v in over.items(): _RR[k] = v
    _RR['SHH'] = shh
    r = None
    for atol in (1e-8, 1e-7, 1e-6, 1e-5):
        try: _RR.integrator.setValue('absolute_tolerance', atol); r = _RR.simulate(0, run, int(run / 6), selections=SEL); break
        except Exception:
            try: _RR.integrator.setValue('maximum_time_step', 5.0)
            except Exception: pass
    if r is None: return (0.0, np.nan, np.nan)
    t = r['time']; m = t >= settle; tt = t[m]; dna = np.asarray(r['Dna'])[m]
    dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]; span = (tt[-1] - tt[0]) / UPH
    return (len(dv) / span if span > 0 else 0.0, np.asarray(r['Cd'])[m].mean(), np.asarray(r['Mk'])[m].mean())


def _threshold(div):
    div = np.asarray(div)
    if div.max() < 1e-4: return np.nan
    lvl = 0.5 * div.max(); above = div >= lvl
    if not above.any(): return np.nan
    i = int(np.argmax(above))
    if i == 0: return float(SHH[0])
    return float(SHH[i-1] + (lvl - div[i-1]) * (SHH[i] - SHH[i-1]) / (div[i] - div[i-1] + 1e-9))


def _fluct(over, T_h=33.0, amp=0.45, shh0=0.35, n_per=6, chunk=6.0):
    _RR.reset()
    for k, v in GNP.items(): _RR[k] = v
    for k, v in over.items(): _RR[k] = v
    _RR['SHH'] = shh0
    for atol in (1e-8, 1e-7, 1e-6):
        try: _RR.integrator.setValue('absolute_tolerance', atol); _RR.simulate(0, 6000, 200, selections=['time']); break
        except Exception: continue
    T = T_h * UPH; fn = lambda t: float(np.clip(shh0 + amp * np.sin(2 * np.pi * (t - 6000) / T), 0.02, 1.1))
    tt, MK, SH = [], [], []
    t = 6000.0; end = 6000 + n_per * T
    while t < end - 1e-6:
        _RR['SHH'] = fn(t); r = None
        for atol in (1e-8, 1e-7, 1e-6, 1e-5):
            try: _RR.integrator.setValue('absolute_tolerance', atol); r = _RR.simulate(t, t + chunk, 3, selections=SEL); break
            except Exception:
                try: _RR.integrator.setValue('maximum_time_step', 3.0)
                except Exception: pass
        if r is None: t += chunk; continue
        MK.append(np.asarray(r['Mk'])[1:].mean()); SH.append(fn(t + chunk / 2)); tt.append(t); t += chunk
    MK = np.array(MK); SH = np.array(SH); w = np.arange(len(MK)) > int(T / chunk)
    mk_cv = MK[w].std() / (MK[w].mean() + 1e-9); sh_cv = SH[w].std() / (SH[w].mean() + 1e-9)
    return mk_cv / (sh_cv + 1e-9)


def _metrics(over):
    ds = [_steady(over, s) for s in SHH]
    div = [d[0] for d in ds]; cd = np.array([d[1] for d in ds]); mk = np.array([d[2] for d in ds])
    band = (SHH >= 0.15) & (SHH <= 0.9)
    gain = float(np.polyfit(np.log(SHH[band]), np.log(cd[band] + 1e-9), 1)[0]) if band.sum() >= 2 else np.nan
    return dict(thr=_threshold(div), gain=gain, mark_slope=float(mk[-3] - mk[2]),
                fluc=_fluct(over), mk_lo=float(mk[2]), mk_hi=float(mk[-3]),
                cd=cd, mk=mk)                         # store full curves -> compute smooth R = Cd_ON/Cd_OFF in analysis


def _task(a):
    pname, val = a
    over = {pname: val} if pname != 'OFF' else {'f0_prc2': 1.0}
    return (pname, val, _metrics(over))


if __name__ == '__main__':
    quick = '--quick' in sys.argv
    # sweep grids (one-at-a-time from BASE)
    def mult(base, facs): return [base * f for f in facs]
    FACS = [0.25, 0.5, 1.0, 2.0, 4.0] if not quick else [0.5, 1.0, 2.0]
    SWEEP = {
        'f0_prc2':   [0.05, 0.10, 0.174826, 0.30, 0.50, 0.80],
        'K_prc2':    [0.0005, 0.001, 0.002, 0.00685695, 0.012, 0.02],
        'n_prc2':    [1.0, 2.0, 3.0, 4.01574, 6.0],
        'a_rw_prc2': mult(BASE['a_rw_prc2'], FACS),
        'del_mk':    mult(BASE['del_mk'], FACS),
        'kEZbas_Cd': [0.0] + mult(BASE['kEZbas_Cd'], [0.5, 1.0, 2.0, 4.0]),
        'kEZE2f':    mult(BASE['kEZE2f'], [0.5, 1.0, 2.0, 4.0]),
        'kDeEZ':     mult(BASE['kDeEZ'], FACS),
    }
    if quick: SWEEP = {'f0_prc2': [0.10, 0.174826, 0.50], 'K_prc2': [0.002, 0.00685695], 'del_mk': mult(BASE['del_mk'], FACS)}
    tasks = [('OFF', 0.0)] + [(p, v) for p, vs in SWEEP.items() for v in vs]
    print(f"Buffering sensitivity: {len(SWEEP)} params, {len(tasks)-1} param-values + OFF ref.", flush=True)
    t0 = time.time(); out = {'BASE': BASE, 'SHH': SHH}
    with mp.Pool(6, initializer=_init, initargs=(None,)) as pool:
        done = 0
        for pname, val, m in pool.imap_unordered(_task, tasks, chunksize=1):
            out.setdefault(pname, []).append((val, m)); done += 1
            if done % 10 == 0: print(f"  {done}/{len(tasks)} ({(time.time()-t0)/60:.1f} min)", flush=True)
    thr_off = out['OFF'][0][1]['thr']; out['thr_off'] = thr_off
    cd_off = out['OFF'][0][1]['cd']; out['cd_off'] = cd_off
    si = int(np.argmin(np.abs(SHH - 0.5)))                # proliferative reference mitogen for R
    # summary. R = Cd_ON/Cd_OFF (repression depth, SMOOTH); buffering = 1-R; static shift ~ 1/R
    print(f"\n  thr_OFF={thr_off:.3f}   Cd_OFF@Shh{SHH[si]:.2f}={cd_off[si]:.2f}   (R=Cd_ON/Cd_OFF; buffering=1-R)\n", flush=True)
    print(f"  {'param':>11} {'value':>10} | {'R':>6} {'buffer':>7} {'1/R':>6} | {'gain':>6} {'fluc_atten':>10} {'mark_slope':>10}", flush=True)
    for pname, vs in SWEEP.items():
        for val, m in sorted(out[pname]):
            R = m['cd'][si] / (cd_off[si] + 1e-9)
            print(f"  {pname:>11} {val:>10.5g} | {R:>6.3f} {1-R:>7.3f} {1/(R+1e-9):>6.2f} | {m['gain']:>6.3f} {m['fluc']:>10.3f} {m['mark_slope']:>+10.3f}", flush=True)
    with open(os.path.join(HERE, 'sim_buffering_parameters_cache.pkl'), 'wb') as fh:
        pickle.dump(out, fh)
    print(f"\n==== DONE in {(time.time()-t0)/60:.1f} min -> sim_buffering_parameters_cache.pkl ====", flush=True)
