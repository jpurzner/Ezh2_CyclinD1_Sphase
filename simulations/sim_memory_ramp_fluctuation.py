"""RE-EXPLORE the memory function on the CURRENT v44 default (Cdh1-EZH2 + graded G1-lengthening + Skp2), with
mitogen RAMP up/down (hysteresis) and FLUCTUATION drives. (JP 2026-07-20.)

Tracks the memory variables through a time-varying Shh: CyclinD1 (fast), EZH2 (baseline+S-phase, slow), H3K27me3
mark Mk, Skp2, and proliferation (divisions). Two experiments:
  A. RAMP Shh up then down at fast/slow rates -> hysteresis LOOPS (up-branch != down-branch = memory). Compares the
     full model to a mark-OFF control (f0_prc2=1: mark exerts no CyclinD1 repression) to isolate the mark's share.
  B. FLUCTUATION: Shh = Shh0 + A*sin(2 pi t/T) at fast/medium/slow periods -> low-pass BUFFERING (CyclinD1 tracks;
     EZH2/mark are attenuated + phase-lagged, more so at fast T). The graded G1-lengthening makes the cycle (and thus
     EZH2's S-phase supply) slow under subthreshold dips, deepening the buffering.

Single-cell deterministic (clearest for the memory function). Piecewise integration with tolerance retry.
Run: ./venv/bin/python simulations/sim_memory_ramp_fluctuation.py   (writes _cache.npz; fig = fig_memory_ramp_fluctuation.py)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

UPH = 62.8
SEL = ['time', 'SHH', 'Cd', 'EZH2', 'Mk', 'Skp2', 'Dna']
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)


def _integrate(rr, shh_fn, t0, t1, chunk=40.0):
    """Piecewise integrate with a time-varying SHH; tolerance retry for stiff low-Shh stretches."""
    T, out = [], {k: [] for k in SEL if k != 'time'}
    t = t0
    while t < t1 - 1e-6:
        rr['SHH'] = float(shh_fn(t))
        r = None
        for atol in (1e-8, 1e-7, 1e-6, 1e-5):
            try:
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + chunk, max(8, int(chunk / 4)), selections=SEL); break
            except Exception:
                rr.integrator.setValue('maximum_time_step', 5.0); continue
        if r is None:
            t += chunk; continue
        T.append(np.asarray(r['time'])[1:])
        for k in out:
            out[k].append(np.asarray(r[k])[1:])
        t += chunk
    T = np.concatenate(T)
    return T, {k: np.concatenate(v) for k, v in out.items()}


def _new(mark_off=False):
    params = {'f0_prc2': 1.0} if mark_off else None      # f0_prc2=1 -> mark exerts ZERO CyclinD1 repression
    rr = te.loada(build_model_v44(params=params))
    rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('maximum_num_steps', 1000000)
    return rr


def _divrate(t, dna, win=80.0):
    """Local proliferation rate: divisions (Dna 1->0) per hour, in a sliding window (per model-time -> /h)."""
    dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    dt = t[dv]
    rate = np.zeros_like(t)
    for i, ti in enumerate(t):
        rate[i] = np.sum(np.abs(dt - ti) < win * UPH / 2) / (win / UPH * 1.0)   # divisions per hour, /model window
    return rate


def ramp(dur_h, shh_lo=0.05, shh_hi=1.0, mark_off=False):
    """Equilibrate at shh_lo, ramp UP to shh_hi over dur_h, then DOWN to shh_lo over dur_h."""
    rr = _new(mark_off); D = dur_h * UPH
    rr['SHH'] = shh_lo
    for atol in (1e-8, 1e-7, 1e-6):
        try: rr.integrator.setValue('absolute_tolerance', atol); rr.simulate(0, 6000, 200, selections=['time']); break
        except Exception: continue
    up = lambda t: shh_lo + (shh_hi - shh_lo) * (t / D)
    dn = lambda t: shh_hi - (shh_hi - shh_lo) * ((t - D) / D)
    tU, oU = _integrate(rr, up, 0, D)
    tD, oD = _integrate(rr, dn, D, 2 * D)
    return (oU['SHH'], oU, _divrate(tU, oU['Dna'])), (oD['SHH'], oD, _divrate(tD, oD['Dna']))


def fluct(period_h, amp=0.45, shh0=0.35, n_per=6):
    rr = _new(False); T = period_h * UPH
    rr['SHH'] = shh0
    for atol in (1e-8, 1e-7, 1e-6):
        try: rr.integrator.setValue('absolute_tolerance', atol); rr.simulate(0, 6000, 200, selections=['time']); break
        except Exception: continue
    fn = lambda t: float(np.clip(shh0 + amp * np.sin(2 * np.pi * (t - 6000) / T), 0.02, 1.1))
    t, o = _integrate(rr, fn, 6000, 6000 + n_per * T, chunk=min(30.0, T / 8))
    return t - 6000, o


if __name__ == '__main__':
    print("Memory function on the current default: mitogen ramp (hysteresis) + fluctuation (buffering).\n")
    res = {}

    # ---- A. ramps at fast (60h) and slow (240h), full model + mark-off control ----
    print("=== A. RAMP up/down: hysteresis loop area (up-branch vs down-branch) ===")
    for dur in (60.0, 240.0):
        for mo, tag in ((False, 'full'), (True, 'markoff')):
            (sU, oU, rU), (sD, oD, rD) = ramp(dur, mark_off=mo)
            res[f'ramp_{int(dur)}_{tag}'] = dict(sU=sU, sD=sD, **{f'U_{k}': oU[k] for k in ('Cd', 'EZH2', 'Mk', 'Skp2')},
                                                 **{f'D_{k}': oD[k] for k in ('Cd', 'EZH2', 'Mk', 'Skp2')}, rU=rU, rD=rD)
            # hysteresis area of CyclinD1 vs SHH (down minus up, normalized)
            grid = np.linspace(0.1, 0.95, 40)
            cdU = np.interp(grid, sU, oU['Cd']); cdD = np.interp(grid, sD[::-1], oD['Cd'][::-1])
            area = np.trapz(cdD - cdU, grid) / (np.trapz(cdU, grid) + 1e-9)
            ezU = np.interp(grid, sU, oU['EZH2']); ezD = np.interp(grid, sD[::-1], oD['EZH2'][::-1])
            ez_area = np.trapz(ezD - ezU, grid) / (np.trapz(ezU, grid) + 1e-9)
            print(f"  ramp {dur:>5.0f}h {tag:>7}: CyclinD1 loop area {area:+.3f} | EZH2 loop area {ez_area:+.3f}")

    # ---- B. fluctuation at fast/medium/slow ----
    print("\n=== B. FLUCTUATION: attenuation + lag of EZH2/mark vs CyclinD1 (low-pass buffering) ===")
    print(f"{'period(h)':>9} | {'Cd CV':>6} {'EZH2 CV':>8} {'Mk CV':>6} | {'EZH2/Cd atten':>13}")
    for T in (14.0, 33.0, 80.0):
        t, o = fluct(T)
        w = t > T * UPH                 # skip first period
        def cv(x): return np.std(x[w]) / (np.mean(x[w]) + 1e-9)
        cd_cv, ez_cv, mk_cv = cv(o['Cd']), cv(o['EZH2']), cv(o['Mk'])
        res[f'fluct_{int(T)}'] = dict(t=t, **{k: o[k] for k in ('SHH', 'Cd', 'EZH2', 'Mk', 'Skp2')})
        print(f"  {T:>7.0f}  | {cd_cv:>6.3f} {ez_cv:>8.3f} {mk_cv:>6.3f} | {ez_cv/(cd_cv+1e-9):>13.3f}")

    np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_memory_ramp_fluctuation_cache.npz'),
             **{f'{k}__{kk}': vv for k, d in res.items() for kk, vv in (d.items() if isinstance(d, dict) else [])})
    import pickle
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sim_memory_ramp_fluctuation_cache.pkl'), 'wb') as fh:
        pickle.dump(res, fh)
    print("\n-> cache saved; figure: fig_memory_ramp_fluctuation.py")
