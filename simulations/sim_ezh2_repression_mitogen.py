"""LARGE overnight study: GNP + MB cells over multiple divisions, with the EZH2 -| CyclinD1 repression ON vs OFF.
(JP 2026-07-20.) Emphasis: when the EZH2->CyclinD1 repression is ON, a cell needs MORE mitogen to proliferate than
when it is OFF -- the repression shifts the mitogen requirement (commitment threshold) UP.

Repression ON  = f0_prc2 at its baked value (mark represses CyclinD1: R = f0 + (1-f0)/(1+(PRC2/K)^n)).
Repression OFF = f0_prc2 = 1.0 (R = 1 for all mark -> ZERO CyclinD1 repression).

Experiments (full v44 default model, multi-division):
  1. DOSE-RESPONSE + mitogen THRESHOLD, GNP & MB, ON vs OFF: sweep Shh, steady-state division rate + CyclinD1 + EZH2
     + mark. Extract the proliferation threshold (Shh at half-max division). => threshold_ON > threshold_OFF.
  2. GNP mitogen RAMP up/down, ON vs OFF, 2 rates: proliferation/CyclinD1 hysteresis + how the repression shifts the
     up- vs down-thresholds (memory).
  3. MULTI-DIVISION trajectories, GNP & MB, ON vs OFF: long runs, division timing + memory variables (EZH2, mark, Skp2).
  4. PARAMETER VARIATION: feedback strength (f0_prc2) x mark de-saturation (a_rw_prc2, del_mk) x mitogen-dose EZH2
     (kEZbas_Cd) -> the mitogen-requirement shift (ON vs OFF) AND whether the mark becomes mitogen-graded (peaks at
     max proliferation, per the ChIP-seq) rather than saturated-flat.

Parallel (ProcessPool) for the dose-response + param sweep; piecewise integration for the ramps.
Run: ./venv/bin/python simulations/sim_ezh2_repression_mitogen.py [--quick]   (writes _cache.pkl; fig later)
"""
import os, sys, pickle, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, tellurium as te
from concurrent.futures import ProcessPoolExecutor, as_completed
from src.build_model_v44_heldt import build_model_v44

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
QUICK = '--quick' in sys.argv
# cell-type conditions (from validate_v44)
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.006)
SEL = ['time', 'Dna', 'Cd', 'EZH2', 'Mk', 'Skp2', 'pRb']


def _rr(overrides):
    rr = te.loada(build_model_v44(params=(overrides or None)))
    rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('maximum_num_steps', 2000000)
    return rr


def _steady(rr, shh, cond, settle=6000, run=13000, npts=None):
    rr.reset()
    for k, v in cond.items(): rr[k] = v
    rr['SHH'] = shh
    npts = npts or int(run / 6)
    r = None
    for atol in (1e-8, 1e-7, 1e-6, 1e-5):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, run, npts, selections=SEL); break
        except Exception:
            try: rr.integrator.setValue('maximum_time_step', 5.0)
            except Exception: pass
    if r is None: return dict(div_h=0.0, Cd=np.nan, EZH2=np.nan, Mk=np.nan, Skp2=np.nan)
    t = r['time']; m = t >= settle; tt = t[m]; dna = np.asarray(r['Dna'])[m]
    dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    span_h = (tt[-1] - tt[0]) / UPH
    return dict(div_h=len(dv) / span_h if span_h > 0 else 0.0,
                Cd=np.asarray(r['Cd'])[m].mean(), EZH2=np.asarray(r['EZH2'])[m].mean(),
                Mk=np.asarray(r['Mk'])[m].mean(), Skp2=np.asarray(r['Skp2'])[m].mean())


def _threshold(shh, div_h, frac=0.5):
    div_h = np.asarray(div_h); shh = np.asarray(shh)
    if div_h.max() < 1e-4: return np.nan
    lvl = frac * div_h.max()
    above = div_h >= lvl
    if not above.any(): return np.nan
    i = np.argmax(above)                       # first Shh where division reaches half-max
    if i == 0: return float(shh[0])
    x0, x1, y0, y1 = shh[i - 1], shh[i], div_h[i - 1], div_h[i]
    return float(x0 + (lvl - y0) * (x1 - x0) / (y1 - y0 + 1e-9))


def dose_task(args):
    """One dose-response curve: (label, celltype-cond, overrides, shh_list)."""
    label, cond, overrides, shh_list = args
    rr = _rr(overrides)
    out = [_steady(rr, s, cond) for s in shh_list]
    div = [o['div_h'] for o in out]
    return label, dict(shh=np.asarray(shh_list), div_h=np.asarray(div),
                       Cd=np.asarray([o['Cd'] for o in out]), EZH2=np.asarray([o['EZH2'] for o in out]),
                       Mk=np.asarray([o['Mk'] for o in out]), Skp2=np.asarray([o['Skp2'] for o in out]),
                       thr=_threshold(shh_list, div))


def ramp(overrides, cond, dur_h, shh_lo=0.05, shh_hi=1.0, chunk=40.0):
    rr = _rr(overrides); rr.reset()
    for k, v in cond.items(): rr[k] = v
    rr['SHH'] = shh_lo
    for atol in (1e-8, 1e-7, 1e-6):
        try: rr.integrator.setValue('absolute_tolerance', atol); rr.simulate(0, 6000, 200, selections=['time']); break
        except Exception: continue
    D = dur_h * UPH
    def run_leg(t0, t1, fn):
        T, SH, DIV, CD, EZ, MK = [], [], [], [], [], []
        t = t0
        while t < t1 - 1e-6:
            rr['SHH'] = float(fn(t)); r = None
            for atol in (1e-8, 1e-7, 1e-6, 1e-5):
                try: rr.integrator.setValue('absolute_tolerance', atol); r = rr.simulate(t, t + chunk, 10, selections=SEL); break
                except Exception:
                    try: rr.integrator.setValue('maximum_time_step', 5.0)
                    except Exception: pass
            if r is None: t += chunk; continue
            tt = np.asarray(r['time'])[1:]; dna = np.asarray(r['Dna'])[1:]
            dv = np.sum((dna[:-1] > 0.9) & (dna[1:] < 0.1))
            T.append(tt.mean()); SH.append(float(fn(t + chunk / 2))); DIV.append(dv / (chunk / UPH))
            CD.append(np.asarray(r['Cd'])[1:].mean()); EZ.append(np.asarray(r['EZH2'])[1:].mean()); MK.append(np.asarray(r['Mk'])[1:].mean())
            t += chunk
        return dict(shh=np.asarray(SH), div_h=np.asarray(DIV), Cd=np.asarray(CD), EZH2=np.asarray(EZ), Mk=np.asarray(MK))
    up = run_leg(0, D, lambda t: shh_lo + (shh_hi - shh_lo) * t / D)
    dn = run_leg(D, 2 * D, lambda t: shh_hi - (shh_hi - shh_lo) * (t - D) / D)
    return up, dn


def trajectory(overrides, cond, shh, T_h=400):
    rr = _rr(overrides); rr.reset()
    for k, v in cond.items(): rr[k] = v
    rr['SHH'] = shh
    for atol in (1e-8, 1e-7, 1e-6):
        try: rr.integrator.setValue('absolute_tolerance', atol); r = rr.simulate(0, T_h * UPH, int(T_h * 4), selections=SEL); break
        except Exception: r = None
    if r is None: return None
    return {k: np.asarray(r[k]) for k in SEL}


if __name__ == '__main__':
    t0 = time.time(); res = {}
    NSHH = 10 if QUICK else 18
    SHH = np.round(np.geomspace(0.04, 1.3, NSHH), 4)
    WK = 6

    # ============ EXP 1: dose-response + threshold, GNP & MB, repression ON vs OFF ============
    print(f"=== EXP 1: dose-response + mitogen threshold (GNP & MB, repression ON vs OFF), {NSHH} Shh pts ===", flush=True)
    tasks = [('GNP_ON', GNP, None, SHH), ('GNP_OFF', GNP, {'f0_prc2': 1.0}, SHH),
             ('MB_ON', MB, None, SHH), ('MB_OFF', MB, {'f0_prc2': 1.0}, SHH)]
    with ProcessPoolExecutor(max_workers=WK) as ex:
        for lab, d in (f.result() for f in as_completed([ex.submit(dose_task, t) for t in tasks])):
            res[f'dose_{lab}'] = d
    for ct in ('GNP', 'MB'):
        on, off = res[f'dose_{ct}_ON']['thr'], res[f'dose_{ct}_OFF']['thr']
        print(f"  {ct}: mitogen threshold ON {on:.3f}  vs OFF {off:.3f}  -> repression needs {on/off if off else float('nan'):.2f}x more mitogen", flush=True)

    # ============ EXP 2: GNP ramp up/down, ON vs OFF, 2 rates ============
    print("\n=== EXP 2: GNP ramp up/down (repression ON vs OFF) ===", flush=True)
    for dur in ([120.0] if QUICK else [60.0, 200.0]):
        for lab, ov in (('ON', None), ('OFF', {'f0_prc2': 1.0})):
            up, dn = ramp(ov, GNP, dur)
            res[f'ramp_{int(dur)}_{lab}'] = dict(up=up, dn=dn)
        print(f"  ramp {dur:.0f}h done (ON & OFF)", flush=True)

    # ============ EXP 3: multi-division trajectories ============
    print("\n=== EXP 3: multi-division trajectories (GNP & MB, ON vs OFF) ===", flush=True)
    for ct, cond, shh in (('GNP', GNP, 0.5), ('MB', MB, 0.5)):
        for lab, ov in (('ON', None), ('OFF', {'f0_prc2': 1.0})):
            tr = trajectory(ov, cond, shh, T_h=(250 if QUICK else 450))
            if tr is not None: res[f'traj_{ct}_{lab}'] = tr
    print("  trajectories done", flush=True)

    # ============ EXP 4: 2D SWEEP -> find a regime with BOTH the mitogen-requirement shift AND the correct
    #      mark direction (peaks at max proliferation, per ChIP-seq). Repression strength (f0_prc2) x mark
    #      de-saturation (a_rw_prc2 & del_mk co-varied). thr_OFF is mark-independent (R=1) -> compute once. ============
    print("\n=== EXP 4: 2D sweep (K_prc2 repression sensitivity  x  mark de-saturation), f0_prc2=0.184 ===", flush=True)
    # KP = PRC2-Hill threshold: LOWER = strong repression at a LOWER mark -> lets a de-saturated (ChIP-seq-correct)
    # mark still repress hard. mark de-saturation (a_rw_prc2, del_mk): idx 0 saturated -> higher = mitogen-graded.
    KP = [0.0074, 0.004, 0.002, 0.001, 0.0005] if not QUICK else [0.0074]
    MARKS = [(0.0023804644163157476, 0.0004037298556555535),   # 0 default (saturated, Mk~0.98)
             (0.0012, 0.0009), (0.0008, 0.0015), (0.0005, 0.0024), (0.0003, 0.004)] if not QUICK else [(0.0023804644163157476, 0.0004037298556555535)]
    off_ref = dose_task(('OFF_ref', GNP, {'f0_prc2': 1.0}, SHH))[1]; res['sweep_OFF_ref'] = off_ref
    combos = []
    for ki, kp in enumerate(KP):
        for mi, (arw, dm) in enumerate(MARKS):
            combos.append((f'k{ki}_m{mi}', GNP, {'f0_prc2': 0.184, 'K_prc2': kp, 'a_rw_prc2': arw, 'del_mk': dm}, SHH))
    with ProcessPoolExecutor(max_workers=WK) as ex:
        for lab, d in (f.result() for f in as_completed([ex.submit(dose_task, t) for t in combos])):
            res[f'sweep_{lab}'] = d
    res['sweep_meta'] = dict(KP=KP, MARKS=MARKS)
    print(f"  {len(combos)} ON combos + 1 OFF ref done. thr_OFF={off_ref['thr']:.3f}", flush=True)
    print(f"\n  {'K_prc2':>7} {'mk_desat':>8} | {'thr_ON':>7} {'shift':>6} | {'Mk_lo':>6} {'Mk_hi':>6} {'Mk_slope':>9} {'BOTH?':>6}")
    both = []
    for ki, kp in enumerate(KP):
        for mi in range(len(MARKS)):
            on = res.get(f'sweep_k{ki}_m{mi}')
            if on is None: continue
            mk_lo, mk_hi = on['Mk'][2], on['Mk'][-3]; slope = mk_hi - mk_lo
            sh = (on['thr'] / off_ref['thr']) if off_ref['thr'] else float('nan')
            ok = (sh > 1.5) and (slope > 0.005)
            if ok: both.append((kp, mi, sh, slope))
            print(f"  {kp:>7.4f} {mi:>8} | {on['thr']:>7.3f} {sh:>6.2f} | {mk_lo:>6.3f} {mk_hi:>6.3f} {slope:>+9.3f} {'YES' if ok else '-':>6}")
    print(f"\n  Regimes with BOTH (mitogen shift >1.5x AND mark peaks-at-proliferation slope>0.005): {len(both)}")
    for kp, mi, sh, sl in both:
        print(f"    K_prc2={kp}, mark-desat idx {mi}: shift {sh:.2f}x, mark slope {sl:+.3f}  <-- 'almost there' candidate")

    with open(os.path.join(HERE, 'sim_ezh2_repression_mitogen_cache.pkl'), 'wb') as fh:
        pickle.dump(res, fh)
    print(f"\n==== DONE in {(time.time()-t0)/60:.1f} min -> sim_ezh2_repression_mitogen_cache.pkl ; figure next ====", flush=True)
