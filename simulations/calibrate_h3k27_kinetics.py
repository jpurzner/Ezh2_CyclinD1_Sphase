"""Calibrate the H3K27me3 de-repression kinetics to JP's target: clear CyclinD1 increase by ~24h,
maximal by ~48h after tazemetostat, with BOTH active demethylation (k_demeth, breaks the arrested-cell
deadlock) and replication-dilution (k_dil, accelerates once cycling resumes).

Stage 1 (this script): sweep k_demeth_cd with k_meth=k_demeth, k_dil=0 (validation-exact: mark tracks
EZH2 at steady state). Reports MB cycling steady-state (must stay ~no-memory) + the de-repression
timecourse (CyclinD1 at 12/24/36/48/72h after EZH2i on a vismo-arrested MB cell).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_memory=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)


def _set_tol(a, r):
    rr.integrator.setValue("absolute_tolerance", a); rr.integrator.setValue("relative_tolerance", r)


def _sim(dur, npts, sel):
    try:
        return rr.simulate(0, dur, npts, selections=sel)
    except Exception:
        _set_tol(1e-6, 1e-4)
        try:
            return rr.simulate(0, dur, npts, selections=sel)
        finally:
            _set_tol(1e-9, 1e-6)


def mb_setup(cdscale=1.5):
    rr.reset()
    rr['SHH'] = 0.5; rr['MYCN_amplification'] = 2.8; rr['Ptch1_copy_number'] = 0.1
    rr['p16'] = 0.306; rr['p18'] = 1.553; rr['kSyP21'] = 0.002
    rr['k_Cd_translation'] = 0.801 * cdscale; rr['P21_div'] = 0.5
    rr['EZH2i'] = 0; rr['GDC0449'] = 0


def evaluate(k_demeth, k_meth=None, k_dil=0.0):
    k_meth = k_demeth if k_meth is None else k_meth
    # ---- steady-state MB cycling (validation proxy) ----
    mb_setup()
    rr['k_demeth_cd'] = k_demeth; rr['k_meth_cd'] = k_meth; rr['k_dil_cd'] = k_dil
    r = _sim(9000, 18000, ['time', 'MPF', 'Cd'])
    m = r['time'] >= 5500
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=400)
    ss_cd, ss_div = r['Cd'][m].mean(), len(pk)
    # ---- de-repression timecourse: vismo-arrest, then EZH2i ----
    mb_setup()
    rr['k_demeth_cd'] = k_demeth; rr['k_meth_cd'] = k_meth; rr['k_dil_cd'] = k_dil
    rr.simulate(0, 6500, 13000)            # equilibrate cycling
    rr['GDC0449'] = 0.85                    # vismo
    _sim(2400, 480, ['time', 'Cd'])        # arrest (40h)
    rr['EZH2i'] = 1                         # add inhibitor
    r = _sim(4500, 900, ['time', 'Cd'])    # 75h de-repression phase
    t = r['time'] / 60.0
    cd_at = {h: r['Cd'][np.argmin(np.abs(t - h))] for h in [0, 12, 24, 36, 48, 72]}
    return ss_cd, ss_div, cd_at


# no-memory MB reference (cdscale=1.5) for the steady-state-preservation check
rr_nomem = te.loada(build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_memory=False))
rr_nomem.integrator.setValue("absolute_tolerance", 1e-9); rr_nomem.integrator.setValue("relative_tolerance", 1e-6)
rr_nomem.reset()
for k, v in dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553,
                 kSyP21=0.002, k_Cd_translation=0.801 * 1.5, P21_div=0.5, EZH2i=0, GDC0449=0).items():
    rr_nomem[k] = v
_rn = rr_nomem.simulate(0, 9000, 18000, selections=['time', 'MPF', 'Cd']); _mn = _rn['time'] >= 5500
print(f"no-memory MB reference: ss Cd = {_rn['Cd'][_mn].mean():.2f}")

print("target: de-repression clear by 24h, max ~48h (repressed~1.4 -> ~3.7); keep MB ss Cd ~ no-mem ref")
print(f"\n{'k_demeth':>8} {'k_meth':>7} {'k_dil':>6} | {'MB ss Cd/div':>12} | CyclinD1 after EZH2i:  0h   12h   24h   36h   48h   72h")
# (k_demeth, k_meth, k_dil): demethylation alone vs both-contribute (slow demeth + replication-dilution)
for kd, km, kdil in [(0.0012, 0.0012, 0.0), (0.0008, 0.0008, 0.008), (0.0008, 0.0008, 0.02),
                     (0.0008, 0.0010, 0.02), (0.0010, 0.0010, 0.01)]:
    ss_cd, ss_div, c = evaluate(kd, km, kdil)
    print(f"{kd:8.4f} {km:7.4f} {kdil:6.3f} | {ss_cd:6.2f} / {ss_div:2d}   | "
          f"                      {c[0]:5.2f} {c[12]:5.2f} {c[24]:5.2f} {c[36]:5.2f} {c[48]:5.2f} {c[72]:5.2f}")
