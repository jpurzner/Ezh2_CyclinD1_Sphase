"""Focused calibration sweep for the two-step Rb: find a regime that (1) cycles, (2) has a real
growth-timed p27-high transient G0, and (3) keeps pRb(hyper) LOW during that G0 (markers agree).

Competing knobs:
  kDsRbmE2f : mono-Rb -> E2f release (needed to break the deadlock; too high -> commits too fast)
  kDeP21Cd  : CyclinD -> p27 clearance, size-gated (sets growth-timed commitment)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

BASE = {'kSyCb': 0.04, 'a25': 0.05}


def run(kds, kdcd, shh=0.5):
    p = dict(BASE); p['kDsRbmE2f'] = kds; p['kDeP21Cd'] = kdcd
    rr = te.loada(build_model_v44(with_two_step_rb=True, params=p))
    rr.integrator.setValue('absolute_tolerance', 1e-9); rr.integrator.setValue('relative_tolerance', 1e-6)
    rr['SHH'] = shh
    try:
        r = rr.simulate(0, 16000, 64000, selections=['time', 'MPF', 'Dna', 'P21', 'pRb', 'aRc', 'mass'])
    except Exception as e:
        return None
    m = r['time'] >= 7000; tt = r['time'][m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.1, distance=int(150 / dt))
    if len(pk) < 2:
        return dict(cyc=False)
    per = float(np.mean(np.diff(tt[pk]))) / 60
    P21 = r['P21'][m]; pRb = r['pRb'][m]; aRc = r['aRc'][m]; Dna = r['Dna'][m]
    inS = (aRc > 0.05) & (Dna < 0.98); inG2 = Dna >= 0.98; preS = ~inS & ~inG2
    g0 = preS & (P21 > 0.1)
    g0_h = per * g0.mean() / ((preS | inS | inG2).mean() or 1)
    prb_in_g0 = float(pRb[g0].mean()) if g0.any() else np.nan      # want LOW
    prb_max = float(pRb.max())
    return dict(cyc=True, per=per, g0=g0_h, prb_g0=prb_in_g0, prb_max=prb_max)


print("two-step Rb sweep (SHH=0.5 GNP); want: cycles, G0(p27) several h, pRb-in-G0 << pRb-max")
print(f"{'kDsRbmE2f':>9} {'kDeP21Cd':>9} | {'cyc':>4} {'period':>7} {'G0(p27)':>8} {'pRb@G0':>7} {'pRbmax':>7}")
for kds in [1.0, 1.5, 2.0, 3.0]:
    for kdcd in [0.04, 0.08, 0.15, 0.30]:
        r = run(kds, kdcd)
        if r is None:
            print(f"{kds:>9.1f} {kdcd:>9.2f} | CRASH"); continue
        if not r['cyc']:
            print(f"{kds:>9.1f} {kdcd:>9.2f} | arrest"); continue
        print(f"{kds:>9.1f} {kdcd:>9.2f} | {'yes':>4} {r['per']:>6.1f}h {r['g0']:>7.1f}h "
              f"{r['prb_g0']:>7.2f} {r['prb_max']:>7.2f}")
