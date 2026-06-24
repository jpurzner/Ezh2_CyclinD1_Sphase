"""Fold-preserving calibration of the replicative-dilution H3K27me3 module: find {k_w,k0,del,K,n}
so GNP & MB both cycle, GNP+HHi/MB+HHi arrest, MB+HHi+EZH2i rescues, and the MB/GNP CyclinD1 fold
(~5.07) + EZH2i fold (~2.2) hold. Evaluates a list of candidate param sets (manual seed search).
Usage: ./venv/bin/python simulations/calibrate_h3k27_dilution.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

_RR = te.loada(build_model_v44(with_hh=True, with_h3k27_dilution=True))
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)

COND = {
    'GNP':          dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0),
    'GNP+HHi':      dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0.95, EZH2i=0),
    'GNP+EZH2i':    dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=1),
    'MB':           dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, HHi=0, EZH2i=0),
    'MB+HHi':       dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, HHi=0.95, EZH2i=0),
    'MB+HHi+EZH2i': dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, HHi=0.95, EZH2i=1),
}
PMK = ['k_w_mk', 'k0_mk', 'del_mk', 'K_mk', 'n_mk']


def _run(pmk, cond):
    for hz, npts, settle in [(11000, 22000, 7000), (8000, 16000, 5000), (5500, 11000, 3500)]:
        try:
            _RR.reset()
            for k, v in zip(PMK, pmk):
                _RR[k] = v
            for k, v in cond.items():
                _RR[k] = v
            _RR['Mk'] = 0.2
            r = _RR.simulate(0, hz, npts, selections=['time', 'MPF', 'Cd', 'Mk'])
            m = r['time'] >= settle
            pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=400)
            per = np.mean(np.diff(r['time'][m][pk])) / 60.0 if len(pk) > 1 else float('nan')
            return len(pk), float(r['Cd'][m].mean()), float(r['Mk'][m].mean()), per
        except Exception:
            continue
    return -1, float('nan'), float('nan'), float('nan')


CANDS = {
    'A K0.30 n6': (0.0023, 0.0002, 0.0015, 0.30, 6),
    'B K0.32 n5': (0.0023, 0.0002, 0.0015, 0.32, 5),
    'C kw0.0018': (0.0018, 0.0002, 0.0015, 0.28, 5),
    'D kw0.0028': (0.0028, 0.0002, 0.0018, 0.32, 6),
}
print("targets: GNP cycle (~3 div), MB cycle, GNP+HHi & MB+HHi arrest (0), MB+HHi+EZH2i rescue (>0),")
print("         MB/GNP CyclinD1 fold ~5.07, EZH2i fold ~2.2\n")
for name, pmk in CANDS.items():
    res = {c: _run(pmk, COND[c]) for c in COND}
    g, mb = res['GNP'], res['MB']
    fold = mb[1] / g[1] if g[1] else float('nan')
    ezfold = res['GNP+EZH2i'][1] / g[1] if g[1] else float('nan')
    print(f"== {name}  {dict(zip(PMK, pmk))}")
    for c in COND:
        d, cd, mk, per = res[c]
        print(f"   {c:13s}: {d:2d} div  Cd {cd:6.2f}  Mk {mk:.3f}  per {per:.1f}h")
    print(f"   --> MB/GNP fold {fold:.2f} (t 5.07) | EZH2i fold {ezfold:.2f} (t 2.2)\n")
