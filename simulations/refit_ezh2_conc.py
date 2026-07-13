"""Phase 0b re-fit harness: with the EZH2 division-halving OFF (concentration convention), EZH2 sits ~47%
higher, over-repressing CyclinD1. Re-fit (kTlEZ, kDeEZ) to restore the EZH2-coupled targets. Reuses
validate_v44's EXACT 28-target scoring by monkeypatching its _MODEL.

Usage:
  ./venv/bin/python simulations/refit_ezh2_conc.py --kTlEZ 0.01027                 # single candidate
  ./venv/bin/python simulations/refit_ezh2_conc.py --grid                          # coarse (kTlEZ,kDeEZ) grid
"""
import os, sys, io, contextlib, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ['EZH2_CONC'] = '1'   # ensure any internal rebuilds are halving-off too
from src.build_model_v44_heldt import build_model_v44
import simulations.validate_v44 as V

KTLEZ0, KDEEZ0 = 0.015098862630302894, 0.0011685177801161655   # current baked (with halving)

def score(kTlEZ, kDeEZ, verbose=False):
    params = dict(V.PARAMS or {}); params['kTlEZ'] = kTlEZ; params['kDeEZ'] = kDeEZ
    V._MODEL = build_model_v44(with_ezh2=True, with_hh=True,
                               with_h3k27_dilution=True, with_h3k27_chain=True,
                               with_two_step_rb=True, with_ezh2_conc=True, params=params)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = V.main()
    txt = buf.getvalue()
    if verbose:
        print(txt)
    return out['passed'], out['total'], txt

if __name__ == '__main__':
    if '--grid' in sys.argv:
        best = None
        for f_tl, f_de in itertools.product([0.68, 0.78, 0.88], [0.65, 0.80, 1.0]):
            kt, kd = KTLEZ0 * f_tl, KDEEZ0 * f_de
            p, t, _ = score(kt, kd)
            tag = f"kTlEZ={kt:.5f}(x{f_tl}) kDeEZ={kd:.5f}(x{f_de})"
            print(f"  {tag:48s} -> {p}/{t}")
            if best is None or p > best[0]:
                best = (p, t, kt, kd, tag)
        print(f"\nBEST: {best[4]} -> {best[0]}/{best[1]}")
    else:
        kt = KTLEZ0; kd = KDEEZ0
        for i, a in enumerate(sys.argv):
            if a == '--kTlEZ': kt = float(sys.argv[i + 1])
            if a == '--kDeEZ': kd = float(sys.argv[i + 1])
        p, t, txt = score(kt, kd, verbose=True)
        print(f"\n>>> kTlEZ={kt:.6f}  kDeEZ={kd:.6f}  ->  {p}/{t}")
