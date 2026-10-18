"""Why doesn't elongated S accumulate EZH2? Track the EZH2-synthesis drivers (E2f, CycE=Ce, CycA=Ca, Cd) and
EZH2/EZH2m in S-phase for MB control vs MB+HU (elongated/arrested S). If CycA (or CycE) COLLAPSES under HU while
E2f stays high, then gating EZH2 synthesis on CycE/CycA (not E2f) is why the boost is missing -- HU-arrested
committed cells keep E2f active (should keep making EZH2) but lose the transient cyclins. Run:
./venv/bin/python simulations/_diag_ezh2_hu.py"""
import os, sys
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase')
import json
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44
# Cdh1-EZH2 fix config (halving convention) if requested: python _diag_ezh2_hu.py fix
if len(sys.argv) > 1 and sys.argv[1] == 'fix':
    CFG = json.load(open(os.path.join(os.path.dirname(__file__), 'optimize_cdh1_ezh2_best.json')))['best']['params']
    M = build_model_v44(params={**CFG, 'k_jmjd3_gli': 0.0}, with_ezh2_conc=False)
    print("=== Cdh1-EZH2 FIX config (halving + Cdh1-gated degradation, eraser off) ===")
else:
    M = build_model_v44()
    print("=== current default (concentration convention) ===")
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=None)
SEL = ['time', 'EZH2', 'EZH2m', 'E2f', 'Ce', 'Ca', 'Cd', 'aRc', 'Dna', 'vfork']


def run(hu):
    rr = te.loada(M); rr.reset(); rr['SHH'] = 0.5
    rr['Ptch1_copy_number'] = 0.1; rr['MYCN_amplification'] = 2.8; rr['p16'] = 0.306; rr['p18'] = 1.553
    rr['HU'] = hu
    r = rr.simulate(0, 12000, 24000, selections=SEL)
    m = r['time'] >= 7000
    d = {k: np.asarray(r[k])[m] for k in SEL}
    syn = d['vfork'] * d['aRc']
    inS = (syn > 0.02) & (d['Dna'] < 0.98)         # BrdU+ S (validate definition)
    return d, inS


for hu, lab in [(0.0, 'control'), (1.0, 'HU')]:
    d, inS = run(hu)
    frac = 100 * inS.mean()
    def sm(k): return float(d[k][inS].mean()) if inS.any() else float('nan')
    print(f"\n=== MB {lab}  (S-phase = {frac:.1f}% of time) ===")
    print(f"  EZH2 {sm('EZH2'):.3f}  EZH2m {sm('EZH2m'):.4f} | drivers in S:  E2f {sm('E2f'):.3f}  Ce {sm('Ce'):.4f}  Ca {sm('Ca'):.4f}  Cd {sm('Cd'):.3f}")
