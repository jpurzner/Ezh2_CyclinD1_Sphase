"""Measure H3K27me3 (Mk) persistence for a given baked config: steady state -> EZH2i on (writing off) ->
time for Mk to fall to half (hours). Eraser OFF (k_jmjd3_gli=0). Reports GNP & MB half-life + steady mark +
ChIP ratio + CyclinD1 fold. Reads the config from optimize_twostep_g0_best.json (or a JSON arg).
Run: ./venv/bin/python simulations/_mark_halflife.py [config.json]"""
import os, sys, json
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase')
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44
UPH = 62.8
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553)

src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), 'optimize_twostep_g0_best.json')
cfg = json.load(open(src))
cfg = cfg.get('best', {}).get('params', cfg) if isinstance(cfg, dict) and 'best' in cfg else cfg
cfg = {**cfg, 'k_jmjd3_gli': 0.0}    # eraser OFF
print(f"config from {os.path.basename(src)}; eraser OFF (k_jmjd3_gli=0)")
M = build_model_v44(params=cfg)


def measure(cond):
    rr = te.loada(M); rr.reset(); rr['SHH'] = 0.5
    for k, v in cond.items(): rr[k] = v
    try:
        r = rr.simulate(0, 9000, 3000, selections=['time', 'Mk', 'Cd'])
    except Exception as e:
        return np.nan, np.nan, np.nan
    m = r['time'] >= 6000
    mk = float(np.mean(np.asarray(r['Mk'])[m])); cd = float(np.mean(np.asarray(r['Cd'])[m]))
    rr['EZH2i'] = 1.0
    try:
        r2 = rr.simulate(9000, 9000 + 6000, 3000, selections=['time', 'Mk'])
    except Exception:
        return mk, cd, np.nan
    mk2 = np.asarray(r2['Mk']); t2 = np.asarray(r2['time']) - 9000
    idx = np.where(mk2 <= mk2[0] / 2)[0]
    thalf = (t2[idx[0]] / UPH) if len(idx) else (6000 / UPH)   # censored at ~95h
    return mk, cd, thalf


mkg, cdg, hg = measure(GNP); mkm, cdm, hm = measure(MB)
print(f"\n{'':>6} {'steady Mk':>10} {'Cd':>8} {'mark t1/2 (h)':>14}")
print(f"{'GNP':>6} {mkg:>10.3f} {cdg:>8.3f} {hg:>14.1f}")
print(f"{'MB':>6} {mkm:>10.3f} {cdm:>8.3f} {hm:>14.1f}")
print(f"\nChIP MB/GNP mark ratio: {mkm/mkg if mkg>1e-6 else float('nan'):.3f}  (eraser off -> ~1.0 expected)")
print(f"CyclinD1 MB/GNP fold  : {cdm/cdg if cdg>1e-6 else float('nan'):.3f}")
print(f"PERSISTENCE: MB mark half-life {hm:.1f}h, GNP {hg:.1f}h (vs ~0.2h in the eraser-ON default)")
