"""Joint calibration for the HU S-entry block: KmHU_fire (model runtime) x SYN_THR (S-gate on the
synthesis flux vfork*aRc). Targets (validation count-fractions): DMSO-S 15.7% [10.2-21.2], HU-S-fold
1.36 [0.82-1.90], HU-G2-fold 0.23 [0.14-0.32], EZH2-in-S boost 1.31 [0.98-1.64]. Uses the branch model
(fire_gate_HU). Single deterministic cell per condition, count-fraction reweight (matches validate_v44)."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MYCN_AMP_MB, PTCH1_MB, P16_MB, P18_MB, KSYP21_MB, P27_THR = 2.8, 0.1, 0.306, 1.553, 0.002, 0.1
SEL = ['time', 'MPF', 'EZH2', 'P21', 'aRc', 'Dna', 'vfork']
M = build_model_v44(with_ezh2=True, with_hh=True)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 400000)
except Exception: pass


def run(hu, kmhu_fire, kdeez=0.00015):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        _R.reset()
        _R['SHH'] = 0.5; _R['Ptch1_copy_number'] = PTCH1_MB; _R['MYCN_amplification'] = MYCN_AMP_MB
        _R['p16'] = P16_MB; _R['p18'] = P18_MB; _R['kSyP21'] = KSYP21_MB
        _R['HU'] = hu; _R['KmHU_fire'] = kmhu_fire; _R['kDeEZ'] = kdeez
        _R.integrator.setValue('absolute_tolerance', atol)
        try: return _R.simulate(0, 12000, 48000, selections=SEL)
        except Exception: pass
    return None


def classify(res, syn_thr, settle=4000):
    t = res['time']; m = t >= settle; tt = t[m]; dt = tt[1] - tt[0]
    P21 = res['P21'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]; MPF = res['MPF'][m]; vf = res['vfork'][m]; ez = res['EZH2'][m]
    syn = vf * aRc                                   # BrdU/EdU analog: active synthesis flux
    in_S = (syn > syn_thr) & (Dna < 0.98)
    in_G2 = (Dna >= 0.98)
    preS = ~in_S & ~in_G2
    sels = {"G0": preS & (P21 > P27_THR), "G1": preS & (P21 <= P27_THR), "S": in_S, "G2": in_G2}
    ezph = {k: (float(ez[v].mean()) if v.any() else np.nan) for k, v in sels.items()}
    pk, _ = find_peaks(MPF, prominence=0.15, distance=int(200 / dt))
    if len(pk) >= 2:
        divt = tt[pk]; lam = np.log(2) / float(np.mean(np.diff(divt)))
        idx = np.searchsorted(divt, tt, side='right') - 1
        age = np.where(idx >= 0, tt - divt[np.clip(idx, 0, None)], np.nan)
        w = np.where(idx >= 0, 2 * lam * np.exp(-lam * np.clip(age, 0, None)), 0.0); W = np.sum(w) or 1.0
        cnt = {k: 100 * float(np.sum(w * v) / W) for k, v in sels.items()}
    else:
        cnt = {k: 100 * float(v.mean()) for k, v in sels.items()}
    sc = sum(cnt.values()) or 1.0; cnt = {k: v / sc * 100 for k, v in cnt.items()}
    return cnt, ezph


TOL = dict(DMSO_S=(15.7, 0.35), HU_S=(1.36, 0.40), HU_G2=(0.23, 0.40), boost=(1.31, 0.25))
def ok(name, v):
    t, tol = TOL[name]; return abs(v - t) <= tol * t if v == v else False


if __name__ == '__main__':
    KM, ST = 0.25, 0.02                                   # strong block + low BrdU thr (G2 depletes, HU-S counted)
    KDEEZ = [0.00015, 0.0003, 0.0005, 0.001, 0.002]       # 3rd knob: EZH2 turnover (default 77h -> faster)
    print(f'3rd-knob check at KmHU_fire={KM}, SYN_THR={ST}: does faster kDeEZ lower the boost & keep G2?')
    print(f'{"kDeEZ":>8} {"t1/2":>6} | {"DMSO_S":>7} {"HU_Sf":>6} {"HU_G2f":>7} {"boost":>6} {"protG2/G0":>9} | pass')
    for kd in KDEEZ:
        mb = run(0.0, KM, kd); hu = run(1.0, KM, kd)
        f0, ez0 = classify(mb, ST); f1, ez1 = classify(hu, ST)
        dS = f0['S']; huS = f1['S']/f0['S'] if f0['S'] else np.nan
        huG2 = f1['G2']/f0['G2'] if f0['G2'] else np.nan
        boost = ez1['S']/ez0['S'] if (ez0['S'] == ez0['S'] and ez0['S']) else np.nan
        pg2 = ez0['G2']/ez0['G0'] if ez0['G0'] else np.nan   # EZH2 protein G2/G0 side-effect (target 1.48)
        passes = sum([ok('DMSO_S', dS), ok('HU_S', huS), ok('HU_G2', huG2), ok('boost', boost)])
        th = np.log(2)/kd/60
        print(f'{kd:8g} {th:5.0f}h | {dS:7.1f} {huS:6.2f} {huG2:7.2f} {boost:6.2f} {pg2:9.2f} | {passes}/4')
