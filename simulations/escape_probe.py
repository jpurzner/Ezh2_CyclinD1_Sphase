"""CDK4/6-alone escape (Yang 2020) phenotype in the integrated model (JP 2026-08-09): does high CyclinD1/CDK6
pull MB cells out of the two-step-Rb permanent lock -> reversible transient G0? Env W_CD_HYPER sets strength.

Reports, at MB fixed CDKI 1x/4x/6x/8x: divisions in 60000 min + mean period (h). OFF -> deep CDKI permanent
(0 div); escape ON -> deep arrest becomes SLOW transient-G0 cycling (finite div, long period = the G0 dwell).
Also baseline cdk6/Cd at 1x (should be ~5.08/7.7, unchanged) + the drive-dependence (deeper arrest = longer dwell).

Run:  W_CD_HYPER=0.05 PYTHONPATH=. ./venv/bin/python simulations/escape_probe.py
"""
import os, json, numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)
WCH = os.environ.get('W_CD_HYPER')
P = dict(BP)
if WCH is not None:
    P['w_cd_hyper'] = float(WCH)
M = build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_cd_hyper_escape=(WCH is not None), params=P)
RR = te.loada(M); RR.integrator.setValue('relative_tolerance', 1e-7); RR.integrator.setValue('absolute_tolerance', 1e-9)
RR.integrator.setValue('maximum_num_steps', 12000000)

def run(mult, T=60000):
    for k, v in P.items():
        try: RR[k] = v
        except Exception: pass
    RR['SHH'] = 0.5; RR['Ptch1_copy_number'] = 0.3
    try: RR['Cd2_expr'] = RR['CD2_EXPR_MB']
    except Exception: pass
    RR['p18'] = 1.73 * mult; RR['p19'] = 0.58 * mult; RR['kSyP21'] = BP['kSyP21'] * mult
    RR.reset()
    r = RR.simulate(0, T, T // 10, selections=['time', 'Dna', 'cdk6', 'Cd'])
    d = np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]
    per = float(np.mean(np.diff(r['time'][d])) / 60) if len(d) >= 2 else None
    mm = r['time'] >= T * 0.7
    return int(len(d)), per, float(r['cdk6'][mm].mean()), float(r['Cd'][mm].mean())

out = {'w_cd_hyper': (float(WCH) if WCH is not None else None)}
for mult in (1.0, 4.0, 6.0, 8.0):
    n, per, cdk6, cd = run(mult)
    out[f'{mult:g}x'] = dict(div=n, period_h=(round(per, 1) if per else None),
                             cdk6=round(cdk6, 2), Cd=round(cd, 2),
                             state=('cycling' if n > 20 else ('transient-G0' if n > 0 else 'permanent-arrest')))
print(json.dumps(out))
