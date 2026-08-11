"""Transient-G0 phenotype probe for the two-compartment PRC2 module (with_proximal_distal).
Given amplifier params (env AMP_PARAMS = JSON, merged on the integrated BP), measures:

  * BASELINE preserved (MB 1x CDKI, SHH0.5): cdk6, Cd, ndiv (should be ~5.08 / ~7.7 / cycling).
  * ARREST DEPTH (MB high CDKI): P_prox, cdk6, Cd in the arrested state (how hard the reservoir represses).
  * TRANSIENT-vs-PERMANENT (hysteresis): arrest at high CDKI, then drop CDKI to 1x mid-run and measure the
    DWELL = time to re-entry (h). None/inf = permanent lock; small = leaky; intermediate = transient dwell.

Run:  AMP_PARAMS='{"a_P":0.3,"f_P":0.5}' PYTHONPATH=. ./venv/bin/python simulations/amplifier_hysteresis.py
"""
import os, json, numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

BP=dict(kPhRbCd=0.1748299602292031,K_CdRb=0.4025007594938155,w_ink4=4.646376316224687,kSyP21=0.0004043793202695936,
        K_cdk6_sink=3.879643001142466,w_p18=1.003950915185932,w_p19=0.6655007434943913,
        k_Cd_tx_Gli_max=0.24197652013556434,k_Cd_tx_MYCN=0.11031706915457935)
AMP=json.loads(os.environ.get('AMP_PARAMS','') or '{}')
P=dict(BP); P.update(AMP)
M=build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_proximal_distal=True, params=P)
RR=te.loada(M); RR.integrator.setValue('relative_tolerance',1e-7); RR.integrator.setValue('absolute_tolerance',1e-9)
RR.integrator.setValue('maximum_num_steps',12000000)

def setp(mult):
    for k,v in P.items():
        try: RR[k]=v
        except: pass
    RR['SHH']=0.5; RR['Ptch1_copy_number']=0.3
    try: RR['Cd2_expr']=RR['CD2_EXPR_MB']
    except: pass
    RR['p18']=1.73*mult; RR['p19']=0.58*mult; RR['kSyP21']=BP['kSyP21']*mult

def divs(dna): return np.where((dna[:-1]>0.9)&(dna[1:]<0.1))[0]

def steady(mult, T=45000):
    setp(mult); RR.reset()
    r=RR.simulate(0,T,T//10,selections=['time','Dna','cdk6','Cd','P_prox'])
    mm=r['time']>=T*0.7; d=divs(r['Dna'])
    return len(d), float(r['cdk6'][mm].mean()), float(r['Cd'][mm].mean()), float(r['P_prox'][mm].mean())

def hysteresis(mult_hi=5.0, t_switch=35000, T=110000):
    setp(mult_hi); RR.reset()
    r1=RR.simulate(0,t_switch,t_switch//10,selections=['time','Dna','cdk6','Cd','P_prox'])
    nb=len(divs(r1['Dna'])); pend=float(r1['P_prox'][-1]); cdk6a=float(r1['cdk6'][-1]); cda=float(r1['Cd'][-1])
    RR['p18']=1.73; RR['p19']=0.58; RR['kSyP21']=BP['kSyP21']
    r2=RR.simulate(t_switch,T,(T-t_switch)//10,selections=['time','Dna'])
    dd=divs(r2['Dna']); na=len(dd)
    dwell=float((r2['time'][dd[0]]-t_switch)/60) if na>0 else None
    return nb, na, dwell, pend, cdk6a, cda

b=steady(1.0); ndH,cdk6H,cdH,pH = steady(6.0)
nb,na,dwell,pend,cdk6a,cda = hysteresis()
out=dict(amp=AMP,
  base_ndiv=b[0], base_cdk6=round(b[1],3), base_Cd=round(b[2],3), base_Pprox=round(b[3],4),
  arrest6x_ndiv=ndH, arrest6x_cdk6=round(cdk6H,3), arrest6x_Cd=round(cdH,3), arrest6x_Pprox=round(pH,3),
  hyst_div_before=nb, hyst_div_after_drop=na, dwell_h=(round(dwell,1) if dwell is not None else None),
  verdict=('permanent-lock' if (nb==0 and na==0) else ('transient-dwell' if (nb==0 and na>0) else 'no-arrest')))
print(json.dumps(out))
