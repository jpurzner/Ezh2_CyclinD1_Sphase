"""Measure the mark-memory CDK6 (dynamic species, reuses PRC2_rep). Targets: MB/GNP ~5x (protein)/~17x (RNA),
vismo ~0.41x (PARTIAL, persistent baseline -> establish-then-withdraw), EZH2-cKO ~2.39x. + divisions.
Run: PYTHONPATH=. ./venv/bin/python simulations/cdk6_markmem_probe.py '<optional json params>'
"""
import sys, json, numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44
P=json.loads(sys.argv[1]) if len(sys.argv)>1 else {}
m=build_model_v44(with_cdk6_gli=True, params=(P or None))

def newrr(shh,pt,p18,p19,cd2):
    rr=te.loada(m); rr.integrator.setValue('relative_tolerance',1e-6); rr.integrator.setValue('maximum_num_steps',1000000); rr.integrator.setValue('absolute_tolerance',1e-8)
    rr['SHH']=shh;rr['Ptch1_copy_number']=pt;rr['p18']=p18;rr['p19']=p19
    if cd2:
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except: pass
    return rr
def cdk6_ndiv(shh,pt,p18,p19,cd2,ezcko=False,vismo_withdraw=False):
    rr=newrr(shh,pt,p18,p19,cd2)
    if ezcko:  # genetic EZH2 cKO: kill EZH2 synthesis + species
        for k in ['kEZbas','kEZbas_Cd','kEZE2f']:
            try: rr[k]=0.0
            except: pass
        try: rr['EZH2']=0.001
        except: pass
    try:
        if vismo_withdraw:  # establish MB (mark eroded), THEN vismo -> persistent baseline
            rr.simulate(0,10000,10000,selections=['time'])
            rr['HHi']=1.0
            r=rr.simulate(0,2000,4000,selections=['time','MPF','cdk6'])  # ~33h withdrawn
        else:
            r=rr.simulate(0,12000,24000,selections=['time','MPF','cdk6'])
    except Exception as e:
        return float('nan'),-1
    tt=r['time']; mm=tt>=(tt[-1]-3000);
    c=float(np.mean(r['cdk6'][mm]))
    T=tt; dt=T[1]-T[0]; pk,_=find_peaks(r['MPF'],prominence=0.15,distance=int(200/dt))
    return c,len(pk)

gnp,ng=cdk6_ndiv(0.5,1.0,0.464,0.36,False)
mb,nm=cdk6_ndiv(0.5,0.3,1.73,0.58,True)
mbv,_=cdk6_ndiv(0.5,0.3,1.73,0.58,True,vismo_withdraw=True)
cko,_=cdk6_ndiv(0.5,1.0,0.464,0.36,False,ezcko=True)
print(f"cdk6:  GNP={gnp:.3f}(ndiv {ng})  MB={mb:.3f}(ndiv {nm})  MB+vismo={mbv:.3f}  GNP+EZH2cKO={cko:.3f}")
print(f"folds: MB/GNP={mb/gnp:.2f} (target ~5-17)   vismo(MBv/MB)={mbv/mb:.2f} (0.41)   EZH2cKO(cko/GNP)={cko/gnp:.2f} (2.39)")
