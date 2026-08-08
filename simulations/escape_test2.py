"""Full transient-G0 toggle: option B (p27 brakes) + slow basal p27 clearance (kDeKPC low -> clearance is
CDK2/Skp2-DEPENDENT, so a high-p27 dwell self-sustains) + CyclinD-hyper-P escape (accumulating CyclinD1 re-enters).
Sweep (kDeKPC, P21_div) for MB; look for a regime where MB CYCLES with a birth-p27-DEPENDENT dwell (period rises
with P21_div = transient G0) and GNP still cycles. Run: PYTHONPATH=. ./venv/bin/python simulations/escape_test2.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

OB=dict(kSeqCd=1.232,kRelCd=0.0141,w_ink4=6.0,kSyP21=0.000707,kPhRbCd=0.339,K_CdRb=0.417)
m=build_model_v44(with_p27_optionB=True, with_cd_hyper_escape=True)

def run(cell,p21div,kdekpc,whyp=0.15,t=18000,npts=36000):
    rr=te.loada(m); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    for k,v in OB.items(): rr[k]=v
    rr['kDeKPC']=kdekpc; rr['w_cd_hyper']=whyp; rr['P21_div']=p21div
    if cell=='GNP': rr['SHH']=0.5; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    else:
        rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            return rr.simulate(0,t,npts,selections=["time","MPF"])
        except Exception: pass
    return None
def ndiv(res):
    if res is None: return (-1,float('nan'))
    t=res['time']; m2=t>=3000; tt=t[m2]; dt=tt[1]-tt[0]
    pk,_=find_peaks(res['MPF'][m2],prominence=0.15,distance=int(200/dt))
    per=np.median(np.diff(tt[pk]/60)) if len(pk)>1 else float('nan')
    return len(pk),per

print("option B + escape: MB period vs (kDeKPC, P21_div). Transient G0 = MB cycles with period RISING vs P21_div.")
print(f"  {'kDeKPC':>8s} | " + " ".join(f"P21d={p}" for p in [0.6,1.5,2.5,3.5]) + " | GNP(0.6)")
for kd in [0.02,0.005,0.001,0.0003]:
    row=[]
    for p in [0.6,1.5,2.5,3.5]:
        n,per=ndiv(run('MB',p,kd)); row.append(f"{n}d/{per:.0f}h" if n>0 else ("arr" if n==0 else "err"))
    ng,perg=ndiv(run('GNP',0.6,kd))
    print(f"  {kd:8.4f} | " + " ".join(f"{x:>9s}" for x in row) + f" | {ng}d/{perg:.0f}h")
print("\n  Want a kDeKPC row where MB period RISES with P21_div (dwell grows) but MB still cycles (n>=3), GNP cycles.")
