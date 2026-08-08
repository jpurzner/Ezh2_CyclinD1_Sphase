"""Test the CyclinD-hyper-P escape: does with_cd_hyper_escape turn the permanent high-p27 arrest into a TRANSIENT G0
(dwell -> re-entry)? Sweep w_cd_hyper x birth-p27 for MB; confirm GNP still cycles. A high-p27 MB birth should now
give ndiv>=3 (cycles) with a LONGER first-cycle (dwell) then normal period -> transient G0. Too-high w_cd_hyper would
let GNP over-commit (breaks Shh-dep). Run: PYTHONPATH=. ./venv/bin/python simulations/escape_test.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def run(model, cell, p21div, whyp, t=16000, npts=32000):
    rr=te.loada(model); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    if cell=='GNP': rr['SHH']=0.5; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    else:
        rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    rr['P21_div']=p21div; rr['w_cd_hyper']=whyp
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            return rr.simulate(0,t,npts,selections=["time","MPF"])
        except Exception as e: last=e
    raise last
def ana(res,settle=0):
    t=res['time']; dt=t[1]-t[0]
    pk,_=find_peaks(res['MPF'],prominence=0.15,distance=int(200/dt)); tpk=t[pk]/60
    first=tpk[0] if len(tpk) else float('nan')
    per=np.median(np.diff(tpk[tpk>tpk[0]+1])) if len(tpk)>2 else float('nan')
    return len(pk), first, per

m=build_model_v44(with_cd_hyper_escape=True)
print("=== escape ON: does high birth-p27 give TRANSIENT G0 (cycles w/ dwell) instead of permanent arrest? ===")
print(f"  {'cell':5s} {'w_cd_hyper':>10s} {'P21_div':>8s} {'ndiv':>5s} {'1st_div(h)':>11s} {'period(h)':>10s}  verdict")
for whyp in [0.10,0.20,0.35]:
    for cell,p21div in [('MB',0.6),('MB',1.8),('MB',2.5),('GNP',0.6)]:
        try:
            n,first,per=ana(run(m,cell,p21div,whyp))
            v = ("TRANSIENT (cycles+dwell)" if (cell=='MB' and p21div>1 and n>=3 and first>35) else
                 "cycles" if n>=3 else "ARREST")
            print(f"  {cell:5s} {whyp:10.2f} {p21div:8.1f} {n:5d} {first:11.1f} {per:10.1f}  {v}")
        except Exception as e:
            print(f"  {cell:5s} {whyp:10.2f} {p21div:8.1f}  ERROR {str(e)[:40]}")
    print()
print("  Want: MB P21_div>1 -> ndiv>=3 with a LONG 1st-div (dwell) then ~normal period = transient G0.")
print("  GNP P21_div=0.6 must still cycle (~16h). If GNP arrests -> escape too weak; if GNP 1st-div too short at")
print("  high w_cd_hyper -> commitment threshold dropped (check Shh-dep in validation).")
