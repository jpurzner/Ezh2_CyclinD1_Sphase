"""Option-B probe: build option A (baked) vs option B (p27-inhibitory), confirm structure + integrate.
GNP and MB, count divisions + report pRb / p27 pools. Option B on option-A params should OVER-brake (needs re-cal).
Run: PYTHONPATH=. ./venv/bin/python simulations/optionB_probe.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def build(optB): return build_model_v44(with_p27_optionB=optB)

# structural check
mA, mB = build(False), build(True)
print("=== structure ===")
print(f"  option A: '+ CdP21)/(K_CdRb' present = {'+ CdP21)/(K_CdRb' in mA}  (expect True)")
print(f"  option B: '+ CdP21)/(K_CdRb' present = {'+ CdP21)/(K_CdRb' in mB}  (expect False)")
print(f"  option B: buffer reaction 'Cd + P21 -> CdP21' present = {'Cd + P21 -> CdP21' in mB}  (expect True -- buffer still consumes Cd)")

def run(model, shh=0.5, ptch1=1.0, p18=0.464, p19=0.36, cd2_mb=False, t=12000, npts=24000):
    rr = te.loada(model)
    rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=ptch1; rr['p18']=p18; rr['p19']=p19
    if cd2_mb:
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            return rr.simulate(0,t,npts,selections=["time","MPF","pRb","Dna","P21","CdP21","Cd","p21a"])
        except Exception as e:
            last=e
    raise last

def divs(res, settle=3000):
    t=res['time']; m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(res['MPF'][m], prominence=0.15, distance=int(200/dt))
    per=np.diff(tt[pk]/60.0) if len(pk)>1 else np.array([])
    return len(pk), (float(np.mean(per)) if len(per) else float('nan'))

def ms(res,sp,settle=4000): t=res['time']; return float(np.mean(res[sp][t>=settle]))

print("\n=== integrate: divisions (period h) | <pRb> <P21> <CdP21> <Cd> ===")
print(f"  {'model/cond':22s} {'ndiv':>5s} {'per(h)':>7s} {'<pRb>':>7s} {'<P21>':>7s} {'<CdP21>':>8s} {'<Cd>':>6s}")
for lab,optB in [('A(baked)',False),('B(p27-inhib)',True)]:
    model=build(optB)
    for cond,(shh,pt,p18,p19,cd2) in [('GNP+SHH',(0.5,1.0,0.464,0.36,False)),('MB',(0.5,0.3,1.73,0.58,True))]:
        try:
            r=run(model,shh=shh,ptch1=pt,p18=p18,p19=p19,cd2_mb=cd2)
            n,per=divs(r)
            print(f"  {lab+' '+cond:22s} {n:5d} {per:7.1f} {ms(r,'pRb'):7.3f} {ms(r,'P21'):7.4f} {ms(r,'CdP21'):8.4f} {ms(r,'Cd'):6.2f}")
        except Exception as e:
            print(f"  {lab+' '+cond:22s}  ERROR: {str(e)[:60]}")
print("\n  (Option B on option-A params is EXPECTED to over-brake -> fewer/no divisions; re-cal needed.)")
