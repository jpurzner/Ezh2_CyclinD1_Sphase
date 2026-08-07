"""Probe the CDK6 drive factor: build with_cdk6_gli, check structure, and report the cdk6 factor value + divisions
at GNP / MB / MB+vismo / GNP+EZH2i / MB+EZH2i. cdk6 should be ~1 at GNP, higher at MB (Gli), drop under vismo,
rise under EZH2i (two-arm governor). Then re-balance kPhRbCd so GNP still cycles ~16h.
Run: PYTHONPATH=. ./venv/bin/python simulations/cdk6_probe.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

mOff = build_model_v44(with_cdk6_gli=False)
mOn  = build_model_v44(with_cdk6_gli=True)
print("=== structure ===")
print(f"  flag off: 'kPhRbCd*cdk6' present = {'kPhRbCd*cdk6' in mOff}  (expect False)")
print(f"  flag on : 'kPhRbCd*cdk6' present = {'kPhRbCd*cdk6' in mOn}   (expect True)")
print(f"  flag on : 'cdk6 :=' present = {'cdk6 :=' in mOn}  (expect True)")

def run(model, shh, ptch1, p18, p19, cd2_mb, hhi=0.0, ezh2i=0.0, params=None, t=11000, npts=22000):
    rr=te.loada(model); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=ptch1; rr['p18']=p18; rr['p19']=p19; rr['HHi']=hhi; rr['EZH2i']=ezh2i
    if cd2_mb:
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    if params:
        for k,v in params.items():
            try: rr[k]=v
            except Exception: pass
    last=None
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            sel=["time","MPF","Cd","EZH2","Gli1"] + (["cdk6"] if 'cdk6 :=' in model else [])
            return rr.simulate(0,t,npts,selections=sel)
        except Exception as e: last=e
    raise last

def ndiv(res, settle=3000):
    t=res['time']; m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(res['MPF'][m],prominence=0.15,distance=int(200/dt))
    per=np.median(np.diff(tt[pk]/60.0)) if len(pk)>1 else float('nan')
    return len(pk),per
def ms(res,sp,settle=4000):
    if sp not in res.colnames: return float('nan')
    t=res['time']; return float(np.mean(res[sp][t>=settle]))

conds=[('GNP',0.5,1.0,0.464,0.36,False,0.0,0.0),
       ('MB',0.5,0.3,1.73,0.58,True,0.0,0.0),
       ('MB+vismo',0.5,0.3,1.73,0.58,True,1.0,0.0),
       ('GNP+EZH2i',0.5,1.0,0.464,0.36,False,0.0,1.0),
       ('MB+EZH2i',0.5,0.3,1.73,0.58,True,0.0,1.0)]

print("\n=== cdk6 factor + divisions (flag ON, PLACEHOLDER params) ===")
print(f"  {'cond':12s} {'cdk6':>6s} {'ndiv':>5s} {'per(h)':>7s} {'<Cd>':>7s} {'<EZH2>':>7s} {'<Gli1>':>7s}")
vals={}
for lab,shh,pt,p18,p19,cd2,hhi,ez in conds:
    r=run(mOn,shh,pt,p18,p19,cd2,hhi=hhi,ezh2i=ez)
    n,per=ndiv(r); c=ms(r,'cdk6'); vals[lab]=c
    print(f"  {lab:12s} {c:6.3f} {n:5d} {per:7.1f} {ms(r,'Cd'):7.2f} {ms(r,'EZH2'):7.3f} {ms(r,'Gli1'):7.3f}")
if vals.get('GNP',0)>0:
    print(f"\n  cdk6 MB/GNP = {vals.get('MB',0)/vals['GNP']:.2f}  (data Cdk6 RNA 17x; functional factor should be modest)")
    print(f"  cdk6 vismo (MB+vismo/MB) = {vals.get('MB+vismo',0)/max(vals.get('MB',1e-9),1e-9):.2f}  (data ~0.41x)")
    print(f"  cdk6 EZH2i (MB+EZH2i/MB) = {vals.get('MB+EZH2i',0)/max(vals.get('MB',1e-9),1e-9):.2f}  (de-repression, expect >1)")
print("\n  NB divisions may be off until kPhRbCd re-balanced (cdk6!=1 at GNP scales the drive).")
