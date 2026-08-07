"""Two-pronged EZH2 governor test: does adding the CDK6 arm strengthen (a) vismo-dependence and (b) the
EZH2i -> vismo-resistance prediction? Compare CDK6 OFF vs ON at MB, MB+vismo, MB+vismo+EZH2i.
The drive proxy = the CyclinD-CDK4/6 -> Rb Vmax term ~ kPhRbCd*cdk6*(Cd+w_Cd2*Cd2). We report divisions + <Cd>
+ <cdk6> (on) so we can see both arms fall under vismo and both de-repress under EZH2i.
Run: PYTHONPATH=. ./venv/bin/python simulations/cdk6_twoarm.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def run(model, hhi=0.0, ezh2i=0.0, t=12000, npts=24000):
    rr=te.loada(model); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58; rr['HHi']=hhi; rr['EZH2i']=ezh2i
    try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
    except Exception: pass
    last=None
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            sel=["time","MPF","Cd"]+(["cdk6"] if 'cdk6 :=' in model else [])
            return rr.simulate(0,t,npts,selections=sel)
        except Exception as e: last=e
    raise last
def ndiv(res,settle=3000):
    t=res['time']; m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(res['MPF'][m],prominence=0.15,distance=int(200/dt)); return len(pk)
def ms(res,sp,settle=4000):
    if sp not in res.colnames: return float('nan')
    t=res['time']; return float(np.mean(res[sp][t>=settle]))

for lab,optb in [('CDK6 OFF (1-arm: CyclinD1 only)',False),('CDK6 ON (2-arm: CyclinD1 + CDK6)',True)]:
    m=build_model_v44(with_cdk6_gli=optb)
    print(f"\n=== {lab} ===")
    print(f"  {'condition':20s} {'ndiv':>5s} {'<Cd>':>7s} {'<cdk6>':>7s}  drive-proxy Cd*cdk6")
    base=None
    for cond,hhi,ez in [('MB',0,0),('MB+vismo',1.0,0),('MB+vismo+EZH2i',1.0,1.0)]:
        r=run(m,hhi=hhi,ezh2i=ez); n=ndiv(r); cd=ms(r,'Cd'); c6=ms(r,'cdk6') if optb else 1.0
        proxy=cd*(c6 if c6==c6 else 1.0)
        print(f"  {cond:20s} {n:5d} {cd:7.2f} {(c6 if c6==c6 else float('nan')):7.3f}  {proxy:8.2f}")
print("\n  Reads:")
print("  - vismo should collapse the drive proxy MORE with CDK6 ON (both Cd AND cdk6 fall) = stronger 2-arm vismo-dependence.")
print("  - MB+vismo+EZH2i: EZH2i de-represses BOTH arms -> more drive restored with CDK6 ON = stronger resistance prediction.")
