"""Coarse option-B regime search: find params where GNP+SHH cycles (~16h), GNP-SHH/HHi arrests (Shh-dependence),
and MB cycles (~24h), under p27-inhibitory CDK4/6 (option B). Levers: kSeqCd (buffer/brake strength), kSyP21 (p27
amount), kPhRbCd (drive), K_CdRb (drive half-max). Reports feasible combos. Run in background.
  PYTHONPATH=. ./venv/bin/python simulations/optionB_recal_coarse.py
"""
import numpy as np, itertools, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def run(model, shh, ptch1, p18, p19, cd2_mb, hhi=0.0, serum=False, t=11000, npts=22000):
    rr = te.loada(model)
    rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=ptch1; rr['p18']=p18; rr['p19']=p19; rr['HHi']=hhi
    if cd2_mb:
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    if serum:
        rr['k_Cd_translation']=0.0
    last=None
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            return rr.simulate(0,t,npts,selections=["time","MPF"])
        except Exception as e: last=e
    raise last

def ndiv(res, settle=3000):
    t=res['time']; m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(res['MPF'][m], prominence=0.15, distance=int(200/dt))
    per=np.diff(tt[pk]/60.0) if len(pk)>1 else np.array([])
    return len(pk), (float(np.median(per)) if len(per) else float('nan'))

# baked option-A values (starting point)
BASE=dict(kSeqCd=2.13243, kRelCd=0.01647, kDeKPC=0.02022, w_ink4=3.6204, kSyP21=0.00116, K_CdRb=0.4032, kPhRbCd=0.1403)

grid=dict(
    kSeqCd=[0.3,0.8,1.5],       # lower -> less p27 buffered on Cd -> weaker option-B brake
    kSyP21=[0.0006,0.00116],    # p27 amount
    kPhRbCd=[0.1403,0.22,0.32], # stronger drive helps GNP cycle
)
keys=list(grid.keys())
print("Coarse option-B search. Want: GNP+SHH div>=3 & period 14-19h; GNP-SHH=0; GNP+HHi=0; serum=0; MB div>=3.")
print(f"{'kSeqCd':>7s} {'kSyP21':>8s} {'kPhRbCd':>8s} | {'GNPd':>4s} {'per':>5s} {'GNP-':>4s} {'HHi':>4s} {'ser':>4s} {'MBd':>4s} {'MBper':>5s} | verdict")
feasible=[]
for combo in itertools.product(*[grid[k] for k in keys]):
    p=dict(BASE);
    for k,v in zip(keys,combo): p[k]=v
    model=build_model_v44(with_p27_optionB=True, params=p)
    try:
        gd,gper=ndiv(run(model,0.5,1.0,0.464,0.36,False))
        gm,_  =ndiv(run(model,0.0,1.0,0.464,0.36,False))          # GNP -SHH
        gh,_  =ndiv(run(model,0.5,1.0,0.464,0.36,False,hhi=1.0))  # GNP +HHi
        gs,_  =ndiv(run(model,0.5,1.0,0.464,0.36,False,serum=True))
        md,mper=ndiv(run(model,0.5,0.3,1.73,0.58,True))
        ok = (gd>=3 and 14<=gper<=19 and gm==0 and gh==0 and gs==0 and md>=3)
        v="FEASIBLE" if ok else ""
        if ok: feasible.append((p,gd,gper,md,mper))
        print(f"{combo[0]:7.2f} {combo[1]:8.5f} {combo[2]:8.3f} | {gd:4d} {gper:5.1f} {gm:4d} {gh:4d} {gs:4d} {md:4d} {mper:5.1f} | {v}")
    except Exception as e:
        print(f"{combo[0]:7.2f} {combo[1]:8.5f} {combo[2]:8.3f} |  ERROR {str(e)[:40]}")
print(f"\n{len(feasible)} feasible combo(s).")
for p,gd,gper,md,mper in feasible[:5]:
    print("  ", {k:round(p[k],5) for k in keys}, f"GNP {gd}d/{gper:.1f}h  MB {md}d/{mper:.1f}h")
