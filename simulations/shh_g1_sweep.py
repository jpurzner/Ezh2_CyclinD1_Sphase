"""Does p27 tune GNP G1 duration NEAR threshold? Fine SHH sweep, GNP period, for:
  (A) option A (baked, p27 inert)
  (B) option B (p27-inhibitory, 31/32 combo)
  (B0) option B with p27 removed (kSyP21~0) — isolates p27's contribution
If option B shows a GRADED lengthening window (period rising before arrest) that option A lacks, and that window
NARROWS/vanishes when p27 is removed, that demonstrates 'p27 tunes G1 near threshold'.
Run: PYTHONPATH=. ./venv/bin/python simulations/shh_g1_sweep.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def model(optB): return build_model_v44(with_p27_optionB=optB,
    params=({'kSeqCd':0.8,'kSyP21':0.0006,'kPhRbCd':0.32} if optB else None))

def period(mdl, shh, ksy=None, t=20000, npts=40000):
    rr=te.loada(mdl); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    if ksy is not None: rr['kSyP21']=ksy
    last=None
    for ms,atol in [(1e9,1e-8),(20.0,1e-7),(5.0,1e-6),(1.0,1e-6)]:
        rr.reset()
        try:
            rr.integrator.setValue("maximum_time_step",ms); rr.integrator.setValue("absolute_tolerance",atol)
            res=rr.simulate(0,t,npts,selections=["time","MPF"])
            tt=res['time']; m=tt>=4000; ttt=tt[m]; dt=ttt[1]-ttt[0]
            pk,_=find_peaks(res['MPF'][m],prominence=0.12,distance=int(180/dt))
            per=np.median(np.diff(ttt[pk]/60.0)) if len(pk)>1 else float('nan')
            return len(pk), per
        except Exception as e: last=e
    return -1, float('nan')

SHH=[0.50,0.45,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.27,0.24,0.20,0.15,0.10]
mA, mB = model(False), model(True)
print("Fine SHH sweep — GNP median cycle period (h); ndiv in parens. Looking for a graded lengthening window before arrest.")
print(f"  {'SHH':>5s} | {'A: baked':>14s} | {'B: p27-inhib':>14s} | {'B0: B, p27~0':>14s}")
for shh in SHH:
    nA,pA=period(mA,shh)
    nB,pB=period(mB,shh)
    nB0,pB0=period(mB,shh,ksy=1e-7)
    def fmt(n,p): return (f"{p:5.1f}h ({n})" if n>1 and p==p else ("arrest" if n==0 else f"({n})"))
    print(f"  {shh:5.2f} | {fmt(nA,pA):>14s} | {fmt(nB,pB):>14s} | {fmt(nB0,pB0):>14s}")
print("\n  Read: if column B shows periods lengthening (e.g. 18->22->28h) across SHH before 'arrest', while A jumps")
print("  straight from ~18h to arrest, option B gives graded G1 lengthening. If B0 (p27 removed) arrests SHARPLY like A,")
print("  the graded window is p27-mediated -> 'p27 tunes GNP G1 near threshold' (Fan-Meyer).")
