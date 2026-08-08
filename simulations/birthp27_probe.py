"""#3 mechanism check: does a cell BORN with high p27 take a TRANSIENT G0 (delayed commitment -> re-entry) in the
full v44 model? Set initial p27 (P21) high with P21_div kept low (a ONE-TIME high-p27 birth), measure time-to-first-
division (the transient-G0 length) + whether it resumes normal cycling. Test option A (default, p27 buffered/inert)
vs option B (p27-inhibitory), GNP vs MB. If high birth-p27 -> long delay then normal period, transient G0 works.
Run: PYTHONPATH=. ./venv/bin/python simulations/birthp27_probe.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

def run(model, cell, birth_p27, t=16000, npts=32000):
    rr=te.loada(model); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    if cell=='GNP': rr['SHH']=0.5; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    else:
        rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass
    rr['P21_div']=0.6                 # subsequent births normal -> isolate the ONE-TIME high-p27 birth
    rr['P21']=birth_p27               # this cell is BORN with high p27
    last=None
    for atol in (1e-8,1e-7,1e-6):
        rr.reset(); rr['P21']=birth_p27
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            return rr.simulate(0,t,npts,selections=["time","MPF","P21","pRb"])
        except Exception as e: last=e
    raise last

def analyze(res, settle_ignore=200):
    t=res['time']; dt=t[1]-t[0]
    pk,_=find_peaks(res['MPF'], prominence=0.15, distance=int(200/dt))
    tpk=t[pk]/60.0
    first=tpk[0] if len(tpk) else float('nan')          # h to 1st division (transient-G0 length if delayed)
    later=np.median(np.diff(tpk)) if len(tpk)>2 else float('nan')  # steady period after re-entry
    return len(pk), first, later

for lab,optb in [('OPTION A (default)',False),('OPTION B (p27-inhibitory)',True)]:
    m=build_model_v44(with_p27_optionB=optb)
    print(f"\n=== {lab} ===")
    print(f"  {'cell':5s} {'birthp27':>8s} {'ndiv':>5s} {'t_1st_div(h)':>12s} {'later_period(h)':>15s}  transient?")
    for cell in ['GNP','MB']:
        base=None
        for bp in [0.6,1.5,2.5,3.5,5.0]:
            try:
                n,first,later=analyze(run(m,cell,bp))
                if base is None and bp==0.6: base=first
                verdict=""
                if n>=2 and first==first and later==later:
                    if base and first>base*1.5: verdict="TRANSIENT G0 (delayed then re-enters)"
                    else: verdict="cycles (little delay)"
                elif n<=1: verdict="ARREST (no re-entry)"
                print(f"  {cell:5s} {bp:8.1f} {n:5d} {first:12.1f} {later:15.1f}  {verdict}")
            except Exception as e:
                print(f"  {cell:5s} {bp:8.1f}  ERROR {str(e)[:40]}")
print("\n  Want: high birth-p27 -> t_1st_div rises (transient G0) then later_period ~normal (re-entry). Then a")
print("  POPULATION with birth-p27 heterogeneity (MB tail) gives a G0 FRACTION. If option A shows no delay, p27 is")
print("  inert there -> the p27-driven G0 needs option B.")
