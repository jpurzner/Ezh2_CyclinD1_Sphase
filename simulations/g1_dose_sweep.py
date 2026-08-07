"""#2 Variable GNP G1 by mitogen dose. Sweep SHH for GNP; report period + G1/S/G2 DURATIONS (h) so we can see
whether lowering mitogen lengthens G1 GRADUALLY while S/G2 stay ~flat (Fan-Meyer shape), vs the current sharp switch.
Mechanism = mu_eff growth-slowdown (mu_min_frac + (1-mu_min_frac)*Cd^n/(K_g1len^n+Cd^n)). Compare param sets.
Usage: PYTHONPATH=. ./venv/bin/python simulations/g1_dose_sweep.py '<label>' '<json params or {}>'
"""
import sys, json, numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

LABEL = sys.argv[1] if len(sys.argv)>1 else 'baseline (latent)'
PARAMS = json.loads(sys.argv[2]) if len(sys.argv)>2 else {}
MODEL = build_model_v44(params=(PARAMS or None))
print(f"===== {LABEL}  params={PARAMS or '(baked)'} =====")

def phase_durations(shh, t=24000, npts=48000):
    """Return (ndiv, period_h, G1_h, S_h, G2M_h, meanCd) for GNP at this SHH."""
    rr=te.loada(MODEL); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    last=None
    for msx,atol in [(1e9,1e-8),(20.0,1e-7),(5.0,1e-6)]:
        rr.reset()
        try:
            rr.integrator.setValue("maximum_time_step",msx); rr.integrator.setValue("absolute_tolerance",atol)
            r=rr.simulate(0,t,npts,selections=["time","MPF","Dna","aRc","vfork","Cd"]); break
        except Exception as e: last=e; r=None
    if r is None: return (-1,float('nan'),float('nan'),float('nan'),float('nan'),float('nan'))
    tt=r['time']; m=tt>=4000; T=tt[m]; dt=T[1]-T[0]
    pk,_=find_peaks(r['MPF'][m],prominence=0.12,distance=int(180/dt))
    per = np.median(np.diff(T[pk]/60.0)) if len(pk)>1 else float('nan')
    # phase time-fractions over the settled window
    syn = r['vfork'][m]*r['aRc'][m]; Dna=r['Dna'][m]
    in_S = (syn>0.02)&(Dna<0.98); in_G2=(Dna>=0.98); preS=~in_S&~in_G2
    fS=in_S.mean(); fG2=in_G2.mean(); fpre=preS.mean(); tot=fS+fG2+fpre or 1
    fS,fG2,fpre=fS/tot,fG2/tot,fpre/tot
    return (len(pk), per, per*fpre, per*fS, per*fG2, float(np.mean(r['Cd'][m])))

print(f"  {'SHH':>5s} {'ndiv':>5s} {'period':>7s} {'G1(h)':>6s} {'S(h)':>6s} {'G2M(h)':>7s} {'<Cd>':>6s}")
for shh in [0.50,0.42,0.35,0.30,0.25,0.20,0.15,0.12,0.10,0.08]:
    n,per,g1,s,g2,cd=phase_durations(shh)
    if n<=0 or per!=per:
        print(f"  {shh:5.2f} {'arrest' if n==0 else n:>5} {'--':>7} {'--':>6} {'--':>6} {'--':>7} {cd:6.2f}")
    else:
        print(f"  {shh:5.2f} {n:5d} {per:7.1f} {g1:6.1f} {s:6.1f} {g2:7.1f} {cd:6.2f}")
print("  Want: G1(h) rises smoothly as SHH falls (graded), S/G2 ~flat, and SHH=0.5 period ~16-18h (validation).")
