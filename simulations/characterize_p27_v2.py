"""Characterize p27 mechanism: option A (baked, p27 inert) vs option B (Fan-Meyer p27-inhibitory).
TEST 1 (the key one for JP's 'p27 tunes G1'): sweep p27 synthesis (kSyP21) at fixed GNP mitogen, measure GNP
  cycle period. Option A -> flat then abrupt arrest (p27 inert until it swamps CDK2). Option B -> GRADED G1
  lengthening as p27 rises (p27 brakes CDK4/6 via the CyclinD1/p27 ratio), then arrest.
TEST 2 (reversibility): CDK4/6i washout, GNP vs MB, robust restart detection (concatenated trace, long window).
Usage:
  PYTHONPATH=. ./venv/bin/python simulations/characterize_p27_v2.py A
  PYTHONPATH=. ./venv/bin/python simulations/characterize_p27_v2.py B '{"kSeqCd":0.8,"kSyP21":0.0006,"kPhRbCd":0.32}'
"""
import sys, json, numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

OPTB = len(sys.argv)>1 and sys.argv[1].upper()=='B'
PARAMS = json.loads(sys.argv[2]) if len(sys.argv)>2 else {}
MODEL = build_model_v44(with_p27_optionB=OPTB, params=(PARAMS or None))
print(f"===== {'OPTION B (p27-inhibitory)' if OPTB else 'OPTION A (baked)'}  params={PARAMS or '(baked)'} =====")
BASE_KSY = PARAMS.get('kSyP21', 0.00116)

def rr_new():
    rr = te.loada(MODEL)
    rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr.integrator.setValue("absolute_tolerance",1e-8)
    return rr
def gnp(rr): rr['SHH']=0.5; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
def mb(rr):
    rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58
    try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
    except Exception: pass
def peaks_h(t,mpf,settle=0,prom=0.12):
    m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(mpf[m],prominence=prom,distance=int(180/dt)); return tt[pk]/60.0

print("\n--- TEST 1: GNP cycle period vs p27 synthesis kSyP21 (fixed mitogen) — does p27 tune G1? ---")
print(f"  (baked/param kSyP21={BASE_KSY:g}; sweeping multiples)")
print(f"  {'kSyP21':>9s} {'xbase':>6s} {'ndiv':>5s} {'median period (h)':>18s}")
for mult in [0.5,1.0,1.5,2.0,3.0,4.0,6.0,8.0]:
    k=BASE_KSY*mult
    rr=rr_new(); gnp(rr); rr['kSyP21']=k
    try:
        res=rr.simulate(0,18000,36000,selections=["time","MPF"])
        pk=peaks_h(res['time'],res['MPF'],settle=3000)
        per=np.median(np.diff(pk)) if len(pk)>1 else float('nan')
        print(f"  {k:9.5f} {mult:6.1f} {len(pk):5d} {per:18.1f}")
    except Exception as e:
        print(f"  {k:9.5f} {mult:6.1f}  ERROR {str(e)[:40]}")

print("\n--- TEST 2: CDK4/6i washout — restart after release (reversibility), robust detection ---")
print("  settle 4000 -> CDK4/6i 3000 min -> release; report restart latency (h) + max MPF/Cd in 200h recovery window")
print(f"  {'cond':6s} {'blkdiv':>6s} {'restart lat(h)':>14s} {'maxMPF_rec':>10s} {'maxCd_rec':>9s}")
for cond,setter in [('GNP',gnp),('MB',mb)]:
    rr=rr_new(); setter(rr); k0=rr['kPhRbCd']
    try:
        s1=rr.simulate(0,4000,8000,selections=["time","MPF","Cd"])
        rr['kPhRbCd']=0.0
        s2=rr.simulate(4000,7000,6000,selections=["time","MPF","Cd"])
        rr['kPhRbCd']=k0
        s3=rr.simulate(7000,19000,24000,selections=["time","MPF","Cd"])
        # concatenate, detect peaks across whole trace, first peak after release (t>7000)
        t=np.concatenate([s1['time'],s2['time'],s3['time']]); mpf=np.concatenate([s1['MPF'],s2['MPF'],s3['MPF']])
        allpk=peaks_h(t,mpf,settle=0)
        blk=[p for p in allpk if 4000/60.0<p<7000/60.0]
        post=[p for p in allpk if p>7000/60.0]
        lat=(post[0]-7000/60.0) if post else float('nan')
        maxmpf=float(np.max(s3['MPF'])); maxcd=float(np.max(s3['Cd']))
        print(f"  {cond:6s} {len(blk):6d} {lat:14.1f} {maxmpf:10.3f} {maxcd:9.2f}")
    except Exception as e:
        print(f"  {cond:6s}  ERROR {str(e)[:60]}")
print("  (Current model has NO differentiation/DREAM arm -> expect BOTH restart, i.e. no GNP/MB reversibility")
print("   asymmetry yet. That asymmetry needs the DREAM/Neurod1 module = design-spec Phase B5.)")
