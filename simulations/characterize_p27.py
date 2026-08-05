"""Characterize the p27 mechanism under option A (baked, inert p27) vs option B (Fan-Meyer p27-inhibitory).
Two decisive tests for JP's goals:
  (1) G1-DURATION vs mitogen (SHH sweep, GNP): does the cycle period lengthen GRADUALLY as mitogen falls, and is
      that lengthening BIGGER under option B (p27 now brakes CDK4/6 via the CyclinD1/p27 ratio)?  -> "p27 tunes GNP G1"
  (2) REVERSIBILITY (CDK4/6i washout, GNP vs MB): block CDK4/6 for a window, restore, measure restart latency.
      Does MB restart FASTER (large re-partitionable CyclinD-CDK4/6 buffer + continuous drive)?  -> MB-specific G0 axis
Usage: PYTHONPATH=. ./venv/bin/python simulations/characterize_p27.py            # option A (baked)
       PYTHONPATH=. ./venv/bin/python simulations/characterize_p27.py B '<json params>'   # option B with recal params
"""
import sys, json, numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

OPTB = len(sys.argv)>1 and sys.argv[1].upper()=='B'
PARAMS = json.loads(sys.argv[2]) if len(sys.argv)>2 else {}
MODEL = build_model_v44(with_p27_optionB=OPTB, params=(PARAMS or None))
LABEL = "OPTION B (p27-inhibitory)" if OPTB else "OPTION A (baked)"
print(f"===== {LABEL}  params={PARAMS if PARAMS else '(baked defaults)'} =====")

def rr_new():
    rr = te.loada(MODEL)
    rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr.integrator.setValue("absolute_tolerance",1e-8)
    return rr

def setcell(rr, cond):
    if cond=='GNP': rr['SHH']=0.5; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    elif cond=='MB':
        rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3; rr['p18']=1.73; rr['p19']=0.58
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except Exception: pass

def peaks_h(t, mpf, settle=0):
    m=t>=settle; tt=t[m]; dt=tt[1]-tt[0]
    pk,_=find_peaks(mpf[m], prominence=0.15, distance=int(200/dt))
    return tt[pk]/60.0

# ---------- Test 1: period vs SHH (GNP) ----------
print("\n--- Test 1: cycle period (h) vs SHH (GNP) — graded G1 lengthening? ---")
print(f"  {'SHH':>5s} {'ndiv':>5s} {'median period (h)':>18s}")
for shh in [0.5,0.35,0.25,0.18,0.12,0.08,0.05]:
    rr=rr_new(); setcell(rr,'GNP'); rr['SHH']=shh
    try:
        res=rr.simulate(0,16000,32000,selections=["time","MPF"])
        pk=peaks_h(res['time'],res['MPF'],settle=3000)
        per=np.median(np.diff(pk)) if len(pk)>1 else float('nan')
        print(f"  {shh:5.2f} {len(pk):5d} {per:18.1f}")
    except Exception as e:
        print(f"  {shh:5.2f}  ERROR {str(e)[:50]}")

# ---------- Test 2: CDK4/6i washout restart latency (GNP vs MB) ----------
print("\n--- Test 2: CDK4/6i washout — restart latency after release (reversibility) ---")
print("  protocol: settle 4000 min -> CDK4/6i (kPhRbCd=0) for 3000 min -> release -> measure time to 1st division")
print(f"  {'cond':6s} {'divs during block':>18s} {'restart latency (h)':>20s} {'<CdP21 during block>':>21s}")
for cond in ['GNP','MB']:
    rr=rr_new(); setcell(rr,cond); k0=rr['kPhRbCd']
    try:
        r1=rr.simulate(0,4000,8000,selections=["time","MPF"])
        rr['kPhRbCd']=0.0
        r2=rr.simulate(4000,7000,6000,selections=["time","MPF","CdP21"])
        rr['kPhRbCd']=k0
        r3=rr.simulate(7000,13000,12000,selections=["time","MPF"])
        blk=peaks_h(r2['time'],r2['MPF'])                       # divisions during the block (want ~0)
        pk3=peaks_h(r3['time'],r3['MPF'])
        latency=(pk3[0]-7000/60.0) if len(pk3) else float('nan') # h from release to first division
        cdp21=float(np.mean(r2['CdP21'])) if 'CdP21' in r2.colnames else float('nan')
        print(f"  {cond:6s} {len(blk):18d} {latency:20.1f} {cdp21:21.4f}")
    except Exception as e:
        print(f"  {cond:6s}  ERROR {str(e)[:60]}")
print("\n  (JP hypothesis: MB restart latency < GNP -> MB is exit-competent / GNP is not.)")
