"""Amplifier probe (JP 2026-08-09 hypothesis): does the mark-accumulation feedback
(out-of-cycle -> less H3K27me3 dilution -> Mk banks -> represses CDK6 + CyclinD1 -> drive erodes)
act as a FATE-changer that lowers the CDKI arrest threshold and produces a compounding transient-G0
dwell -- while keeping baseline CDK6/CyclinD1 (hence the 31/32 calibration) preserved?

Integrated governor: with_p27_optionB + with_cdk6_gli, best params (integrated_optimize_best.json).
Layer candidate amplifier params on top via env AMP_PARAMS (JSON). Prints a JSON metric bundle.

REFERENCE baseline (base K, 1x CDKI, SHH 0.5, MB): cdk6~5.084, Cd~7.71, MB proliferates,
arrest threshold ~5x (SHH0.5). The amplifier is currently SATURATED (cdk6 flat, Cd -3% at 3x).
Goal of a candidate: baseline cdk6/Cd MATCHED, big dynamic gap (cdk6/Cd drop at 3x), arrest
threshold LOWERED (toward ~2-3x reachable by CDKI variance), GNP still Shh-dependent (SHH0 -> 0 div).

Run:  AMP_PARAMS='{"K_prc2_cdk6":0.035,"n_prc2_cdk6":6}' PYTHONPATH=. ./venv/bin/python simulations/amplifier_probe.py
"""
import os, json, numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

BP=dict(kPhRbCd=0.1748299602292031,K_CdRb=0.4025007594938155,w_ink4=4.646376316224687,kSyP21=0.0004043793202695936,
        K_cdk6_sink=3.879643001142466,w_p18=1.003950915185932,w_p19=0.6655007434943913,
        k_Cd_tx_Gli_max=0.24197652013556434,k_Cd_tx_MYCN=0.11031706915457935)
# baked defaults the auto-match scales (from _ts_bake / builder antimony)
CDK6_BAS0, CDK6_GLI0 = 0.0002, 0.0188
CDTX_BASAL0 = 0.000529
TARGET_CDK6, TARGET_CD = 5.084, 7.707     # reference baseline (base K, 1x CDKI, SHH 0.5, MB)
AUTO_MATCH = os.environ.get('AUTO_MATCH','1')=='1'
AMP=json.loads(os.environ.get('AMP_PARAMS','') or '{}')
P=dict(BP); P.update(AMP)

# Build & COMPILE the model ONCE; all amplifier params (K/n) + matched scales + conditions are set at
# runtime via rr[k]=v (global params are runtime-settable), avoiding a te.loada recompile per sim (~1-2s each).
M=build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, params=BP)
RR=te.loada(M)
RR.integrator.setValue('relative_tolerance',1e-7); RR.integrator.setValue('absolute_tolerance',1e-9)
RR.integrator.setValue('maximum_num_steps',8000000)
try: CD2_DEFAULT=RR['Cd2_expr']
except Exception: CD2_DEFAULT=None

def sim(mult, shh, gnp=False, T=30000, sel=('time','Dna','Mk','cdk6','Cd'), model=None):
    rr=RR
    npts=max(3000, T//10)   # ~10-min resolution: enough for Dna-based division detection, ~4x faster than T points
    for k,v in P.items():
        try: rr[k]=v
        except: pass
    if gnp:
        rr['SHH']=shh; rr['Ptch1_copy_number']=1.0; p18,p19=0.464,0.36
        if CD2_DEFAULT is not None: rr['Cd2_expr']=CD2_DEFAULT   # restore GNP D2 (reset() keeps params; avoid MB leak)
    else:
        rr['SHH']=shh; rr['Ptch1_copy_number']=0.3; p18,p19=1.73,0.58
        try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
        except: pass
    rr['p18']=p18*mult; rr['p19']=p19*mult; rr['kSyP21']=BP['kSyP21']*mult
    rr.reset()
    r=rr.simulate(0,T,npts,selections=list(sel))
    return r

def divisions(r):
    dna=r['Dna']; return np.where((dna[:-1]>0.9)&(dna[1:]<0.1))[0]

def steady(mult, shh, gnp=False):
    r=sim(mult,shh,gnp); d=divisions(r)
    if len(d)>=3:
        t=r['time']; return len(d),(t[d[-1]]-t[d[-2]])/60, r['Mk'][d[-1]], r['cdk6'][d[-1]], r['Cd'][d[-1]]
    return len(d), float('nan'), float('nan'), float('nan'), float('nan')

def arrest_threshold(shh):
    for mx in (2,3,4,5,6):
        if len(divisions(sim(mx,shh)))==0: return mx
    return 99

def compounding(mult, shh):
    """cycle-length progression: ratio of 5th to 2nd inter-division interval (>1 = deepening)."""
    r=sim(mult,shh); d=divisions(r); t=r['time']
    if len(d)>=6:
        iv=np.diff(t[d])/60
        return float(iv[4]/iv[1]) if iv[1]>0 else float('nan'), float(iv[-1]/iv[0]) if iv[0]>0 else float('nan')
    return float('nan'), float('nan')

# --- auto-match baseline cdk6 & Cd (scale synthesis/transcription) so a candidate's higher K/n
#     does not confound the threshold via a shifted baseline. Fixed-point (coupled via PRC2_rep). ---
if AUTO_MATCH:
    P.setdefault('k_cdk6_bas',CDK6_BAS0); P.setdefault('k_cdk6_Gli',CDK6_GLI0)
    P.setdefault('k_Cd_tx_Gli_max',BP['k_Cd_tx_Gli_max']); P.setdefault('k_Cd_tx_MYCN',BP['k_Cd_tx_MYCN']); P.setdefault('k_Cd_tx_basal',CDTX_BASAL0)
    for _ in range(4):
        b=steady(1.0,0.5)
        if b[3]==b[3] and b[3]>0:
            f=TARGET_CDK6/b[3]; P['k_cdk6_bas']*=f; P['k_cdk6_Gli']*=f
        if b[4]==b[4] and b[4]>0:
            g=TARGET_CD/b[4]; P['k_Cd_tx_Gli_max']*=g; P['k_Cd_tx_MYCN']*=g; P['k_Cd_tx_basal']*=g

# --- baseline (1x) and dynamic (3x) at SHH 0.5 ---
b1=steady(1.0,0.5); b3=steady(3.0,0.5)
cdk6_drop = 100*(b1[3]-b3[3])/b1[3] if b1[3]==b1[3] and b1[3] else float('nan')
cd_drop   = 100*(b1[4]-b3[4])/b1[4] if b1[4]==b1[4] and b1[4] else float('nan')
comp5,compL = compounding(3.0,0.5)

out=dict(
  amp_params=AMP,
  baseline_cdk6=round(b1[3],3), baseline_Cd=round(b1[4],3), baseline_cyc_h=round(b1[1],2), baseline_ndiv=int(b1[0]),
  slow_cdk6=round(b3[3],3), slow_Cd=round(b3[4],3), slow_cyc_h=round(b3[1],2),
  cdk6_drop_pct=round(cdk6_drop,1), Cd_drop_pct=round(cd_drop,1),
  arrest_thr_shh050=arrest_threshold(0.5), arrest_thr_shh025=arrest_threshold(0.25),
  compound_5v2=round(comp5,3) if comp5==comp5 else None, compound_last_v_first=round(compL,3) if compL==compL else None,
  gnp_ndiv_shh05=int(len(divisions(sim(1.0,0.5,gnp=True)))),
  gnp_ndiv_shh0=int(len(divisions(sim(1.0,0.0,gnp=True)))),
  resolved_params={k:(round(v,8) if isinstance(v,float) else v) for k,v in P.items()},
)
print(json.dumps(out))
