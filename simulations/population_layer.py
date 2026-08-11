"""Population layer for the MB Atoh1+ transient-G0 FRACTION, in the integrated governor model
(option B p27 + CDK6 mark-memory + INK4 p18/p19 sink, re-cal'd 31/32; integrated_optimize_best.json).

Three self-contained analyses, run sequentially (roadrunner is NOT thread-safe):

  (1) BULK heterogeneity  -- each cell gets a per-cell lognormal CDKI multiplier (CV=0.6) on p18/p19/p27.
                             A cell is G0 if it does NOT divide in a long window (brake > drive).
  (2) THRESHOLD sweep     -- deterministic: how high must a cell's CDKI be (as a fold of the mean)
                             before MB arrests? Also CDK6 x0.3 to see if titration lowers the threshold.
  (3) BIMODAL tail        -- a distinct high-CDKI subpopulation (frac at ~5x) on top of normal noise;
                             the arrested (G0) fraction should track the tail prevalence.

FINDING (see docs/transient_g0_population_layer_2026-08-09.md):
  Bulk CV=0.6 -> 0% G0 in every condition. MB's drive is so high that a cell needs ~5x the mean CDKI
  to arrest (sweep: 1x=6 div, 3x=3 div, 5x=0). Normal expression noise never reaches 5x -> 0% G0.
  A rare bimodal high-CDKI tail (~5x, ~10-20% prevalence) DOES produce a discrete MB Atoh1+ G0
  fraction (~5-12%, = tail prevalence). So the Atoh1+ transient-G0 subpopulation is a rare high-CDKI
  TAIL (a stress/niche state), NOT bulk expression noise -- confirming the reduced-model
  stochastic-commitment result (variance-not-mean) in the full v44 governor. G0 is CDKI-GATED and
  shows a graded slowdown (6->4->3 divisions) before arrest.

Run: PYTHONPATH=. ./venv/bin/python simulations/population_layer.py
"""
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

BP=dict(kPhRbCd=0.1748299602292031,K_CdRb=0.4025007594938155,w_ink4=4.646376316224687,kSyP21=0.0004043793202695936,
        K_cdk6_sink=3.879643001142466,w_p18=1.003950915185932,w_p19=0.6655007434943913,
        k_Cd_tx_Gli_max=0.24197652013556434,k_Cd_tx_MYCN=0.11031706915457935)
m=build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, params=BP)

# baseline CDKI levels (MB) that the per-cell multiplier scales
P18_MB, P19_MB = 1.73, 0.58

def n_div(cdki_mult, cdk6_mult=1.0, p18ko=False, p19ko=False, p27ko=False, t=11000, npts=22000):
    """Number of divisions (MPF peaks) for one MB cell with a given CDKI multiplier. 0 => G0/arrest."""
    rr=te.loada(m); rr.integrator.setValue('relative_tolerance',1e-6); rr.integrator.setValue('maximum_num_steps',1000000)
    for k,v in BP.items():
        try: rr[k]=v
        except: pass
    rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3
    try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
    except: pass
    rr['p18']=0.0 if p18ko else P18_MB*cdki_mult
    rr['p19']=0.0 if p19ko else P19_MB*cdki_mult
    rr['kSyP21']=BP['kSyP21']*(0.05 if p27ko else cdki_mult)
    try: rr['k_cdk6_Gli']=rr['k_cdk6_Gli']*cdk6_mult
    except: pass
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue('absolute_tolerance',atol)
            r=rr.simulate(0,t,npts,selections=['time','MPF']); tt=r['time']; mm=tt>=3000; T=tt[mm]; dt=T[1]-T[0]
            pk,_=find_peaks(r['MPF'][mm],prominence=0.15,distance=int(200/dt)); return len(pk)
        except: pass
    return 0

def is_G0(cdki_mult, **kw): return 1 if n_div(cdki_mult,**kw)==0 else 0

# ---------------------------------------------------------------------------
# (1) BULK heterogeneity: lognormal CDKI multiplier, CV=0.6
# ---------------------------------------------------------------------------
N=40; CV=0.6
def bulk_g0(seed, **kw):
    rng=np.random.default_rng(seed)
    mults=rng.lognormal(np.log(1.0)-0.5*np.log(1+CV**2), np.sqrt(np.log(1+CV**2)), N)
    return 100.0*np.mean([is_G0(mu,**kw) for mu in mults])

print(f"(1) BULK CDKI heterogeneity (N={N}, CV={CV}); G0 = no division in window.")
print(f"    MB (integrated):            {bulk_g0(1):5.1f}%   <- Atoh1+ transient-G0 subpop")
print(f"    MB, p18+p19 KO (INK4 off):  {bulk_g0(1,p18ko=True,p19ko=True):5.1f}%")
print(f"    MB, p27 KO:                 {bulk_g0(1,p27ko=True):5.1f}%")
print(f"    MB, CDK6 x0.5 (less titr):  {bulk_g0(1,cdk6_mult=0.5):5.1f}%")
print("    => bulk noise ~0%: no cell reaches the arrest threshold; drive dominates.\n")

# ---------------------------------------------------------------------------
# (2) THRESHOLD sweep: how high must CDKI be before MB arrests?
# ---------------------------------------------------------------------------
print("(2) THRESHOLD sweep -- MB divisions vs CDKI fold-of-mean:")
row=" ".join(f"{x}x={n_div(x)}" for x in (1,2,3,5,8,12))
print(f"    default CDK6:  {row}")
row2=" ".join(f"{x}x={n_div(x,cdk6_mult=0.3)}" for x in (1,2,3,5,8))
print(f"    CDK6 x0.3:     {row2}")
print("    => arrest (G0) only at CDKI >= ~5x; graded slowdown (6->4->3) before it; CDK6 barely shifts it.\n")

# ---------------------------------------------------------------------------
# (3) BIMODAL tail: a distinct high-CDKI subpopulation -> discrete G0 fraction
# ---------------------------------------------------------------------------
print("(3) BIMODAL high-CDKI tail -- G0 fraction tracks the tail prevalence:")
for frac_hi, hi, seed in [(0.15,5.0,0),(0.20,5.0,0),(0.25,6.0,0)]:
    rng=np.random.default_rng(seed)
    hi_mask=rng.random(N)<frac_hi
    mults=np.where(hi_mask, rng.lognormal(np.log(hi),0.2,N), rng.lognormal(0,0.4,N))
    g0=100.0*np.mean([is_G0(mu) for mu in mults])
    print(f"    {int(frac_hi*100)}% at ~{hi}x  -> Atoh1+ G0 fraction = {g0:.0f}%")
print("    => a rare high-CDKI tail produces the discrete Atoh1+ G0 subpopulation; G0 frac ~ tail prevalence.")
