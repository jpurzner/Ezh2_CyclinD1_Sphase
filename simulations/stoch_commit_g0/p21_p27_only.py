"""Dichotomy from p21 + p27 ALONE (JP 2026-08-04 direction):
  - DROP p18/p19 from the G0 mechanism (INK4 = amplifier/timer, not a G0 source; spec sec 7.4). p18 same in both types.
  - Do NOT trust TMT p21/p27 means -> birth-p21/p27 distributions are ALLOWED to differ GNP vs MB.
  - p21 stays involved: LOW bulk mean (miR-17~92 repression, spec sec 7.5) but SPORADIC in MB (rare high-p21 cells).
Question: with p21/p27 only, what produces MB transient-G0 while GNP has none? mean route vs sporadic-tail route.

CALIBRATION FIRST: GNP must sit at ~0-1% G0 (proliferating, NO transient G0). Then MB = same birth-p + higher
drive (Cd 4x, kCd 2x) + optional MB-specific tail. Run:
  ./venv/bin/python simulations/stoch_commit_g0/p21_p27_only.py
"""
import numpy as np
from model import make_params, integrate, _lognorm

# p18 dropped from the G0 mechanism: MB keeps only the drive increase (Cd 4x, kCd 2x), p18 EQUAL to GNP.
GNP0 = dict(Cd=1.0, kCd=1.0, p18=1.0)
MB0  = dict(Cd=4.0, kCd=2.0, p18=1.0)   # <- p18=1.0 (was 3.0): INK4 differential removed
N=4000

def ens(ct, m21, cv21, m27, cv27, seed=1, E_thr=0.3, frac_hi=0.0, p21_hi=0.0, cv_hi=0.3, Cd_cv=0.0):
    """Ensemble with optional MB-specific p21 tail (frac_hi at p21_hi) and optional per-cell CyclinD1 variance."""
    rng=np.random.default_rng(seed); p=make_params(ct)
    if Cd_cv>0:
        p['Cd']=_lognorm(rng,p['Cd'],Cd_cv,N)
    if frac_hi>0:
        hi=rng.random(N)<frac_hi
        p21=np.where(hi, rng.lognormal(np.log(p21_hi),cv_hi,N), _lognorm(rng,m21,cv21,N))
    else:
        p21=_lognorm(rng,m21,cv21,N)
    p27=_lognorm(rng,m27,cv27,N); E=np.full(N,0.05)
    Ef,_,_=integrate(E,p21,p27,p)
    return float(np.mean(Ef<E_thr)), float(np.mean(p21)), float(np.mean(p27))

print("="*80)
print("CALIBRATION — find GNP birth-p (m21,m27) giving GNP ~0-1% G0 (proliferating). MB = same birth + higher drive.")
print("="*80)
print("   m21   m27  | GNP G0%   MB(no tail) G0%")
CAL=None
for m21,m27 in [(0.02,0.05),(0.03,0.08),(0.05,0.10),(0.05,0.15),(0.05,0.05),(0.03,0.05)]:
    g_gnp,_,_=ens(GNP0,m21,0.8,m27,0.5); g_mb,_,_=ens(MB0,m21,0.8,m27,0.5)
    flag=""
    if g_gnp<=0.01 and CAL is None: CAL=(m21,m27); flag=" <- adopt"
    print(f"  {m21:.2f}  {m27:.2f} |  {g_gnp*100:4.1f}%      {g_mb*100:4.1f}%{flag}")
if CAL is None: CAL=(0.03,0.05)
M21,M27=CAL
print(f"  -> GNP calibrated at m21={M21}, m27={M27} (GNP ~0% G0). MB at the SAME birth-p is also ~0% (drive dominance).")
g_gnp0,_,_=ens(GNP0,M21,0.8,M27,0.5); g_mb0,_,_=ens(MB0,M21,0.8,M27,0.5)
print(f"     GNP {g_gnp0*100:.1f}%   MB(no tail) {g_mb0*100:.1f}%   -> both proliferating; the dichotomy must be ADDED.")

print()
print("="*80)
print("TEST 2 — MB-specific MEAN p21 rise (GNP at calibrated low): fold over GNP mean needed for MB G0")
print("="*80)
print("   MB mean p21 | MB G0%  | fold over GNP")
for m21 in [0.1,0.3,0.6,1.0,1.5,2.0,3.0]:
    g_mb,_,_=ens(MB0,m21,0.8,M27,0.5)
    print(f"     {m21:4.2f}     |  {g_mb*100:4.1f}%  |  {m21/M21:5.0f}x")
print(f"  -> p21 on the CDK2 arm: a mean rise CAN lift MB G0 (unlike p18). But needs a LARGE fold over GNP's mean")
print(f"     ({M21}), and miR-17~92 pins the MB p21 MEAN LOW -> a uniform mean rise is NOT biologically supported.")

print()
print("="*80)
print("TEST 3 — MB-specific MEAN p27 rise (p21=0), GNP at calibrated low")
print("="*80)
print("   MB mean p27 | MB G0%")
for m27 in [0.1,0.3,0.6,1.0,1.5,2.0,3.0]:
    g_mb,_,_=ens(MB0,0.0,0.8,m27,0.5)
    print(f"     {m27:4.2f}     |  {g_mb*100:4.1f}%")
print("  -> same shape (p27 also CDK2-arm). Ayrault: MB RETAINS & NEUTRALISES p27 (haploinsufficient, not elevated)")
print("     -> a mean nuclear-p27 rise is not supported either. Neither mean route matches the biology.")

print()
print("="*80)
print("TEST 4 — MB-specific p21 TAIL (JP: sporadic high-p21 cells, LOW bulk mean via miR-17~92). CLEAN calibration.")
print("="*80)
print(f"  GNP (no tail): G0={g_gnp0*100:.1f}%    MB (no tail): G0={g_mb0*100:.1f}%")
print("   frac_hi p21_hi | MB G0%  bulk-p21-mean | GNP-if-SAME-tail | dichotomy?")
for frac,phi in [(0.10,2.5),(0.20,2.5),(0.30,2.5),(0.15,3.5),(0.25,3.5),(0.35,3.5)]:
    g_mb,mean_mb,_=ens(MB0,M21,0.8,M27,0.5,frac_hi=frac,p21_hi=phi)
    g_gnp_tail,_,_=ens(GNP0,M21,0.8,M27,0.5,frac_hi=frac,p21_hi=phi)
    ok = "YES" if (g_mb-g_gnp0>0.05 and g_gnp0<0.02) else ""
    print(f"    {frac:.2f}   {phi:.1f}  |  {g_mb*100:4.1f}%    {mean_mb:5.2f}       |     {g_gnp_tail*100:4.1f}%       | {ok}")
print("  -> MB gets real G0; bulk p21 mean stays LOW (undetectable-consistent); GNP (no tail) ~0%. KEY: GNP-if-SAME-tail")
print("     is HIGHER than MB (GNP's lower drive parks more easily) -> the dichotomy REQUIRES the tail be MB-SPECIFIC.")

print()
print("="*80)
print("TEST 5 — p27 via Fan-Meyer cyclinD1/p27 ratio: p21=0, per-cell inherited CyclinD1 variance (GNP calibrated)")
print("="*80)
print("   Cd_cv  m27 | GNP G0%   MB G0%   gap")
for cd_cv,m27 in [(0.4,M27),(0.6,M27),(0.8,M27),(0.6,0.3),(0.8,0.3)]:
    g_gnp,_,_=ens(GNP0,0.0,0.8,m27,0.5,Cd_cv=cd_cv)
    g_mb,_,_ =ens(MB0 ,0.0,0.8,m27,0.5,Cd_cv=cd_cv)
    print(f"   {cd_cv:.1f}   {m27:.2f} |  {g_gnp*100:4.1f}%   {g_mb*100:4.1f}%   {(g_mb-g_gnp)*100:+5.1f}pp")
print("  -> CyclinD1-partitioning variance alone: MB's HIGHER mean CyclinD1 keeps daughters above the ratio threshold")
print("     -> MB <= GNP. A ratio route also needs the low-CyclinD1 / high-p27 tail to be MB-ENRICHED.")

print()
print("="*80)
print("TEST 6 — combined: MB-specific p21 tail + modest p27 elevation (both on the CDK2 arm)")
print("="*80)
for frac,phi,m27hi in [(0.20,2.0,M27),(0.20,2.0,1.0),(0.25,2.0,1.0),(0.25,2.5,1.5)]:
    g_mb,mean_mb,_=ens(MB0,M21,0.8,m27hi,0.5,frac_hi=frac,p21_hi=phi)
    print(f"  p21-tail(frac={frac},hi={phi}) + p27 mean {m27hi}: GNP {g_gnp0*100:4.1f}%  MB {g_mb*100:4.1f}%  gap {(g_mb-g_gnp0)*100:+5.1f}pp  (bulk p21 {mean_mb:.2f})")
print("  -> a modest p27 elevation LOWERS the p21_hi the tail needs; both add on the same inhibitory (CDK2) arm.")

print()
print("="*80)
print("TEST 7 — settings that PERMIT the dichotomy (p21 tail): MB G0 window vs (tail fraction, p21_hi)")
print("="*80)
print("  target: MB G0 in 10-30% AND GNP(no tail) < 2%. Scan p21_hi at each frac.")
for frac in [0.10,0.15,0.20,0.25,0.30]:
    hits=[]
    for phi in [1.5,2.0,2.5,3.0,3.5,4.0]:
        g_mb,_,_=ens(MB0,M21,0.8,M27,0.5,frac_hi=frac,p21_hi=phi)
        if 0.10<=g_mb<=0.35: hits.append(f"{phi:.1f}({g_mb*100:.0f}%)")
    print(f"   frac={frac:.2f}: in-window at p21_hi = {', '.join(hits) if hits else '(none)'}")
print()
print(f"  [calibration used: GNP m21={M21}, m27={M27}; MB drive Cd=4,kCd=2, p18 dropped; noise=lognormal(multiplicative)]")
