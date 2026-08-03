"""Can the GNP-no-G0 / MB-transient-G0 dichotomy arise from the CKI-vs-CyclinD1 competition?
Constraint (proteomics): GNP->MB differs ONLY in Cd (~4x), CDK4/6 (~2x, folded into kCd), p18 (~3x).
p27 is FLAT (same birth in GNP and MB); p21 is undetectable (low, same). So GNP and MB share birth-p.
"""
import numpy as np
from model import BASE, GNP, MB, make_params, integrate, run_ensemble, _lognorm

# ---- 1. bistability: separatrix scan (E_final vs birth-p27 at fixed low p21), GNP vs MB ----
print("=== 1. Separatrix scan: E_final vs birth-p27 (p21_birth=0.05). Sharp step = bistable. ===")
print("  p27_b :  GNP E_fin   MB E_fin")
for p27b in [0.0,0.1,0.2,0.3,0.4,0.5,0.8,1.2,2.0]:
    eg=integrate(np.array([0.05]),np.array([0.05]),np.array([p27b]),make_params(GNP))[0][0]
    em=integrate(np.array([0.05]),np.array([0.05]),np.array([p27b]),make_params(MB))[0][0]
    print(f"  {p27b:5.2f} :   {eg:6.3f}     {em:6.3f}   {'(GNP G0)' if eg<0.3 else ''}{' (MB G0)' if em<0.3 else ''}")
print("  -> the birth-p27 SEPARATRIX (commit->G0) is HIGHER for MB than GNP  =  MB needs MORE p27 to arrest.")

# ---- 2. calibrate a SHARED flat birth-p so GNP ~0% G0; then MB with the SAME birth-p ----
print("\n=== 2. Dichotomy under the proteomics constraint (SHARED flat birth-p; only Cd/kCd/p18 differ) ===")
BIRTH=dict(m21=0.05,cv21=0.8,m27=0.12,cv27=0.5)   # calibrated so GNP ~0% G0
for lab,ct in [('GNP',GNP),('MB',MB)]:
    r=run_ensemble(ct,**BIRTH)
    print(f"  {lab}: G0 = {r['g0_frac']*100:5.1f}%")
print("  -> with p27/p21 FLAT (as measured), MB does NOT gain G0. Direction of the CyclinD1 effect:")

# ---- 3. THE core result: G0 fraction vs CyclinD1 (drive). Does high CyclinD1 suppress G0? ----
print("\n=== 3. G0 fraction vs CyclinD1 level (Cd), p18=1, shared birth-p (raised m27=0.3 to expose the trend) ===")
B2=dict(m21=0.05,cv21=0.8,m27=0.30,cv27=0.5)
for cd in [0.5,1.0,2.0,3.0,4.0,6.0]:
    r=run_ensemble(dict(Cd=cd,kCd=1.0,p18=1.0),**B2)
    print(f"  Cd={cd:.1f}: G0 = {r['g0_frac']*100:5.1f}%")
print("  -> higher CyclinD1 => LESS G0 (drive pushes cells over the R-point). MB's Cd-up moves the WRONG way.")

# ---- 4. JP's question: can elevated p18 (INK4) rescue G0 at MB's high Cd? how much is needed? ----
print("\n=== 4. G0 vs p18 at MB drive (Cd=4, kCd=2). How much INK4 to overcome high CyclinD1? ===")
for p18 in [1,3,6,10,20,40]:
    r=run_ensemble(dict(Cd=4.0,kCd=2.0,p18=p18),**B2)
    print(f"  p18={p18:2d}x: G0 = {r['g0_frac']*100:5.1f}%")
print("  -> p18 acts on the SATURATED CDK4/6 arm; the measured 3x is far from enough at Cd=4/kCd=2.")

# ---- 5. what WOULD give MB>GNP G0: MB-elevated CDK2-arm inhibitor (p27/p21) -- but NOT observed ----
print("\n=== 5. G0 vs birth-p27 at MB drive (Cd=4,kCd=2,p18=3) -- the route the data DON'T support ===")
for m27 in [0.12,0.3,0.6,1.0,1.5,2.5]:
    r=run_ensemble(MB,m21=0.05,cv21=0.8,m27=m27,cv27=0.5)
    print(f"  MB birth-p27 mean={m27:.2f}: G0 = {r['g0_frac']*100:5.1f}%   (proteomics: p27 FLAT vs GNP -> no elevation)")

# ---- 6. birth-p21 heterogeneity (Spencer route) at MB drive ----
print("\n=== 6. G0 vs birth-p21 mean & CV at MB drive (the Spencer birth-p21 route) ===")
for m21 in [0.05,0.3,0.8,1.5]:
    r=run_ensemble(MB,m21=m21,cv21=1.0,m27=0.12,cv27=0.5)
    print(f"  MB birth-p21 mean={m21:.2f} (cv1.0): G0 = {r['g0_frac']*100:5.1f}%   (proteomics: p21 UNDETECTABLE -> mean low)")
