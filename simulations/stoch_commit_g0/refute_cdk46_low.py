"""REFUTATION ANGLE: cdk46-low-exit (Yang-2020 p21-independent route).

Claim under test: the data-grounded MB (Cd x4, kCd x2, p18 x3; SHARED birth-p with GNP)
CANNOT have G0 fraction >= GNP + 5pp, and in fact commits harder.

This angle: G0 = CDK4/6 activity LOW at mitotic exit (Yang et al 2020, Nat Cell Biol),
a SEPARATE bistable axis from the CDK2 (E2F) bifurcation. The distinguishing physics is
how INK4 (p18) inhibits CDK4/6:

  - Baseline model (analysis.py): Michaelis brake  Dd = kCd*Cd/(Kd*(1+w18*p18)+Cd).
    Here cyclin D (Cd) is the DRIVE numerator, so Cd x4 * kCd x2 = 8x drive swamps the
    additive p18 brake. MB Dd is HIGHER than GNP. Claim wins trivially.

  - This angle's biology: INK4 proteins are STOICHIOMETRIC inhibitors -- p18 binds and
    sequesters CDK4/6 monomers, blocking cyclin-D loading. So the correct competition is
    p18 (abundance, up 3x) vs CDK4/6 (abundance, up 2x), NOT p18 vs cyclinD. Under
    sequestration the drive is CDK4/6 (2x) and the brake is p18 (3x): brake outpaces drive.
    Cyclin D only matters if free CDK4/6 remains.

We test, honestly, whether stoichiometric sequestration + a CDK4/6-low G0 definition can
make the DATA-GROUNDED MB have G0 >= GNP + 5pp with SHARED birth-p, and whether the
required regime is biologically defensible.
"""
import numpy as np
from model import BASE, GNP, MB, make_params, integrate, run_ensemble, _lognorm

rng_master = np.random.default_rng(0)

# --------------------------------------------------------------------------------------
# Active CDK4/6-cyclinD complex under STOICHIOMETRIC INK4 sequestration.
#   CDK4/6 monomer abundance  ~ kCd   (GNP 1, MB 2)
#   INK4 (p18) abundance      ~ p18   (GNP 1, MB 3)  -- sequesters CDK4/6 1:1 (scaled by w18)
#   cyclinD abundance         ~ Cd    (GNP 1, MB 4)
# free CDK4/6 available for cyclinD = max(0, kCd - w18*p18)
# active complex = cyclinD-loaded free CDK4/6 = free * Cd/(Km+Cd)   (cyclinD loading, saturable)
# --------------------------------------------------------------------------------------
def active_cdk46(p, w18, Km=0.5):
    free = np.maximum(0.0, p['kCd'] - w18*p['p18'])
    load = p['Cd']/(Km + p['Cd'])
    return free*load

def report(label, gnp_val, mb_val, thr_desc=""):
    gap = mb_val - gnp_val
    verdict = "COUNTEREXAMPLE" if gap >= 5.0 else "no"
    print(f"  {label:38s} GNP={gnp_val:5.1f}%  MB={mb_val:5.1f}%  gap={gap:+5.1f}pp  [{verdict}] {thr_desc}")
    return gap

print("="*90)
print("STEP 0. Deterministic active CDK4/6 under stoichiometric sequestration (no noise), vs w18")
print("        free = max(0, kCd - w18*p18); GNP(kCd1,p18 1) vs MB(kCd2,p18 3). Cd loading included.")
print("="*90)
print(f"  {'w18':>5}  {'GNP free':>9} {'MB free':>9} {'GNP active':>11} {'MB active':>10}  note")
for w18 in [0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.8,1.0]:
    g = active_cdk46(make_params(GNP), w18)
    m = active_cdk46(make_params(MB), w18)
    gf = max(0.0, 1-w18*1); mf = max(0.0, 2-w18*3)
    note = "MB<GNP active" if m < g else ""
    print(f"  {w18:5.2f}  {gf:9.3f} {mf:9.3f} {g:11.3f} {m:10.3f}  {note}")

print("\n  -> Cd x4 nearly saturates the loading term (Cd/(0.5+Cd)=0.89 vs GNP 0.67), so it does NOT")
print("     rescue MB when free CDK4/6 collapses. Whether MB active < GNP active depends on w18.")
print("     For MB free < GNP free we need 2-3*w18 < 1-w18  =>  w18 > 0.5.")

# --------------------------------------------------------------------------------------
# STEP 1. Stochastic CDK4/6-low ensemble.  Cell-to-cell variability (birth) on the
# abundances Cd, kCd(CDK4/6), p18 -- SHARED CV across GNP and MB (only the MEANS differ,
# per proteomics). G0 = active CDK4/6 below threshold at mitotic exit (Yang-2020 axis).
# We calibrate the threshold so GNP sits at ~0-6% G0, then read MB at the SAME threshold.
# --------------------------------------------------------------------------------------
def cdk46_ensemble(ct, w18, N=8000, cv=0.3, thr=None, seed=1, Km=0.5, gnp_target_q=0.06):
    rng = np.random.default_rng(seed)
    p = make_params(ct)
    Cd  = _lognorm(rng, p['Cd'],  cv, N)
    kCd = _lognorm(rng, p['kCd'], cv, N)
    p18 = _lognorm(rng, p['p18'], cv, N)
    free = np.maximum(0.0, kCd - w18*p18)
    active = free * (Cd/(Km+Cd))
    return active, p

def frac_below(active, thr):
    return 100.0*float(np.mean(active < thr))

print("\n" + "="*90)
print("STEP 1. Stochastic CDK4/6-low G0 (shared CV on abundances; threshold calibrated to GNP)")
print("="*90)
for w18 in [0.4, 0.5, 0.6, 0.65, 0.7, 0.8]:
    ga,_ = cdk46_ensemble(GNP, w18)
    ma,_ = cdk46_ensemble(MB,  w18, seed=2)
    # calibrate threshold so GNP G0 = 6% (its calibrated value); read MB at same thr
    thr = np.quantile(ga, 0.06)
    g = frac_below(ga, thr); m = frac_below(ma, thr)
    report(f"w18={w18:.2f} (thr@GNP6%={thr:.3f})", g, m)
print("  NOTE: GNP is PINNED to 6% by construction (threshold = its 6th percentile). The test is")
print("        whether MB exceeds 11% at that SAME threshold.")

print("\n" + "="*90)
print("STEP 2. But is GNP biologically allowed to sit near the sequestration cliff? Check GNP's")
print("        OWN median free-CDK4/6 headroom. GNPs are Shh-high, INK4-LOW, highly proliferative.")
print("="*90)
print(f"  {'w18':>5} {'GNP median free':>16} {'GNP frac fully-seq':>20} {'MB median free':>16}")
for w18 in [0.4,0.5,0.6,0.65,0.7,0.8]:
    rng=np.random.default_rng(5)
    gk=_lognorm(rng,1.0,0.3,20000); gp=_lognorm(rng,1.0,0.3,20000)
    gfree=np.maximum(0.0,gk-w18*gp)
    rng=np.random.default_rng(6)
    mk=_lognorm(rng,2.0,0.3,20000); mp=_lognorm(rng,3.0,0.3,20000)
    mfree=np.maximum(0.0,mk-w18*mp)
    print(f"  {w18:5.2f} {np.median(gfree):16.3f} {100*np.mean(gfree<=1e-6):19.1f}% {np.median(mfree):16.3f}")
print("  -> the w18 window that makes MB<GNP (w18>=0.6) forces GNP to have LITTLE free CDK4/6 too")
print("     (>50% sequestered), contradicting the biology of a proliferating INK4-low GNP.")

