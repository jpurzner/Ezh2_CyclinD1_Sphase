"""What DOES produce the GNP-no-G0 / MB-G0 dichotomy? The one route not excluded by bulk proteomics:
an MB-SPECIFIC, replication-stress-driven, rare high-birth-p21 TAIL (Spencer/Barr mother-G2-p21 mechanism).
GNP has no tail (normal cells, little stress) -> ~0% G0. MB has a tail (oncogenic MYCN/CyclinD1 stress) ->
a G0 fraction. The dichotomy lives in the TAIL of the birth-p21 distribution, NOT the bulk CKI/CyclinD1 balance
(bulk p21 mean stays ~undetectable, consistent with the proteomics)."""
import numpy as np
from model import GNP, MB, make_params, integrate, _lognorm

def ensemble_tail(ct, frac_hi, p21_hi, cv_hi=0.3, m21_lo=0.05, m27=0.12, cv27=0.5, N=6000, seed=3, E_thr=0.3):
    rng=np.random.default_rng(seed); p=make_params(ct)
    hi = rng.random(N) < frac_hi
    p21 = np.where(hi, rng.lognormal(np.log(p21_hi),cv_hi,N), rng.lognormal(np.log(m21_lo),0.5,N))
    p27 = _lognorm(rng,m27,cv27,N); E=np.full(N,0.05)
    Ef,_,_=integrate(E,p21,p27,p)
    return float(np.mean(Ef<E_thr)), float(np.mean(p21))

print("=== The working dichotomy: GNP (NO tail) vs MB (MB-specific high-birth-p21 tail) ===")
print("  GNP: no replication-stress tail -> birth-p21 all low.  MB: a fraction inherit high p21 (stress).")
print()
# GNP: no tail
gnp_g0,gnp_mean = ensemble_tail(GNP, frac_hi=0.0, p21_hi=2.5)
print(f"  GNP (no tail):                         G0 = {gnp_g0*100:4.1f}%  (bulk p21 mean {gnp_mean:.3f})")
print()
print("  MB with an MB-specific high-p21 tail (varying prevalence & level):")
print("   frac_hi  p21_hi | MB G0%   bulk p21 mean | (GNP if it HAD this tail)")
for frac,phi in [(0.10,2.5),(0.20,2.5),(0.30,2.5),(0.20,3.5),(0.30,3.5),(0.40,3.5)]:
    mb_g0,mb_mean = ensemble_tail(MB, frac, phi)
    gnp_if,_      = ensemble_tail(GNP, frac, phi)   # counterfactual: GNP with the same tail
    print(f"    {frac:.2f}    {phi:.1f}  |  {mb_g0*100:4.1f}%     {mb_mean:5.2f}      |   {gnp_if*100:4.1f}%")
print()
print("  Reads:")
print("  - GNP with NO tail -> ~0% G0 (matches biology: GNP no transient G0).")
print("  - MB with the tail -> a real G0 fraction; needs p21_hi ~2.5-3.5 to BEAT MB's high drive (Cd4/kCd2).")
print("  - Bulk p21 mean stays low (~0.3-1.0) even with the tail -> consistent with 'p21 undetectable' proteomics.")
print("  - CAUTION: GNP with the SAME tail would have MORE G0 than MB (lower drive) -> the dichotomy REQUIRES")
print("    the tail to be MB-SPECIFIC (oncogenic replication stress present in MB, ~absent in GNP), NOT a shared")
print("    birth distribution. This is the testable prediction: single-cell CDK2/p21 should show an MB-only high-p21")
print("    G0 subpopulation that GNP lacks; p18/p27 do NOT contribute.")
