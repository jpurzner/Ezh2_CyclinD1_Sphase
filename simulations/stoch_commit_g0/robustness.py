"""Robustness of the 'stoichiometry can't give MB>GNP G0' conclusion + the one viable route (rare birth-p21 tail)."""
import numpy as np
from model import BASE, GNP, MB, make_params, integrate, run_ensemble, _lognorm

# ---- A. Is it 'saturation' or just drive>>brake? Decompose the MB drive/brake ratio (both drive forms) ----
print("=== A. Why MB commits harder: drive vs brake (Dd = kCd*Cd/(Kd*(1+p18)+Cd)) ===")
for lab,ct in [('GNP',GNP),('MB',MB)]:
    p=make_params(ct)
    Dd_sat = p['kCd']*p['Cd']/(p['Kd']*(1+p['w18']*p['p18'])+p['Cd'])
    Dd_lin = p['kCd']*p['Cd']/(1+p['w18']*p['p18'])           # non-saturating comparison
    print(f"  {lab}: drive kCd*Cd={p['kCd']*p['Cd']:.1f}  brake(1+p18)={1+p['p18']:.0f}  ->  Dd_sat={Dd_sat:.2f}  Dd_linear={Dd_lin:.2f}")
print("  -> MB drive (kCd*Cd) up 8x, brake (1+p18) up only 2x => net drive HIGHER in MB under EITHER form.")
print("     So it's not a saturation artifact: the CyclinD1/CDK4-6 increase simply dwarfs the p18 increase.")

# ---- B. How much p18 (INK4) is actually needed to match GNP's G0 at MB drive? (precise threshold) ----
print("\n=== B. p18 needed at MB drive (Cd4,kCd2) to reach GNP-like ~6% G0 ===")
B2=dict(m21=0.05,cv21=0.8,m27=0.12,cv27=0.5)
for p18 in [3,8,12,15,18,22]:
    r=run_ensemble(dict(Cd=4.0,kCd=2.0,p18=p18),**B2)
    print(f"  p18={p18:2d}x (measured=3x): G0={r['g0_frac']*100:5.1f}%")

# ---- C. THE viable route: a rare high-birth-p21 TAIL (Spencer). Bulk mean stays low (undetectable) but a
#         small fraction of cells inherit high p21 (replication stress). Can that give MB G0 at high drive? ----
print("\n=== C. Rare high-birth-p21 subpopulation at MB drive (mean stays ~undetectable; heavy tail) ===")
def mixture_ensemble(ct, frac_hi, p21_hi, N=4000, seed=1, E_thr=0.3):
    rng=np.random.default_rng(seed); p=make_params(ct)
    hi = rng.random(N) < frac_hi
    p21 = np.where(hi, rng.lognormal(np.log(p21_hi), 0.3, N), rng.lognormal(np.log(0.05), 0.5, N))
    p27 = _lognorm(rng,0.12,0.5,N); E=np.full(N,0.05)
    Ef,_,_=integrate(E,p21,p27,p)
    return float(np.mean(Ef<E_thr)), float(np.mean(p21))
for frac,p21hi in [(0.10,1.5),(0.20,1.5),(0.20,2.5),(0.10,3.0)]:
    g0_mb,meanp = mixture_ensemble(MB,frac,p21hi)
    g0_gnp,_    = mixture_ensemble(GNP,frac,p21hi)
    print(f"  {int(frac*100)}% cells @ p21~{p21hi}: MB G0={g0_mb*100:4.1f}%  GNP G0={g0_gnp*100:4.1f}%  (bulk mean p21~{meanp:.2f})")
print("  -> a rare high-p21 tail CAN make G0 cells even at MB drive; but note it also makes GNP G0 unless the")
print("     TAIL itself is MB-specific (more replication stress in MB). The dichotomy then lives in the TAIL, not the mean.")

# ---- D. So what single factor flips GNP(no G0)->MB(G0)? Test each MB change in isolation at a shared birth-p ----
print("\n=== D. Each MB change in isolation (shared birth-p, m27=0.30 to sit near separatrix): does it ADD G0? ===")
B3=dict(m21=0.05,cv21=0.8,m27=0.30,cv27=0.5)
base=run_ensemble(dict(Cd=1.0,kCd=1.0,p18=1.0),**B3)['g0_frac']
print(f"  GNP baseline: G0={base*100:.1f}%")
for lab,ct in [('+Cd4x only',dict(Cd=4.0,kCd=1.0,p18=1.0)),
               ('+kCd2x only',dict(Cd=1.0,kCd=2.0,p18=1.0)),
               ('+p18 3x only',dict(Cd=1.0,kCd=1.0,p18=3.0)),
               ('full MB',MB)]:
    g=run_ensemble(ct,**B3)['g0_frac']
    print(f"  {lab:14s}: G0={g*100:5.1f}%  (Delta {(g-base)*100:+.1f})")
print("  -> Cd and kCd REMOVE G0; only p18 adds it, but 3x is swamped by the 8x drive in full MB.")
