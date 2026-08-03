"""REFUTATION ANGLE: CdP27 sequestration buffer (v44 CdP21 mechanism).

Biology (Sherr & Roberts / Cheng): CyclinD1-CDK4/6 titrates CIP/KIP (p27) into an
inactive CyclinD1-p27 complex, removing FREE p27 from CDK2 inhibition. MB has high Cd
=> sequesters more p27. Question: does this open a CDK2-arm G0 route in MB, or (expected)
relieve CDK2 and commit MORE?

Equilibrium buffer, recomputed every RHS eval (Cd constant per cell, p27 total is the species):
    Cd_free + p27_free <-> CdP27 ,   Kb = Cd_free*p27_free / CdP27
    complex C = 0.5[(Cdt+p27t+Kb) - sqrt((Cdt+p27t+Kb)^2 - 4 Cdt p27t)]
    free p27 = p27t - C   (used in the CDK2 inhibition term 'inh')

I test many wirings honestly:
  V1  free p27 in inh; degradation acts on TOTAL p27
  V2  free p27 in inh; degradation acts on FREE p27 (Skp2 targets free/CDK2-bound)
  V3  p27 is a REQUIRED assembly factor for CDK4/6 (activity gated on bound complex) + free p27 -> CDK2
All with SHARED birth-p (GNP-calibrated) and data-grounded MB (Cd4,kCd2,p18 3).
A valid counterexample = MB G0 >= GNP G0 + 5pp.
"""
import numpy as np
from model import BASE, GNP, MB, make_params, _lognorm

def _free_p27(Cdt, p27t, Kb):
    s = Cdt + p27t + Kb
    C = 0.5*(s - np.sqrt(np.maximum(s*s - 4.0*Cdt*p27t, 0.0)))
    return np.maximum(p27t - C, 0.0), C   # free p27, complex

def integrate_buf(E, p21, p27, p, Kb=0.3, variant='V1', T=120.0, dt=0.05):
    Cdt = p['Cd']   # total CyclinD1 available to buffer p27 (per cell, constant)
    def rhs(E, p21, p27):
        E=np.maximum(E,0); p21=np.maximum(p21,0); p27=np.maximum(p27,0)
        fp27, C = _free_p27(Cdt, p27, Kb)
        # CyclinD-CDK4/6 drive
        if variant=='V3':
            # p27 as assembly factor: CDK4/6 activity scales with amount of p27 loaded onto CyclinD (complex C)
            assemble = C/(0.2 + C)                      # saturating assembly benefit
            Dd = p['kCd']*p['Cd']*(0.4+0.6*assemble)/(p['Kd']*(1+p['w18']*p['p18'])+p['Cd'])
        else:
            Dd = p['kCd']*p['Cd']/(p['Kd']*(1+p['w18']*p['p18'])+p['Cd'])
        inh = 1.0/(1.0+(p21+fp27)/p['Kp'])              # only FREE p27 inhibits CDK2
        sa = p['Erb']*E**p['n']/(p['Krb']**p['n']+E**p['n'])
        dE = p['kb'] + p['km']*(Dd+sa)*inh - p['kdE']*E
        deg = p['d0'] + p['dS']*E**p['md']/(p['Kd2']**p['md']+E**p['md'])
        dp21 = p['s21'] - deg*p21
        if variant=='V2':
            dp27 = p['s27'] - deg*fp27                  # Skp2 degrades only free p27
        else:
            dp27 = p['s27'] - deg*p27
        return dE, dp21, dp27
    n=int(T/dt)
    for _ in range(n):
        a1=rhs(E,p21,p27)
        a2=rhs(E+0.5*dt*a1[0],p21+0.5*dt*a1[1],p27+0.5*dt*a1[2])
        a3=rhs(E+0.5*dt*a2[0],p21+0.5*dt*a2[1],p27+0.5*dt*a2[2])
        a4=rhs(E+dt*a3[0],p21+dt*a3[1],p27+dt*a3[2])
        E  =np.maximum(E  +(dt/6)*(a1[0]+2*a2[0]+2*a3[0]+a4[0]),0)
        p21=np.maximum(p21+(dt/6)*(a1[1]+2*a2[1]+2*a3[1]+a4[1]),0)
        p27=np.maximum(p27+(dt/6)*(a1[2]+2*a2[2]+2*a3[2]+a4[2]),0)
    return E,p21,p27

def ens_buf(ct, Kb, variant, N=3000, m21=0.05, cv21=0.8, m27=0.12, cv27=0.5, E_thr=0.3, seed=0):
    rng=np.random.default_rng(seed); p=make_params(ct)
    p21_0=_lognorm(rng,m21,cv21,N); p27_0=_lognorm(rng,m27,cv27,N); E=np.full(N,0.05)
    Ef,_,_=integrate_buf(E,p21_0,p27_0,p,Kb=Kb,variant=variant)
    return float(np.mean(Ef<E_thr))

# Shared, GNP-calibrated birth-p. p27 FLAT / p21 low, IDENTICAL for GNP and MB.
BIRTHS = [
    dict(m21=0.05,cv21=0.8,m27=0.12,cv27=0.5),
    dict(m21=0.05,cv21=0.8,m27=0.30,cv27=0.5),   # sit nearer separatrix
    dict(m21=0.05,cv21=1.0,m27=0.50,cv27=0.6),   # higher shared p27 + heavier tail
]

print("=== CdP27 sequestration buffer: data-grounded GNP vs MB, SHARED birth-p ===")
print("    (free p27 = total - CyclinD1*p27 complex; only free p27 inhibits CDK2)\n")
best = (-1e9, None)
for variant in ['V1','V2','V3']:
    for Kb in [0.05, 0.2, 0.5, 1.0, 2.0]:
        for bi,B in enumerate(BIRTHS):
            g_gnp = ens_buf(GNP, Kb, variant, **B)*100
            g_mb  = ens_buf(MB,  Kb, variant, **B)*100
            d = g_mb - g_gnp
            flag = '  <== COUNTEREXAMPLE' if d>=5.0 else ''
            print(f"  {variant} Kb={Kb:4.2f} birth{bi} (m27={B['m27']:.2f}): "
                  f"GNP={g_gnp:5.1f}%  MB={g_mb:5.1f}%  (MB-GNP={d:+5.1f}pp){flag}")
            if d>best[0]:
                best=(d,f"{variant} Kb={Kb} birth{bi} m27={B['m27']} -> GNP {g_gnp:.1f}% MB {g_mb:.1f}%")
        print()

print("="*70)
print(f"BEST MB-GNP margin achieved: {best[0]:+.1f} pp")
print(f"  at: {best[1]}")
print(f"  VALID COUNTEREXAMPLE (>=+5pp): {'YES' if best[0]>=5.0 else 'NO'}")
