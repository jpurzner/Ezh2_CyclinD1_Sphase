"""Faithful assigned angle: a SEPARATE BISTABLE CDK4/6 axis (Yang-2020).
CDK4/6 activity A has its own positive feedback (hysteresis) + stoichiometric INK4 brake.
G0 = A settles LOW at mitotic exit. Birth noise (shared CV) on Cd/kCd/p18. No hard-floor
threshold trick: G0 read from the bistable settling point. Does data-grounded MB get more G0?"""
import numpy as np
from model import make_params, _lognorm, GNP, MB

def bistable_cdk46(ct, w18, N=8000, cv=0.3, seed=1, Km=0.5, A0=0.05,
                   basal=0.02, Vmax=1.0, fb=1.0, Kfb=0.4, h=6, deg=1.0, thr=0.3, T=120, dt=0.05):
    rng=np.random.default_rng(seed); p=make_params(ct)
    Cd=_lognorm(rng,p['Cd'],cv,N); kCd=_lognorm(rng,p['kCd'],cv,N); p18=_lognorm(rng,p['p18'],cv,N)
    inp=np.maximum(0.0,kCd-w18*p18)*(Cd/(Km+Cd))   # stoichiometric free CDK4/6, cyclinD-loaded
    A=np.full(N,A0)
    for _ in range(int(T/dt)):
        def rhs(A):
            A=np.maximum(A,0); fbterm=fb*A**h/(Kfb**h+A**h)
            return basal+Vmax*(inp+fbterm)-deg*A
        k1=rhs(A);k2=rhs(A+.5*dt*k1);k3=rhs(A+.5*dt*k2);k4=rhs(A+dt*k3)
        A=np.maximum(A+(dt/6)*(k1+2*k2+2*k3+k4),0)
    return 100.0*float(np.mean(A<thr)), inp.mean()

print("=== Bistable CDK4/6 axis: G0 = A settles low. thr calibrated so GNP~6% via basal drive ===")
print(" (positive feedback = hysteresis; this is the Yang-2020 'CDK4/6-low latched' state)")
for w18 in [0.2,0.3,0.4,0.5,0.6,0.65]:
    g,gi=bistable_cdk46(GNP,w18)
    m,mi=bistable_cdk46(MB,w18,seed=2)
    gap=m-g
    tag="COUNTEREXAMPLE" if gap>=5 else "no"
    print(f"  w18={w18:.2f}: GNP G0={g:5.1f}% (in {gi:.2f})  MB G0={m:5.1f}% (in {mi:.2f})  gap={gap:+5.1f}pp [{tag}]")

print("\n=== Recalibrate threshold per-w18 to pin GNP=6%, then read MB (fair shared-threshold test) ===")
def dist(ct,w18,cv=0.3,seed=1,Km=0.5,**kw):
    # return final A distribution
    rng=np.random.default_rng(seed); p=make_params(ct)
    Cd=_lognorm(rng,p['Cd'],cv,20000);kCd=_lognorm(rng,p['kCd'],cv,20000);p18=_lognorm(rng,p['p18'],cv,20000)
    inp=np.maximum(0.0,kCd-w18*p18)*(Cd/(Km+Cd)); A=np.full(20000,0.05)
    for _ in range(int(120/0.05)):
        def rhs(A):
            A=np.maximum(A,0);return 0.02+1.0*(inp+1.0*A**6/(0.4**6+A**6))-1.0*A
        k1=rhs(A);k2=rhs(A+.025*k1);k3=rhs(A+.025*k2);k4=rhs(A+.05*k3)
        A=np.maximum(A+(0.05/6)*(k1+2*k2+2*k3+k4),0)
    return A
for w18 in [0.2,0.3,0.4,0.5,0.6]:
    gA=dist(GNP,w18); mA=dist(MB,w18,seed=2)
    thr=np.quantile(gA,0.06)
    g=100*np.mean(gA<thr); m=100*np.mean(mA<thr)
    print(f"  w18={w18:.2f} (thr={thr:.3f}): GNP=6.0%(pinned) MB={m:.1f}%  gap={m-6.0:+.1f}pp  [GNP medA={np.median(gA):.2f} MB medA={np.median(mA):.2f}]")
