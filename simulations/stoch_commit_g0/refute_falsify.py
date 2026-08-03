"""Falsify the two load-bearing assumptions of the bistable-CDK4/6-low 'counterexample':
 (E1) cyclinD (Cd x4) is functionally INERT because CDK4/6 is strictly limiting (loading saturates).
 (E2) multiplicative shared-CV noise (so MB's 3x p18 has 3x absolute variance -> fat lower tail).
If the MB>GNP G0 gap collapses when EITHER assumption is relaxed to an equally/more defensible
alternative, the counterexample is not robust and the claim stands."""
import numpy as np
from model import make_params, _lognorm, GNP, MB

def bistable(ct, w18, active_fn, noise='mult', N=20000, cv=0.3, seed=1, thr_q=0.06, gnpA=None):
    rng=np.random.default_rng(seed); p=make_params(ct)
    if noise=='mult':
        Cd=_lognorm(rng,p['Cd'],cv,N); kCd=_lognorm(rng,p['kCd'],cv,N); p18=_lognorm(rng,p['p18'],cv,N)
    else:
        Cd=np.maximum(1e-6,rng.normal(p['Cd'],cv,N)); kCd=np.maximum(1e-6,rng.normal(p['kCd'],cv,N)); p18=np.maximum(0,rng.normal(p['p18'],cv,N))
    inp=active_fn(Cd,kCd,p18,w18); A=np.full(N,0.05)
    for _ in range(int(120/0.05)):
        def rhs(A):
            A=np.maximum(A,0); return 0.02+1.0*(inp+1.0*A**6/(0.4**6+A**6))-1.0*A
        k1=rhs(A);k2=rhs(A+.025*k1);k3=rhs(A+.025*k2);k4=rhs(A+.05*k3)
        A=np.maximum(A+(0.05/6)*(k1+2*k2+2*k3+k4),0)
    return A

# active-complex forms:
def f_cdklim(Cd,kCd,p18,w18):   # CDK4/6 strictly limiting, cyclinD in excess (the counterexample's form)
    return np.maximum(0.0,kCd-w18*p18)*(Cd/(0.5+Cd))
def f_colim(Cd,kCd,p18,w18):    # cyclinD CO-limiting: complex = min(free CDK4/6, cyclinD)  (mass-action-ish)
    return np.minimum(np.maximum(0.0,kCd-w18*p18), Cd)
def f_massaction(Cd,kCd,p18,w18):  # reversible: free CDK4/6 competes between cyclinD binding and INK4 seq
    freeCDK=np.maximum(0.0,kCd-w18*p18)   # INK4 sequesters first
    return freeCDK*Cd/(freeCDK+Cd+1e-9)   # then cyclinD-CDK4/6 complex (both co-limiting, symmetric)

def run(active_fn, noise, label):
    print(f"\n--- {label} ---")
    for w18 in [0.3,0.4,0.5,0.6]:
        gA=bistable(GNP,w18,active_fn,noise=noise,seed=1)
        mA=bistable(MB,w18,active_fn,noise=noise,seed=2)
        thr=np.quantile(gA,0.06)
        g=100*np.mean(gA<thr); m=100*np.mean(mA<thr)
        tag="COUNTEREX" if (m-6.0)>=5 else "no"
        print(f"  w18={w18:.2f} thr={thr:.3f}: GNP=6.0%(pin) MB={m:5.1f}%  gap={m-6.0:+5.1f}pp [{tag}]  medA G/M={np.median(gA):.2f}/{np.median(mA):.2f}")

print("="*80)
print("E1 TEST: relax 'cyclinD inert'. Compare active-complex forms (multiplicative noise).")
print("="*80)
run(f_cdklim, 'mult', "CDK4/6-limiting, cyclinD in excess  [the counterexample's form]")
run(f_colim,  'mult', "cyclinD CO-limiting: complex = min(free CDK4/6, cyclinD)")
run(f_massaction,'mult',"cyclinD-CDK4/6 mass-action (both co-limiting) -- Cd x4 now boosts MB")

print("\n"+"="*80)
print("E2 TEST: relax 'multiplicative noise' -> additive (equal absolute variance G/M).")
print("="*80)
run(f_cdklim, 'add', "CDK4/6-limiting form, ADDITIVE noise (kills MB's inflated variance)")
