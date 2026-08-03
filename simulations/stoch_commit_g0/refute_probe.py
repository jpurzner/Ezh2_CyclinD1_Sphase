"""Probe the STEP-1 surprise: at w18=0.4-0.5, MB has HIGHER MEAN active CDK4/6 yet MORE G0.
Is this a genuine, defensible mechanism (heteroscedastic INK4 tail) or a fragile artifact?"""
import numpy as np
from model import make_params, integrate, _lognorm, GNP, MB

def sample(ct, w18, N=200000, cv=0.3, seed=1, Km=0.5, noise='mult'):
    rng=np.random.default_rng(seed); p=make_params(ct)
    if noise=='mult':
        Cd=_lognorm(rng,p['Cd'],cv,N); kCd=_lognorm(rng,p['kCd'],cv,N); p18=_lognorm(rng,p['p18'],cv,N)
    else:  # additive: SAME absolute sd for GNP and MB (sd = cv on the GNP mean of each species)
        sdC,sdK,sdP = cv*1.0, cv*1.0, cv*1.0
        Cd=np.maximum(1e-6,rng.normal(p['Cd'],sdC,N)); kCd=np.maximum(1e-6,rng.normal(p['kCd'],sdK,N)); p18=np.maximum(0,rng.normal(p['p18'],sdP,N))
    free=np.maximum(0.0,kCd-w18*p18); active=free*(Cd/(Km+Cd))
    return active

print("=== A. Distribution stats (multiplicative shared-CV noise) ===")
for w18 in [0.4,0.5]:
    g=sample(GNP,w18); m=sample(MB,w18,seed=2)
    thr=np.quantile(g,0.06)
    print(f" w18={w18}: GNP mean={g.mean():.3f} sd={g.std():.3f} | MB mean={m.mean():.3f} sd={m.std():.3f} | thr(GNP6%)={thr:.3f}")
    print(f"          frac<thr: GNP={100*np.mean(g<thr):.1f}%  MB={100*np.mean(m<thr):.1f}%   MB/GNP sd ratio={m.std()/g.std():.2f}")
    print(f"          frac at floor(active==0): GNP={100*np.mean(g<=1e-9):.2f}%  MB={100*np.mean(m<=1e-9):.2f}%")

print("\n=== B. Same test but ADDITIVE noise (equal absolute variance GNP vs MB) -> kills heteroscedasticity ===")
for w18 in [0.4,0.5]:
    g=sample(GNP,w18,noise='add'); m=sample(MB,w18,seed=2,noise='add')
    thr=np.quantile(g,0.06)
    print(f" w18={w18}: frac<thr GNP={100*np.mean(g<thr):.1f}%  MB={100*np.mean(m<thr):.1f}%  gap={100*(np.mean(m<thr)-np.mean(g<thr)):+.1f}pp")

print("\n=== C. Sensitivity of the counterexample to CV (multiplicative). Effect should GROW with CV ===")
for cv in [0.1,0.2,0.3,0.5]:
    g=sample(GNP,0.4,cv=cv); m=sample(MB,0.4,cv=cv,seed=2)
    thr=np.quantile(g,0.06)
    print(f" cv={cv}: GNP=6.0%(pinned) MB={100*np.mean(m<thr):.1f}%  gap={100*(np.mean(m<thr)-0.06):+.1f}pp")

print("\n=== D. Does the CDK4/6-low state actually LATCH to G0 when coupled to E2F self-activation? ===")
# Feed the noisy active CDK4/6 (Dd) into the E-integrator as the drive, birth-E low, SHARED birth-p.
def coupled(ct, w18, N=8000, cv=0.3, seed=1, Km=0.5, E_thr=0.3, m27=0.12):
    rng=np.random.default_rng(seed); p=make_params(ct)
    Cd=_lognorm(rng,p['Cd'],cv,N); kCd=_lognorm(rng,p['kCd'],cv,N); p18=_lognorm(rng,p['p18'],cv,N)
    Ddix=np.maximum(0.0,kCd-w18*p18)*(Cd/(Km+Cd))
    p21=_lognorm(rng,0.05,0.8,N); p27=_lognorm(rng,m27,0.5,N); E=np.full(N,0.05)
    # custom integrate: use Ddix as the (fixed, cell-specific) CDK4/6 drive instead of the Michaelis form
    dt=0.05; nT=int(120/dt)
    for _ in range(nT):
        def rhs(E,p21,p27):
            E=np.maximum(E,0);p21=np.maximum(p21,0);p27=np.maximum(p27,0)
            inh=1.0/(1.0+(p21+p27)/p['Kp']); sa=p['Erb']*E**p['n']/(p['Krb']**p['n']+E**p['n'])
            dE=p['kb']+p['km']*(Ddix+sa)*inh-p['kdE']*E
            deg=p['d0']+p['dS']*E**p['md']/(p['Kd2']**p['md']+E**p['md'])
            return dE,p['s21']-deg*p21,p['s27']-deg*p27
        a1=rhs(E,p21,p27);a2=rhs(E+.5*dt*a1[0],p21+.5*dt*a1[1],p27+.5*dt*a1[2])
        a3=rhs(E+.5*dt*a2[0],p21+.5*dt*a2[1],p27+.5*dt*a2[2]);a4=rhs(E+dt*a3[0],p21+dt*a3[1],p27+dt*a3[2])
        E=np.maximum(E+(dt/6)*(a1[0]+2*a2[0]+2*a3[0]+a4[0]),0)
        p21=np.maximum(p21+(dt/6)*(a1[1]+2*a2[1]+2*a3[1]+a4[1]),0)
        p27=np.maximum(p27+(dt/6)*(a1[2]+2*a2[2]+2*a3[2]+a4[2]),0)
    return 100.0*float(np.mean(E<E_thr))
for w18 in [0.3,0.4,0.5,0.6]:
    g=coupled(GNP,w18); m=coupled(MB,w18,seed=2)
    print(f" w18={w18}: coupled E-low G0  GNP={g:.1f}%  MB={m:.1f}%  gap={m-g:+.1f}pp")
