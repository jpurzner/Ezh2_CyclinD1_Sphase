"""Reduced stochastic-commitment model (p21 / p18 / p27) to test the GNP-vs-MB transient-G0 dichotomy.
VECTORIZED ensemble (all cells integrated simultaneously as numpy arrays).

Question (JP): can transient G0 appear in MB but NOT GNP via increased CKI expression, even with MB's high
CyclinD1? Which settings permit it; if not, why?

Structure (Yao Rb-E2F bistable switch + Overton p21/p27<->CDK2 double-negative + INK4/CIP brake on a
SATURATING CyclinD-CDK4/6 drive):
  Dd  = kCd*Cd / (Kd*(1 + w18*p18) + Cd)                 # CyclinD-CDK4/6 (mono-Rb); saturates in Cd; INK4(p18) brake
  inh = 1 / (1 + (p21+p27)/Kp)                            # CIP/KIP inhibition of the CyclinE-CDK2 (E2F self-act) arm
  dE  = kb + km*(Dd + Erb*E^n/(Krb^n+E^n))*inh - kdE*E    # E2F: mono-Rb drive + CDK2 self-activation, CIP-inhibited
  dp21= s21 - (d0 + dS*E)*p21 ;  dp27= s27 - (d0 + dS*E)*p27   # CDK2/Skp2 degradation (double-negative)
Fate after settling: E>E_thr = committed (CDK2-inc); E<E_thr = transient G0 (CDK2-low). Ensemble over birth
p21/p27 (lognormal, inherited from mother G2 = noise). p18 = cell-type INK4 constant.
"""
import numpy as np

BASE = dict(kCd=1.0, Kd=0.5, w18=1.0, km=1.0, Erb=1.0, Krb=0.4, n=6, kb=0.01, kdE=1.0,
            Kp=0.3, s21=0.03, s27=0.02, d0=0.006, dS=2.0, Kd2=0.5, md=4)
# d0 small + Skp2 degradation ULTRASENSITIVE in E (Kd2,md): p21/p27 cleared only when COMMITTED (E high, Skp2=E2F
# target) and MAINTAINED at low E -> genuine bistable double-negative. This is the bifurcation the literature needs.
GNP = dict(Cd=1.0, kCd=1.0, p18=1.0)
MB  = dict(Cd=4.0, kCd=2.0, p18=3.0)   # CyclinD1 ~4x, CDK4/6 ~2x, p18 ~3x; p27 flat, p21 low (birth)

def make_params(ct, extra=None):
    p = dict(BASE); p.update(ct)
    if extra: p.update(extra)
    return p

def integrate(E, p21, p27, p, T=120.0, dt=0.05):
    """Vectorized RK4 over arrays E,p21,p27 (shape N)."""
    def rhs(E, p21, p27):
        E=np.maximum(E,0); p21=np.maximum(p21,0); p27=np.maximum(p27,0)
        Dd = p['kCd']*p['Cd']/(p['Kd']*(1+p['w18']*p['p18'])+p['Cd'])
        inh = 1.0/(1.0+(p21+p27)/p['Kp'])
        sa = p['Erb']*E**p['n']/(p['Krb']**p['n']+E**p['n'])
        dE = p['kb'] + p['km']*(Dd+sa)*inh - p['kdE']*E
        deg = p['d0'] + p['dS']*E**p['md']/(p['Kd2']**p['md']+E**p['md'])   # Skp2 (E2F target): ultrasensitive in E
        return dE, p['s21']-deg*p21, p['s27']-deg*p27
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

def _lognorm(rng,mean,cv,N):
    if mean<=0: return np.zeros(N)
    sig=np.sqrt(np.log(1+cv**2)); mu=np.log(mean)-0.5*sig**2
    return rng.lognormal(mu,sig,N)

def run_ensemble(ct, N=2000, m21=0.05, cv21=0.8, m27=0.5, cv27=0.4, E_thr=0.3, seed=0, extra=None, E0=0.05):
    rng=np.random.default_rng(seed); p=make_params(ct,extra)
    p21_0=_lognorm(rng,m21,cv21,N); p27_0=_lognorm(rng,m27,cv27,N); E=np.full(N,E0)
    Ef,_,_=integrate(E,p21_0,p27_0,p)
    g0=float(np.mean(Ef<E_thr))
    return dict(g0_frac=g0,E=Ef,p21_birth=p21_0,p27_birth=p27_0,
                E_commit=float(np.mean(Ef[Ef>=E_thr])) if (Ef>=E_thr).any() else np.nan,
                E_g0=float(np.mean(Ef[Ef<E_thr])) if (Ef<E_thr).any() else np.nan)

if __name__=='__main__':
    print("=== bistability check (GNP params, vary birth p21; hold p27=0.5) ===")
    p=make_params(GNP)
    for x in [0.0,0.2,0.5,1.0,2.0,4.0]:
        Ef,pf,_=integrate(np.array([0.05]),np.array([x]),np.array([0.5]),p)
        print(f"  p21_birth={x:.1f}: E_final={Ef[0]:.3f}  -> {'G0' if Ef[0]<0.3 else 'commit'}")
    # hysteresis: same params, high-E vs low-E start (fixed low p births)
    lo=integrate(np.array([0.05]),np.array([0.1]),np.array([0.1]),p)[0][0]
    hi=integrate(np.array([2.0]),np.array([0.1]),np.array([0.1]),p)[0][0]
    print(f"  hysteresis @ low birth-p: E from low-start={lo:.3f}, from high-start={hi:.3f}  ({'BISTABLE' if abs(lo-hi)>0.2 else 'monostable'})")
    print()
    print("=== baseline dichotomy: GNP vs MB (Cd4x,kCd2x,p18 3x; p27 flat; p21 low birth cv0.8) ===")
    for lab,ct in [('GNP',GNP),('MB',MB)]:
        r=run_ensemble(ct)
        print(f"  {lab}: G0 = {r['g0_frac']*100:5.1f}%   E_commit~{r['E_commit']:.2f} E_G0~{r['E_g0']:.2f}")
