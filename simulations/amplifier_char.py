"""Characterize a validated amplifier regime (JP 2026-08-09 H3K27me3-amplifier hypothesis).
Given the fully-resolved parameter set (env RP_PARAMS = JSON, from amplifier_probe's resolved_params),
measure the three things that decide whether the amplifier gives a *transient* G0 population:

  (1) HYSTERESIS / transient-vs-permanent -- run MB at high CDKI until it arrests, then drop CDKI back to
      1x mid-simulation and see whether it RE-ENTERS the cycle. Re-entry => reversible/transient G0
      (good); stays arrested => permanent lock (bad). Also the reverse control (stays cycling at 1x).
  (2) ENTRY COMPOUNDING -- at a near-threshold CDKI, the progression of successive cycle lengths and the
      banking of Mk / erosion of cdk6+Cd over the first several divisions ("over a couple of divisions").
  (3) POPULATION G0 fraction vs CDKI variance (CV) -- lognormal per-cell CDKI on p18/p19/p27; fraction
      that fail to divide in the window = transient-G0 subpopulation. Reports at CV 0.6 / 0.9 / 1.2.

Run:  RP_PARAMS="$(cat resolved.json)" PYTHONPATH=. ./venv/bin/python simulations/amplifier_char.py
"""
import os, json, numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

RP=json.loads(os.environ.get('RP_PARAMS','') or '{}')
BASE=dict(kPhRbCd=0.1748299602292031,K_CdRb=0.4025007594938155,w_ink4=4.646376316224687,kSyP21=0.0004043793202695936,
          K_cdk6_sink=3.879643001142466,w_p18=1.003950915185932,w_p19=0.6655007434943913,
          k_Cd_tx_Gli_max=0.24197652013556434,k_Cd_tx_MYCN=0.11031706915457935)
P=dict(BASE); P.update(RP)
KSYP21_BASE=P.get('kSyP21',BASE['kSyP21'])
M=build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, params=P)
RR=te.loada(M); RR.integrator.setValue('relative_tolerance',1e-7); RR.integrator.setValue('absolute_tolerance',1e-9)
RR.integrator.setValue('maximum_num_steps',8000000)

def set_params(rr, mult):
    for k,v in P.items():
        try: rr[k]=v
        except: pass
    rr['SHH']=0.5; rr['Ptch1_copy_number']=0.3
    try: rr['Cd2_expr']=rr['CD2_EXPR_MB']
    except: pass
    rr['p18']=1.73*mult; rr['p19']=0.58*mult; rr['kSyP21']=KSYP21_BASE*mult

def divs(dna): return np.where((dna[:-1]>0.9)&(dna[1:]<0.1))[0]

def ndiv(mult, T=30000):
    set_params(RR, mult); RR.reset()
    r=RR.simulate(0,T,max(3000,T//10),selections=['time','Dna'])
    return len(divs(r['Dna']))

def hysteresis(mult_hi=3.0, t_switch=20000, T=45000):
    """arrest at mult_hi, then drop to 1x at t_switch; count divisions before/after the switch."""
    set_params(RR, mult_hi); RR.reset()
    npts1=max(2000,t_switch//10)
    r1=RR.simulate(0,t_switch,npts1,selections=['time','Dna'])
    n_before=len(divs(r1['Dna']))
    # drop CDKI to 1x (keep species state = continue from arrested state)
    RR['p18']=1.73*1.0; RR['p19']=0.58*1.0; RR['kSyP21']=KSYP21_BASE*1.0
    npts2=max(2000,(T-t_switch)//10)
    r2=RR.simulate(t_switch,T,npts2,selections=['time','Dna'])
    n_after=len(divs(r2['Dna']))
    return n_before, n_after

def entry(mult, T=40000):
    set_params(RR, mult); RR.reset()
    r=RR.simulate(0,T,max(4000,T//8),selections=['time','Dna','Mk','cdk6','Cd']); d=divs(r['Dna']); t=r['time']
    if len(d)>=4:
        iv=np.diff(t[d])/60
        return [round(x,1) for x in iv[:6]], round(float(r['Mk'][d[0]]),4), round(float(r['Mk'][d[-1]]),4), \
               round(float(r['cdk6'][d[0]]),3), round(float(r['cdk6'][d[-1]]),3), round(float(r['Cd'][d[0]]),3), round(float(r['Cd'][d[-1]]),3)
    return [], None,None,None,None,None,None

def population(cv, N=40, seed=0):
    rng=np.random.default_rng(seed)
    mu=np.log(1.0)-0.5*np.log(1+cv**2); sd=np.sqrt(np.log(1+cv**2))
    mults=rng.lognormal(mu,sd,N)
    g0=sum(1 for m in mults if ndiv(m)==0)
    return round(100.0*g0/N,1)

out=dict(rp_keys=sorted([k for k in RP.keys()]))
nb,na=hysteresis()
out['hysteresis']=dict(div_before_at3x=nb, div_after_drop_to_1x=na,
                       verdict=('transient/reversible' if (nb==0 and na>0) else ('permanent-lock' if (nb==0 and na==0) else 'did-not-arrest-at-3x')))
iv,mk0,mkL,c60,c6L,cd0,cdL=entry(2.0)
out['entry_at_2x']=dict(cycle_intervals_h=iv, Mk_first=mk0, Mk_last=mkL, cdk6_first=c60, cdk6_last=c6L, Cd_first=cd0, Cd_last=cdL)
out['arrest_thr']=next((mx for mx in (2,3,4,5,6) if ndiv(mx)==0), 99)
out['population_G0_pct']={f'CV{cv}': population(cv) for cv in (0.6,0.9,1.2)}
print(json.dumps(out))
