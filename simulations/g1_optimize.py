"""#2 re-balance: find params giving a WIDE graded GNP G1 dose-response AND validation ~30/32.
The graded window needs a lower commitment threshold (kPhRbCd up) which erodes the SHH=0 arrest -> restore it via
lower k_Cd_tx_basal; restore periods via mu / k_mu_cki. Objective = validation passed (32) + a GRADED-G1 bonus so the
re-balance can't just disengage K_g1len. Drives validate_v44.py (subprocess) + 2 inline GNP sims (SHH 0.5, 0.25).
Run: ./venv/bin/python simulations/g1_optimize.py [--rand=80] [--workers=6] [--climb=3]
"""
import sys, os, json, subprocess, tempfile, random, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from concurrent.futures import ThreadPoolExecutor
from src.build_model_v44_heldt import build_model_v44

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT,'venv','bin','python'); VALIDATE=os.path.join(ROOT,'simulations','validate_v44.py')
TMP=tempfile.mkdtemp(prefix='g1_'); BEST=os.path.join(ROOT,'simulations','g1_optimize_best.json'); LOG=os.path.join(ROOT,'simulations','g1_optimize_log.jsonl')

SPACE=[('mu',0.0007,0.0014,True),('k_mu_cki',0.30,0.75,True),('kPhRbCd',0.15,0.28,True),('K_CdRb',0.30,0.50,True),
       ('k_Cd_tx_basal',0.00015,0.00055,True),('K_g1len',0.9,1.8,True),('mu_min_frac',0.15,0.40,True),('n_g1len',2,3,False)]
SEED={'mu':0.0011,'k_mu_cki':0.55,'kPhRbCd':0.20,'K_CdRb':0.40,'k_Cd_tx_basal':0.00030,'K_g1len':1.3,'mu_min_frac':0.2,'n_g1len':2}

def sample(rng): return {n:(10**rng.uniform(math.log10(lo),math.log10(hi)) if lg else round(rng.uniform(lo,hi))) for n,lo,hi,lg in SPACE}

def gnp_g1(model, shh):
    """GNP G1 duration (h) at this SHH; (ndiv, G1_h)."""
    rr=te.loada(model); rr.integrator.setValue("relative_tolerance",1e-6); rr.integrator.setValue("maximum_num_steps",1000000)
    rr['SHH']=shh; rr['Ptch1_copy_number']=1.0; rr['p18']=0.464; rr['p19']=0.36
    for atol in (1e-8,1e-7,1e-6):
        rr.reset()
        try:
            rr.integrator.setValue("absolute_tolerance",atol)
            r=rr.simulate(0,20000,40000,selections=["time","MPF","Dna","aRc","vfork"]); break
        except Exception: r=None
    if r is None: return (0,float('nan'))
    tt=r['time']; m=tt>=4000; T=tt[m]; dt=T[1]-T[0]
    pk,_=find_peaks(r['MPF'][m],prominence=0.12,distance=int(180/dt))
    if len(pk)<2: return (len(pk),float('nan'))
    per=np.median(np.diff(T[pk]/60.0))
    syn=r['vfork'][m]*r['aRc'][m]; Dna=r['Dna'][m]
    preS=~((syn>0.02)&(Dna<0.98))&~(Dna>=0.98)
    return (len(pk), per*preS.mean())

def evaluate(params):
    # SUBPROCESS validation ONLY here (roadrunner is NOT thread-safe -> no inline sims in the pool). K_g1len engagement
    # is forced by the search bounds; the graded-G1 metric is computed single-threaded post-hoc (add_grade below).
    rp=os.path.join(TMP,f'r_{abs(hash(json.dumps(params,sort_keys=True)))%10**9}.json')
    env=dict(os.environ,CDKI_SPECIES='1',P21_PIP_DEGRON='1',PYTHONPATH=ROOT,VALIDATE_PARAMS=json.dumps(params),VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY,VALIDATE],env=env,cwd=ROOT,capture_output=True,timeout=600); out=json.load(open(rp))
    except Exception as e: return dict(params=params,passed=0,grade=0,err=str(e)[:40])
    passed=int(out.get('passed',0)); cd=0.0
    for c in out.get('checks',[]):
        if c.get('name')=='CyclinD1 MB/GNP': cd=float(c.get('actual',0))
    return dict(params=params,passed=passed,cd=cd,grade=None)

def add_grade(r):
    # single-threaded graded-G1 metric for a candidate (safe: no threads)
    try:
        model=build_model_v44(params=r['params'])
        n5,g5=gnp_g1(model,0.5); n25,g25=gnp_g1(model,0.25)
        r['g1_hi']=g5; r['g1_lo']=g25
        r['grade']=float(g25/g5) if (n5>=2 and n25>=2 and g5==g5 and g25==g25 and g5>0) else 0.0
    except Exception:
        r['grade']=0.0
    return r

def score(r):
    g=r.get('grade')
    gbonus = (min(g,2.5) if g is not None else 0.0)
    gate = (1.0 if (g is not None and g>=1.4) else 0.0)
    return (r['passed'] + gate, gbonus, -abs(r.get('cd',0)-5.07))
def _arg(f,d):
    for a in sys.argv:
        if a.startswith(f+'='): return int(a.split('=')[1])
    return d

if __name__=='__main__':
    N=_arg('--rand',80);W=_arg('--workers',6);C=_arg('--climb',3); rng=random.Random(13); logfh=open(LOG,'w')
    def batch(cs,tag):
        res=[]
        with ThreadPoolExecutor(max_workers=W) as pool:
            for r in pool.map(evaluate,cs):
                res.append(r); logfh.write(json.dumps({'tag':tag,'passed':r['passed'],'grade':round(r.get('grade') or 0,2),'g1_hi':r.get('g1_hi'),'g1_lo':r.get('g1_lo'),'cd':r.get('cd'),'params':r['params']})+'\n'); logfh.flush()
        return res
    def grade_top(allr, k=12):
        # add the single-threaded graded-G1 metric to the top-k (by passed) that lack it, then re-sort by full score
        allr.sort(key=lambda r:(r['passed'], -abs(r.get('cd',0)-5.07)), reverse=True)
        for r in allr[:k]:
            if r.get('grade') is None: add_grade(r)
        return sorted(allr, key=score, reverse=True)
    def rpt(tag,b): print(f"{tag}: passed {b['passed']}/32 grade {(b.get('grade') or 0):.2f} (G1 {b.get('g1_hi',0) or 0:.0f}->{b.get('g1_lo',0) or 0:.0f}h) CyclinD1 {b.get('cd'):.2f}",flush=True)
    allr=grade_top(batch([SEED]+[sample(rng) for _ in range(N)],'rand'))
    rpt("RANDOM best", allr[0])
    for rd in range(C):
        base=allr[0]['params']; jit=[]
        for _ in range(W*4):
            p=dict(base)
            for n,lo,hi,lg in SPACE:
                if rng.random()<0.5: p[n]=(min(hi,max(lo,p[n]*rng.uniform(0.8,1.25))) if lg else round(min(hi,max(lo,p[n]+rng.choice([-1,0,1])))))
            jit.append(p)
        allr=grade_top(allr+batch(jit,f'climb{rd}'))
        rpt(f"CLIMB {rd}", allr[0])
    b=allr[0]; json.dump({'passed':b['passed'],'grade':b.get('grade'),'g1_hi':b.get('g1_hi'),'g1_lo':b.get('g1_lo'),'cd':b.get('cd'),'params':b['params']},open(BEST,'w'),indent=2)
    print(f"\n=== BEST: passed {b['passed']}/32, graded-G1 ratio {b.get('grade'):.2f} (G1 {b.get('g1_hi',0):.1f}h -> {b.get('g1_lo',0):.1f}h) ===\n{json.dumps(b['params'],indent=2)}")
