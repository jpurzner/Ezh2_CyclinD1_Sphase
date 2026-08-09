"""Joint re-calibration of the INTEGRATED model: option B (p27-inhibitory) + CDK6 mark-memory + INK4(p18/p19) sink.
Stacking the three separately-calibrated pieces drops validation to 28/32 (CyclinD1 folds slip). This re-fits the
drive/brake balance to recover ~30/32. Drives validate_v44.py with P27_OPTIONB=1 CDK6_GLI=1. Max passed; tiebreak
CyclinD1 MB/GNP -> 5.07. (Pure subprocess validation -> thread-safe; no inline roadrunner.)
Run: ./venv/bin/python simulations/integrated_optimize.py [--rand=90] [--workers=6] [--climb=3]
"""
import sys, os, json, subprocess, tempfile, random, math
from concurrent.futures import ThreadPoolExecutor
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT,'venv','bin','python'); VALIDATE=os.path.join(ROOT,'simulations','validate_v44.py')
TMP=tempfile.mkdtemp(prefix='integ_'); BEST=os.path.join(ROOT,'simulations','integrated_optimize_best.json'); LOG=os.path.join(ROOT,'simulations','integrated_optimize_log.jsonl')

SPACE=[  # option B + CDK6 sink/INK4 + CyclinD1-fold drive
    ('kPhRbCd',        0.15, 0.40, True),
    ('K_CdRb',         0.30, 0.55, True),
    ('w_ink4',         2.0,  8.0,  True),
    ('kSyP21',         0.0004,0.0016,True),
    ('K_cdk6_sink',    1.0,  6.0,  True),   # higher = less INK4 sink -> more brake
    ('w_p18',          0.6,  3.0,  True),
    ('w_p19',          0.6,  3.0,  True),
    ('k_Cd_tx_Gli_max',0.12, 0.30, True),   # CyclinD1 Gli drive (fold)
    ('k_Cd_tx_MYCN',   0.03, 0.12, True),    # MB-specific CyclinD1 (fold)
]
SEED={'kPhRbCd':0.339,'K_CdRb':0.417,'w_ink4':6.0,'kSyP21':0.000707,'K_cdk6_sink':2.0,'w_p18':1.0,'w_p19':1.0,'k_Cd_tx_Gli_max':0.17494,'k_Cd_tx_MYCN':0.062}

def sample(rng): return {n:(10**rng.uniform(math.log10(lo),math.log10(hi)) if lg else rng.uniform(lo,hi)) for n,lo,hi,lg in SPACE}
def evaluate(params):
    rp=os.path.join(TMP,f'r_{abs(hash(json.dumps(params,sort_keys=True)))%10**9}.json')
    env=dict(os.environ,P27_OPTIONB='1',CDK6_GLI='1',CDKI_SPECIES='1',P21_PIP_DEGRON='1',PYTHONPATH=ROOT,
             VALIDATE_PARAMS=json.dumps(params),VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY,VALIDATE],env=env,cwd=ROOT,capture_output=True,timeout=600); out=json.load(open(rp))
    except Exception as e: return dict(params=params,passed=0,cd=0.0,err=str(e)[:50])
    cd=0.0
    for c in out.get('checks',[]):
        if c.get('name')=='CyclinD1 MB/GNP': cd=float(c.get('actual',0))
    return dict(params=params,passed=int(out.get('passed',0)),cd=cd)
def score(r): return (r['passed'], -abs(r.get('cd',0)-5.07))
def _arg(f,d):
    for a in sys.argv:
        if a.startswith(f+'='): return int(a.split('=')[1])
    return d
if __name__=='__main__':
    N=_arg('--rand',90);W=_arg('--workers',6);C=_arg('--climb',3); rng=random.Random(19); logfh=open(LOG,'w')
    def batch(cs,tag):
        res=[]
        with ThreadPoolExecutor(max_workers=W) as pool:
            for r in pool.map(evaluate,cs):
                res.append(r); logfh.write(json.dumps({'tag':tag,'passed':r['passed'],'cd':r.get('cd'),'params':r['params']})+'\n'); logfh.flush()
        return res
    allr=sorted(batch([SEED]+[sample(rng) for _ in range(N)],'rand'),key=score,reverse=True)
    b=allr[0]; print(f"RANDOM best: {b['passed']}/32 CyclinD1 {b.get('cd'):.2f}",flush=True)
    for rd in range(C):
        base=allr[0]['params']; jit=[]
        for _ in range(W*4):
            p=dict(base)
            for n,lo,hi,lg in SPACE:
                if rng.random()<0.5: p[n]=min(hi,max(lo,p[n]*rng.uniform(0.8,1.25)))
            jit.append(p)
        allr=sorted(allr+batch(jit,f'climb{rd}'),key=score,reverse=True)
        b=allr[0]; print(f"CLIMB {rd}: {b['passed']}/32 CyclinD1 {b.get('cd'):.2f}",flush=True)
    b=allr[0]; json.dump({'passed':b['passed'],'cd':b.get('cd'),'params':b['params']},open(BEST,'w'),indent=2)
    print(f"\n=== BEST integrated: {b['passed']}/32  CyclinD1 {b.get('cd'):.2f} ===\n{json.dumps(b['params'],indent=2)}")
