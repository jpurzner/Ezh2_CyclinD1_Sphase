"""Focused re-calibration of the option-B (Fan-Meyer p27-inhibitory CDK4/6) model toward max validation.
Drives validate_v44.py via subprocess with P27_OPTIONB=1. Coarse 3-param combo already gives 31/32; this refines
kSeqCd/kRelCd/kDeKPC/w_ink4/kSyP21/kPhRbCd/K_CdRb to (a) maximize passed, (b) pull CyclinD1 MB/GNP -> 5.07.
Run: ./venv/bin/python simulations/optionB_optimize.py [--rand=80] [--workers=6] [--climb=3]
"""
import sys, os, json, subprocess, tempfile, random
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='optB_')
BEST = os.path.join(ROOT, 'simulations', 'optionB_optimize_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optionB_optimize_log.jsonl')

# (name, lo, hi, log-scale) — around the feasible option-B region (coarse combo kSeqCd0.8/kSyP21.0006/kPhRbCd0.32=31/32)
SPACE = [
    ('kSeqCd',  0.2,  1.5, True),
    ('kRelCd',  0.008,0.05,True),
    ('kDeKPC',  0.004,0.03,True),
    ('w_ink4',  1.5,  6.0, True),
    ('kSyP21',  0.0004,0.0016,True),
    ('kPhRbCd', 0.20, 0.40, True),
    ('K_CdRb',  0.30, 0.50, True),
]
SEED = {'kSeqCd':0.8,'kRelCd':0.01647,'kDeKPC':0.02022,'w_ink4':3.6204,'kSyP21':0.0006,'kPhRbCd':0.32,'K_CdRb':0.4032}

def sample(rng):
    p={}
    for n,lo,hi,lg in SPACE:
        p[n]= (10**rng.uniform(__import__('math').log10(lo),__import__('math').log10(hi))) if lg else rng.uniform(lo,hi)
    return p

def evaluate(params):
    rp = os.path.join(TMP, f'r_{abs(hash(json.dumps(params,sort_keys=True)))%10**9}.json')
    env = dict(os.environ, P27_OPTIONB='1', CDKI_SPECIES='1', P21_PIP_DEGRON='1',
               VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=600)
        out = json.load(open(rp))
    except Exception as e:
        return dict(params=params, passed=0, cd=0.0, err=str(e)[:60])
    passed = int(out.get('passed',0))
    checks = {c['name']: c for c in out.get('checks',[])} if 'checks' in out else {}
    cd = 0.0
    for c in out.get('checks', []):
        if c.get('name')=='CyclinD1 MB/GNP': cd=float(c.get('actual',0))
    return dict(params=params, passed=passed, cd=cd)

def score(r):  # maximize passed, then CyclinD1 fold closeness to 5.07
    return (r['passed'], -abs(r.get('cd',0)-5.07))

def _arg(f,d):
    for a in sys.argv:
        if a.startswith(f+'='): return int(a.split('=')[1])
    return d

if __name__=='__main__':
    N=_arg('--rand',80); WORKERS=_arg('--workers',6); CLIMB=_arg('--climb',3)
    rng=random.Random(7)
    cands=[SEED]+[sample(rng) for _ in range(N)]
    logfh=open(LOG,'w')
    def run_batch(cs,tag):
        res=[]
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for r in pool.map(evaluate, cs):
                res.append(r); logfh.write(json.dumps({'tag':tag,'passed':r['passed'],'cd':r.get('cd'),'params':r['params']})+'\n'); logfh.flush()
        return res
    allr=run_batch(cands,'rand')
    allr.sort(key=score, reverse=True)
    b=allr[0]; print(f"RANDOM best: passed {b['passed']}/32  CyclinD1 {b.get('cd'):.2f}", flush=True)
    for rd in range(CLIMB):
        base=allr[0]['params']; jit=[]
        for _ in range(WORKERS*4):
            p=dict(base)
            for n,lo,hi,lg in SPACE:
                if rng.random()<0.5:
                    p[n]=min(hi,max(lo,p[n]*rng.uniform(0.8,1.25)))
            jit.append(p)
        allr=sorted(allr+run_batch(jit,f'climb{rd}'), key=score, reverse=True)
        b=allr[0]; print(f"CLIMB {rd}: passed {b['passed']}/32  CyclinD1 {b.get('cd'):.2f}", flush=True)
    b=allr[0]
    json.dump({'passed':b['passed'],'cd':b.get('cd'),'params':b['params']}, open(BEST,'w'), indent=2)
    print(f"\n=== BEST option-B: {b['passed']}/32  CyclinD1 MB/GNP {b.get('cd'):.2f} ===")
    print(json.dumps(b['params'], indent=2))
