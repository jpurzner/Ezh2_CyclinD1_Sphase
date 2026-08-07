"""Re-balance the CDK6-arm model (with_cdk6_gli) to hold ~30/32. Drives validate_v44.py with CDK6_GLI=1.
Co-fits the drive-balance (kPhRbCd, K_CdRb, w_ink4) + cdk6 factor params (cdk6_basal, cdk6_Gli, K_Gli_cdk6,
K_EZH2_cdk6) so the CDK6 second arm is engaged but validation is preserved. Maximize passed; tiebreak CyclinD1->5.07.
Run: ./venv/bin/python simulations/cdk6_optimize.py [--rand=80] [--workers=6] [--climb=3]
"""
import sys, os, json, subprocess, tempfile, random, math
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='cdk6_')
BEST = os.path.join(ROOT, 'simulations', 'cdk6_optimize_best.json')
LOG = os.path.join(ROOT, 'simulations', 'cdk6_optimize_log.jsonl')

SPACE = [
    ('kPhRbCd',    0.10, 0.30, True),   # baked 0.1403; cdk6~1 at GNP so drive ~preserved, allow re-balance
    ('K_CdRb',     0.30, 0.55, True),
    ('w_ink4',     2.5,  6.0,  True),
    ('cdk6_basal', 0.15, 0.45, True),
    ('cdk6_Gli',   1.5,  4.0,  True),
    ('K_Gli_cdk6', 0.20, 0.45, True),
    ('K_EZH2_cdk6',4.0,  14.0, True),
]
SEED = {'kPhRbCd':0.1403,'K_CdRb':0.4032,'w_ink4':3.6204,
        'cdk6_basal':0.3,'cdk6_Gli':2.7,'K_Gli_cdk6':0.3,'K_EZH2_cdk6':8.0}

def sample(rng):
    return {n:(10**rng.uniform(math.log10(lo),math.log10(hi)) if lg else rng.uniform(lo,hi)) for n,lo,hi,lg in SPACE}

def evaluate(params):
    rp = os.path.join(TMP, f'r_{abs(hash(json.dumps(params,sort_keys=True)))%10**9}.json')
    env = dict(os.environ, CDK6_GLI='1', CDKI_SPECIES='1', P21_PIP_DEGRON='1',
               VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=600)
        out = json.load(open(rp))
    except Exception as e:
        return dict(params=params, passed=0, cd=0.0, err=str(e)[:50])
    cd = 0.0
    for c in out.get('checks', []):
        if c.get('name')=='CyclinD1 MB/GNP': cd=float(c.get('actual',0))
    return dict(params=params, passed=int(out.get('passed',0)), cd=cd)

def score(r): return (r['passed'], -abs(r.get('cd',0)-5.07))
def _arg(f,d):
    for a in sys.argv:
        if a.startswith(f+'='): return int(a.split('=')[1])
    return d

if __name__=='__main__':
    N=_arg('--rand',80); WORKERS=_arg('--workers',6); CLIMB=_arg('--climb',3)
    rng=random.Random(11); logfh=open(LOG,'w')
    def run_batch(cs,tag):
        res=[]
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for r in pool.map(evaluate, cs):
                res.append(r); logfh.write(json.dumps({'tag':tag,'passed':r['passed'],'cd':r.get('cd'),'params':r['params']})+'\n'); logfh.flush()
        return res
    allr=run_batch([SEED]+[sample(rng) for _ in range(N)],'rand'); allr.sort(key=score,reverse=True)
    b=allr[0]; print(f"RANDOM best: {b['passed']}/32 CyclinD1 {b.get('cd'):.2f}", flush=True)
    for rd in range(CLIMB):
        base=allr[0]['params']; jit=[]
        for _ in range(WORKERS*4):
            p=dict(base)
            for n,lo,hi,lg in SPACE:
                if rng.random()<0.5: p[n]=min(hi,max(lo,p[n]*rng.uniform(0.8,1.25)))
            jit.append(p)
        allr=sorted(allr+run_batch(jit,f'climb{rd}'),key=score,reverse=True)
        b=allr[0]; print(f"CLIMB {rd}: {b['passed']}/32 CyclinD1 {b.get('cd'):.2f}", flush=True)
    b=allr[0]; json.dump({'passed':b['passed'],'cd':b.get('cd'),'params':b['params']}, open(BEST,'w'), indent=2)
    print(f"\n=== BEST CDK6-arm: {b['passed']}/32  CyclinD1 {b.get('cd'):.2f} ===\n{json.dumps(b['params'],indent=2)}")
