"""D (CyclinD1 mRNA scale) WITHOUT the Cd_mRNA initial condition (unreachable via the override regex).
If the degeneracy is exact, the only residual should be the decaying initial transient."""
import os, json, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor
ROOT='/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase'
PY=os.path.join(ROOT,'venv/bin/python'); VAL=os.path.join(ROOT,'simulations/validate_v44.py')
V=dict(k_Cd_tx_basal=0.000529,k_Cd_tx_Gli_max=0.152118,k_Cd_tx_MYCN=0.031157,
       k_Cd_translation=0.801,K_tx_mk=5.119476334025431)
def validate(p):
    tf=tempfile.NamedTemporaryFile('w',suffix='.json',delete=False); tf.close()
    env=dict(os.environ,VALIDATE_PARAMS=json.dumps(p),VALIDATE_RESULTS=tf.name,
             TWO_STEP_RB='1',H3K27_CHAIN='1',EZH2_CONC='0')
    try:
        subprocess.run([PY,VAL],cwd=ROOT,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=900)
        return {c['name']:c['actual'] for c in json.load(open(tf.name))['checks']}
    except Exception as e: return {'__error__':repr(e)}
    finally:
        try: os.unlink(tf.name)
        except OSError: pass
def make(L):
    p={k:V[k]*L for k in ('k_Cd_tx_basal','k_Cd_tx_Gli_max','k_Cd_tx_MYCN')}
    p['k_Cd_translation']=V['k_Cd_translation']/L; p['K_tx_mk']=V['K_tx_mk']*L
    return p
base=validate({})
with ThreadPoolExecutor(max_workers=4) as ex: outs=list(ex.map(lambda L: validate(make(L)),(2.0,0.5)))
for L,o in zip((2.0,0.5),outs):
    if '__error__' in o: print('ERROR',o); continue
    d={t:abs(o.get(t,base[t])-base[t])/max(abs(base[t]),1e-12) for t in base}
    w=max(d,key=d.get)
    print(f'D(L={L}) maxreldiff {d[w]:.3e}  n>1e-6 {sum(1 for v in d.values() if v>1e-6)}  worst {w}')
