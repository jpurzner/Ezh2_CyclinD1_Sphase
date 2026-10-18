"""Sensitivity of the 29 validation targets to each load-bearing parameter (perturb +25%, re-validate,
elasticity = dlog(target_actual)/dlog(param)). Baked-default model (cell-type split). Parallel subprocess."""
import os, sys, json, subprocess, tempfile, numpy as np, time
from concurrent.futures import ThreadPoolExecutor
ROOT='/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase'
PY=os.path.join(ROOT,'venv/bin/python'); VAL=os.path.join(ROOT,'simulations/validate_v44.py')
OUT='/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/81291f88-7d19-4149-be27-53b1162ee52d/scratchpad'
# load-bearing params -> landed value (from _ts_bake + the split)
LANDED=dict(mu=0.000769,k_mu_cki=0.307,kPhRbCd=0.1433,K_CdRb=0.3284,k_Cd_tx_Gli_max=0.143566,
 k_Cd_tx_basal=0.000236,k_Cd_tx_MYCN=0.031157,M_commit=2.476413,kDeP21Cd=0.630028,kDsRbmE2f=1.882106,
 f_commit_carry=0.374502,kSyDna=0.046863,M_size=3.140623,kEZbas=0.000622,kEZbas_Cd=0.000318,
 kEZE2f=0.012264,Kez_cd=9.175748,kDeEZ=0.000277,kDeEZ_C1=0.000386,K_E2f_EZ=0.449773,k_jmjd3_gli=0.005,
 n_rw=2,k_dil_S=0.693,k_h33=0.003,a0_prc2=0.002875,a_rw_prc2=0.002266,g_prc2=0.070491,K_prc2=0.006857,
 n_prc2=4.015745,f0_prc2=0.174826,del_mk=0.00044,KmHU_fire=0.279729,kSySkp2_Cd=0.244571,
 K_Skp2_Cd=10.296078,mu_min_frac=0.4,K_g1len=0.6,k_Cd_mRNA_deg=0.001)
def validate(params):
    tf=tempfile.NamedTemporaryFile('w',suffix='.json',delete=False); tf.close()
    env=dict(os.environ,VALIDATE_PARAMS=json.dumps(params),VALIDATE_RESULTS=tf.name,
             TWO_STEP_RB='1',H3K27_CHAIN='1',EZH2_CONC='0')  # DECOUPLE_COMMIT default '1' now
    try:
        subprocess.run([PY,VAL],cwd=ROOT,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=400)
        R=json.load(open(tf.name)); return {c['name']:c['actual'] for c in R['checks']}
    except Exception: return None
    finally:
        try: os.unlink(tf.name)
        except: pass
def perturb(name):
    p=dict(LANDED); p[name]=LANDED[name]*1.25
    return name, validate({name:p[name]})
if __name__=='__main__':
    t0=time.time()
    base=validate({}); assert base, 'baseline failed'
    names=list(LANDED); rows={}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for i,(name,act) in enumerate(ex.map(perturb,names)):
            rows[name]=act
            if i%8==0: print(f' [{i}/{len(names)}] {name} ({time.time()-t0:.0f}s)',flush=True)
    # elasticity: (actual_pert - actual_base)/|actual_base| / 0.25 ; binaries -> flip (sign of change, cap 1)
    tgts=list(base)
    E={}
    for name in names:
        act=rows[name]
        if act is None: E[name]={t:0.0 for t in tgts}; continue
        e={}
        for t in tgts:
            b=base[t]; a=act.get(t,b)
            if abs(b)<1e-9: e[t]=float(np.clip(a-b,-1,1))   # binary/zero target -> raw change
            else: e[t]=float((a-b)/abs(b)/0.25)
        E[name]=e
    json.dump({'landed':LANDED,'base':base,'elasticity':E,'targets':tgts},open(OUT+'/sensitivity.json','w'),indent=0)
    print('saved sensitivity.json; base passed=',sum(1 for t in tgts if t)); print('done',f'{time.time()-t0:.0f}s')
