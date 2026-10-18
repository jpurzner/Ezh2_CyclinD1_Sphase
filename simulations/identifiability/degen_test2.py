"""Re-test degeneracy claims B and C WITH the initial conditions scaled (omitted in v1).

B (EZH2 protein scale):  kTlEZ x L, EZH2(0) x L, {a0_prc2, a_rw_prc2} / L
C (EZH2 mRNA scale):     {kEZbas, kEZbas_Cd, kEZE2f} x L, EZH2m(0) x L, kTlEZ / L
D (CyclinD1 mRNA scale): {k_Cd_tx_basal, k_Cd_tx_Gli_max, k_Cd_tx_MYCN} x L, Cd_mRNA(0) x L,
                         k_Cd_translation / L, K_tx_mk x L
"""
import os, json, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase'
PY = os.path.join(ROOT, 'venv/bin/python')
VAL = os.path.join(ROOT, 'simulations/validate_v44.py')

V = dict(a0_prc2=0.0028749753507882796, a_rw_prc2=0.0022657109194688194,
         kTlEZ=0.015098862630302894, EZH2=0.5, EZH2m=0.1,
         kEZbas=0.0006220377408563062, kEZbas_Cd=0.00031843672436329825,
         kEZE2f=0.01226393864518667,
         k_Cd_tx_basal=0.000529, k_Cd_tx_Gli_max=0.152118, k_Cd_tx_MYCN=0.031157,
         Cd_mRNA=2.6, k_Cd_translation=0.801, K_tx_mk=5.119476334025431)


def validate(params):
    tf = tempfile.NamedTemporaryFile('w', suffix='.json', delete=False); tf.close()
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=tf.name,
               TWO_STEP_RB='1', H3K27_CHAIN='1', EZH2_CONC='0')
    try:
        subprocess.run([PY, VAL], cwd=ROOT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)
        return {c['name']: c['actual'] for c in json.load(open(tf.name))['checks']}
    except Exception as e:
        return {'__error__': repr(e)}
    finally:
        try: os.unlink(tf.name)
        except OSError: pass


def make(kind, L):
    p = {}
    if kind == 'B':
        p['kTlEZ'] = V['kTlEZ'] * L; p['EZH2'] = V['EZH2'] * L
        for k in ('a0_prc2', 'a_rw_prc2'): p[k] = V[k] / L
    elif kind == 'C':
        for k in ('kEZbas', 'kEZbas_Cd', 'kEZE2f'): p[k] = V[k] * L
        p['EZH2m'] = V['EZH2m'] * L; p['kTlEZ'] = V['kTlEZ'] / L
    elif kind == 'D':
        for k in ('k_Cd_tx_basal', 'k_Cd_tx_Gli_max', 'k_Cd_tx_MYCN'): p[k] = V[k] * L
        p['Cd_mRNA'] = V['Cd_mRNA'] * L
        p['k_Cd_translation'] = V['k_Cd_translation'] / L
        p['K_tx_mk'] = V['K_tx_mk'] * L
    return p


CASES = [(k, L) for k in ('B', 'C', 'D') for L in (2.0, 0.5)]

if __name__ == '__main__':
    base = validate({})
    with ThreadPoolExecutor(max_workers=6) as ex:
        outs = list(ex.map(lambda c: validate(make(*c)), CASES))
    print(f"{'case':12s} {'maxreldiff':>12s} {'n>1e-6':>7s}  worst")
    rep = {}
    for (k, L), o in zip(CASES, outs):
        if '__error__' in o:
            print(f'{k}(L={L}) ERROR {o["__error__"]}'); continue
        d = {t: abs(o.get(t, base[t]) - base[t]) / max(abs(base[t]), 1e-12) for t in base}
        w = max(d, key=d.get)
        print(f'{k}(L={L:<4})  {d[w]:12.3e} {sum(1 for v in d.values() if v>1e-6):7d}  {w}')
        rep[f'{k}_{L}'] = {'max_rel_diff': d[w], 'worst': w}
    json.dump(rep, open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     'degen_test2.json'), 'w'), indent=1)
    print('saved degen_test2.json')
