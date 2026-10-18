"""Verify the claimed EXACT scaling degeneracies in the PRC2/EZH2 sector.

Claim A (PRC2 amplitude):  {a0_prc2, a_rw_prc2} x L,  K_prc2 x L,  {kme1,kme2,kme3} / L
Claim B (EZH2 protein):    kTlEZ x L,  {a0_prc2, a_rw_prc2} / L
Claim C (EZH2 mRNA):       {kEZbas, kEZbas_Cd, kEZE2f} x L,  kTlEZ / L
Control:                   a0_prc2 x L alone (should NOT be invariant)

If exact, every one of the 32 validation outputs is bit-comparable to baseline.
"""
import os, sys, json, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase'
PY = os.path.join(ROOT, 'venv/bin/python')
VAL = os.path.join(ROOT, 'simulations/validate_v44.py')

BAKE = dict(a0_prc2=0.0028749753507882796, a_rw_prc2=0.0022657109194688194,
            K_prc2=0.006856946735209077, kme1=9.706031386488556,
            kme2=8.266881635062989, kme3=4.555295001018562,
            kTlEZ=0.015098862630302894, kEZbas=0.0006220377408563062,
            kEZbas_Cd=0.00031843672436329825, kEZE2f=0.01226393864518667)


def validate(params):
    tf = tempfile.NamedTemporaryFile('w', suffix='.json', delete=False); tf.close()
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=tf.name,
               TWO_STEP_RB='1', H3K27_CHAIN='1', EZH2_CONC='0')
    try:
        subprocess.run([PY, VAL], cwd=ROOT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)
        R = json.load(open(tf.name))
        return {c['name']: c['actual'] for c in R['checks']}
    except Exception as e:
        return {'__error__': repr(e)}
    finally:
        try: os.unlink(tf.name)
        except OSError: pass


def make(kind, L):
    p = {}
    if kind == 'A':
        for k in ('a0_prc2', 'a_rw_prc2', 'K_prc2'): p[k] = BAKE[k] * L
        for k in ('kme1', 'kme2', 'kme3'):           p[k] = BAKE[k] / L
    elif kind == 'B':
        p['kTlEZ'] = BAKE['kTlEZ'] * L
        for k in ('a0_prc2', 'a_rw_prc2'):           p[k] = BAKE[k] / L
    elif kind == 'C':
        for k in ('kEZbas', 'kEZbas_Cd', 'kEZE2f'):  p[k] = BAKE[k] * L
        p['kTlEZ'] = BAKE['kTlEZ'] / L
    elif kind == 'ctrl':
        p['a0_prc2'] = BAKE['a0_prc2'] * L
    return p


CASES = [(k, L) for k in ('A', 'B', 'C', 'ctrl') for L in (2.0, 0.5)]

if __name__ == '__main__':
    base = validate({})
    with ThreadPoolExecutor(max_workers=8) as ex:
        outs = list(ex.map(lambda c: validate(make(*c)), CASES))
    rep = {}
    print(f"{'case':10s} {'maxreldiff':>12s}  {'n>1e-6':>7s}  worst target")
    for (k, L), o in zip(CASES, outs):
        if '__error__' in o:
            print(f'{k}(L={L}) ERROR {o["__error__"]}'); continue
        d = {t: abs(o.get(t, base[t]) - base[t]) / max(abs(base[t]), 1e-12) for t in base}
        worst = max(d, key=d.get)
        n = sum(1 for v in d.values() if v > 1e-6)
        print(f'{k}(L={L:<4}) {d[worst]:12.3e}  {n:7d}  {worst}')
        rep[f'{k}_{L}'] = {'max_rel_diff': d[worst], 'worst': worst, 'n_changed': n}
    json.dump({'base': base, 'report': rep}, open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'degen_test.json'), 'w'), indent=1)
    print('\nsaved degen_test.json')
