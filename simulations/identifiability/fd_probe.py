"""Finite-difference step-size convergence probe.

Question: are the validation outputs smooth enough in log-parameter for central differences
to resolve sensitivities, or does event-driven / division-counting non-smoothness dominate?

Perturb a few representative parameters by relative steps h = 0.2, 0.1, 0.05, 0.02, 0.01, 0.005
(central), and report the elasticity estimate at each h. A smooth output gives a plateau;
a non-smooth one wanders or blows up as h shrinks.
"""
import os, sys, json, subprocess, tempfile, time
from concurrent.futures import ThreadPoolExecutor

ROOT = '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase'
PY = os.path.join(ROOT, 'venv/bin/python')
VAL = os.path.join(ROOT, 'simulations/validate_v44.py')
OUT = os.path.dirname(os.path.abspath(__file__))

# landed values, from the prior session's sensitivity.py (still the baked defaults)
LANDED = dict(kSyDna=0.046863, kPhRbCd=0.1433, k_Cd_tx_Gli_max=0.143566,
              kDeEZ=0.000277, del_mk=0.00044, K_prc2=0.006857)

STEPS = [0.20, 0.10, 0.05, 0.02, 0.01, 0.005]


def validate(params):
    tf = tempfile.NamedTemporaryFile('w', suffix='.json', delete=False); tf.close()
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=tf.name,
               TWO_STEP_RB='1', H3K27_CHAIN='1', EZH2_CONC='0')
    try:
        subprocess.run([PY, VAL], cwd=ROOT, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=600)
        R = json.load(open(tf.name))
        return {c['name']: c['actual'] for c in R['checks']}
    except Exception as e:
        print('  FAILED', params, e, flush=True)
        return None
    finally:
        try: os.unlink(tf.name)
        except OSError: pass


def job(spec):
    name, h, sign = spec
    p = {name: LANDED[name] * (1.0 + sign * h)}
    return spec, validate(p)


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(n, h, s) for n in LANDED for h in STEPS for s in (+1, -1)]
    print(f'{len(jobs)} runs', flush=True)
    res = {}
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, (spec, act) in enumerate(ex.map(job, jobs)):
            res[spec] = act
            if i % 10 == 0:
                print(f'  [{i}/{len(jobs)}] {time.time()-t0:.0f}s', flush=True)
    base = json.load(open(os.path.join(OUT, 'val_base.json')))
    base = {c['name']: c['actual'] for c in base['checks']}
    out = {}
    for n in LANDED:
        out[n] = {}
        for h in STEPS:
            up, dn = res.get((n, h, +1)), res.get((n, h, -1))
            if up is None or dn is None:
                continue
            e = {}
            for t in base:
                b = base[t]
                if abs(b) < 1e-9:
                    continue
                # central log-elasticity: dlog(y)/dlog(p) ~ (y+ - y-)/(2 h b)
                e[t] = (up.get(t, b) - dn.get(t, b)) / (2 * h * b)
            out[n][h] = e
    json.dump({'base': base, 'central': out, 'landed': LANDED, 'steps': STEPS},
              open(os.path.join(OUT, 'fd_probe.json'), 'w'), indent=0)
    print('saved fd_probe.json', f'{time.time()-t0:.0f}s')
