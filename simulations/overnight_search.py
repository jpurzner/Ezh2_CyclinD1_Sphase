"""Overnight WIDE parameter search over EZH2 + H3K27me3-mark parameters, scored against the EXACT
current validate_v44 (all 27 targets, via VALIDATE_PARAMS env override). Goal: push past 24/27 --
recover CyclinD1 MB+HHi/MB WITHOUT losing EZH2 transcript S/G0 or any other passing target.

Hill-climb around the best-so-far + periodic wide random restarts; parallel subprocess evals; writes
best.json (best-so-far) + jsonl (every eval) so progress survives a crash and is reviewable any time.

Run (background):  ./venv/bin/python simulations/overnight_search.py [hours] [workers]
Review:            cat simulations/overnight_search_best.json
"""
import sys, os, json, time, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor
import re
import numpy as np

HOURS = float(sys.argv[1]) if len(sys.argv) > 1 else 7.5
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 7
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(REPO, 'venv', 'bin', 'python')
VAL = os.path.join('simulations', 'validate_v44.py')
BEST_F = os.path.join(REPO, 'simulations', 'overnight_search_best.json')
LOG_F = os.path.join(REPO, 'simulations', 'overnight_search.jsonl')

# param -> (seed, lo, hi) ; log-uniform search. Seed = current committed 24/27 model.
PARAMS = {
    'kEZbas':            (0.00041,  0.0001, 0.005),
    'kEZbas_Cd':         (0.0014,   0.0002, 0.012),
    'kEZE2f':            (0.0059,   0.002,  0.05),
    'kDeEZ':             (0.00084,  0.0002, 0.0015),
    'kTlEZ':             (0.0083,   0.002,  0.025),
    'kDeEZm':            (0.0098,   0.003,  0.03),
    'Kez_cd':            (2.94,     0.5,    8.0),
    'K_EZH2_repression': (0.539,    0.1,    2.0),
    'K_mk':              (0.305,    0.05,   0.9),
    'n_mk':              (4.15,     1.5,    8.0),
    'f0_mk':             (0.233,    0.05,   0.65),
    'k_w_mk':            (0.0016,   0.0003, 0.012),
    'k0_mk':             (0.000258, 5e-5,   0.002),
    'g_mk':              (0.3,      0.05,   0.85),
    'K_tx_mk':           (5.0,      1.0,    20.0),
    'p_tx_mk':           (2.0,      1.0,    5.0),
    'del_mk':            (0.0007,   0.0001, 0.006),
}
NAMES = list(PARAMS)
SEED = {k: v[0] for k, v in PARAMS.items()}
SUM_RE = re.compile(r'VALIDATION SUMMARY:\s+(\d+)/(\d+)')
TGT_RE = re.compile(r'\[[+-]\]\s+(.+?)\s*:\s+([-\d.]+)\s+\(target\s+([-\d.]+)')


def evaluate(params):
    """Run validate_v44 with these params; return (npass, relerr, breakdown) or (-1, 1e9, {})."""
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), H3K27_DILUTION='1')
    tf = tempfile.NamedTemporaryFile(prefix='vr_', suffix='.json', delete=False, dir='/tmp')
    tf.close(); env['VALIDATE_RESULTS'] = tf.name
    try:
        r = subprocess.run([PY, VAL], cwd=REPO, env=env, capture_output=True, text=True, timeout=240)
        out = r.stdout
    except Exception:
        out = ''
    finally:
        try: os.unlink(tf.name)
        except Exception: pass
    m = SUM_RE.search(out)
    if not m:
        return (-1, 1e9, {})
    npass = int(m.group(1))
    relerr = 0.0; brk = {}
    for name, val, tgt in TGT_RE.findall(out):
        try:
            v = float(val); t = float(tgt); brk[name.strip()] = (v, t)
            if abs(t) > 1e-6:
                relerr += abs(v - t) / abs(t)
        except Exception:
            pass
    return (npass, relerr, brk)


def perturb(base, rng, pmut=0.4, sigma=0.3):
    p = dict(base)
    for k in NAMES:
        if rng.random() < pmut:
            _, lo, hi = PARAMS[k]
            p[k] = float(np.clip(p[k] * np.exp(rng.normal(0, sigma)), lo, hi))
    return p


def restart(rng):
    return {k: float(np.exp(rng.uniform(np.log(PARAMS[k][1]), np.log(PARAMS[k][2])))) for k in NAMES}


def save_best(best):
    tmp = BEST_F + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(best, f, indent=2)
    os.replace(tmp, BEST_F)


if __name__ == '__main__':
    rng = np.random.default_rng(20260708)
    t0 = time.time(); deadline = t0 + HOURS * 3600
    open(LOG_F, 'w').close()
    print(f'overnight search: {HOURS}h, {WORKERS} workers, {len(NAMES)} params', flush=True)

    np0, re0, brk0 = evaluate(SEED)
    best = dict(npass=np0, relerr=re0, params=SEED, breakdown=brk0, evals=1, found_at=0.0)
    save_best(best)
    print(f'seed: {np0}/27  relerr={re0:.3f}', flush=True)

    n_eval = 1; rounds = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        while time.time() < deadline:
            rounds += 1
            cands = []
            for _ in range(WORKERS):
                cands.append(restart(rng) if rng.random() < 0.15 else perturb(best['params'], rng))
            for params, (npass, relerr, brk) in zip(cands, pool.map(evaluate, cands)):
                n_eval += 1
                with open(LOG_F, 'a') as f:
                    f.write(json.dumps(dict(npass=npass, relerr=round(relerr, 4), params=params)) + '\n')
                better = (npass > best['npass']) or (npass == best['npass'] and relerr < best['relerr'] - 1e-9)
                if npass >= 0 and better:
                    best = dict(npass=npass, relerr=relerr, params=params, breakdown=brk,
                                evals=n_eval, found_at=round((time.time() - t0) / 3600, 3))
                    save_best(best)
                    fails = [n for n, (v, t) in brk.items() if abs(v - t) > 0.001 and (t == 0 or abs(v - t) / abs(t) > 1e-6)]
                    print(f'[{best["found_at"]:.2f}h #{n_eval}] NEW BEST {npass}/27 relerr={relerr:.3f}  fails={fails}', flush=True)
            if rounds % 20 == 0:
                print(f'[{(time.time()-t0)/3600:.2f}h] {n_eval} evals, best {best["npass"]}/27 (relerr {best["relerr"]:.3f})', flush=True)

    print(f'DONE: {n_eval} evals in {(time.time()-t0)/3600:.2f}h. Best {best["npass"]}/27 relerr={best["relerr"]:.3f}', flush=True)
    print('best params:', json.dumps(best['params']), flush=True)
