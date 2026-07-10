"""Overnight optimizer for the Gli->Jmjd3 eraser reconception (JP, 2026-07-09).

New biology (Shi 2014 + JP): the CyclinD1-promoter H3K27me3 is a WRITER(EZH2, cycle-gated) vs ERASER
(Gli->Jmjd3/Kdm6b, mitogen-gated) race. High-Gli MB actively ERASES the mark -> ChIP shows MB has ~HALF
the GNP H3K27me3 on CyclinD1 (opposite the old model). Two new model params (default 0, validation-preserving):
  k_jmjd3_gli  -- Gli-driven demethylase rate (makes MB mark drop)
  w_ezdir      -- weight of EZH2-DIRECT vs mark-Hill CyclinD1 repression (moves the 5.07 fold off the mark,
                  so the mark can be LOW in MB without breaking the dose fold)

Goal: find params where (i) MB/GNP promoter H3K27me3 ~ 0.5 (the ChIP) AND (ii) all non-soft validation
targets still pass (CyclinD1 MB/GNP 5.07, EZH2i rescue, transcript S/G0, arrests, folds...).

Minimizes  loss = 10*n_hardfail + 8*|chip - 0.5|   over the mark+eraser+direct-blend params.
Subprocess-drives the exact validate_v44 (VALIDATE_PARAMS/VALIDATE_RESULTS). Random sample -> hill-climb.

Run:  ./venv/bin/python simulations/optimize_jmjd3_chip.py [--rand=320] [--workers=8] [--climb=5]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='jmjd3opt_')
BEST = os.path.join(ROOT, 'simulations', 'optimize_jmjd3_chip_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optimize_jmjd3_chip_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
CHIP_TARGET = 0.5                       # measured MB/GNP promoter H3K27me3 (~half)

SPACE = [
    ('k_jmjd3_gli',       0.010, 0.150, True),
    ('w_ezdir',           0.00,  1.00,  False),
    ('del_mk',            0.0012, 0.006, True),
    ('k_w_mk',            0.002, 0.020, True),
    ('k0_mk',             0.0002, 0.0015, True),
    ('g_mk',              0.05,  0.50,  True),
    ('K_mk',              0.12,  0.55,  False),
    ('n_mk',              1.20,  3.50,  False),
    ('f0_mk',             0.05,  0.30,  False),
    ('K_EZH2_repression', 0.20,  1.20,  True),
]
ANCHOR = {'k_jmjd3_gli': 0.06, 'w_ezdir': 0.7, 'del_mk': 0.00328, 'k_w_mk': 0.0045, 'k0_mk': 0.00054,
          'g_mk': 0.1055, 'K_mk': 0.303, 'n_mk': 1.676, 'f0_mk': 0.145, 'K_EZH2_repression': 0.5254}
_ctr = [0]


def sample(rng):
    p = {}
    for name, lo, hi, lg in SPACE:
        p[name] = float(np.exp(rng.uniform(np.log(lo), np.log(hi)))) if lg else float(rng.uniform(lo, hi))
    return p


def perturb(p, rng, scale):
    q = dict(p)
    for name, lo, hi, lg in SPACE:
        if lg:
            q[name] = float(np.clip(p[name] * np.exp(rng.normal(0, scale)), lo, hi))
        else:
            q[name] = float(np.clip(p[name] + rng.normal(0, scale * (hi - lo)), lo, hi))
    return q


def evaluate(params):
    i = _ctr[0]; _ctr[0] += 1
    rp = os.path.join(TMP, f'r_{i % 4096}.json')
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=420)
        with open(rp) as fh:
            out = json.load(fh)
    except Exception:
        return None
    checks = {c['name']: c for c in out.get('checks', [])}
    marks = out.get('marks', {})
    if 'CyclinD1 MB/GNP' not in checks or 'MB' not in marks or 'GNP + SHH' not in marks:
        return None
    gnp_mk = float(marks['GNP + SHH']); mb_mk = float(marks['MB'])
    if gnp_mk <= 1e-4:
        return None
    chip = mb_mk / gnp_mk
    n_hardfail = sum(1 for n, c in checks.items() if n not in EXCLUDE and not c['pass'])
    # primary: pass validation + hit ChIP ratio. tiebreak (0.3*w_ezdir): prefer solutions that keep the
    # mark FUNCTIONAL (lower w_ezdir = mark still represses CyclinD1 / gates threshold) over pure-readout.
    loss = 10.0 * n_hardfail + 8.0 * abs(chip - CHIP_TARGET) + 0.3 * float(params.get('w_ezdir', 0.0))
    good = (n_hardfail == 0 and 0.40 <= chip <= 0.60)
    return dict(params=params, loss=float(loss), n_hardfail=int(n_hardfail), chip=float(chip),
                good=bool(good), passed=int(out.get('passed', 0)),
                mbgnp_cd=float(checks['CyclinD1 MB/GNP']['actual']),
                ezi=float(checks.get('EZH2i CycD1 fold (GNP)', {}).get('actual', 0)),
                transcript=float(checks.get('EZH2 transcript S/G0 (1.8-2.5)', {}).get('actual', 0)),
                gnp_mk=gnp_mk, mb_mk=mb_mk)


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='):
            return type(d)(a.split('=', 1)[1])
    return d


def batch(cands, pool, logfh, tag):
    res = []
    futs = {pool.submit(evaluate, c): c for c in cands}
    done = 0
    for fut in as_completed(futs):
        r = fut.result(); done += 1
        if r is None:
            continue
        logfh.write(json.dumps({'tag': tag, **{k: r[k] for k in ('loss', 'n_hardfail', 'chip', 'good', 'passed', 'mbgnp_cd', 'params')}}) + '\n'); logfh.flush()
        res.append(r)
        if done % 25 == 0:
            good = [x for x in res if x['good']]
            bl = min((x['loss'] for x in res), default=99)
            print(f'  [{tag}] {done}/{len(cands)} eval, {len(good)} GOOD, best loss {bl:.3f}', flush=True)
    return res


if __name__ == '__main__':
    N = _arg('--rand', 320); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 5)
    rng = np.random.default_rng(20260710)
    print(f'JMJD3/ChIP optimizer: {N} random + {CLIMB} climb, {WORKERS} workers. target MB/GNP mark={CHIP_TARGET}', flush=True)
    print(f'  loss = 10*n_hardfail + 8*|chip-0.5| ; soft/excluded: {sorted(EXCLUDE)}', flush=True)
    logfh = open(LOG, 'w')
    allr = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        cands = [ANCHOR] + [sample(rng) for _ in range(N)]
        allr = batch(cands, pool, logfh, 'rand')
        allr.sort(key=lambda x: x['loss'])
        if allr:
            b = allr[0]
            print(f"RANDOM best: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} MB/GNP_cd {b['mbgnp_cd']:.2f} passed {b['passed']}", flush=True)
        for rd in range(CLIMB):
            seeds = allr[:6] if allr else [{'params': ANCHOR}]
            scale = 0.28 * (0.6 ** rd)
            cands = [perturb(s['params'], rng, scale) for s in seeds for _ in range(8)]
            res = batch(cands, pool, logfh, f'climb{rd}')
            allr = sorted(allr + res, key=lambda x: x['loss'])
            b = allr[0]
            ng = sum(1 for x in allr if x['good'])
            print(f"CLIMB {rd}: best loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} MB/GNP_cd {b['mbgnp_cd']:.2f} | {ng} GOOD so far", flush=True)
    logfh.close()
    if not allr:
        print('NO RESULTS'); sys.exit(0)
    with open(BEST, 'w') as fh:
        json.dump(dict(best=allr[0], top12=allr[:12]), fh, indent=2)
    b = allr[0]
    print('\n=== BEST (min loss) ===')
    print(f"  loss {b['loss']:.3f} | n_hardfail {b['n_hardfail']} | passed {b['passed']}/27 | good={b['good']}")
    print(f"  ChIP MB/GNP mark {b['chip']:.3f} (target 0.5; GNP {b['gnp_mk']:.3f} MB {b['mb_mk']:.3f})")
    print(f"  CyclinD1 MB/GNP {b['mbgnp_cd']:.2f} (5.07) | EZH2i fold {b['ezi']:.2f} (2.2) | transcript {b['transcript']:.2f}")
    print(f"  NEW params: k_jmjd3_gli {b['params']['k_jmjd3_gli']:.4f} (eraser) | w_ezdir {b['params']['w_ezdir']:.3f} (0=mark-only,1=EZH2-direct)")
    print(f"  full params = {json.dumps(b['params'])}")
    gd = [x for x in allr if x['good']]
    if gd:
        wz = sorted(x['params']['w_ezdir'] for x in gd)
        print(f"  {len(gd)} GOOD solutions; w_ezdir range [{wz[0]:.2f}, {wz[-1]:.2f}] median {wz[len(wz)//2]:.2f} (lower = mark keeps CyclinD1/threshold function)")
    print(f"  -> {BEST}")
