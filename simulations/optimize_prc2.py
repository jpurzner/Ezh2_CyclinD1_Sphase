"""Optimizer for the FORMAL PRC2-COMPLEX repression mechanism (JP 2026-07-10, replaces w_ezdir).

Repression is carried by PRC2 occupancy at Ccnd1: PRC2 = EZH2*(a0_prc2 + a_rw_prc2*Mk)*(1 - eviction).
The SAME PRC2 writes H3K27me3 (methylation = PRC2*(1-Mk)) and represses CyclinD1 (Hill on PRC2). The mark
is a read-write amplifier; CyclinD1 stays EZH2-responsive because all repression is PRC2-mediated; MB keeps
occupancy via high EZH2 (5.07 dose) while Gli/Jmjd3 keeps the mark low (ChIP MB/GNP~0.5).

Free: k_jmjd3_gli, a0_prc2, a_rw_prc2, g_prc2, K_prc2, n_prc2, f0_prc2, del_mk.
Loss = 10*n_hardfail + 6*|chip-0.5| + 2*|CyclinD1 MB/GNP - 5.07|/5.07 + 2*|EZH2i - 2.2|/2.2.
Subprocess-drives validate_v44 (now builds with_prc2=True by default). Random -> hill-climb.

Run:  ./venv/bin/python simulations/optimize_prc2.py [--rand=320] [--workers=8] [--climb=5]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='prc2opt_')
BEST = os.path.join(ROOT, 'simulations', 'optimize_prc2_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optimize_prc2_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
CHIP_TARGET = 0.5

SPACE = [
    ('k_jmjd3_gli', 0.02,   0.15,  True),
    ('a0_prc2',     0.0002, 0.004, True),
    ('a_rw_prc2',   0.001,  0.012, True),
    ('g_prc2',      0.02,   0.40,  True),
    ('K_prc2',      0.0005, 0.006, True),
    ('n_prc2',      1.50,   4.00,  False),
    ('f0_prc2',     0.05,   0.30,  False),
    ('del_mk',      0.0015, 0.006, True),
    # Feature B: EZH2 transcription -- shift toward more E2f-gating (writer OFF in arrest, Palbo drop ~56%)
    ('kEZbas',      0.0001, 0.0015, True),
    ('kEZbas_Cd',   0.0002, 0.0025, True),
    ('kEZE2f',      0.004,  0.025,  True),
    ('Kez_cd',      2.0,    10.0,   True),
    ('K_E2f_EZ',    0.15,   0.50,   False),
]
ANCHOR = {'k_jmjd3_gli': 0.0984, 'a0_prc2': 0.000900, 'a_rw_prc2': 0.001457, 'g_prc2': 0.2366,
          'K_prc2': 0.000868, 'n_prc2': 3.144, 'f0_prc2': 0.2302, 'del_mk': 0.005152,
          'kEZbas': 0.000517, 'kEZbas_Cd': 0.001525, 'kEZE2f': 0.006195, 'Kez_cd': 5.605, 'K_E2f_EZ': 0.3}
_ctr = [0]


def sample(rng):
    return {n: (float(np.exp(rng.uniform(np.log(lo), np.log(hi)))) if lg else float(rng.uniform(lo, hi)))
            for n, lo, hi, lg in SPACE}


def perturb(p, rng, scale):
    q = dict(p)
    for n, lo, hi, lg in SPACE:
        q[n] = float(np.clip(p[n] * np.exp(rng.normal(0, scale)), lo, hi)) if lg \
            else float(np.clip(p[n] + rng.normal(0, scale * (hi - lo)), lo, hi))
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
    cd = float(checks['CyclinD1 MB/GNP']['actual'])
    ezi = float(checks.get('EZH2i CycD1 fold (GNP)', {}).get('actual', 0))
    palbo = float(checks.get('EZH2 Palbo mRNA drop (~0.44)', {}).get('actual', 0.44))
    n_hardfail = sum(1 for n, c in checks.items() if n not in EXCLUDE and not c['pass'])
    loss = (10.0 * n_hardfail + 6.0 * abs(chip - CHIP_TARGET) + 2.0 * abs(cd - 5.07) / 5.07
            + 2.0 * abs(ezi - 2.2) / 2.2 + 1.5 * abs(palbo - 0.44) / 0.44)
    good = (n_hardfail == 0 and 0.40 <= chip <= 0.60)
    return dict(params=params, loss=float(loss), n_hardfail=int(n_hardfail), chip=float(chip), good=bool(good),
                passed=int(out.get('passed', 0)), mbgnp_cd=cd, ezi=ezi, palbo=float(palbo),
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
        logfh.write(json.dumps({'tag': tag, **{k: r[k] for k in ('loss', 'n_hardfail', 'chip', 'good', 'passed', 'mbgnp_cd', 'ezi', 'params')}}) + '\n'); logfh.flush()
        res.append(r)
        if done % 25 == 0:
            print(f"  [{tag}] {done}/{len(cands)} eval, {sum(1 for x in res if x['good'])} GOOD, best loss {min((x['loss'] for x in res), default=99):.3f}", flush=True)
    return res


if __name__ == '__main__':
    N = _arg('--rand', 320); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 5)
    rng = np.random.default_rng(20260710)
    print(f'PRC2 optimizer: {N} random + {CLIMB} climb, {WORKERS} workers. ChIP target {CHIP_TARGET}', flush=True)
    logfh = open(LOG, 'w')
    allr = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        allr = batch([ANCHOR] + [sample(rng) for _ in range(N)], pool, logfh, 'rand')
        allr.sort(key=lambda x: x['loss'])
        if allr:
            b = allr[0]; print(f"RANDOM best: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} EZi {b['ezi']:.2f}", flush=True)
        for rd in range(CLIMB):
            seeds = allr[:6] if allr else [{'params': ANCHOR}]
            cands = [perturb(s['params'], rng, 0.28 * (0.6 ** rd)) for s in seeds for _ in range(8)]
            allr = sorted(allr + batch(cands, pool, logfh, f'climb{rd}'), key=lambda x: x['loss'])
            b = allr[0]
            print(f"CLIMB {rd}: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} EZi {b['ezi']:.2f} | {sum(1 for x in allr if x['good'])} GOOD", flush=True)
    logfh.close()
    if not allr:
        print('NO RESULTS'); sys.exit(0)
    with open(BEST, 'w') as fh:
        json.dump(dict(best=allr[0], top12=allr[:12]), fh, indent=2)
    b = allr[0]
    print('\n=== BEST (min loss) ===')
    print(f"  loss {b['loss']:.3f} | n_hardfail {b['n_hardfail']} | passed {b['passed']}/27 | good={b['good']}")
    print(f"  ChIP MB/GNP mark {b['chip']:.3f} | CyclinD1 MB/GNP {b['mbgnp_cd']:.2f} (5.07) | EZH2i {b['ezi']:.2f} (2.2) | transcript {b['transcript']:.2f} | Palbo drop {b['palbo']:.2f} (0.44)")
    print(f"  params = {json.dumps(b['params'])}")
    print(f"  -> {BEST}")
