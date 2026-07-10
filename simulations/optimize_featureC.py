"""Feature C joint optimizer (2026-07-10): commitment carryover + cell-cycle + EZH2 gating + PRC2.

Goal: with the Spencer commitment carryover ENGAGED (f_commit_carry >= 0.3, so committed cells retain
E2f in G1), co-tune the cell-cycle (kSyDna, M_commit), EZH2 E2f-gating (kEZbas/kEZbas_Cd/kEZE2f/Kez_cd/
K_E2f_EZ), and PRC2 (a0/a_rw/g/K/n/f0_prc2, del_mk, k_jmjd3_gli) params to land ALL of:
  - Point 3: CyclinD1-controlled commitment (carryover on)
  - Feature B: EZH2 Palbo mRNA drop ~0.44 (now decouple-able: cycling-G0 keeps E2f, only sustained arrest->0)
  - the existing anchors: ChIP MB/GNP mark 0.5, CyclinD1 MB/GNP 5.07, EZH2i 2.2, transcript S/G0 ~2, phase fractions

Loss = 10*n_hardfail + 6*|chip-.5| + 2*|CD-5.07|/5.07 + 2*|EZi-2.2|/2.2 + 2*|palbo-.44|/.44 + 1.5*|transcript-2|/2.
Subprocess-drives validate_v44 (28 targets). Random -> hill-climb.

Run:  ./venv/bin/python simulations/optimize_featureC.py [--rand=500] [--workers=8] [--climb=6]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='featC_')
BEST = os.path.join(ROOT, 'simulations', 'optimize_featureC_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optimize_featureC_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
CHIP_TARGET = 0.5

SPACE = [
    ('f_commit_carry', 0.30, 0.80, False),   # Feature C: commitment carryover (ENGAGED)
    ('kSyDna',   0.042, 0.095, True),        # S-phase duration FLOORED higher (BrdU Ts~3h): low kSyDna lengthens S -> inflates MB S count%
    ('M_commit', 0.80,  1.50,  False),       # G0/G1 commitment size threshold
    ('kEZbas',   0.0001, 0.0015, True),      # EZH2 transcription: shift toward E2f-gating for the Palbo drop
    ('kEZbas_Cd', 0.0002, 0.0025, True),
    ('kEZE2f',   0.004,  0.030,  True),
    ('Kez_cd',   2.0,    10.0,   True),
    ('K_E2f_EZ', 0.15,   0.50,   False),
    ('k_jmjd3_gli', 0.02,  0.15,  True),     # PRC2 / eraser
    ('a0_prc2',  0.0001, 0.0025, True),      # accessory recruitment (mark-INDEPENDENT) -- push DOWN
    ('a_rw_prc2', 0.002, 0.018,  True),      # EED READ-WRITE amplification -- push UP (mark carries the repression)
    ('g_prc2',   0.02,   0.40,   True),
    ('K_prc2',   0.0003, 0.006,  True),
    ('n_prc2',   1.50,   4.00,   False),
    ('f0_prc2',  0.05,   0.30,   False),
    ('del_mk',   0.0015, 0.006,  True),
]
ANCHOR = {'f_commit_carry': 0.372, 'kSyDna': 0.052, 'M_commit': 1.117,
          'kEZbas': 0.000402, 'kEZbas_Cd': 0.000382, 'kEZE2f': 0.00575, 'Kez_cd': 3.535, 'K_E2f_EZ': 0.283,
          'k_jmjd3_gli': 0.1094, 'a0_prc2': 0.0004, 'a_rw_prc2': 0.007, 'g_prc2': 0.0919,
          'K_prc2': 0.001060, 'n_prc2': 1.643, 'f0_prc2': 0.1935, 'del_mk': 0.002063}
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
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=480)
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
    transcript = float(checks.get('EZH2 transcript S/G0 (1.8-2.5)', {}).get('actual', 0) or 0)
    # read-write DOMINANCE: fraction of GNP PRC2 occupancy carried by the read-write (a_rw*Mk) vs
    # accessory (a0) arm. High rw_frac = the H3K27me3 mark is a functional AMPLIFIER of EZH2->CyclinD1
    # (sharper gain, raises the entry threshold), not a passive readout. Reward rw_frac -> >=0.55.
    a0v = float(params['a0_prc2']); arwv = float(params['a_rw_prc2'])
    rw_frac = arwv * gnp_mk / (a0v + arwv * gnp_mk) if (a0v + arwv * gnp_mk) > 0 else 0.0
    mbs = float(checks.get('MB S count% (flow; BrdU Ts~3h)', {}).get('actual', 0))
    n_hardfail = sum(1 for n, c in checks.items() if n not in EXCLUDE and not c['pass'])
    loss = (10.0 * n_hardfail + 6.0 * abs(chip - CHIP_TARGET) + 2.0 * abs(cd - 5.07) / 5.07
            + 2.0 * abs(ezi - 2.2) / 2.2 + 2.0 * abs(palbo - 0.44) / 0.44
            + 1.5 * abs(transcript - 2.0) / 2.0 + 4.0 * max(0.0, 0.55 - rw_frac)
            + 1.5 * abs(mbs - 15.7) / 15.7)
    good = (n_hardfail == 0 and 0.40 <= chip <= 0.60 and rw_frac >= 0.45)
    return dict(params=params, loss=float(loss), n_hardfail=int(n_hardfail), chip=float(chip), good=bool(good),
                passed=int(out.get('passed', 0)), mbgnp_cd=cd, ezi=ezi, palbo=palbo, transcript=transcript,
                rw_frac=float(rw_frac), mbs=mbs,
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
        logfh.write(json.dumps({'tag': tag, **{k: r[k] for k in ('loss', 'n_hardfail', 'chip', 'good', 'passed', 'mbgnp_cd', 'ezi', 'palbo', 'transcript', 'rw_frac', 'mbs', 'params')}}) + '\n'); logfh.flush()
        res.append(r)
        if done % 30 == 0:
            print(f"  [{tag}] {done}/{len(cands)} eval, {sum(1 for x in res if x['good'])} GOOD, best loss {min((x['loss'] for x in res), default=99):.3f}", flush=True)
    return res


if __name__ == '__main__':
    N = _arg('--rand', 500); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 6)
    rng = np.random.default_rng(20260710)
    print(f'Feature C joint optimizer: {N} random + {CLIMB} climb, {WORKERS} workers ({len(SPACE)} params)', flush=True)
    logfh = open(LOG, 'w')
    allr = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        allr = batch([ANCHOR] + [sample(rng) for _ in range(N)], pool, logfh, 'rand')
        allr.sort(key=lambda x: x['loss'])
        if allr:
            b = allr[0]; print(f"RANDOM best: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} EZi {b['ezi']:.2f} palbo {b['palbo']:.2f} txn {b['transcript']:.2f} fc {b['params']['f_commit_carry']:.2f}", flush=True)
        for rd in range(CLIMB):
            seeds = allr[:6] if allr else [{'params': ANCHOR}]
            cands = [perturb(s['params'], rng, 0.28 * (0.6 ** rd)) for s in seeds for _ in range(8)]
            allr = sorted(allr + batch(cands, pool, logfh, f'climb{rd}'), key=lambda x: x['loss'])
            b = allr[0]
            print(f"CLIMB {rd}: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} EZi {b['ezi']:.2f} palbo {b['palbo']:.2f} txn {b['transcript']:.2f} S% {b['mbs']:.1f} fc {b['params']['f_commit_carry']:.2f} | {sum(1 for x in allr if x['good'])} GOOD", flush=True)
    logfh.close()
    if not allr:
        print('NO RESULTS'); sys.exit(0)
    with open(BEST, 'w') as fh:
        json.dump(dict(best=allr[0], top12=allr[:12]), fh, indent=2)
    b = allr[0]
    print('\n=== BEST (min loss) ===')
    print(f"  loss {b['loss']:.3f} | n_hardfail {b['n_hardfail']} | passed {b['passed']}/28 | good={b['good']}")
    print(f"  ChIP {b['chip']:.3f} | CyclinD1 {b['mbgnp_cd']:.2f} | EZH2i {b['ezi']:.2f} | Palbo {b['palbo']:.2f} | transcript {b['transcript']:.2f} | MB S% {b['mbs']:.1f} | f_commit_carry {b['params']['f_commit_carry']:.2f}")
    print(f"  read-write fraction {b['rw_frac']:.2f} (mark's share of PRC2 occupancy; >0.5 = mark amplifies EZH2's impact) | a0 {b['params']['a0_prc2']:.5f} a_rw {b['params']['a_rw_prc2']:.5f}")
    print(f"  params = {json.dumps(b['params'])}")
    print(f"  -> {BEST}")
