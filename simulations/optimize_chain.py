"""Joint optimizer for the SERIAL-METHYLATION chain model (with_h3k27_chain=True).

Goal: a FASTER me2->me3 (kme3 bounded to a fast range -> me3 lag ~9-14 h, not ~18 h) that STILL hits
all 28 validation targets, by co-tuning the chain rates + Gli/Jmjd3 eraser + PRC2 repression + the
cell-cycle/EZH2 params. Same subprocess-driven, band-penalty design as optimize_featureC, but
validate_v44 is run with H3K27_CHAIN=1 so it builds the explicit me1/me2/me3 chain.

Run:  ./venv/bin/python simulations/optimize_chain.py [--rand=400] [--workers=8] [--climb=7]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='chain_')
BEST = os.path.join(ROOT, 'simulations', 'optimize_chain_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optimize_chain_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
CHIP_TARGET = 0.5

SPACE = [
    ('f_commit_carry', 0.30, 0.80, False),
    ('kSyDna',   0.042, 0.095, True),
    ('M_commit', 0.80,  1.50,  False),
    ('kEZbas',   0.0001, 0.0015, True),
    ('kEZbas_Cd', 0.0002, 0.0025, True),
    ('kEZE2f',   0.004,  0.030,  True),
    ('Kez_cd',   2.0,    10.0,   True),
    ('K_E2f_EZ', 0.15,   0.50,   False),
    ('k_jmjd3_gli', 0.02,  0.30,  True),     # eraser -- WIDENED up (faster kme3 needs a stronger eraser for ChIP)
    ('a0_prc2',  0.0001, 0.0025, True),
    ('a_rw_prc2', 0.002, 0.018,  True),
    ('g_prc2',   0.02,   0.40,   True),
    ('K_prc2',   0.0003, 0.008,  True),      # repression K -- widened up (re-balance CyclinD1 fold)
    ('n_prc2',   1.50,   4.00,   False),
    ('f0_prc2',  0.05,   0.30,   False),
    ('del_mk',   0.0015, 0.006,  True),
    # --- serial-methylation chain ---
    ('kme1',     4.0,    20.0,   True),      # FAST me0->me1
    ('kme2',     4.0,    20.0,   True),      # FAST me1->me2
    ('kme3',     4.0,    11.0,   True),      # me2->me3: BOUNDED FAST (lag ~9-14 h) -- the point of this run
    ('d_me',     0.0008, 0.004,  True),      # precursor turnover
]
ANCHOR = {'f_commit_carry': 0.5672456517307928, 'kSyDna': 0.04265121846938828, 'M_commit': 1.2843242130809425,
          'kEZbas': 0.00043040856433872047, 'kEZbas_Cd': 0.0013508978696408793, 'kEZE2f': 0.021216027072906932,
          'Kez_cd': 7.830608487520997, 'K_E2f_EZ': 0.5, 'k_jmjd3_gli': 0.11,
          'a0_prc2': 0.00037763461927116913, 'a_rw_prc2': 0.004297010112265872, 'g_prc2': 0.04244445931486211,
          'K_prc2': 0.0034813553293576824, 'n_prc2': 3.385026804758643, 'f0_prc2': 0.05054636367014487,
          'del_mk': 0.0015, 'kme1': 10.0, 'kme2': 10.0, 'kme3': 6.0, 'd_me': 0.0015}
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
    env = dict(os.environ, H3K27_CHAIN='1', VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=rp)
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
    a0v = float(params['a0_prc2']); arwv = float(params['a_rw_prc2'])
    rw_frac = arwv * gnp_mk / (a0v + arwv * gnp_mk) if (a0v + arwv * gnp_mk) > 0 else 0.0
    mbs = float(checks.get('MB S count% (flow; BrdU Ts~3h)', {}).get('actual', 0))
    n_hardfail = sum(1 for n, c in checks.items() if n not in EXCLUDE and not c['pass'])
    bnd = lambda x, t, b: max(0.0, abs(x - t) / t - b)
    loss = (10.0 * n_hardfail
            + 3.0 * bnd(chip, CHIP_TARGET, 0.15) + 1.5 * bnd(cd, 5.07, 0.15)
            + 1.5 * bnd(ezi, 2.2, 0.20) + 1.5 * bnd(palbo, 0.44, 0.25)
            + 1.0 * bnd(transcript, 2.0, 0.25) + 1.0 * bnd(mbs, 15.7, 0.30)
            + 5.0 * max(0.0, 0.55 - rw_frac))
    good = (n_hardfail == 0 and 0.40 <= chip <= 0.60 and rw_frac >= 0.45)
    return dict(params=params, loss=float(loss), n_hardfail=int(n_hardfail), chip=float(chip), good=bool(good),
                passed=int(out.get('passed', 0)), mbgnp_cd=cd, ezi=ezi, palbo=palbo, transcript=transcript,
                rw_frac=float(rw_frac), mbs=mbs, gnp_mk=gnp_mk, mb_mk=mb_mk)


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
    N = _arg('--rand', 400); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 7)
    rng = np.random.default_rng(20260713)
    print(f'CHAIN joint optimizer: {N} random + {CLIMB} climb, {WORKERS} workers ({len(SPACE)} params); kme3 bounded FAST', flush=True)
    logfh = open(LOG, 'w')
    allr = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        allr = batch([ANCHOR] + [sample(rng) for _ in range(N)], pool, logfh, 'rand')
        allr.sort(key=lambda x: x['loss'])
        if allr:
            b = allr[0]; print(f"RANDOM best: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} kme3 {b['params']['kme3']:.1f}", flush=True)
        for rd in range(CLIMB):
            seeds = allr[:6] if allr else [{'params': ANCHOR}]
            cands = [perturb(s['params'], rng, 0.28 * (0.6 ** rd)) for s in seeds for _ in range(8)]
            allr = sorted(allr + batch(cands, pool, logfh, f'climb{rd}'), key=lambda x: x['loss'])
            b = allr[0]
            print(f"CLIMB {rd}: loss {b['loss']:.3f} nfail {b['n_hardfail']} chip {b['chip']:.3f} CD {b['mbgnp_cd']:.2f} EZi {b['ezi']:.2f} palbo {b['palbo']:.2f} S% {b['mbs']:.1f} kme3 {b['params']['kme3']:.1f} | {sum(1 for x in allr if x['good'])} GOOD", flush=True)
    logfh.close()
    if not allr:
        print('NO RESULTS'); sys.exit(0)
    with open(BEST, 'w') as fh:
        json.dump(dict(best=allr[0], top12=allr[:12]), fh, indent=2)
    b = allr[0]
    print('\n=== BEST (min loss) ===')
    print(f"  loss {b['loss']:.3f} | n_hardfail {b['n_hardfail']} | passed {b['passed']}/28 | good={b['good']}")
    print(f"  ChIP {b['chip']:.3f} | CyclinD1 {b['mbgnp_cd']:.2f} | EZH2i {b['ezi']:.2f} | Palbo {b['palbo']:.2f} | transcript {b['transcript']:.2f} | MB S% {b['mbs']:.1f}")
    print(f"  chain: kme1 {b['params']['kme1']:.1f} kme2 {b['params']['kme2']:.1f} kme3 {b['params']['kme3']:.2f} d_me {b['params']['d_me']:.4f} | rw_frac {b['rw_frac']:.2f}")
    print(f"  params = {json.dumps(b['params'])}")
    print(f"  -> {BEST}")
