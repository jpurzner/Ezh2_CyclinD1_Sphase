"""CO-OPTIMIZE the H3K27me3 mark to MAXIMIZE withdrawal memory while keeping every non-soft validation
target on-data -- in particular PINNING CyclinD1 MB/GNP to the measured 5.07 (the naive strong-memory bake
drifted it to 6.75). Free params: k_w_mk, del_mk, K_mk, g_mk (transcription->PRC2 eviction, the key
decoupler that can keep the PROLIFERATING folds calibrated while the mark accumulates only in withdrawal),
n_mk, f0_mk.

Objective  : memory = mark level in the withdrawn state (MB+HHi, the de-weighted condition) + contrast bonus.
Hard constr: MB/GNP in [4.75,5.45] AND every validate_v44 check passes EXCEPT the 3 de-weighted
             (MB+HHi/MB soft per JP; MYCN GNP+HHi and MB G2+M duration = pre-existing ignore).
Method     : random sample -> keep feasible -> hill-climb around the top-K. Subprocess-drives validate_v44
             via VALIDATE_PARAMS/VALIDATE_RESULTS (exact 27-target scoring). Checkpoints best + jsonl.

Run:  ./venv/bin/python simulations/sim_mark_memory_cooptimize.py [--rand=260] [--workers=8] [--climb=4]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='cooptmk_')
BEST = os.path.join(ROOT, 'simulations', 'sim_mark_memory_cooptimize_best.json')
LOG = os.path.join(ROOT, 'simulations', 'sim_mark_memory_cooptimize_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
MBGNP_LO, MBGNP_HI = 4.75, 5.45

# free params: (name, lo, hi, log?)
SPACE = [
    ('k_w_mk', 0.004, 0.08, True),
    ('del_mk', 0.0004, 0.004, True),
    ('K_mk',   0.25, 0.95, False),
    ('g_mk',   0.05, 0.60, True),
    ('n_mk',   1.20, 4.00, False),
    ('f0_mk',  0.10, 0.30, False),
]
_lock_counter = [0]


def sample(rng):
    p = {}
    for name, lo, hi, lg in SPACE:
        p[name] = float(np.exp(rng.uniform(np.log(lo), np.log(hi)))) if lg else float(rng.uniform(lo, hi))
    return p


def perturb(p, rng, scale=0.25):
    q = dict(p)
    for name, lo, hi, lg in SPACE:
        if lg:
            q[name] = float(np.clip(p[name] * np.exp(rng.normal(0, scale)), lo, hi))
        else:
            q[name] = float(np.clip(p[name] + rng.normal(0, scale * (hi - lo)), lo, hi))
    return q


def evaluate(params):
    i = _lock_counter[0]; _lock_counter[0] += 1
    rp = os.path.join(TMP, f'res_{i % 4096}.json')
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(params), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=360)
        with open(rp) as fh:
            out = json.load(fh)
    except Exception:
        return None
    checks = {c['name']: c for c in out.get('checks', [])}
    if 'CyclinD1 MB/GNP' not in checks:
        return None
    marks = out.get('marks', {})
    mbgnp = checks['CyclinD1 MB/GNP']['actual']
    feasible = (MBGNP_LO <= mbgnp <= MBGNP_HI)
    n_hardfail = 0
    for name, c in checks.items():
        if name in EXCLUDE:
            continue
        if not c['pass']:
            feasible = False; n_hardfail += 1
    m_mbhhi = float(marks.get('MB + HHi', 0.0) or 0.0)
    m_mb = float(marks.get('MB', 0.0) or 0.0)
    m_gnphhi = float(marks.get('GNP + HHi', 0.0) or 0.0)
    memory = m_mbhhi + 0.3 * max(0.0, m_mbhhi - m_mb)              # withdrawal mark + rise-on-withdrawal bonus
    return dict(params=params, feasible=bool(feasible), memory=float(memory),
                passed=int(out.get('passed', 0)), n_hardfail=int(n_hardfail), mbgnp=float(mbgnp),
                mark_mbhhi=m_mbhhi, mark_mb=m_mb, mark_gnphhi=m_gnphhi,
                ezi=float(checks['EZH2i CycD1 fold (GNP)']['actual']),
                transcript=float(checks.get('EZH2 transcript S/G0 (1.8-2.5)', {}).get('actual', 0)))


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='):
            return type(d)(a.split('=', 1)[1])
    return d


def run_batch(cands, pool, logfh, tag):
    results = []
    futs = {pool.submit(evaluate, c): c for c in cands}
    done = 0
    for fut in as_completed(futs):
        r = fut.result(); done += 1
        if r is None:
            continue
        logfh.write(json.dumps({'tag': tag, **{k: r[k] for k in ('feasible', 'memory', 'passed', 'mbgnp', 'mark_mbhhi', 'ezi', 'params')}}) + '\n'); logfh.flush()
        results.append(r)
        if done % 20 == 0:
            fe = [x for x in results if x['feasible']]
            bm = max([x['memory'] for x in fe], default=0)
            print(f'  [{tag}] {done}/{len(cands)} eval, {len(fe)} feasible, best memory {bm:.3f}', flush=True)
    return results


if __name__ == '__main__':
    N_RAND = _arg('--rand', 260); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 4)
    rng = np.random.default_rng(20260709)
    print(f'CO-OPTIMIZE mark memory: {N_RAND} random + {CLIMB} climb rounds, {WORKERS} workers', flush=True)
    print(f'  constraints: MB/GNP in [{MBGNP_LO},{MBGNP_HI}] + all non-soft targets pass (soft: {sorted(EXCLUDE)})', flush=True)
    logfh = open(LOG, 'w')
    all_feasible = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        # baseline anchor (current baked) + random sample
        anchor = {'k_w_mk': 0.025, 'del_mk': 0.001, 'K_mk': 0.72, 'g_mk': 0.10554, 'n_mk': 1.6755, 'f0_mk': 0.14497}
        cands = [anchor] + [sample(rng) for _ in range(N_RAND)]
        res = run_batch(cands, pool, logfh, 'rand')
        all_feasible = [r for r in res if r['feasible']]
        all_feasible.sort(key=lambda x: -x['memory'])
        print(f'RANDOM done: {len(all_feasible)} feasible / {len(res)} evaluated', flush=True)
        if all_feasible:
            print(f'  top random memory {all_feasible[0]["memory"]:.3f}  (MB/GNP {all_feasible[0]["mbgnp"]:.2f}, mark_MBHHi {all_feasible[0]["mark_mbhhi"]:.3f})', flush=True)

        # hill-climb around the top-K feasible
        for rd in range(CLIMB):
            seeds = all_feasible[:6] if all_feasible else [{'params': anchor}]
            scale = 0.30 * (0.6 ** rd)
            cands = []
            for s in seeds:
                cands += [perturb(s['params'], rng, scale) for _ in range(8)]
            res = run_batch(cands, pool, logfh, f'climb{rd}')
            fe = [r for r in res if r['feasible']]
            all_feasible = sorted(all_feasible + fe, key=lambda x: -x['memory'])
            bm = all_feasible[0]['memory'] if all_feasible else 0
            print(f'CLIMB {rd}: +{len(fe)} feasible, best memory {bm:.3f}', flush=True)

    logfh.close()
    if not all_feasible:
        print('NO FEASIBLE POINT FOUND (constraints too tight?)'); sys.exit(0)
    best = all_feasible[0]
    with open(BEST, 'w') as fh:
        json.dump(dict(best=best, top10=all_feasible[:10]), fh, indent=2)
    print('\n=== BEST (max withdrawal memory, all non-soft targets on-data) ===')
    print(f'  memory={best["memory"]:.3f}  mark_MBHHi={best["mark_mbhhi"]:.3f}  mark_GNPHHi={best["mark_gnphhi"]:.3f}  mark_MB(cyc)={best["mark_mb"]:.3f}')
    print(f'  MB/GNP={best["mbgnp"]:.2f}  EZH2i_fold={best["ezi"]:.2f}  transcript_S/G0={best["transcript"]:.2f}  passed={best["passed"]}/27')
    print(f'  params={json.dumps(best["params"])}')
    print(f'  -> {BEST}')
