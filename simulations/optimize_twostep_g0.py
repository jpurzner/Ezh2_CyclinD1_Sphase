"""FOCUSED joint re-optimizer for the two-step Rb default (root-cause fix for the critical review).

Goal: genuinely land the CyclinD1 MB/GNP fold at ~5.07 (currently 7.1, passing only on the wide 0.45
band) AND keep a real MB p27-high G0 dwell (~15%) -- WITHOUT the inflated birth-p27=1.8 crutch. The
decoupler (confirmed by _diag_p27_toggle): the MITOGEN-DOSE EZH2 term (kEZbas_Cd, Kez_cd) keeps EZH2
high during G0 (MB Cd stays high in G0), so G0 no longer suppresses whole-cycle EZH2 -> no longer
inflates Cd. Also targets the two standing failures (GNP+HHi de-repression; MB HU G2 fold).

Searches ONLY the EZH2 / G0 / HU knobs; leaves the validated two-step commitment + PRC2-chain machinery
fixed. P21_div_MB is set MB-specifically via the P21_DIV_MB env var. MB G0 read from validate's
phase_dmso_duration['G0']. Same subprocess band-penalty design as optimize_twostep.

Run:  ./venv/bin/python simulations/optimize_twostep_g0.py [--rand=250] [--workers=8] [--climb=6]
"""
import sys, os, json, subprocess, tempfile
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
TMP = tempfile.mkdtemp(prefix='twostepg0_')
BEST = os.path.join(ROOT, 'simulations', 'optimize_twostep_g0_best.json')
LOG = os.path.join(ROOT, 'simulations', 'optimize_twostep_g0_log.jsonl')
EXCLUDE = {'CyclinD1 MB+HHi/MB', 'MYCN GNP+HHi/GNP', 'MB G2+M duration ~2.5h (direct)'}
CHIP_TARGET = 0.5

# model params -> VALIDATE_PARAMS ; P21_div_MB is special (-> P21_DIV_MB env, MB-specific).
# NB the chain-default CyclinD1 repression is the PRC2-occupancy Hill (f0_prc2 + (1-f0)/(1+(PRC2/K_prc2)^n));
# K_EZH2_repression is DEAD in this path. EZH2i de-repression is set by K_prc2 (lower=deeper) / a_rw_prc2 (more PRC2).
SPACE = [
    ('kEZbas',      0.0003, 0.0025, True),   # cycle-flat EZH2 baseline
    ('kEZbas_Cd',   0.0003, 0.0060, True),   # MITOGEN-DOSE cycle-flat EZH2 (fold-vs-G0 decoupler)
    ('kEZE2f',      0.0040, 0.0250, True),   # cycle-dependent EZH2 (suppressed in G0)
    ('Kez_cd',      2.5,    15.0,   True),   # EZH2 Cd-saturation -> MB/GNP EZH2 contrast
    ('k_jmjd3_gli', 0.05,   0.30,   True),   # Gli->Jmjd3 eraser (ChIP MB<GNP; GNP+HHi; fold)
    ('K_prc2',      0.0030, 0.0100, True),   # *** PRC2-Hill threshold: lower = deeper repression -> bigger EZH2i ***
    ('a_rw_prc2',   0.0020, 0.0080, True),   # *** PRC2 read-write deposition: more PRC2 -> bigger EZH2i ***
    ('a0_prc2',     0.0005, 0.0030, True),   # PRC2 nucleation
    ('f0_prc2',     0.05,   0.20,   False),  # readout leaky floor
    ('n_prc2',      2.5,    4.2,    False),  # readout Hill steepness
    ('del_mk',      0.0015, 0.0050, True),   # H3K27 mark turnover (EZH2i de-repression speed)
    ('KmHU_fire',   0.30,   3.0,    True),   # HU origin-firing block (HU S vs G2 fold)
    ('kSyDna',      0.042,  0.090,  True),   # S-phase duration
    ('kDeP21Cd',    0.30,   1.50,   True),   # CDK4/6-ACTIVITY-gated p27 clearance (2026-07-15 reframe): sets G0/cycling balance + CDK4/6i fill
]
# NB: the birth-p27 crutch (P21_div_MB) is RETIRED (transient-G0 reframe) -- fixed at the GNP baseline 0.6,
# NOT searched. MB baseline G0 now emerges from the CDK4/6-INK4 balance (population tail); CDK4/6i fills it.
# shipped chain/PRC2 defaults (extracted from the built model / _ts_bake)
DEFAULTS = {'kEZbas': 0.0006988044688820413, 'kEZbas_Cd': 0.00032921999678280794,
            'kEZE2f': 0.010236776694678679, 'Kez_cd': 9.353863909522834, 'k_jmjd3_gli': 0.15038550273688858,
            'K_prc2': 0.007468156574911385, 'a_rw_prc2': 0.0038491897545513795, 'a0_prc2': 0.0017696627119243948,
            'f0_prc2': 0.1335879266811463, 'n_prc2': 3.8097549179278962, 'del_mk': 0.0022862766051927256,
            'KmHU_fire': 0.3, 'kSyDna': 0.04286412792304176, 'kDeP21Cd': 0.705}
ANCHOR = {**DEFAULTS}
# SEED3 = current 28/28-ish basin post transient-G0 reframe (activity-gated clearance, crutch retired, pRb marker)
SEED3 = {**DEFAULTS, 'kSyDna': 0.062}
# SEED2 = deeper S-phase / stronger clearance variant
SEED2 = {**DEFAULTS, 'kSyDna': 0.068, 'kDeP21Cd': 0.9}
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
    env = dict(os.environ, TWO_STEP_RB='1', H3K27_CHAIN='1',   # P21_DIV_MB stays at its default (0.6, crutch retired)
               VALIDATE_PARAMS=json.dumps(dict(params)), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=600)
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
    ezfold = float(checks.get('EZH2 MB/GNP', {}).get('actual', 0))
    ezi = float(checks.get('EZH2i CycD1 fold (GNP)', {}).get('actual', 0))
    hhi = float(checks.get('CyclinD1 GNP+HHi/GNP', {}).get('actual', 1.0))
    hu_s = float(checks.get('MB HU S fold', {}).get('actual', 1.36))
    hu_g2 = float(checks.get('MB HU G2 fold', {}).get('actual', 0.23))
    mbs = float(checks.get('MB S count% (flow; BrdU Ts~3h)', {}).get('actual', 0))
    g2g0 = float(checks.get('EZH2 protein G2/G0 (1.48)', {}).get('actual', 0))
    mb_g0 = float(out.get('phase_dmso_duration', {}).get('G0', 0.0))   # % duration-fraction
    n_hardfail = sum(1 for n, c in checks.items() if n not in EXCLUDE and not c['pass'])
    bnd = lambda x, t, b: max(0.0, abs(x - t) / t - b)
    kd21 = float(params.get('kDeP21Cd', 0.705))
    loss = (10.0 * n_hardfail
            + 2.5 * bnd(cd, 5.07, 0.12)          # *** fold: land it TIGHT ***
            + 2.0 * bnd(ezfold, 2.05, 0.18)      # EZH2 MB/GNP toward 2.05 (same knob as fold)
            + 2.0 * bnd(hhi, 0.157, 0.30)        # GNP+HHi de-repression
            + 2.0 * bnd(hu_g2, 0.23, 0.40)       # MB HU G2 fold
            + 3.0 * bnd(chip, CHIP_TARGET, 0.15)
            + 2.0 * bnd(ezi, 1.7, 0.25)          # EZH2i de-repression: catalytic-inhibitor ceiling ~1.7 (NOT 2.2)
            + 2.0 * bnd(mbs, 15.7, 0.20)         # *** PROTECT MB S% (don't inflate) ***
            + 1.0 * bnd(g2g0, 1.48, 0.30)        # EZH2 cycle-phase gradient (structurally hard; low weight)
            + 1.0 * bnd(hu_s, 1.36, 0.30)
            + 1.0 * max(0.0, mb_g0 - 8.0))       # deterministic MB G0 should be ~0 (crutch retired; MB 20% is population-level)
    good = (n_hardfail == 0 and abs(cd - 5.07) / 5.07 <= 0.15
            and abs(mbs - 15.7) / 15.7 <= 0.25)
    return dict(params=params, loss=float(loss), n_hardfail=int(n_hardfail), chip=float(chip),
                good=bool(good), passed=int(out.get('passed', 0)), mbgnp_cd=cd, ezfold=ezfold,
                ezi=ezi, hhi=hhi, hu_g2=hu_g2, mb_g0=mb_g0, mbs=mbs, g2g0=g2g0, p21mb=kd21)


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
        logfh.write(json.dumps({'tag': tag, **{k: r[k] for k in ('loss', 'n_hardfail', 'passed', 'mbgnp_cd', 'ezfold', 'hhi', 'hu_g2', 'mb_g0', 'chip', 'p21mb', 'good', 'params')}}) + '\n'); logfh.flush()
        res.append(r)
        if done % 25 == 0:
            print(f"  [{tag}] {done}/{len(cands)}, {sum(1 for x in res if x['good'])} GOOD, best {min((x['loss'] for x in res), default=99):.3f}", flush=True)
    return res


if __name__ == '__main__':
    N = _arg('--rand', 250); WORKERS = _arg('--workers', 8); CLIMB = _arg('--climb', 6)
    rng = np.random.default_rng(20260714)
    print(f'FOCUSED two-step G0 optimizer: {N} random + {CLIMB} climb, {WORKERS} workers ({len(SPACE)} knobs)', flush=True)
    logfh = open(LOG, 'w')
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        allr = batch([ANCHOR, SEED2, SEED3] + [sample(rng) for _ in range(N)], pool, logfh, 'rand')
        allr.sort(key=lambda x: x['loss'])
        if allr:
            b = allr[0]; print(f"RANDOM best: loss {b['loss']:.3f} nfail {b['n_hardfail']} passed {b['passed']}/28 CD {b['mbgnp_cd']:.2f} G0 {b['mb_g0']:.1f} ezf {b['ezfold']:.2f}", flush=True)
        for rd in range(CLIMB):
            seeds = allr[:6] if allr else [{'params': ANCHOR}]
            cands = [perturb(s['params'], rng, 0.25 * (0.6 ** rd)) for s in seeds for _ in range(8)]
            allr = sorted(allr + batch(cands, pool, logfh, f'climb{rd}'), key=lambda x: x['loss'])
            b = allr[0]
            print(f"CLIMB {rd}: loss {b['loss']:.3f} nfail {b['n_hardfail']} passed {b['passed']}/28 CD {b['mbgnp_cd']:.2f} ezf {b['ezfold']:.2f} G0 {b['mb_g0']:.1f} HHi {b['hhi']:.3f} HUg2 {b['hu_g2']:.3f} p27 {b['p21mb']:.2f} | {sum(1 for x in allr if x['good'])} GOOD", flush=True)
    logfh.close()
    if not allr:
        print('NO RESULTS'); sys.exit(0)
    with open(BEST, 'w') as fh:
        json.dump(dict(best=allr[0], top12=allr[:12]), fh, indent=2)
    b = allr[0]
    print('\n=== BEST (min loss) ===')
    print(f"  loss {b['loss']:.3f} | n_hardfail {b['n_hardfail']} | passed {b['passed']}/28 | good={b['good']}")
    print(f"  CyclinD1 fold {b['mbgnp_cd']:.2f} | EZH2 MB/GNP {b['ezfold']:.2f} | MB G0 {b['mb_g0']:.1f}% | GNP+HHi {b['hhi']:.3f} | HU G2 {b['hu_g2']:.3f} | ChIP {b['chip']:.3f} | birth-p27 {b['p21mb']:.2f}")
    print(f"  params = {json.dumps(b['params'])}")
    print(f"  -> {BEST}")
