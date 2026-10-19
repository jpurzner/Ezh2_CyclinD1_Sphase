"""Overnight master (JP 2026-10-18, run while JP sleeps). Runs, in sequence, with heavy statistics:

  PHASE A  two-step model (baked default, model-tuned protein t1/2):
           - fine sweep for the config giving MB snapshot-G0 ~= 34% with low permanent + GNP ~0
           - LOCK that config: high-N (300) MB + GNP reference run, saved for the paper figure
  PHASE B  measured-protein INK4 variant (p18 23.9h, p21 0.98h; p27/Skp2 keep regulated clearance):
           - validate (subprocess) -> pass count
           - 2x2 attribution at N=150: {model-tuned, measured-protein} x {burst p21+p18, burst p21-only}
             -> does stable p18 drop out as a driver, leaving p21 alone? + GNP controls
  PHASE C  BAKE confirmation: full validation with the new defaults (with_cdki_mrna ON) -> expect 30/32,
           and legacy (CDKI_MRNA=0) -> expect the prior 30/32.

Everything is logged with timestamps to overnight/progress.log; all metrics -> overnight/*.json; a human
SUMMARY.md is written at the end. Each phase is wrapped so one failure cannot abort the rest.
"""
import sys, os, json, time, subprocess, datetime, numpy as np, tellurium as te
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from src.build_model_v44_heldt import build_model_v44

OUT = os.path.join(HERE, 'overnight'); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, 'progress.log')
DT = 30.0; WARM = 7000.0; PRB_G0 = 1.5; SIG_FROZ = 0.35
VENV_PY = os.path.join(ROOT, 'venv', 'bin', 'python')
QUICK = os.environ.get('QUICK') == '1'   # smoke-test: tiny N, 1-config sweep, short REC, skip long validations
# per-phase sizes (QUICK shrinks everything so the whole pipeline runs in ~3 min)
N_SWEEP, N_LOCK, N_ATTR, N_GNP = (8, 10, 8, 8) if QUICK else (60, 300, 150, 120)
REC_SWEEP, REC_LOCK, REC_ATTR = (1500.0, 1800.0, 1500.0) if QUICK else (5000.0, 7000.0, 6000.0)


def log(msg):
    line = f"[{datetime.datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with open(LOG, 'a') as f:
        f.write(line + "\n")


def make_rr(measured_protein=False):
    rr = te.loada(build_model_v44(with_cdki_mrna=True, cdki_measured_protein=measured_protein))
    rr.integrator.setValue('maximum_num_steps', 20000000)
    return rr


def _sim(rr, t, sel):
    for pts, rt, at in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6), (120, 1e-3, 1e-5)]:
        try:
            rr.integrator.setValue('relative_tolerance', rt); rr.integrator.setValue('absolute_tolerance', at)
            return rr.simulate(t, t + DT, pts, selections=sel)
        except Exception:
            continue
    return None


def run_cell(rr, ctx, cfg, froz_i, warm_i, rng, rec_sel=('pRb',)):
    """One cell. cfg: a_on, a_off, a_on18, duty, tau, p21_lvl, gnp_p21_fac, p21_only, rec."""
    rr.reset(); rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3
    KS0 = rr['kSyP21']
    if ctx == 'GNP':
        rr['SHH'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        rr['p18'] = 0.464 * froz_i; rr['p19'] = 0.36 * froz_i
        rr['kSyp21a'] = cfg['p21_lvl'] * cfg['gnp_p21_fac'] * froz_i
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = 1.73 * froz_i; rr['p19'] = 0.58 * froz_i
        rr['kSyp21a'] = cfg['p21_lvl'] * froz_i
    rr['kSyP21'] = KS0 * froz_i
    duty, tau = cfg['duty'], cfg['tau']; k_on = duty / tau; k_off = (1 - duty) / tau
    g21 = rng.random() < duty; g18 = rng.random() < duty
    sel = list(rec_sel)

    def step(t0, t1, rec):
        nonlocal g21, g18
        tr = {k: [] for k in sel}; t = t0
        while t < t1:
            if g21 and rng.random() < k_off * DT: g21 = False
            elif (not g21) and rng.random() < k_on * DT: g21 = True
            if g18 and rng.random() < k_off * DT: g18 = False
            elif (not g18) and rng.random() < k_on * DT: g18 = True
            rr['p21a_tx'] = cfg['a_on'] if g21 else cfg['a_off']
            if cfg['p21_only']:
                rr['p18_tx'] = 1.0
            else:
                rr['p18_tx'] = cfg['a_on18'] if g18 else 0.5
            r = _sim(rr, t, sel)
            if r is None: return None
            if rec:
                for k in sel: tr[k].append(float(r[k][-1]))
            t += DT
        return tr
    if step(0, warm_i, False) is None: return None
    tr = step(0, cfg['rec'], True)
    if tr is None or len(tr['pRb']) < cfg['rec'] / DT * 0.8: return None
    out = {k: np.array(v) for k, v in tr.items()}
    return out if np.all(np.isfinite(out['pRb'])) else None


def summarize(rr, ctx, cfg, N, seed=11, rec_sel=('pRb',)):
    rng = np.random.default_rng(seed)
    froz = np.clip(rng.lognormal(0, SIG_FROZ, N), 0.5, 2.0)
    warm = rng.uniform(600.0, WARM, N)
    cells = [c for c in (run_cell(rr, ctx, cfg, froz[i], warm[i], rng, rec_sel) for i in range(N)) if c is not None]
    n = len(cells)
    if n == 0: return None
    allpr = np.concatenate([c['pRb'] for c in cells])
    occ = float((allpr < PRB_G0).mean())
    trans = float(np.mean([(c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any() for c in cells]))
    perm = float(np.mean([(c['pRb'] < PRB_G0).all() for c in cells]))
    cyc = float(np.mean([(c['pRb'] >= PRB_G0).all() for c in cells]))
    res = {'ctx': ctx, 'n': n, 'N': N, 'occ': occ, 'transient': trans, 'perm': perm, 'cyc': cyc,
           'mortality': float(1 - n / N), 'cfg': {k: cfg[k] for k in cfg}}
    if 'EZH2' in rec_sel:
        ez = np.concatenate([c['EZH2'] for c in cells]); q = allpr < PRB_G0
        res['ezh2_q'] = float(np.median(ez[q])) if q.any() else None
        res['ezh2_c'] = float(np.median(ez[~q])) if (~q).any() else None
    for key in rec_sel:
        if key == 'pRb': continue
        res[f'cv_{key}'] = float(np.median([c[key].std() / (c[key].mean() + 1e-9) for c in cells]))
    # keep a few example transient trajectories (subsampled) if rich selection
    if len(rec_sel) > 1:
        ex = [i for i, c in enumerate(cells) if (c['pRb'] < PRB_G0).any() and (c['pRb'] >= PRB_G0).any()][:3]
        sub = slice(None, None, 2)
        res['examples'] = [{k: cells[i][k][sub].round(4).tolist() for k in rec_sel} for i in ex]
        res['prb_pool'] = allpr[::9].round(3).tolist()
    return res


def base_cfg(**kw):
    c = dict(a_on=3.0, a_off=0.15, a_on18=1.6, duty=0.30, tau=700.0,
             p21_lvl=0.03, gnp_p21_fac=0.1, p21_only=False, rec=6000.0)
    c.update(kw); return c


_VAL_DONE = [0]
def validate(env_extra, label):
    if QUICK and _VAL_DONE[0] >= 1:
        log(f"  [QUICK] skipping validation [{label}]")
        return {'label': label, 'summary': '(QUICK skip)', 'fails': []}
    _VAL_DONE[0] += 1
    env = dict(os.environ); env.pop('QUICK', None); env.update(env_extra)
    log(f"  validating [{label}] env={env_extra} ...")
    try:
        r = subprocess.run([VENV_PY, os.path.join(HERE, 'validate_v44.py')], cwd=ROOT, env=env,
                           capture_output=True, text=True, timeout=3600)
        line = [l for l in r.stdout.splitlines() if 'VALIDATION SUMMARY' in l]
        fails = [l.strip() for l in r.stdout.splitlines() if l.strip().startswith('[-]')]
        summ = line[0] if line else '(no summary line)'
        log(f"  [{label}] {summ}")
        for fl in fails: log(f"      {fl}")
        return {'label': label, 'summary': summ, 'fails': fails}
    except Exception as e:
        log(f"  [{label}] validation ERROR: {e}")
        return {'label': label, 'summary': f'ERROR {e}', 'fails': []}


RESULTS = {'started': str(datetime.datetime.now())}

# ============================== PHASE A ==============================
try:
    log("=== PHASE A: two-step (model-tuned protein) — fine sweep for MB ~34% ===")
    t0 = time.time()
    rrA = make_rr(measured_protein=False)
    sweep = []
    _grid = [(2.6, 0.30)] if QUICK else [(a, d) for a in [2.0, 2.3, 2.6, 2.9] for d in [0.28, 0.33]]
    for a_on, duty in _grid:
            cfg = base_cfg(a_on=a_on, duty=duty, tau=700, rec=REC_SWEEP)
            r = summarize(rrA, 'MB', cfg, N=N_SWEEP, seed=11)
            if r:
                sweep.append(r)
                log(f"  MB a_on={a_on} duty={duty}: occ={r['occ']:.1%} trans={r['transient']:.0%} perm={r['perm']:.0%} n={r['n']}")
    RESULTS['phaseA_sweep'] = sweep
    # pick config closest to 0.34 with perm < 0.06
    elig = [r for r in sweep if r['perm'] < 0.06] or sweep
    best = min(elig, key=lambda r: abs(r['occ'] - 0.34))
    bc = best['cfg']
    log(f"  -> locked config: a_on={bc['a_on']} duty={bc['duty']} tau={bc['tau']} (sweep occ {best['occ']:.1%})")
    # LOCK: high-N reference, rich recording, MB + GNP
    lock_cfg = base_cfg(a_on=bc['a_on'], duty=bc['duty'], tau=bc['tau'], rec=REC_LOCK)
    rec_rich = ('pRb', 'EZH2', 'p18_prot', 'p21a', 'p21a_m', 'P21', 'Dna')
    log(f"  locking MB (N={N_LOCK}, REC={REC_LOCK:.0f}) ...")
    RESULTS['phaseA_lock_MB'] = summarize(rrA, 'MB', lock_cfg, N=N_LOCK, seed=21, rec_sel=rec_rich)
    log(f"    MB lock: occ={RESULTS['phaseA_lock_MB']['occ']:.1%} trans={RESULTS['phaseA_lock_MB']['transient']:.0%} perm={RESULTS['phaseA_lock_MB']['perm']:.0%} n={RESULTS['phaseA_lock_MB']['n']}")
    log(f"  locking GNP (N={N_LOCK}, REC={REC_LOCK:.0f}) ...")
    RESULTS['phaseA_lock_GNP'] = summarize(rrA, 'GNP', lock_cfg, N=N_LOCK, seed=22, rec_sel=rec_rich)
    log(f"    GNP lock: occ={RESULTS['phaseA_lock_GNP']['occ']:.1%} cyc={RESULTS['phaseA_lock_GNP']['cyc']:.0%} n={RESULTS['phaseA_lock_GNP']['n']}")
    json.dump(RESULTS, open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
    log(f"  PHASE A done in {(time.time()-t0)/60:.1f} min")
except Exception as e:
    import traceback; log(f"PHASE A FAILED: {e}\n{traceback.format_exc()}")

# ============================== PHASE B ==============================
try:
    log("=== PHASE B: measured-protein INK4 variant — validation + 2x2 attribution ===")
    t0 = time.time()
    RESULTS['phaseB_validate'] = validate({'CDKI_MRNA': '1', 'CDKI_MEASURED_PROTEIN': '1'}, 'two-step + measured-protein')
    bc = RESULTS.get('phaseA_lock_MB', {}).get('cfg') or base_cfg()
    attrib = {}
    rrT = make_rr(measured_protein=False)   # model-tuned
    rrM = make_rr(measured_protein=True)    # measured INK4
    for mlab, rrx in [('tuned', rrT), ('measured', rrM)]:
        for blab, p21only in [('p21+p18', False), ('p21only', True)]:
            cfg = base_cfg(a_on=bc['a_on'], duty=bc['duty'], tau=bc['tau'], rec=REC_ATTR, p21_only=p21only)
            r = summarize(rrx, 'MB', cfg, N=N_ATTR, seed=31, rec_sel=('pRb', 'p18_prot', 'p21a', 'p21a_m'))
            attrib[f'MB_{mlab}_{blab}'] = r
            if r: log(f"  MB [{mlab}/{blab}]: occ={r['occ']:.1%} trans={r['transient']:.0%} perm={r['perm']:.0%} "
                      f"CV(p18={r.get('cv_p18_prot',0):.2f} p21={r.get('cv_p21a',0):.2f}) n={r['n']}")
    # GNP controls (should be ~0 in both)
    for mlab, rrx in [('tuned', rrT), ('measured', rrM)]:
        cfg = base_cfg(a_on=bc['a_on'], duty=bc['duty'], tau=bc['tau'], rec=REC_ATTR)
        r = summarize(rrx, 'GNP', cfg, N=N_GNP, seed=32)
        attrib[f'GNP_{mlab}'] = r
        if r: log(f"  GNP [{mlab}]: occ={r['occ']:.1%} cyc={r['cyc']:.0%} n={r['n']}")
    RESULTS['phaseB_attribution'] = attrib
    json.dump(RESULTS, open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
    log(f"  PHASE B done in {(time.time()-t0)/60:.1f} min")
except Exception as e:
    import traceback; log(f"PHASE B FAILED: {e}\n{traceback.format_exc()}")

# ============================== PHASE C ==============================
try:
    log("=== PHASE C: BAKE confirmation ===")
    t0 = time.time()
    RESULTS['phaseC_baked_default'] = validate({}, 'baked default (with_cdki_mrna ON)')
    RESULTS['phaseC_legacy'] = validate({'CDKI_MRNA': '0'}, 'legacy single-step (CDKI_MRNA=0)')
    json.dump(RESULTS, open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
    log(f"  PHASE C done in {(time.time()-t0)/60:.1f} min")
except Exception as e:
    import traceback; log(f"PHASE C FAILED: {e}\n{traceback.format_exc()}")

# ============================== SUMMARY ==============================
RESULTS['finished'] = str(datetime.datetime.now())
json.dump(RESULTS, open(os.path.join(OUT, 'results.json'), 'w'), indent=1)
try:
    def g(path, k, d='?'):
        o = RESULTS.get(path) or {}
        return o.get(k, d) if isinstance(o, dict) else d
    mb = RESULTS.get('phaseA_lock_MB') or {}; gnp = RESULTS.get('phaseA_lock_GNP') or {}
    at = RESULTS.get('phaseB_attribution') or {}
    def occ(d): return f"{d['occ']:.0%}" if isinstance(d, dict) and d else 'n/a'
    lines = []
    lines.append("# Overnight run — CDKI two-step mRNA/protein + transient-G0\n")
    lines.append(f"Started {RESULTS['started']} · finished {RESULTS['finished']}\n")
    lines.append("## Phase A — locked ~34% reference (two-step, model-tuned protein)")
    if mb:
        lines.append(f"- **MB**: snapshot-G0 **{mb['occ']:.0%}**, transient {mb['transient']:.0%}, permanent {mb['perm']:.0%} (n={mb['n']}/{mb['N']}, mortality {mb['mortality']:.0%})")
        lines.append(f"  - config: A_ON={mb['cfg']['a_on']}, duty={mb['cfg']['duty']}, tau={mb['cfg']['tau']:.0f} min; p21 mRNA CV {mb.get('cv_p21a_m','?'):.2f} -> protein CV {mb.get('cv_p21a','?'):.2f}" if isinstance(mb.get('cv_p21a_m'), float) else f"  - config: A_ON={mb['cfg']['a_on']}, duty={mb['cfg']['duty']}, tau={mb['cfg']['tau']:.0f}")
        if mb.get('ezh2_q') and mb.get('ezh2_c'):
            lines.append(f"  - EZH2 quiescent {mb['ezh2_q']:.2f} vs cycling {mb['ezh2_c']:.2f} (fold {mb['ezh2_q']/mb['ezh2_c']:.2f})")
    if gnp:
        lines.append(f"- **GNP**: snapshot-G0 **{gnp['occ']:.0%}**, pure-cycling {gnp['cyc']:.0%} (n={gnp['n']}/{gnp['N']}, mortality {gnp['mortality']:.0%})")
    lines.append("\n## Phase B — measured-protein INK4 (p18 23.9h, p21 0.98h)")
    vb = RESULTS.get('phaseB_validate') or {}
    lines.append(f"- validation: {vb.get('summary','?')}")
    lines.append("- 2x2 attribution (MB snapshot-G0):")
    lines.append(f"  - model-tuned protein:  p21+p18 = {occ(at.get('MB_tuned_p21+p18'))}, p21-only = {occ(at.get('MB_tuned_p21only'))}")
    lines.append(f"  - measured protein:     p21+p18 = {occ(at.get('MB_measured_p21+p18'))}, p21-only = {occ(at.get('MB_measured_p21only'))}")
    lines.append(f"  - GNP controls: tuned {occ(at.get('GNP_tuned'))}, measured {occ(at.get('GNP_measured'))}")
    lines.append("  - read: if measured p21+p18 ~= measured p21-only, stable p18 has dropped out as a driver -> p21 alone carries it.")
    lines.append("\n## Phase C — bake confirmation")
    lines.append(f"- baked default (two-step ON): {g('phaseC_baked_default','summary')}")
    lines.append(f"- legacy single-step (CDKI_MRNA=0): {g('phaseC_legacy','summary')}")
    lines.append("\n## Files\n- overnight/results.json (all metrics) · overnight/progress.log (timeline)")
    open(os.path.join(OUT, 'SUMMARY.md'), 'w').write("\n".join(lines) + "\n")
    log("SUMMARY.md written.")
except Exception as e:
    import traceback; log(f"SUMMARY FAILED: {e}\n{traceback.format_exc()}")
log("=== OVERNIGHT MASTER COMPLETE ===")
