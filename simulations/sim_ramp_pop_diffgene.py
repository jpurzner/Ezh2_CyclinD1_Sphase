"""Enhanced redo (JP 2026-07-21, overnight): the GNP mitogen ramp AND the distribution-population, now WITH the
poised differentiation gene (with_diff_gene=True) so the premature-activation tripwire is tracked through both.
  A. GNP SHH ramp up/down: CyclinD1 governed + division rate + CyclinD1-mark + diff-gene mark(me3_diff) + Diff.
     Physiological ramp -> the diff gene stays POISED (robust). B. an OVERDRIVE ramp (growth-rate mu ramped up/down)
     pushes division fast -> me3_diff collapses -> Diff fires (premature activation = the speed limit) in a ramp.
  C/D. POPULATION (cells sampling the calibrated distributions) across a growth-drive sweep: the FRACTION that
     prematurely activates the diff gene, and the per-cell me3_diff vs division rate (the fastest cells trip).
Run: ./venv/bin/python simulations/sim_ramp_pop_diffgene.py   (writes _cache.pkl; fig = fig_ramp_pop_diffgene.py)"""
import os, sys, pickle, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, multiprocessing as mp, tellurium as te
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
KTL0, KTLEZ0, MU0 = 0.801, 0.0150988626303029, 0.0005
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
SEL = ['time', 'SHH', 'Cd', 'Mk', 'Mk_diff', 'Diff', 'Dna']


def _divrate(t, dna, win=100.0):
    dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]; dt = t[dv]
    return np.array([np.sum(np.abs(dt - ti) < win * UPH / 2) / (win / UPH) for ti in t])


def _integrate(rr, drive_fn, t0, t1, chunk=40.0):
    T, out = [], {k: [] for k in SEL if k != 'time'}
    t = t0
    while t < t1 - 1e-6:
        for k, v in drive_fn(t).items(): rr[k] = v
        r = None
        for atol in (1e-8, 1e-7, 1e-6, 1e-5):
            try:
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + chunk, max(8, int(chunk / 4)), selections=SEL); break
            except Exception:
                rr.integrator.setValue('maximum_time_step', 5.0); continue
        if r is None: t += chunk; continue
        T.append(np.asarray(r['time'])[1:])
        for k in out: out[k].append(np.asarray(r[k])[1:])
        t += chunk
    T = np.concatenate(T)
    return T, {k: np.concatenate(v) for k, v in out.items()}


def _new(M):
    rr = te.loada(M); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('maximum_num_steps', 1000000)
    return rr


def ramp(M, key, lo, hi, dur_h=140.0, others=None):
    rr = _new(M); D = dur_h * UPH
    rr[key] = lo
    for k, v in (others or {}).items(): rr[k] = v
    for atol in (1e-8, 1e-7, 1e-6):
        try: rr.integrator.setValue('absolute_tolerance', atol); rr.simulate(0, 6000, 100, selections=['time']); break
        except Exception: continue
    up = lambda t: {key: lo + (hi - lo) * (t / D)}
    dn = lambda t: {key: hi - (hi - lo) * ((t - D) / D)}
    tU, oU = _integrate(rr, up, 0, D); tD, oD = _integrate(rr, dn, D, 2 * D)
    t = np.concatenate([tU, tD]); o = {k: np.concatenate([oU[k], oD[k]]) for k in oU}
    th = (t - t[0]) / UPH
    drive = np.where(th < dur_h, lo + (hi - lo) * th / dur_h, hi - (hi - lo) * (th - dur_h) / dur_h)
    drive = np.clip(drive, min(lo, hi), max(lo, hi))
    return dict(t=th, drive=drive, Cd=o['Cd'], Mk=o['Mk'], Mk_diff=o['Mk_diff'],
                Diff=o['Diff'], divrate=_divrate(t, o['Dna']), split=len(tU))


# ---- population: per-cell premature activation across a growth-drive (mu) sweep ----
_RR = None
def _init(M):
    global _RR
    _RR = te.loada(M); _RR.integrator.setValue('relative_tolerance', 1e-6); _RR.integrator.setValue('maximum_num_steps', 2000000)

def _cell(a):
    ci, cellp, mus = a
    out = []
    for mu in mus:
        _RR.reset()
        for k, v in GNP.items(): _RR[k] = v
        for k, v in cellp.items(): _RR[k] = v
        _RR['SHH'] = 0.6; _RR['mu'] = mu
        r = None
        for atol in (1e-8, 1e-7, 1e-6, 1e-5):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(0, 20000, 3000, selections=['time', 'Dna', 'Mk_diff', 'Diff']); break
            except Exception:
                try: _RR.integrator.setValue('maximum_time_step', 5.0)
                except Exception: pass
        if r is None: out.append((np.nan, np.nan, np.nan)); continue
        t = np.asarray(r['time']); m = t > 11000; dna = np.asarray(r['Dna'])[m]
        dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]; span = (t[m][-1] - t[m][0]) / UPH
        out.append((len(dv) / span * 24 if span > 0 else 0.0, float(np.asarray(r['Mk_diff'])[m].mean()), float(np.asarray(r['Diff'])[m].mean())))
    return ci, np.array(out)


if __name__ == '__main__':
    M = build_model_v44(with_diff_gene=True)
    res = {}
    print("A/B. ramps with the differentiation tripwire...", flush=True)
    res['ramp_shh'] = ramp(M, 'SHH', 0.2, 1.0, dur_h=140.0)           # physiological mitogen ramp (start above commitment; no re-entry burst)
    res['ramp_mu'] = ramp(M, 'mu', 0.0004, 0.0016, dur_h=140.0, others={'SHH': 0.6})  # OVERDRIVE (growth-rate) ramp
    print(f"  SHH ramp: max Diff={res['ramp_shh']['Diff'].max():.2f} (poised if ~0); "
          f"mu-overdrive ramp: max Diff={res['ramp_mu']['Diff'].max():.2f} (fires if ~1)", flush=True)

    print("C/D. population premature-activation across a growth-drive sweep...", flush=True)
    N = 80; rng = np.random.default_rng(2027)
    cd_scale = np.clip(np.exp(rng.normal(0, 0.70, N)), 0.2, 5.0)
    ez_scale = np.clip(np.exp(rng.normal(0, 0.51, N)), 0.3, 3.0)
    p21 = np.clip(np.exp(rng.normal(np.log(0.45), 0.55, N)), 0.1, 2.5)
    mus = np.round(np.linspace(0.0005, 0.0016, 9), 5)
    cellps = [dict(k_Cd_translation=KTL0 * cd_scale[i], kTlEZ=KTLEZ0 * ez_scale[i], P21_div=float(p21[i])) for i in range(N)]
    t0 = time.time()
    with mp.Pool(6, initializer=_init, initargs=(M,)) as pool:
        arrs = {}
        for ci, arr in pool.imap_unordered(_cell, [(i, cellps[i], mus) for i in range(N)], chunksize=2):
            arrs[ci] = arr
    A = np.stack([arrs[i] for i in range(N)])       # [N, len(mus), 3 (divday, me3_diff, Diff)]
    res['pop'] = dict(mu=mus, cd_scale=cd_scale, A=A)
    frac = (A[:, :, 2] > 0.5).mean(axis=0)
    print(f"  premature-activation fraction vs mu: " + " ".join(f"{m*1e3:.2f}:{f*100:.0f}%" for m, f in zip(mus, frac)), flush=True)
    print(f"  ({(time.time()-t0)/60:.1f} min)", flush=True)
    pickle.dump(res, open(os.path.join(HERE, 'sim_ramp_pop_diffgene_cache.pkl'), 'wb'))
    print("-> cached; figure: fig_ramp_pop_diffgene.py", flush=True)
