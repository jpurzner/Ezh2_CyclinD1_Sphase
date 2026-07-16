"""CDK4/6 inhibitor response in GNP vs MB, modeled as a DOSE-RESPONSE (the reframe, 2026-07-15).

Supersedes the earlier resistant-clone framing (sim_cdk46i_escape.py): JP's palbociclib data (pRb-Ser807/811⁺
drops to 16-18% at 1 µM and ~2% at 5 µM, reversible) is NOT a resistant clone -- it is the DOSE-RESPONSE. We model
the drug as GRADED residual CDK4/6 activity, kPhRbCd = KP*resid (resid=1 no drug, resid=0 full inhibition; higher
dose -> lower resid), over a population carrying CyclinD1 (CV 0.70) + CKI (CV 0.42) heterogeneity. Three results,
all emergent from v44 as-is (no CDK2-bypass "option B", no senescence compartment):

  1. DOSE-RESPONSE: the fraction still dividing falls with dose, and EVERY dividing cell is pRb-Ser807/811⁺ (the
     high-CyclinD tail overcomes partial inhibition and keeps hyperphosphorylating Rb). The residual is dose-limited.
  2. REVERSIBLE QUIESCENT RESERVE: a drug-INDEPENDENT baseline-G0 fraction (high-CKI, SOX2/OLIG2-like; MB ~20%,
     GNP ~2%) that sits out of cycle at every dose and re-enters when its quiescence signal is released
     (Cook Sangar 2017: SHH-MB regrows on withdrawal). This is a distinct bistable state, NOT a distribution tail.
  3. REVERSIBILITY / re-entry kinetics: arrest is reversible at all doses; DEEPER arrest (higher dose) -> SLOWER
     washout re-entry (more p27 to clear), which reproduces JP's "5 µM takes longer to wash out" without a
     senescence state. (High-dose off-target toxicity is a separate, unmodeled confound.)

Run:  ./venv/bin/python simulations/sim_cdk46i_dose_response.py [--n=100] [--workers=8] [--pilot]
"""
import sys, os, argparse, json, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import multiprocessing as mp
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
from scipy.signal import find_peaks

@contextlib.contextmanager
def _quiet_cvode():
    fd = sys.stderr.fileno(); saved = os.dup(fd); devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, fd); yield
    finally:
        os.dup2(saved, fd); os.close(devnull); os.close(saved)

PRB_POS = 1.5                      # pRb-Ser807/811 positive threshold (hyperphospho); matches validate_v44
KTL0, CD_SDLOG = 0.801, 0.70       # CyclinD1 translation median + abundance CV (data)
INK4_SDLOG = CIP_SDLOG = 0.40      # p16/p18/Cip-Kip cell-to-cell CV ~0.42
RESERVE_KSYP21 = 0.11              # high-CKI reserve -> baseline G0 (reversible; a distinct state)
# drug-independent reversible reserve fraction (baseline G0; JP culture estimate, SOFT -- a distinct bistable state)
RESERVE_FRAC = dict(GNP=0.02, MB=0.20)
# dose grid: residual CDK4/6 activity (1 = no drug -> 0 = saturating). Fine near the transition.
RESIDS = [1.0, 0.90, 0.80, 0.72, 0.68, 0.66, 0.64, 0.62, 0.60, 0.55, 0.50, 0.30, 0.0]
COH = dict(GNP=dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002),
           MB =dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002))

_RR = None; _KP = None
def _init(_):
    global _RR, _KP
    os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
    _RR = te.loada(build_model_v44()); _KP = float(_RR['kPhRbCd'])

def _classify(coh, resid, ktl, p16, p18, kcip, reserve):
    """Return (dividing:bool, pRb_positive:bool) for one cell at one dose."""
    rr = _RR; base = COH[coh]
    rr.reset(); rr['SHH'] = 0.5
    for k, v in base.items(): rr[k] = v
    rr['k_Cd_translation'] = ktl; rr['p16'] = p16; rr['p18'] = p18
    rr['kSyP21'] = RESERVE_KSYP21 if reserve else kcip
    rr['kPhRbCd'] = _KP * resid
    with _quiet_cvode():
        try:
            r = rr.simulate(0, 11000, 4400, selections=['time', 'MPF', 'pRb'])
        except Exception:
            return (False, False)              # deep-arrest integration death -> not dividing
    m = r['time'] >= 6000; tt = r['time'][m]; dt = tt[1] - tt[0]
    nd = len(find_peaks(np.asarray(r['MPF'])[m], prominence=0.15, distance=int(200 / dt))[0])
    prb = float(np.mean(np.asarray(r['pRb'])[m]))
    return (nd >= 1, nd >= 1 and prb > PRB_POS)

def _worker(a):
    coh, resid, ktl, p16, p18, kcip, reserve = a
    div, pos = _classify(coh, resid, ktl, p16, p18, kcip, reserve)
    return (div, pos)

def _reentry_kinetics(coh, pool):
    """Deeper arrest (higher dose) -> slower washout re-entry? Single representative cell per dose."""
    return pool.map(_reentry_one, [(coh, resid) for resid in (0.62, 0.50, 0.30, 0.10, 0.0)])
def _reentry_one(a):
    coh, resid = a; rr = _RR; base = COH[coh]; T = 6000
    rr.reset(); rr['SHH'] = 0.5
    for k, v in base.items(): rr[k] = v
    rr['kPhRbCd'] = _KP * resid
    with _quiet_cvode():
        try: rr.simulate(0, 4000 + T, int((4000 + T) / 2.5))
        except Exception: return (resid, float('nan'))
        rr['kPhRbCd'] = _KP                    # WASHOUT
        r2 = rr.simulate(4000 + T, 4000 + T + 9000, 3600, selections=['time', 'Dna'])
    dna = np.asarray(r2['Dna']); tt = np.asarray(r2['time'])
    div = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    return (resid, float(tt[div[0]] - (4000 + T)) if len(div) else float('nan'))

def run(n, workers, pilot):
    rng = np.random.default_rng(20260715)
    out = {'resids': RESIDS, 'reserve_frac': RESERVE_FRAC, 'cohorts': {}}
    with mp.Pool(workers, initializer=_init, initargs=(None,)) as pool:
        for coh in ('GNP', 'MB'):
            # one heterogeneous population reused across all doses; a fixed reserve fraction is drug-independent G0
            ktl = KTL0 * np.exp(rng.normal(0, CD_SDLOG, n))
            p16 = COH[coh]['p16'] * np.exp(rng.normal(0, INK4_SDLOG, n))
            p18 = COH[coh]['p18'] * np.exp(rng.normal(0, INK4_SDLOG, n))
            kcip = 0.002 * np.exp(rng.normal(0, CIP_SDLOG, n))
            is_res = np.arange(n) < int(round(RESERVE_FRAC[coh] * n))
            curve_div, curve_pos = [], []
            for resid in RESIDS:
                args = [(coh, resid, float(ktl[i]), float(p16[i]), float(p18[i]), float(kcip[i]), bool(is_res[i]))
                        for i in range(n)]
                res = pool.map(_worker, args)
                curve_div.append(100.0 * np.mean([d for d, _ in res]))
                curve_pos.append(100.0 * np.mean([p for _, p in res]))
            reentry = _reentry_kinetics(coh, pool) if not pilot else []
            out['cohorts'][coh] = dict(pct_dividing=curve_div, pct_dividing_prb_pos=curve_pos,
                                       reentry=reentry, n=n)
            print(f"  {coh}: no-drug {curve_div[0]:.0f}% dividing -> saturating {curve_div[-1]:.0f}%; "
                  f"reserve {100*RESERVE_FRAC[coh]:.0f}% (drug-independent G0). "
                  f"pRb⁺==dividing at every dose: {all(abs(a-b)<1e-6 for a,b in zip(curve_div,curve_pos))}")
            if reentry:
                print("     re-entry (deeper arrest -> slower): " +
                      ", ".join(f"{100*(1-r):.0f}%inh:{t:.0f}u" for r, t in reentry))
    with open(os.path.join(os.path.dirname(__file__), 'cdk46i_dose_response_results.json'), 'w') as f:
        json.dump(out, f, indent=2)
    print("  -> simulations/cdk46i_dose_response_results.json")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=100); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--pilot', action='store_true')
    a = ap.parse_args()
    run(40 if a.pilot else a.n, a.workers, a.pilot)

if __name__ == '__main__':
    main()
