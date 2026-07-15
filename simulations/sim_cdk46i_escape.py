"""CDK4/6i therapy-escape: a MIXTURE model of GNP vs MB under CDK4/6 inhibition.

WHAT THIS IS (read before citing numbers). This is a two-subpopulation MIXTURE model, NOT an emergent-fraction
prediction. It does two separable things:

  (1) MECHANISM VERIFICATION (this is the model's contribution -- these are ODE outputs):
      * a HIGH-CKI cell (kSyP21~0.11) deterministically arrests in a reversible, pRb-hypophospho, p27-high G0;
      * a LOW total-Rb cell (Rb pool ~0.5) deterministically keeps DIVIDING under CDK4/6i (kPhRbCd=0), because
        with little Rb to hold E2F it crosses the restriction point independent of cyclinD-CDK4/6;
      * a NORMAL cell cycles at baseline and ARRESTS under CDK4/6i (p27 accumulates, pRb->hypophospho), and
        RE-ENTERS the cycle when the drug is washed out (reversible quiescence -- the Chahin/palbociclib result).
      All four are genuine ODE behaviours, each verified by deterministic sweeps (see docs section 7); MB's
      Hh-autonomy is why lowering SHH/MYCN or raising INK4 does NOT arrest it (only high CKI does).

  (2) COMPOSITION at ASSUMED fractions (these are INPUTS, not predictions):
      The subpopulation SIZES -- MB 65% normal / 20% high-CKI reserve / 15% low-Rb resistant, GNP 100% normal --
      are JP's cell-culture ESTIMATES (~20% baseline G0; ~10-20% still dividing under CDK4/6i), a SOFT target
      (docs section 3/6: "an informed guess, NOT a fit"). They are HAND-SET here. The reported baseline-G0 and
      still-dividing percentages therefore EQUAL those input fractions -- the CyclinD1-abundance heterogeneity
      (CD_SDLOG=0.70) is included for realism but does NOT move the classification (the state overrides dominate),
      so it does not make the magnitudes emergent. What the model asserts is that these mechanisms, at these
      assumed sizes, ARE SUFFICIENT to reproduce the observed GNP-clean-arrest vs MB-two-route-escape contrast --
      and that GNP lacks both subpopulations by biological assumption (normal tissue: no RB1 loss, no tumour
      quiescent reserve), not by derivation.

FALSIFIABLE (separates the two escape routes experimentally): reserve = p27-high / pRb-hypophospho / SOX2+ (or
OLIG2+), Ki67-, re-enters on release of the quiescence signal; resistant = Rb-low / RB1-loss, Ki67+ DURING CDK4/6i
(their pRb stays low even while cycling, so Ki67 -- not pRb -- marks them).

Conditions per cell: baseline and CDK4/6i (kPhRbCd=0). Classification is from the trajectory (divisions + pRb),
recorded from t=0 so early-arresting cells are scored on real data, not an integration-death fallback. REVERSIBILITY
is verified separately (--reversibility): a clean single-integration arrest-then-release for the reserve (release
its CKI signal) and the CDK4/6i-arrested normal cell (wash out the drug), each measuring re-entry (re-division).
This is kept OUT of the population loop because the chunked, tight-tolerance integration kills stiff arrested cells
before a mid-run washout is reached (an artifact, not a biological 0% re-entry).

Run:  ./venv/bin/python simulations/sim_cdk46i_escape.py [--n=200] [--workers=8] [--pilot] [--reversibility]
"""
import sys, os, argparse, json, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import multiprocessing as mp
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

@contextlib.contextmanager
def _quiet_cvode():
    """Silence CVODE's non-fatal C++ stderr flood (stiff arrested-cell integration -- the trajectory recorded
    up to the failure is genuinely arrested, and classification uses it directly). Redirects fd 2 to /dev/null."""
    fd = sys.stderr.fileno(); saved = os.dup(fd); devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, fd); yield
    finally:
        os.dup2(saved, fd); os.close(devnull); os.close(saved)

SEL = ['time', 'MPF', 'pRb', 'Dna']
SETTLE, MEAS, CHUNK, NP = 4000.0, 9000.0, 60.0, 9
PRB_G0 = 1.5                       # pRb-hypophospho G0 marker (matches validate_v44 PRB_G0_THR)
CD_SDLOG = 0.70                    # CyclinD1 abundance CV 0.70 (data); realism only -- does NOT move classification
KTL0 = 0.801                       # k_Cd_translation median
# Per-cell CKI heterogeneity (2026-07-15, JP): p16/p18/Cip-Kip are now DISTRIBUTIONS, not unitary values --
# each cell draws a multiplier ~ LogN(0, SDLOG) on its cohort CKI level. NB this adds realistic cell-to-cell
# variation in commitment/arrest PROPENSITY, but does NOT by itself create baseline G0 or CDK4/6i escape: the
# palbo block is UPSTREAM of the CKIs (kPhRbCd=0 zeros Rb mono-phosphorylation regardless of CKI level), and the
# baseline-G0 arrest threshold sits ~40x above the Cip/Kip median (a distinct high-CKI state, not a tail). The
# reserve/resist subpopulations therefore remain distinct states; the CKI draws vary everything else.
INK4_SDLOG = 0.40                  # p16/p18 cell-to-cell CV ~0.42 (INK4 protein heterogeneity)
CIP_SDLOG  = 0.40                  # Cip/Kip (p27/p21) synthesis CV ~0.42

# grounded subpopulation overrides (verified deterministic probes -- see docstring)
RESERVE = dict(kSyP21=0.11)                    # high CKI -> reversible baseline G0 (SOX2/OLIG2 reserve)
RESIST_SPECIES = dict(Rb=0.5, pRb=0.0)         # low total Rb pool -> RB1-loss / low-Rb CDK4/6i resistance

# cohorts: identity params + subpopulation mixture (SIZES ARE HAND-SET -- JP culture estimates, a SOFT target)
GNP = dict(name='GNP',
           params=dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002),
           mix=[('normal', 1.00), ('reserve', 0.00), ('resist', 0.00)])
MB  = dict(name='MB',
           params=dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002),
           mix=[('normal', 0.65), ('reserve', 0.20), ('resist', 0.15)])

_RR = None
_KPHRB0 = None
def _init(_):
    global _RR, _KPHRB0
    os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
    _RR = te.loada(build_model_v44())
    _KPHRB0 = float(_RR['kPhRbCd'])            # baked CDK4/6 rate to restore for non-CDK4/6i cells (reset() gotcha)

def _simulate(base_params, state, kcd, drug_schedule):
    """Integrate one cell; return (DNA, PRB, times) recorded from t=0. drug_schedule(t)->True means CDK4/6i active."""
    rr = _RR
    for k, v in base_params.items(): rr[k] = v         # set ALL varying params every call (reset() leaks params)
    if state == 'reserve':
        for k, v in RESERVE.items(): rr[k] = v
    rr['k_Cd_translation'] = float(kcd)
    rr['kPhRbCd'] = _KPHRB0
    rr.reset(); rr['SHH'] = 0.5
    if state == 'resist':
        for k, v in RESIST_SPECIES.items(): rr[k] = v  # low Rb pool -- species, set AFTER reset
    DNA, PRB, TIME = [], [], []
    t = 0.0; TEND = SETTLE + MEAS
    with _quiet_cvode():
        while t < TEND:
            rr['kPhRbCd'] = 0.0 if drug_schedule(t) else _KPHRB0
            r = None
            for atol in (1e-9, 1e-8, 1e-7):
                try:
                    rr.integrator.setValue('absolute_tolerance', atol)
                    r = rr.simulate(t, t + CHUNK, NP, selections=SEL); break
                except Exception:
                    r = None
            if r is None: break                        # integration death -> classify on what we have (real data)
            DNA.append(np.asarray(r['Dna'])[1:]); PRB.append(np.asarray(r['pRb'])[1:]); TIME.append(np.asarray(r['time'])[1:])
            t += CHUNK
    if not DNA: return np.array([]), np.array([]), np.array([])
    return np.concatenate(DNA), np.concatenate(TIME) * 0 + 1, np.concatenate(PRB)  # placeholder time not needed

def _ndiv(DNA, mask=None):
    d = DNA if mask is None else DNA[mask]
    if d.size < 2: return 0
    return int(np.sum((d[:-1] > 0.9) & (d[1:] < 0.1)))

def _classify_traj(DNA, PRB):
    """cycle if the recorded trajectory divides; else g0 iff pRb-hypophospho. Uses real recorded data (records
    from t=0), so an early-arresting cell that later fails integration is still scored on its genuine arrest."""
    if DNA.size == 0: return 'g0'                       # nothing integrated at all -> deep arrest
    if _ndiv(DNA) >= 1: return 'cycle'
    return 'g0' if float(np.mean(PRB)) < PRB_G0 else 'cycle'

def _worker(args):
    idx, base_params, state, kcd = args
    n_none = lambda t: False; n_all = lambda t: True
    dna_b, _, prb_b = _simulate(base_params, state, kcd, n_none)   # baseline
    dna_i, _, prb_i = _simulate(base_params, state, kcd, n_all)    # CDK4/6i
    return (state, _classify_traj(dna_b, prb_b), _classify_traj(dna_i, prb_i))

def run_cohort(cohort, n, pool, rng):
    states = []
    for name, frac in cohort['mix']:
        states += [name] * int(round(frac * n))
    states = (states + ['normal'] * n)[:n]
    kcd = KTL0 * np.exp(rng.normal(0, CD_SDLOG, size=n))
    # per-cell CKI DISTRIBUTIONS: multiply each cohort CKI by a lognormal draw (0 stays 0, e.g. GNP p16)
    p16m = np.exp(rng.normal(0, INK4_SDLOG, size=n))
    p18m = np.exp(rng.normal(0, INK4_SDLOG, size=n))
    cipm = np.exp(rng.normal(0, CIP_SDLOG, size=n))
    def cell_params(i):
        p = dict(cohort['params'])
        p['p16'] = p['p16'] * float(p16m[i]); p['p18'] = p['p18'] * float(p18m[i])
        p['kSyP21'] = p['kSyP21'] * float(cipm[i])
        return p
    args = [(i, cell_params(i), states[i], float(kcd[i])) for i in range(n)]
    res = pool.map(_worker, args)
    base = [b for (_, b, _) in res]; cdki = [c for (_, _, c) in res]
    g0_base = 100.0 * np.mean([b == 'g0' for b in base])
    g0_cdki = 100.0 * np.mean([c == 'g0' for c in cdki])
    div_cdki = 100.0 * np.mean([c == 'cycle' for c in cdki])
    esc_states = [st for (st, _, c) in res if c == 'cycle']
    return dict(name=cohort['name'], n=n, g0_base=g0_base, g0_cdki=g0_cdki, div_cdki=div_cdki,
                esc_normal=esc_states.count('normal'), esc_resist=esc_states.count('resist'),
                esc_reserve=esc_states.count('reserve'))

def verify_reversibility():
    """Clean single-integration arrest-then-release, verifying re-entry (the 'reversible' claim). NOT chunked/tight
    -- that is what let stiff arrested cells die mid-run in the population loop. Returns re-division counts."""
    _init(None)
    M = build_model_v44(); KP = float(te.loada(M)['kPhRbCd'])
    MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553)
    def arrest_then_release(setter_on, setter_off):
        rr = te.loada(M); rr['SHH'] = 0.5
        for k, v in MB.items(): rr[k] = v
        setter_on(rr); rr.reset(); rr['SHH'] = 0.5
        for k, v in MB.items(): rr[k] = v
        setter_on(rr)
        r1 = rr.simulate(0, 8000, 3200, selections=['time', 'MPF', 'Dna'])
        n_arrest = _ndiv(np.asarray(r1['Dna']))
        setter_off(rr)                                    # RELEASE
        r2 = rr.simulate(8000, 18000, 4000, selections=['time', 'MPF', 'Dna'])
        return n_arrest, _ndiv(np.asarray(r2['Dna']))
    # reserve: arrest by high CKI, release the CKI signal
    na, nr = arrest_then_release(lambda rr: rr.__setitem__('kSyP21', 0.11),
                                 lambda rr: rr.__setitem__('kSyP21', 0.002))
    print(f"  RESERVE (high-CKI arrest -> release CKI):  divisions during arrest={na}, after release={nr}  "
          f"-> {'RE-ENTERS (reversible)' if nr >= 1 else 'stays arrested'}")
    # CDK4/6i-arrested normal: arrest by drug, wash out
    na2, nr2 = arrest_then_release(lambda rr: rr.__setitem__('kPhRbCd', 0.0),
                                   lambda rr: rr.__setitem__('kPhRbCd', KP))
    print(f"  CDK4/6i-arrested normal (drug -> washout):  divisions during drug={na2}, after washout={nr2}  "
          f"-> {'RE-ENTERS (reversible)' if nr2 >= 1 else 'stays arrested'}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--pilot', action='store_true')
    ap.add_argument('--reversibility', action='store_true')
    a = ap.parse_args()
    if a.reversibility:
        verify_reversibility(); return
    n = 40 if a.pilot else a.n
    rng = np.random.default_rng(20260715)
    with mp.Pool(a.workers, initializer=_init, initargs=(None,)) as pool:
        out = {}
        for cohort in (GNP, MB):
            r = run_cohort(cohort, n, pool, rng)
            out[cohort['name']] = r
            print(f"  {r['name']:>3} (n={r['n']}): baseline G0 {r['g0_base']:5.1f}%  |  CDK4/6i still-dividing "
                  f"{r['div_cdki']:5.1f}% (all {r['esc_resist']} resistant, {r['esc_normal']} normal)")
    print("  NOTE: baseline-G0 and still-dividing % EQUAL the hand-set mixture (JP culture estimate, soft target);"
          " the model verifies the MECHANISMS + reversibility, not the fractions.")
    with open(os.path.join(os.path.dirname(__file__), 'cdk46i_escape_results.json'), 'w') as f:
        json.dump(out, f, indent=2)
    print("  -> simulations/cdk46i_escape_results.json")

if __name__ == '__main__':
    main()
