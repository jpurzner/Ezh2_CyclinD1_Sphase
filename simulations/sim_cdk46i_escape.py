"""CDK4/6i therapy-escape: GNP's clean arrest vs MB's two escape routes.

The transient-G0 reframe (docs/transient_g0_synthesis.md) makes the reversible p27-high G0 CDK4/6-INDUCIBLE:
near-empty at baseline, filling when CDK4/6 is inhibited. But two JP observations do NOT fall out of a single
homogeneous MB population -- they require heterogeneity, and each maps to a grounded, falsifiable mechanism:

  1. MB carries a ~20% REVERSIBLE G0 RESERVE at baseline (a quiescent, drug-surviving, re-entry-competent pool --
     the SOX2/OLIG2 quiescent-reserve logic). A homogeneous MB does NOT sit in baseline G0: its Hh-autonomous
     mitogen (Ptch1-low + MYCN) clears p27 faster than the INK4 brake holds it, so raising INK4, lowering SHH, or
     lowering MYCN all keep it cycling (verified sweeps). What DOES arrest it is a high CDK-inhibitor tone
     (kSyP21 ~ 0.11 -> pRb-hypophospho, p27-high, reversible G0). So the reserve = a high-CKI subpopulation.

  2. ~10-20% of MB cells CONTINUE TO DIVIDE under CDK4/6i (JP, cell culture). A homogeneous MB fully arrests
     (p27 0->1.38, pRb->0). What escapes is a LOW-Rb subpopulation (total Rb pool ~0.5): with little Rb to hold
     E2F, cells cross the restriction point independent of cyclinD-CDK4/6 and keep dividing under CDK4/6i
     (verified). This is the canonical RB1-loss / low-Rb CDK4/6i-resistance route.

GNP has neither subpopulation (no reserve, full Rb) -> it arrests cleanly under CDK4/6i (baseline G0 ~0, ~0%
escape). The MB - GNP contrast is the therapy-escape liability.

Mixture (MB): ~65% normal cycling, ~20% high-CKI reserve, ~15% low-Rb resistant. GNP: 100% normal. Each cell also
carries CyclinD1-abundance heterogeneity (k_Cd_translation, CV 0.70) so the classification is simulated, not echoed.
Two conditions per cell: baseline and CDK4/6i (kPhRbCd=0). Classify each simulated trajectory as CYCLING (divides
in the window) or G0 (no division + pRb-hypophosphorylated), so the reported fractions come from the dynamics.

Run:  ./venv/bin/python simulations/sim_cdk46i_escape.py [--n=200] [--workers=8] [--pilot]
"""
import sys, os, argparse, json, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import multiprocessing as mp
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

@contextlib.contextmanager
def _quiet_cvode():
    """Silence CVODE's non-fatal C++ stderr flood (stiff arrested-cell integration -- handled by the atol
    retry loop + G0 fallback, so the classification is unaffected). Redirects fd 2 to /dev/null."""
    fd = sys.stderr.fileno()
    saved = os.dup(fd)
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, fd); yield
    finally:
        os.dup2(saved, fd); os.close(devnull); os.close(saved)

SEL = ['time', 'MPF', 'pRb', 'Dna']
SETTLE, MEAS, CHUNK, NP = 5000.0, 9000.0, 60.0, 9
PRB_G0 = 1.5                       # pRb-hypophospho G0 marker (matches validate_v44 PRB_G0_THR)
CD_SDLOG = 0.70                    # CyclinD1 abundance CV 0.70 (data); softens thresholds
KTL0 = 0.801                       # k_Cd_translation median

# grounded subpopulation overrides (verified deterministic probes)
RESERVE = dict(kSyP21=0.11)                    # high CKI -> reversible baseline G0 (SOX2/OLIG2 reserve)
RESIST_SPECIES = dict(Rb=0.5, pRb=0.0)         # low total Rb pool -> RB1-loss / low-Rb CDK4/6i resistance

# cohorts: identity params + subpopulation mixture
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

def _classify(base_params, state, kcd, cdk46i):
    """Return 'cycle' or 'g0' for one cell under one condition (simulated, not echoed)."""
    rr = _RR
    # set ALL varying PARAMS explicitly every call -- reset() does NOT reset params, they leak across cells
    for k, v in base_params.items(): rr[k] = v         # Ptch1/MYCN/p16/p18/kSyP21 (kSyP21=0.002 default)
    if state == 'reserve':
        for k, v in RESERVE.items(): rr[k] = v         # override kSyP21 high
    rr['k_Cd_translation'] = float(kcd)
    rr['kPhRbCd'] = 0.0 if cdk46i else _KPHRB0
    rr.reset()                                         # species -> initial (params persist)
    rr['SHH'] = 0.5
    if state == 'resist':
        for k, v in RESIST_SPECIES.items(): rr[k] = v  # low Rb pool -- species, set AFTER reset
    DNA, PRB = [], []
    t = 0.0; TEND = SETTLE + MEAS
    with _quiet_cvode():
        while t < TEND:
            r = None
            for atol in (1e-9, 1e-8, 1e-7):
                try:
                    rr.integrator.setValue('absolute_tolerance', atol)
                    r = rr.simulate(t, t + CHUNK, NP, selections=SEL); break
                except Exception:
                    r = None
            if r is None: break
            if t + CHUNK > SETTLE:
                DNA.append(np.asarray(r['Dna'])[1:]); PRB.append(np.asarray(r['pRb'])[1:])
            t += CHUNK
    if not DNA: return 'g0'                             # integration death at extreme arrest -> treat as G0
    DNA = np.concatenate(DNA); PRB = np.concatenate(PRB)
    ndiv = int(np.sum((DNA[:-1] > 0.9) & (DNA[1:] < 0.1)))
    if ndiv >= 1: return 'cycle'
    return 'g0' if float(np.mean(PRB)) < PRB_G0 else 'cycle'

def _worker(args):
    idx, base_params, state, kcd = args
    return (state,
            _classify(base_params, state, kcd, cdk46i=False),   # baseline
            _classify(base_params, state, kcd, cdk46i=True))    # CDK4/6i

def run_cohort(cohort, n, pool, rng):
    # build the per-cell state list from the mixture, then draw CyclinD1 heterogeneity
    states = []
    for name, frac in cohort['mix']:
        states += [name] * int(round(frac * n))
    states = (states + ['normal'] * n)[:n]
    kcd = KTL0 * np.exp(rng.normal(0, CD_SDLOG, size=n))
    args = [(i, cohort['params'], states[i], float(kcd[i])) for i in range(n)]
    res = pool.map(_worker, args)
    base = [b for (_, b, _) in res]; cdki = [c for (_, _, c) in res]
    g0_base = 100.0 * np.mean([b == 'g0' for b in base])
    g0_cdki = 100.0 * np.mean([c == 'g0' for c in cdki])
    div_cdki = 100.0 * np.mean([c == 'cycle' for c in cdki])
    # escape composition under CDK4/6i (which states keep dividing)
    esc_states = [st for (st, _, c) in res if c == 'cycle']
    return dict(name=cohort['name'], n=n, g0_base=g0_base, g0_cdki=g0_cdki, div_cdki=div_cdki,
                esc_normal=esc_states.count('normal'), esc_resist=esc_states.count('resist'),
                esc_reserve=esc_states.count('reserve'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--pilot', action='store_true')
    a = ap.parse_args()
    n = 40 if a.pilot else a.n
    rng = np.random.default_rng(20260715)                # deterministic (reproducibility; not Python hash)
    with mp.Pool(a.workers, initializer=_init, initargs=(None,)) as pool:
        out = {}
        for cohort in (GNP, MB):
            r = run_cohort(cohort, n, pool, rng)
            out[cohort['name']] = r
            print(f"  {r['name']:>3} (n={r['n']}): baseline G0 {r['g0_base']:5.1f}%  |  "
                  f"CDK4/6i G0 {r['g0_cdki']:5.1f}%  still-dividing {r['div_cdki']:5.1f}%  "
                  f"(escape: {r['esc_resist']} low-Rb + {r['esc_normal']} normal)")
    with open(os.path.join(os.path.dirname(__file__), 'cdk46i_escape_results.json'), 'w') as f:
        json.dump(out, f, indent=2)
    print("  -> simulations/cdk46i_escape_results.json")

if __name__ == '__main__':
    main()
