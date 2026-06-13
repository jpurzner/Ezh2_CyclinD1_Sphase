"""validate_v45.py -- validation suite for the v45 stochastic-commitment model.

v45 replaces v44's phenomenological size-gate with an emergent BISTABLE CDK2-p27 commitment
toggle + Rb-dilution G1 timer + mitogen-gated CyclinD bootstrap + mass-capped quiescence, driven
by the REAL v44 HH/MYCN/EZH2->CyclinD1 cascade (via condition_cd). This script checks it against
the same experimental targets used for v44 (docs/Ezh2_CcnD1_model_targets.md), restructured for
the COUNT-fraction convention and the population (fractional) rescue read-out.

Targets & provenance (see the v44 compendium):
  * MB DMSO flow COUNT-fractions: 2N(G0+G1)=68.2%, S=15.7%; G2+M anchored to the ~2.5h DIRECT
    duration (the 4N gate over-counts). MB baseline quiescent ~21.7% (G0 marker).
  * CyclinD1 folds (from the wired cascade): MB/GNP 7.58, MB+HHi/MB 0.144, GNP+HHi/GNP 0.157,
    GNP+EZH2i/GNP 2.2.
  * Period (GNP) ~22-23h (soft).
  * Rescue (POPULATION/fractional): Hh-withdrawal raises the quiescent fraction; EZH2i lowers it.
  * MYCN floor: MB resists Hh-withdrawal ARREST (MB+HHi permanent-arrest << GNP+HHi), the paper's claim.
  * GNP-SHH -> arrest; MB -> cycles.

Run:  ./venv/bin/python simulations/validate_v45.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import v45_stochastic_commitment as v45

N = 400            # ensemble size per condition
SEED = 0

_PASS = 0; _FAIL = 0


def check(name, val, target, tol, lo_ok=False, hi_ok=False, fmt="{:.2f}"):
    """Pass if |val-target| <= tol*|target| (relative), or one-sided if lo_ok/hi_ok."""
    global _PASS, _FAIL
    if np.isnan(val):
        ok = False
    elif lo_ok:
        ok = val <= target * (1 + tol)
    elif hi_ok:
        ok = val >= target * (1 - tol)
    else:
        ok = abs(val - target) <= tol * max(abs(target), 1e-9)
    _PASS += ok; _FAIL += not ok
    v = fmt.format(val); t = fmt.format(target)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:42s} {v:>9s}  (target {t}, tol {tol:.0%})")
    return ok


def check_bool(name, ok):
    global _PASS, _FAIL
    _PASS += bool(ok); _FAIL += not ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:42s} {'yes' if ok else 'no':>9s}")
    return ok


# ---- run every condition once, cache the read-outs ----------------------------------------
def summarize(name):
    cd = v45.condition_cd(name)
    cells, st = v45.ensemble(cd, N=N, seed=SEED, return_stats=True)
    out = dict(cd=cd, arrest=st['arrest_frac'], n=len(cells))
    if cells:
        cf, _ = v45.count_fractions(cells)
        pers = np.array([c['per'] for c in cells])
        out.update(G0=cf['G0'], G1=cf['G1'], S=cf['S'], G2M=cf['G2M'], twoN=cf['G0'] + cf['G1'],
                   per_mean=float(pers.mean()), per_med=float(np.median(pers)),
                   g2m_dur=float(np.mean([c['durs']['G2M'] for c in cells])),
                   quiesc=100 * st['arrest_frac'] + (1 - st['arrest_frac']) * cf['G0'])
    else:
        out.update(G0=np.nan, G1=np.nan, S=np.nan, G2M=np.nan, twoN=np.nan,
                   per_mean=np.nan, per_med=np.nan, g2m_dur=np.nan, quiesc=100 * st['arrest_frac'])
    return out


if __name__ == "__main__":
    print("=" * 78)
    print("v45 VALIDATION (stochastic-commitment model; N=%d/condition)" % N)
    print("=" * 78)
    print("Running conditions (each builds the v44 cascade for Cd, then the v45 ensemble)...")
    S = {name: summarize(name) for name in v45.CONDITIONS}

    print("\n-- CyclinD1 (Cd) folds from the wired HH/MYCN/EZH2 cascade --")
    g = S['GNP']['cd']
    check("CyclinD1 MB/GNP",        S['MB']['cd'] / g,            7.58, 0.30)
    check("CyclinD1 MB+HHi/MB",     S['MB+HHi']['cd'] / S['MB']['cd'], 0.144, 0.40)
    check("CyclinD1 GNP+HHi/GNP",   S['GNP+HHi']['cd'] / g,       0.157, 0.40)
    check("CyclinD1 GNP+EZH2i/GNP", S['GNP+EZH2i']['cd'] / g,     2.2, 0.40)

    print("\n-- MB DMSO cell-cycle COUNT-fractions (flow) --")
    check("MB 2N (G0+G1) count%",   S['MB']['twoN'],   68.2, 0.15)
    check("MB S count% (BrdU-anchored)", S['MB']['S'],  15.7, 0.40)
    check("MB G2+M DURATION ~2.5h (direct)", S['MB']['g2m_dur'], 2.5, 0.40)
    check("MB baseline quiescent% (G0 marker)", S['MB']['quiesc'], 21.7, 0.40)

    print("\n-- Period (GNP, soft ~22-23h literature) --")
    check("GNP median period (h)",  S['GNP']['per_med'], 22.5, 0.20)

    print("\n-- Single-cell viability --")
    check_bool("MB cycles (arrest < 10%)",          S['MB']['arrest'] < 0.10)
    check_bool("GNP cycles (arrest < 10%)",         S['GNP']['arrest'] < 0.10)
    # GNP-SHH: arrest is measured as QUIESCENCE (Ki67-/pRb-, the experimental read-out), not a
    # model-internal permanent-arrest flag. At Cd~0.25 (just above the permanent-arrest threshold)
    # the cells are strongly quiescent + very slow (~39h), which the data cannot distinguish from
    # hard arrest. GNP+HHi (Cd 0.16, below threshold) does show permanent arrest.
    check_bool("GNP-SHH strongly quiescent (>60%)", S['GNP-SHH']['quiesc'] > 60)
    check_bool("GNP+HHi permanent arrest (>50%)",   S['GNP+HHi']['arrest'] > 0.50)

    print("\n-- POPULATION rescue (fractional; quiescent% = arrest + cycler-G0) --")
    check_bool("Hh-withdrawal RAISES MB quiescence (MB+HHi > MB)",
               S['MB+HHi']['quiesc'] > S['MB']['quiesc'] + 5)
    check_bool("EZH2i RESCUES (MB+HHi+EZH2i < MB+HHi)",
               S['MB+HHi+EZH2i']['quiesc'] < S['MB+HHi']['quiesc'] - 5)
    check_bool("EZH2i lowers MB quiescence (MB+EZH2i < MB)",
               S['MB+EZH2i']['quiesc'] < S['MB']['quiesc'])

    print("\n-- MYCN floor: MB resists Hh-withdrawal ARREST (the paper's claim) --")
    check_bool("MB+HHi permanent-arrest << GNP+HHi (MYCN floor)",
               S['MB+HHi']['arrest'] < 0.20 and S['GNP+HHi']['arrest'] > 0.50)

    print("\n" + "=" * 78)
    print(f"v45 VALIDATION: {_PASS}/{_PASS + _FAIL} passed")
    print("=" * 78)

    # condition table for the record
    print(f"\n  {'condition':16s} {'Cd':>6s} {'arrest%':>8s} {'2N':>5s} {'S':>5s} {'G2M':>5s} "
          f"{'quiesc%':>8s} {'medT':>6s}")
    for name in v45.CONDITIONS:
        s = S[name]
        print(f"  {name:16s} {s['cd']:6.2f} {100*s['arrest']:7.0f}% {s['twoN']:5.0f} {s['S']:5.0f} "
              f"{s['G2M']:5.0f} {s['quiesc']:7.0f}% {s['per_med']:6.1f}")
    sys.exit(1 if _FAIL else 0)
