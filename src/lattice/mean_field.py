"""Mean-field reduction of the H3K27 nucleosome lattice.

Two jobs:
  1. Give the CORRECT steady state of the 5-state single-site chain, so that
     kappa-style reasoning uses the real flux into me3 (= P*k_me3*[me2]) rather
     than pretending me2 is saturated.
  2. Serve as the candidate ANTIMONY SURROGATE. Everything here is a small set
     of ODEs / algebraic fixed points with no arrays, so it is expressible as
     Antimony reactions; the lattice is not. Comparing the two is the test of
     whether an Antimony approximation is good enough.

Single-site continuous-time Markov chain over {u, ac, me1, me2, me3}:

    u  --P*k_me1-->  me1 --P*k_me2--> me2 --P*k_me3--> me3
    u  <--d--------  me1 <--erase---- me2 <--erase---- me3
    u  --k_ac*act--> ac   --g_ac-->   u
    (exchange sends every state -> u at rate `exch`)

with P = EZH2*(1-EZH2i)*(a0 + nucl + a_rw*m^n_rw)*(1 - w_inh*active), i.e. the
read-write feedback closes through m, the local me3 density. Self-consistency
m = pi_me3(P(m)) is where bistability lives.
"""
import numpy as np
from promoter_lattice import Params, Inputs, U, AC, ME1, ME2, ME3, NSTATES


def site_generator(P, par: Params, erase, exch, ac_rate):
    """Q matrix (NSTATES x NSTATES) for one nucleosome at PRC2 activity P."""
    Q = np.zeros((NSTATES, NSTATES))
    Q[U, ME1] = P * par.k_me1
    Q[U, AC] = ac_rate
    Q[AC, U] = par.g_ac + exch
    Q[ME1, ME2] = P * par.k_me2
    Q[ME1, U] = par.del_mk + exch
    Q[ME2, ME3] = P * par.k_me3
    Q[ME2, ME1] = erase
    Q[ME2, U] = exch
    Q[ME3, ME2] = erase
    Q[ME3, U] = exch
    np.fill_diagonal(Q, -Q.sum(1))
    return Q


def stationary(P, par: Params, erase, exch, ac_rate):
    """Stationary distribution over the five states."""
    Q = site_generator(P, par, erase, exch, ac_rate)
    A = np.vstack([Q.T[:-1], np.ones(NSTATES)])
    b = np.zeros(NSTATES); b[-1] = 1.0
    return np.linalg.lstsq(A, b, rcond=None)[0]


def prc2(m, par: Params, inp: Inputs, nucl=0.0, active=0.0):
    boost = par.a0_prc2 + par.nucleation_gain * nucl + par.a_rw_prc2 * m ** par.n_rw
    return inp.EZH2 * (1 - inp.EZH2i) * boost * max(0.0, 1 - par.w_inh_active * active)


def me3_response(m, par: Params, inp: Inputs, nucl=0.0, active=0.0,
                 psi_kdm=0.0, body=1.0, psi_p300=0.0):
    """pi_me3 given that the local me3 density is m — the RHS of the fixed point."""
    erase = par.del_mk + par.k_kdm6b * inp.Gli1 * psi_kdm
    exch = par.k_h33 * inp.tx * body
    ac_rate = par.k_ac * inp.activity * psi_p300
    return stationary(prc2(m, par, inp, nucl, active), par, erase, exch, ac_rate)[ME3]


def fixed_points(par: Params, inp: Inputs, grid=2001, **kw):
    """All self-consistent m = pi_me3(P(m)), with stability.

    Returns list of (m, stable). Two stable roots = bistable region.
    """
    ms = np.linspace(0, 1, grid)
    g = np.array([me3_response(m, par, inp, **kw) for m in ms]) - ms
    roots = []
    for i in range(len(ms) - 1):
        if g[i] == 0 or g[i] * g[i + 1] < 0:
            lo, hi = ms[i], ms[i + 1]
            for _ in range(60):                       # bisection
                mid = 0.5 * (lo + hi)
                gm = me3_response(mid, par, inp, **kw) - mid
                if (me3_response(lo, par, inp, **kw) - lo) * gm <= 0:
                    hi = mid
                else:
                    lo = mid
            r = 0.5 * (lo + hi)
            eps = 1e-4
            slope = ((me3_response(r + eps, par, inp, **kw) - (r + eps))
                     - (me3_response(r - eps, par, inp, **kw) - (r - eps))) / (2 * eps)
            roots.append((float(r), bool(slope < 0)))
    return roots


def kappa_param(par: Params, inp: Inputs, psi_kdm=0.0, body=1.0, nucl=0.0, active=0.0):
    """Hodges & Crabtree's kappa as a PARAMETER ratio (not a realised flux).

    Propagation rate = 1 / (transit time u -> me3 through the chain) at a fully
    marked neighbourhood; turnover = erase + exchange. At a fixed point realised
    gain equals realised loss by definition, so a flux ratio is tautological --
    this is the non-tautological version, and it is the quantity H&C report as
    optimal in [1.0, 1.5].
    """
    P = prc2(1.0, par, inp, nucl, active)
    if P <= 0:
        return 0.0
    transit = 1.0 / (P * par.k_me1) + 1.0 / (P * par.k_me2) + 1.0 / (P * par.k_me3)
    loss = par.del_mk + par.k_kdm6b * inp.Gli1 * psi_kdm + par.k_h33 * inp.tx * body
    return float((1.0 / transit) / loss)


def bifurcation(par: Params, inp: Inputs, arw_grid, **kw):
    """Map stable/unstable fixed points vs read-write strength a_rw_prc2.

    Returns dict with, per a_rw: the low stable branch, the unstable branch and
    the high stable branch (nan where absent). A region with three roots is
    bistable.
    """
    from dataclasses import replace
    lo, mid, hi, kap = [], [], [], []
    for arw in arw_grid:
        p2 = replace(par, a_rw_prc2=float(arw))
        fps = fixed_points(p2, inp, grid=801, **kw)
        st = sorted([m for m, s in fps if s])
        un = sorted([m for m, s in fps if not s])
        lo.append(st[0] if st else np.nan)
        hi.append(st[-1] if len(st) > 1 else (st[0] if st else np.nan))
        mid.append(un[0] if un else np.nan)
        kap.append(kappa_param(p2, inp, **{k: v for k, v in kw.items()
                                           if k in ('psi_kdm', 'body', 'nucl', 'active')}))
    return dict(a_rw=np.asarray(arw_grid), low=np.array(lo), unstable=np.array(mid),
                high=np.array(hi), kappa=np.array(kap),
                bistable=np.array([not np.isnan(m) for m in mid]))


def solve_for_me3(target_me3: float, par: Params, inp: Inputs,
                  lo=0.05, hi=40.0, **kw):
    """Pick a_rw_prc2 so the flank's stable me3 lands on `target_me3`."""
    import copy
    def top(arw):
        p2 = copy.replace(par, a_rw_prc2=arw) if hasattr(copy, 'replace') else _with(par, arw)
        fps = [m for m, st in fixed_points(p2, inp, grid=401, **kw) if st]
        return max(fps) if fps else 0.0
    if top(hi) < target_me3:
        raise ValueError(f'target me3={target_me3} unreachable even at a_rw={hi}')
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if top(mid) < target_me3:
            lo = mid
        else:
            hi = mid
    arw = 0.5 * (lo + hi)
    got = top(arw)
    if abs(got - target_me3) > 0.05:
        raise ValueError(
            f'target me3={target_me3:.2f} lies in the UNSTABLE gap of a bistable '
            f'system (nearest attainable {got:.3f} at a_rw={arw:.2f}); pick a '
            f'target on the low or high branch')
    return arw


def _with(par: Params, arw):
    from dataclasses import replace
    return replace(par, a_rw_prc2=arw)
