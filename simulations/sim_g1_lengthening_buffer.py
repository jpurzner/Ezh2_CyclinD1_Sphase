"""Graded G1 / cell-cycle LENGTHENING (no transient G0) as the buffering mechanism (JP 2026-07-20).

JP's mechanism + research grounding (Contestabile 2009 GNP: Tc 16.25h, G1 7.63h, S 7.11h, G2 1.18h, M 0.33h):
CyclinD1 responds to mitogen FAST; EZH2 = a mitogen-independent BASELINE + an E2F/S-phase-expressed component, so
EZH2 responds SLOW -- it falls only when the cell spends LESS time in S, i.e. when the CYCLE LENGTHENS. There are
NO transient G0s in GNP, so the slowdown is GRADED G1 lengthening, not a discrete arrest. => EZH2 (and the
H3K27me3 mark) LOW-PASS-FILTER mitogen; the filter time-constant IS the cell-cycle-length response.

The current v44 model does NOT do this: Tc is ~constant (~20-22h) across the proliferative range then abruptly
arrests below SHH~0.08 (measured). This prototype shows what graded G1 lengthening buys.

Reduced model (time in hours, grounded in the GNP numbers):
  D (CyclinD1): fast tracker of mitogen M.        tauD ~ 2 h
  S-phase frequency fS(D) = 1/Tc(D): Tc = Ts + TG1(D), Ts fixed ~7 h, TG1 LENGTHENS as D falls (graded; arrest only
     as D->0). So low mitogen -> long G1 -> low S-frequency, continuously (NO discrete G0).
  EZH2 E = BASELINE + S-phase component:  dE/dt = (Ebas + kEs*fS - E)/tauE,  tauE ~ one cycle (slow), Ebas > 0.
  H3K27me3 Mk: written by E, decayed by turnover + replication dilution (~ fS):  dMk/dt = kw*E*(1-Mk) - (del+mu*fS)*Mk
Run: ./venv/bin/python simulations/sim_g1_lengthening_buffer.py   (writes _cache.npz; fig = fig_g1_lengthening_buffer.py)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# GNP cycle grounding (Contestabile 2009 control, hours)
TS = 7.1          # S duration (fixed; GF-unbiased, trustworthy)
TG2M = 1.5        # G2+M (short in GNP)
TC_MIN = 16.25    # min cycle length at saturating mitogen
TG1_MIN = TC_MIN - TS - TG2M   # ~7.65 h, matches Contestabile G1 7.63 h


def Tc_of_D(D, TG1max=60.0, Kg=0.9, hg=3.0):
    """Cell-cycle length: G1 LENGTHENS as CyclinD1 falls (graded). TG1 from TG1_MIN up to TG1max; arrest ~ D->0."""
    frac = D**hg / (Kg**hg + D**hg)                 # 0..1 proliferative competence
    TG1 = TG1_MIN + (TG1max - TG1_MIN) * (1 - frac)
    return TS + TG2M + TG1                          # hours


def simulate(M_t, T=1200.0, dt=0.05, kD=3.0, b=0.05, tauD=2.0, Ebas=1.5, kEs=45.0, tauE=18.0,
             kw=0.5, del_=0.004, mu=8.0, lengthen=True, hard_G0=False, D0=1.5, E0=2.6, Mk0=0.7):
    n = int(T / dt); t = np.arange(n) * dt
    D = np.empty(n); E = np.empty(n); Mk = np.empty(n); fS = np.empty(n); Tc = np.empty(n)
    d, e, mk = D0, E0, Mk0
    for i in range(n):
        M = M_t(t[i])
        d = max(0.0, d + ((kD * (b + M)) - d) / tauD * dt)
        if hard_G0:                                  # current-model style: cycle at fixed Tc above threshold, else ARREST
            f = (1.0 / TC_MIN) if d > 0.75 else 0.0
        elif lengthen:                               # graded G1 lengthening (JP)
            f = 1.0 / Tc_of_D(d)
        else:                                        # no lengthening: fixed Tc whenever above a low floor
            f = (1.0 / TC_MIN) if d > 0.15 else 0.0
        e = max(0.0, e + (Ebas + kEs * f - e) / tauE * dt)
        mk = min(1.0, max(0.0, mk + (kw * e * (1 - mk) - (del_ + mu * f) * mk) * dt))
        D[i] = d; E[i] = e; Mk[i] = mk; fS[i] = f; Tc[i] = (1.0 / f) if f > 1e-6 else np.nan
    return t, D, E, Mk, fS, Tc


def steady(M, **kw):
    t, D, E, Mk, fS, Tc = simulate(lambda _t: M, T=3000, **kw)
    s = slice(-2000, None)
    return D[s].mean(), E[s].mean(), Mk[s].mean(), np.nanmean(Tc[s]), fS[s].mean()


if __name__ == '__main__':
    print("Graded G1-lengthening buffering prototype (GNP-grounded).\n")
    res = {}

    # ---- A. steady-state vs mitogen: does the cycle LENGTHEN gradually? EZH2 baseline floor? ----
    Ms = np.array([0.02, 0.05, 0.1, 0.2, 0.35, 0.55, 0.8, 1.2])
    print("=== A. steady state vs mitogen (graded G1 lengthening) ===")
    print(f"{'M':>5} | {'CyclinD1':>8} {'Tc(h)':>7} {'EZH2':>6} {'H3K27me3':>9}")
    dr = []
    for M in Ms:
        d, e, mk, tc, f = steady(M)
        dr.append((M, d, e, mk, tc, f)); print(f"{M:>5.2f} | {d:>8.2f} {tc:>7.1f} {e:>6.2f} {mk:>9.3f}")
    res['dr'] = np.array(dr)

    # ---- B. WITHDRAWAL: CyclinD1 fast, EZH2/mark slow (via cycle lengthening) = buffering ----
    print("\n=== B. mitogen withdrawal (M 0.9 -> 0.12 subthreshold at t=200h): timescale separation ===")
    st = 200.0
    def M_wd(t): return 0.9 if t < st else 0.12
    t, D, E, Mk, fS, Tc = simulate(M_wd, T=st + 500)
    res['wd'] = np.array([t, D, E, Mk, fS, Tc]); res['wd_settle'] = np.array([st])
    def thalf(x, base_hi):
        idx = np.where((t >= st) & (x <= 0.5 * base_hi + 0.5 * x[-1]))[0]
        return (t[idx[0]] - st) if len(idx) else np.nan
    D_hi = D[(t > st - 20) & (t < st)].mean(); E_hi = E[(t > st - 20) & (t < st)].mean()
    Mk_hi = Mk[(t > st - 20) & (t < st)].mean()
    print(f"  CyclinD1 half-fall {thalf(D, D_hi):.1f} h | EZH2 half-fall {thalf(E, E_hi):.1f} h | "
          f"H3K27me3 half-fall {thalf(Mk, Mk_hi):.1f} h  -> CyclinD1 fast, EZH2/mark slow (buffered)")

    # ---- C. contrast: graded lengthening vs abrupt G0 vs no-lengthening, under a transient DIP ----
    print("\n=== C. transient subthreshold DIP (M 0.9->0.12 for 60 h, then back to 0.9): EZH2/mark retained? ===")
    def M_dip(t): return 0.12 if (st <= t < st + 60) else 0.9
    modes = [('graded G1 lengthening', dict(lengthen=True, hard_G0=False)),
             ('abrupt G0 (current model)', dict(lengthen=False, hard_G0=True)),
             ('no lengthening (fixed Tc)', dict(lengthen=False, hard_G0=False))]
    print(f"{'mode':>26} | {'EZH2 min':>8} {'mark min':>8} {'mark @ end':>10}")
    for lab, kw in modes:
        t2, D2, E2, Mk2, fS2, Tc2 = simulate(M_dip, T=st + 300, **kw)
        w = (t2 >= st) & (t2 < st + 120)
        print(f"{lab:>26} | {E2[w].min():>8.2f} {Mk2[w].min():>8.3f} {Mk2[-1]:>10.3f}")
        res[f'dip_{lab.split()[0]}'] = np.array([t2, D2, E2, Mk2, fS2, Tc2])
    res['dip_settle'] = np.array([st])

    np.savez(os.path.join(HERE, 'sim_g1_lengthening_buffer_cache.npz'), **res)
    print("\n-> figure: fig_g1_lengthening_buffer.py")
