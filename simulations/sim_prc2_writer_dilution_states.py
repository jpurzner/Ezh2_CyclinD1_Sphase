"""Two proliferation STATES from a writing(EZH2)-vs-dilution(division) balance on H3K27me3 (JP, 2026-07-19).

JP's picture (GNP, no transient G0): H3K27me3 on differentiation genes is set by EZH2 WRITING vs replication
DILUTION, and BOTH scale with proliferation:
  - fast dividing: EZH2 very high (E2F-driven) + heavy per-cycle dilution -> writing OUTPACES dilution -> mark HIGH
    during division (silences differentiation genes -> proliferation sustained);
  - division stops: dilution falls away BUT EZH2 collapses -> no replacement -> mark lost -> differentiation genes
    de-repress and transcribe.
=> two states: fast-dividing/high-replacement (mark ON) vs arrested/no-replacement (mark OFF).

The key NEW element vs sim_prc2_buffering_memory: DILUTION is DIVISION-COUPLED (mu*p), so when a cell slows,
dilution drops and the mark is PRESERVED just as EZH2 falls -- a positive feedback that could make the ON state
STICKY (hysteresis/memory). Model (dimensionless, p = proliferation index in [0,1]):
  drive = b + M + s * Mk^n/(Ks^n + Mk^n)                 # CyclinD1 ~ basal + mitogen + mark-SUSTAIN (diff-off)
  p     = drive^hp/(Kp^hp + drive^hp)                    # proliferation index (division rate)
  dE/dt = (kE*p - E)/tauE                                # EZH2 tracks proliferation (E2F), lag tauE
  dMk/dt= kw*E*(1-Mk) - (del + mu*p)*Mk                  # written by EZH2; decay = turnover del + DILUTION mu*p
Mark maximal-at-division needs writing to outpace dilution: kw*kE > mu (tested). Run:
  ./venv/bin/python simulations/sim_prc2_writer_dilution_states.py   (writes _cache.npz; fig = fig_prc2_writer_dilution_states.py)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P = dict(b=0.05, s=0.9, n=4.0, Ks=0.45, Kp=0.9, hp=4.0, kE=1.0, tauE=0.6, kw=1.4, del_=0.02, mu=0.9)


def deriv(E, Mk, M, p):
    drive = P['b'] + M + P['s'] * Mk**P['n'] / (P['Ks']**P['n'] + Mk**P['n'])
    p_t = drive**P['hp'] / (P['Kp']**P['hp'] + drive**P['hp'])
    dE = (P['kE'] * p_t - E) / P['tauE']
    dMk = P['kw'] * E * (1 - Mk) - (P['del_'] + P['mu'] * p_t) * Mk
    return dE, dMk, p_t


def run(M_t, T=800, dt=0.01, E0=0.5, Mk0=0.5, mu=None):
    old = P['mu']
    if mu is not None: P['mu'] = mu
    n = int(T / dt); t = np.arange(n) * dt
    E = np.empty(n); Mk = np.empty(n); p = np.empty(n)
    e, mk = E0, Mk0
    for i in range(n):
        de, dmk, pt = deriv(e, mk, M_t(t[i]), None)
        e = max(0.0, e + de * dt); mk = min(1.0, max(0.0, mk + dmk * dt))
        E[i] = e; Mk[i] = mk; p[i] = pt
    P['mu'] = old
    return t, E, Mk, p


def steady(M, E0=0.5, Mk0=0.5, mu=None):
    t, E, Mk, p = run(lambda _t: M, T=600, E0=E0, Mk0=Mk0, mu=mu)
    return E[-1000:].mean(), Mk[-1000:].mean(), p[-1000:].mean()


if __name__ == '__main__':
    print("Writing(EZH2)-vs-dilution(division) two-state model.\n")
    res = {}

    # ---- A. dose-response: is the mark maximal at high division? does EZH2 track proliferation? ----
    Ms = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.65, 0.9, 1.2])
    dr = np.array([steady(M, E0=0.9, Mk0=0.9) for M in Ms])       # from the ON branch
    res['dr_M'] = Ms; res['dr'] = dr
    print("=== A. dose-response (ON branch): M | EZH2 | H3K27me3 | proliferation ===")
    for M, (e, mk, p) in zip(Ms, dr):
        print(f"  M={M:>4.2f} | EZH2 {e:>5.2f}  Mk {mk:>5.2f}  p {p:>5.2f}")
    print(f"  writing-outpaces-dilution check: kw*kE = {P['kw']*P['kE']:.2f} vs mu = {P['mu']:.2f} "
          f"-> {'mark rises/stays with division (max at division) ✓' if P['kw']*P['kE'] > P['mu'] else 'dilution wins'}")

    # ---- B. HYSTERESIS: ramp mitogen up then down; bistable? Compare division-coupled dilution vs constant turnover ----
    print("\n=== B. mitogen up/down ramp -> hysteresis (bistability = memory) ===")
    def hysteresis(mu):
        up = np.linspace(0.0, 0.8, 40); dn = up[::-1]
        pu, pd, mu_up, mu_dn = [], [], [], []
        e, mk = 0.05, 0.05                                        # start OFF (arrested/differentiated)
        for M in up:
            t, E, Mk, p = run(lambda _t: M, T=200, E0=e, Mk0=mk, mu=mu); e, mk = E[-1], Mk[-1]
            pu.append(p[-1]); mu_up.append(mk)
        for M in dn:
            t, E, Mk, p = run(lambda _t: M, T=200, E0=e, Mk0=mk, mu=mu); e, mk = E[-1], Mk[-1]
            pd.append(p[-1]); mu_dn.append(mk)
        return up, np.array(pu), dn, np.array(pd), np.array(mu_up), np.array(mu_dn[::-1])
    for mu, lab in [(P['mu'], 'division-coupled dilution (mu=0.9)'), (0.0, 'constant turnover (mu=0)')]:
        up, pu, dn, pd, mku, mkd = hysteresis(mu)
        # hysteresis width: max gap between down-ramp and up-ramp proliferation at matched M
        gap = float(np.max(pd[::-1] - pu))
        # switch-off mitogen (down ramp) vs switch-on mitogen (up ramp)
        m_on = up[np.argmax(pu > 0.5)] if (pu > 0.5).any() else np.nan
        m_off = dn[np.argmax(pd < 0.5)] if (pd < 0.5).any() else np.nan
        print(f"  {lab}: switch ON at M={m_on:.2f}, switch OFF at M={m_off:.2f}, hysteresis width Δp={gap:.2f}")
        res[f'hys_{int(mu*10)}'] = np.array([up, pu, dn, pd, mku, mkd], dtype=object)

    # ---- C. WITHDRAWAL memory window: from fast-dividing, M->0; how long do Mk (mark) & p (division) persist? ----
    print("\n=== C. withdrawal from fast division (M 0.6 -> 0): memory window ===")
    st = 200
    def M_wd(t): return 0.6 if t < st else 0.0
    for mu, lab in [(P['mu'], 'mu=0.9 (dilution stops when division stops)'), (0.0, 'mu=0 (constant turnover)')]:
        t, E, Mk, p = run(M_wd, T=st + 600, E0=0.9, Mk0=0.9, mu=mu)
        after = t >= st
        # time for the MARK to fall to half, and whether proliferation LATCHES (p>0.3 at end)
        mk_hi = Mk[(t > st - 20) & (t < st)].mean()
        idx = np.where(after & (Mk <= 0.5 * mk_hi))[0]
        thalf = (t[idx[0]] - st) if len(idx) else 600
        banked = np.trapz(p[after], t[after])
        print(f"  {lab}: mark half-life {thalf:>5.1f}, divisions banked after withdrawal {banked:>6.1f}, "
              f"final p {p[-1]:.2f} ({'LATCH (persists)' if p[-1] > 0.3 else 'arrests'})")
        res[f'wd_{int(mu*10)}'] = np.array([t, E, Mk, p]); res['wd_settle'] = np.array([st])

    # ---- D. the TWO STATES side by side (ON vs OFF attractor at the same low mitogen) ----
    print("\n=== D. two coexisting states at M=0.15 (bistable region) ===")
    for E0, Mk0, lab in [(0.9, 0.9, 'from ON (fast-dividing)'), (0.05, 0.05, 'from OFF (arrested)')]:
        e, mk, p = steady(0.15, E0=E0, Mk0=Mk0)
        print(f"  {lab}: EZH2 {e:.2f}  H3K27me3 {mk:.2f}  proliferation {p:.2f}  dilution(mu*p) {P['mu']*p:.2f}")

    np.savez(os.path.join(HERE, 'sim_prc2_writer_dilution_states_cache.npz'), **res)
    print("\n-> figure: fig_prc2_writer_dilution_states.py")
