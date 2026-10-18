"""JP's SPEED-LIMITER framework in the mark equation (2026-07-21): H3K27me3 is a graded
writing-vs-dilution rate-sensor; the writer (EZH2) is PROLIFERATION-COUPLED so dividing drives
writing; writing SATURATES so there is a MAX SUSTAINABLE division rate (past it the mark collapses
-> differentiation genes de-repress). CyclinD1-as-a-target = a NEGATIVE-FEEDBACK GOVERNOR that holds
division at a sustainable rate. Bistability is in the DIVISION decision (go/no-go), NOT the mark.

Self-consistent coupled model (all from the equation):
    Cd    = Cd_drive(mitogen) * R(P)                 R = mark repression of CyclinD1 (governor)
    r     = rate(Cd)                                  commitment-thresholded division rate (1/Tc)
    w     = W0*(WBAS + r/(KR+r))                      writer EZH2 = proliferation-coupled + SATURATING
    P     = min(1, w/(DEL + ln2*r))                   mark = writing / (turnover + replicative dilution)
Solve the fixed point in P for each mitogen, WITH vs WITHOUT the R(P) governor.
Run: ./venv/bin/python simulations/sim_speed_limiter.py"""
import numpy as np
LN2 = np.log(2)
DEL = 0.020            # mark turnover (1/h)
CD_MAX, KG = 3.7, 0.42
F0, KP, NP = 0.18, 0.50, 4.0      # R(P): high mark -> strong CyclinD1 repression
RATE_MAX, KC, H = 0.180, 1.00, 6.0  # rate(Cd): commitment threshold at Cd~KC; max 1/Tc (fast enough to outrun the writer if ungoverned)
W0, WBAS, KR = 0.045, 0.311, 0.0217  # writer: baseline (WBAS fraction) -> proliferation-coupled -> SATURATES at W0
THETA_DIFF = 0.50      # mark below this -> differentiation genes de-repress

def cd_drive(m): return CD_MAX * m ** 2 / (KG ** 2 + m ** 2)
def R(P): return F0 + (1 - F0) / (1 + (P / KP) ** NP)
def rate(Cd): return RATE_MAX * Cd ** H / (KC ** H + Cd ** H)
def writer(r): return W0 * (WBAS + (1 - WBAS) * r / (KR + r))   # EZH2: baseline + proliferation-coupled, saturating
def mark_of_rate(r): return min(1.0, writer(r) / (DEL + LN2 * r))

def solve(m, feedback=True, P0=1.0, drive=1.0):
    """drive = oncogenic CyclinD1 multiplier (GNP=1; MB>1 for constitutive Gli/MYCN)."""
    P = P0
    for _ in range(400):
        Cd = drive * cd_drive(m) * (R(P) if feedback else 1.0)
        r = rate(Cd)
        P = 0.6 * P + 0.4 * mark_of_rate(r)
    Cd = drive * cd_drive(m) * (R(P) if feedback else 1.0); r = rate(Cd)
    return dict(P=P, Cd=Cd, r=r, Tc=(1 / r if r > 1e-4 else np.inf))

if __name__ == '__main__':
    print("mark maintenance P(r) vs division rate r (writer saturates -> mark collapses when too fast):")
    for r in (0.0, 0.02, 0.04, 0.0625, 0.09, 0.13, 0.18):
        tc = 1 / r if r > 0 else np.inf
        print(f"  r={r:.3f} (Tc={tc:5.1f}h)  w={writer(r):.4f}  dilution={LN2*r:.4f}  P={mark_of_rate(r):.3f}"
              f"  {'<-- below diff threshold (DIFFERENTIATION)' if mark_of_rate(r) < THETA_DIFF else ''}")
    print(f"\n  max sustainable rate (P=θ_diff={THETA_DIFF}): ", end='')
    rr = np.linspace(0.001, 0.3, 3000); Pr = np.array([mark_of_rate(x) for x in rr])
    sus = rr[Pr >= THETA_DIFF]
    print(f"r_max={sus.max():.3f}  (Tc_min={1/sus.max():.1f}h)")

    print("\nGovernor + BISTABLE division test (solve from HIGH-mark vs LOW-mark init):")
    print(f"  {'m':>5} | {'r(hi)':>7} {'r(lo)':>7} {'split?':>6} | {'P(gov)':>7} Tc(gov) | {'r(no-gov)':>9} {'P(nog)':>7}")
    for m in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.2):
        gh = solve(m, True, P0=1.0); gl = solve(m, True, P0=0.0); n = solve(m, False)
        split = 'YES' if abs(gh['r'] - gl['r']) > 0.01 else ''
        print(f"  {m:>5.2f} | {gh['r']:>7.4f} {gl['r']:>7.4f} {split:>6} | {gh['P']:>7.3f} {gh['Tc']:>6.1f} | {n['r']:>9.4f} {n['P']:>7.3f}")
