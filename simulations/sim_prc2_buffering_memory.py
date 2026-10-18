"""PRC2 / H3K27me3 as a BUFFERING MEMORY for proliferation: under what conditions does the mark let cell
division PERSIST through mitogen fluctuations + withdrawal? (JP overnight, 2026-07-19.)

Constraints (JP, new): NO transient G0 in GNP; the mark ACCUMULATES DURING MAXIMAL DIVISION and DECLINES with
differentiation (ChIP-seq). So the WRITER is proliferation-coupled (W ~ proliferation) in BOTH topologies below,
and we report the mark-timing to confirm it peaks at high division.

TWO topologies for "PRC2 acts on the CyclinD1 axis", scale-free (time in CyclinD1-decay units):
  A ("REPRESS"):  mark REPRESSES CyclinD1 directly (JP's literal description).
        drive = (b + M);  rep(Mk) = 1 - g*Mk^n/(Krep^n+Mk^n);  dD/dt = kD*drive*rep - D
  C ("SUSTAIN"):  mark represses an ANTI-proliferative brake -> NET SUSTAINS the proliferative drive
        (double-negative; = PRC2 silencing differentiation genes to hold the proliferative state).
        dD/dt = kD*(b + M + s*Mk^n/(Ks^n+Mk^n)) - D          # mark ADDS a mitogen-INDEPENDENT drive
  common:  pi(D)=D^hp/(Kp^hp+D^hp) (division rate);  W=pi (proliferation-coupled writer);
           dMk/dt = kw*W*(a0 + a_rw*Mk)*(1-Mk) - del*Mk       # cooperative read-write + turnover (latch if a_rw>>a0, n>=2)

Feedback OFF (EZH2i) = no mark action (g=0 / s=0). The MEMORY is the slow Mk. KEY question the sweeps answer:
  - Topology A buffers the DOSE-RESPONSE (flattens D vs M) but CANNOT sustain division at M=0 (rep only scales the
    mitogen term -> D->kD*b*rep <= unrepressed basal). So A gives robustness to SMALL/slow dips, NOT withdrawal memory.
  - Topology C SUSTAINS division after withdrawal for ~1/del (mark-decay-limited), strongly if the read-write LATCHES
    (a_rw>>a0, n>=2): once proliferating (mark ON), it holds through withdrawal until the mark decays -> a genuine
    proliferative MEMORY, and it decays with differentiation (mark lost) -> matches the ChIP-seq.

Run: ./venv/bin/python simulations/sim_prc2_buffering_memory.py   (writes _cache.npz; figure = fig_prc2_buffering_memory.py)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def simulate(M_t, topo, g, deln, n, s=0.0, T=400.0, dt=0.02, kD=6.0, b=0.06, Kp=1.0, hp=4.0,
             Krep=0.5, Ks=0.5, kw=1.0, a0=0.05, a_rw=1.0, D0=0.5, Mk0=0.3):
    nstep = int(T / dt); t = np.arange(nstep) * dt
    D = np.empty(nstep); Mk = np.empty(nstep); pi = np.empty(nstep)
    d, mk = D0, Mk0
    for i in range(nstep):
        M = M_t(t[i])
        if topo == 'A':
            drive = (b + M) * (1.0 - g * mk**n / (Krep**n + mk**n))
        else:  # 'C' sustain
            drive = b + M + s * mk**n / (Ks**n + mk**n)
        p = d**hp / (Kp**hp + d**hp)
        d = max(0.0, d + (kD * drive - d) * dt)
        mk = min(1.0, max(0.0, mk + (kw * p * (a0 + a_rw * mk) * (1.0 - mk) - deln * mk) * dt))
        D[i] = d; Mk[i] = mk; pi[i] = p
    return t, D, Mk, pi


def steady(M, **kw):
    _, D, Mk, pi = simulate(lambda _t: M, T=500, **kw)
    return D[-800:].mean(), Mk[-800:].mean(), pi[-800:].mean()


def dose_response(topo, g, deln, n, s, a_rw=1.0):
    Ms = np.array([0.03, 0.06, 0.1, 0.15, 0.22, 0.32, 0.45, 0.65, 0.9, 1.2])
    out = np.array([steady(M, topo=topo, g=g, deln=deln, n=n, s=s, a_rw=a_rw) for M in Ms])
    D, Mk = out[:, 0], out[:, 1]
    band = (Ms >= 0.15) & (Ms <= 0.9)
    gain = np.polyfit(np.log(Ms[band]), np.log(D[band] + 1e-9), 1)[0]
    mark_up = Mk[band][-1] - Mk[band][0]
    return Ms, D, Mk, gain, mark_up


def persistence(topo, g, deln, n, s, a_rw=1.0, M_hi=0.9, M_lo=0.0, settle=250, post=200):
    def M_t(t): return M_hi if t < settle else M_lo
    t, D, Mk, pi = simulate(M_t, topo=topo, g=g, deln=deln, n=n, s=s, a_rw=a_rw, T=settle + post)
    after = t >= settle
    banked = np.trapz(pi[after], t[after])
    pi_hi = pi[(t > settle - 20) & (t < settle)].mean()
    idx = np.where((t >= settle) & (pi <= 0.5 * pi_hi))[0]
    thalf = (t[idx[0]] - settle) if len(idx) else post
    return banked, thalf, pi_hi, (t, D, Mk, pi, settle)


def fluct_cv(topo, g, deln, n, s, M0=0.35, amp=0.6, period=40.0, settle=250, cycles=8):
    def M_t(t): return M0 * (1.0 + amp * np.sin(2 * np.pi * (t - settle) / period)) if t >= settle else M0
    t, D, Mk, pi = simulate(M_t, topo=topo, g=g, deln=deln, n=n, s=s, T=settle + cycles * period)
    w = t >= settle + period
    return pi[w].std() / (pi[w].mean() + 1e-9)


if __name__ == '__main__':
    print("PRC2 buffering-memory exploration (reduced). Writer = proliferation-coupled (mark peaks at division).\n")
    res = {}

    # ============ A. topology A (mark REPRESSES CyclinD1): dose-response buffering, but NO withdrawal memory ============
    print("=== Topology A (mark REPRESSES CyclinD1) : dose-response buffering + withdrawal (fair, baseline-normalized) ===")
    print(f"{'g':>5} {'1/del':>6} {'n':>3} | {'DRgain':>6} {'markUp':>6} | {'thalf_ON':>8} {'thalf_OFF':>9} {'basePi_ON':>9}")
    A = []
    for n in (1.0, 2.0, 4.0):
        for deln in (0.15, 0.05, 0.015):
            for g in (0.0, 0.6, 0.85, 0.97):
                Ms, D, Mk, gain, mark_up = dose_response('A', g, deln, n, 0.0)
                b_on, th_on, pih_on, _ = persistence('A', g, deln, n, 0.0)
                A.append((g, deln, n, gain, mark_up, b_on, th_on, pih_on))
                if g in (0.0, 0.85):
                    th_off = A[[i for i, r in enumerate(A) if r[0] == 0.0 and r[1] == deln and r[2] == n][0]][6]
                    print(f"{g:>5.2f} {1/deln:>6.1f} {n:>3.0f} | {gain:>6.3f} {mark_up:>6.3f} | {th_on:>8.1f} {th_off:>9.1f} {pih_on:>9.3f}")
    res['A'] = np.array(A)

    # dose-response curves for the figure (representative: n=2, slow mark), g on/off
    Ms, D0, Mk0, _, _ = dose_response('A', 0.0, 0.05, 2.0, 0.0)
    _, Dg, Mkg, _, _ = dose_response('A', 0.85, 0.05, 2.0, 0.0)
    res['A_dr'] = np.array([Ms, D0, Dg, Mk0, Mkg])

    # ============ C. topology C (mark SUSTAINS proliferation): find the GRADED-MEMORY sweet spot ============
    # We want: finite-but-LONG persistence memory (bank_extra>0, thalf>>off) AND the mark still DECLINES with
    # differentiation (mark_up>0, i.e. NOT a permanent latch that stays ON at low mitogen). Sweep sustain s x
    # mark timescale 1/del x read-write cooperativity (a_rw: low=monostable/graded, high=bistable/latch).
    print("\n=== Topology C (mark SUSTAINS CyclinD1): MONOSTABLE (low read-write) vs BISTABLE (high read-write) ===")
    print("  want: finite-but-long memory (thalf>>off, not censored) AND mark DECLINES after withdrawal (fits ChIP-seq)")
    print(f"{'a_rw':>5} {'s':>5} {'1/del':>6} | {'bank_ex':>7} {'thalf':>6} {'markDrop':>8} {'markUp':>6} {'latch?':>6} {'MEM?':>5}")
    C = []; POST = 600
    for a_rw in (0.05, 2.0):                                 # 0.05 = monostable (mark tracks proliferation); 2.0 = bistable latch
        for deln in (0.03, 0.012):
            for s in (0.0, 0.2, 0.4, 0.7, 1.0, 1.4, 2.0):
                b_c, th_c, pih, trc = persistence('C', 0.0, deln, 4.0, s, a_rw=a_rw, post=POST)
                b_off, th_off, _, _ = persistence('C', 0.0, deln, 4.0, 0.0, a_rw=a_rw, post=POST)
                t, D, Mk, pi, st = trc
                hi = Mk[(t > st - 20) & (t < st)].mean(); lo = Mk[-1]
                mark_drop = (hi - lo) / (hi + 1e-9)
                latched = pi[-1] > 0.3
                _, _, _, gx, mup = dose_response('C', 0.0, deln, 4.0, s, a_rw=a_rw)
                memory = (not latched) and (th_c > 2 * th_off) and (mark_drop > 0.4)
                C.append((a_rw, s, deln, b_c - b_off, th_c, mark_drop, mup, int(latched), int(memory)))
                if deln == 0.012:
                    print(f"{a_rw:>5.2f} {s:>5.1f} {1/deln:>6.1f} | {b_c-b_off:>7.2f} {th_c:>6.1f} {mark_drop:>8.2f} {mup:>6.3f} "
                          f"{'YES' if latched else 'no':>6} {'yes' if memory else '-':>5}")
    res['C'] = np.array(C)

    # representative WITHDRAWAL traces for the figure: A (repress, no memory) vs feedback-off vs C-latch (memory)
    def pack(tr): t, D, Mk, pi, st = tr; return np.array([t, D, Mk, pi])
    _, _, _, trA = persistence('A', 0.85, 0.03, 2.0, 0.0, post=250)
    _, _, _, trOff = persistence('A', 0.0, 0.03, 2.0, 0.0, post=250)
    _, _, _, trC = persistence('C', 0.0, 0.02, 4.0, 0.7, post=250)     # supercritical sustain -> latch = memory
    res['trace_A'] = pack(trA); res['trace_off'] = pack(trOff); res['trace_C'] = pack(trC)
    res['settle'] = np.array([trA[4]])

    # ---- DEMO 1: latch persists through a transient mitogen DIP (fluctuation) but a committed differentiation
    #      erasure RELEASES it. Shows the two-input logic: mitogen withdrawal != differentiation. ----
    st = 250              # match the withdrawal-trace settle so the figure aligns (t - settle)
    def M_dip(t):        # equilibrate high, a deep transient dip, then recover to a SUBTHRESHOLD level (0.12)
        if t < st: return 0.9
        if t < st + 40: return 0.0
        return 0.12       # below the no-memory proliferation threshold -> only the latch sustains division here
    t1, D1, Mk1, pi1 = simulate(M_dip, topo='C', g=0.0, deln=0.02, n=4.0, s=0.7, a_rw=1.0, T=st + 300)
    # same protocol WITHOUT the sustain memory (s=0): division collapses during the dip
    t1b, D1b, Mk1b, pi1b = simulate(M_dip, topo='C', g=0.0, deln=0.02, n=4.0, s=0.0, a_rw=1.0, T=st + 300)
    res['dip_mem'] = np.array([t1, [M_dip(x) for x in t1], D1, pi1, Mk1])
    res['dip_nomem'] = np.array([t1b, [M_dip(x) for x in t1b], D1b, pi1b, Mk1b])

    # DEMO 2: latch (memory) vs a DIFFERENTIATION erasure pulse (active demethylase: transient boost of del at t=350)
    def M_wd(t): return 0.9 if t < st else 0.0
    def run_erase(erase):
        nstep = int((st + 250) / 0.02); tt = np.arange(nstep) * 0.02
        d, mk = 0.5, 0.3; Dr = np.empty(nstep); Mr = np.empty(nstep); Pr = np.empty(nstep)
        for i in range(nstep):
            M = M_wd(tt[i])
            deln = 0.02 + (0.6 if (erase and st + 150 <= tt[i] < st + 180) else 0.0)   # differentiation demethylase burst (t-settle in [150,180])
            drive = 0.06 + M + 0.7 * mk**4 / (0.5**4 + mk**4)
            p = d**4 / (1.0 + d**4)
            d = max(0.0, d + (6.0 * drive - d) * 0.02)
            mk = min(1.0, max(0.0, mk + (1.0 * p * (0.05 + 1.0 * mk) * (1 - mk) - deln * mk) * 0.02))
            Dr[i] = d; Mr[i] = mk; Pr[i] = p
        return tt, Dr, Mr, Pr
    te0, De0, Me0, Pe0 = run_erase(False)     # withdrawal only -> latch persists (memory)
    te1, De1, Me1, Pe1 = run_erase(True)      # withdrawal + differentiation erasure -> latch released
    res['erase_off'] = np.array([te0, De0, Me0, Pe0]); res['erase_on'] = np.array([te1, De1, Me1, Pe1])
    print(f"\n  DEMO latch through a mitogen DIP: with sustain-memory divisions persist (final pi {pi1[-1]:.2f}) "
          f"vs collapse without (final pi {pi1b[-1]:.2f}).")
    print(f"  DEMO differentiation erasure: withdrawal-only keeps latch (final pi {Pe0[-1]:.2f}); "
          f"withdrawal+erasure releases it (final pi {Pe1[-1]:.2f}).")

    # ============ fluctuation buffering (topology A), frequency-resolved ============
    print("\n=== Fluctuation buffering (topology A), frequency-resolved: CV_ON/CV_OFF (<1 => buffered) ===")
    periods = np.array([8, 15, 30, 60, 120.0])
    print(f"{'period':>7} | " + " ".join(f"g={g}" for g in (0.0, 0.6, 0.85)))
    F = []
    for T in periods:
        row = [fluct_cv('A', g, 0.05, 2.0, 0.0, period=T) for g in (0.0, 0.6, 0.85)]
        F.append([T] + row)
        print(f"{T:>7.0f} | " + " ".join(f"{v:>5.3f}" for v in row))
    res['F'] = np.array(F)

    np.savez(os.path.join(HERE, 'sim_prc2_buffering_memory_cache.npz'), **res)
    # ---- verdict ----
    A = res['A']; th_on = A[:, 6]
    # persistence gain of topology A vs its own feedback-off (paired by del,n)
    gains_A = []
    for (g, deln, n, *_rest), th in zip(A, th_on):
        if g == 0.0: continue
        th_off = A[(A[:, 0] == 0) & (A[:, 1] == deln) & (A[:, 2] == n)][0, 6]
        gains_A.append(th - th_off)
    C = res['C']
    print(f"\nVERDICT:")
    print(f"  Topology A (repress CyclinD1): max withdrawal-persistence gain vs feedback-off = {max(gains_A):+.2f} "
          f"(median {np.median(gains_A):+.2f})  -> A does NOT create withdrawal memory.")
    n_mem = int(C[:, 8].sum()); n_latch = int(C[:, 7].sum())     # memory=col8, latched=col7
    print(f"  Topology C (sustain): GRADED-memory regimes (finite persistence + mark declines, fits ChIP-seq): {n_mem}/{len(C)}; "
          f"permanent latches (mark stays ON, violates ChIP-seq): {n_latch}/{len(C)}.")
    mem = C[C[:, 8] == 1]
    if len(mem):
        mono = mem[mem[:, 0] < 0.5]
        print(f"  Graded-memory band: a_rw in {{{','.join(sorted(set(f'{x:.2f}' for x in mem[:,0])))}}}, s in "
              f"[{mem[:,1].min():.2f}, {mem[:,1].max():.2f}]; thalf up to {mem[:,4].max():.0f} (vs off ~2), "
              f"mark_drop {mem[:,5].min():.2f}-{mem[:,5].max():.2f} (mark declines with differentiation). "
              f"{'MONOSTABLE (low read-write) supplies the graded band.' if len(mono) else ''}")
    print("  -> figure: fig_prc2_buffering_memory.py")
