"""Prototype (reduced 1-locus stochastic model): does PUNCTATE transcription resolve the tension --
MB<GNP (ChIP) + a LONGER-persisting H3K27me3 mark, at a WEAK Gli eraser -- where continuous eviction cannot?

Mark M in [0,1] at the Ccnd1 locus:
  dM/dt = write - erase - evict*M
  write = kw*EZH2*(a0 + a_rw*M)*(1-M)        # PRC2 read-write (cooperative), blocked when EZH2i
  erase = (del + ke*Gli)*M                     # turnover + Gli->Jmjd3/Kdm6b eraser
  evict = g*B(t)  [PUNCTATE, B in {0,1} bursts]   OR   g*<B>=g*lam/(lam+mu)  [CONTINUOUS, matched mean]
Transcription = Poisson bursts: on-rate lam (proportional to transcription, HIGHER in MB), off-rate mu.
MB: EZH2 high, Gli high, lam high (Cd_mRNA~15).  GNP: EZH2 lower, Gli low, lam low (Cd_mRNA~3).
Compare mean steady-state M (MB vs GNP) and mark half-life after EZH2i, punctate vs continuous, weak vs strong eraser."""
import numpy as np
KW, A0, A_RW, DEL, G, MU = 0.020, 0.05, 0.45, 0.0010, 0.030, 0.10
DT, T = 0.5, 260000.0

def sim(EZH2, Gli, lam, punctate, ke, seed, write_on=True, M0=0.6):
    rng = np.random.default_rng(seed); M = M0; B = 0
    mean_evict = G * lam / (lam + MU)
    n = int(T / DT); burn = n // 2; acc = []
    for i in range(n):
        if punctate:
            if B == 0:
                if rng.random() < lam * DT: B = 1
            else:
                if rng.random() < MU * DT: B = 0
            evict = G * B
        else:
            evict = mean_evict
        w = (KW * EZH2 * (A0 + A_RW * M) * (1 - M)) if write_on else 0.0
        M += (w - (DEL + ke * Gli) * M - evict * M) * DT
        M = min(1.0, max(0.0, M))
        if i >= burn: acc.append(M)
    return float(np.mean(acc))

def halflife(EZH2, Gli, lam, punctate, ke, seed):
    """settle with writing, then EZH2i (write off); return time (in DT-units*DT) for M to halve."""
    rng = np.random.default_rng(seed); M = 0.7; B = 0
    mean_evict = G * lam / (lam + MU)
    # settle
    for i in range(int(120000 / DT)):
        if punctate:
            if B == 0:
                if rng.random() < lam * DT: B = 1
            elif rng.random() < MU * DT: B = 0
            evict = G * B
        else: evict = mean_evict
        M += (KW * EZH2 * (A0 + A_RW * M) * (1 - M) - (DEL + ke * Gli) * M - evict * M) * DT
        M = min(1.0, max(0.0, M))
    M0 = M; t = 0.0
    for i in range(int(300000 / DT)):
        if punctate:
            if B == 0:
                if rng.random() < lam * DT: B = 1
            elif rng.random() < MU * DT: B = 0
            evict = G * B
        else: evict = mean_evict
        M += (0.0 - (DEL + ke * Gli) * M - evict * M) * DT   # writing OFF (EZH2i)
        M = min(1.0, max(0.0, M)); t += DT
        if M <= M0 / 2: return t
    return t

MB  = dict(EZH2=3.6, Gli=0.19, lam=0.10)   # high transcription -> frequent bursts
GNP = dict(EZH2=2.0, Gli=0.027, lam=0.02)  # low transcription -> rare bursts

for ke, ke_name in [(0.22, 'STRONG eraser (current)'), (0.03, 'WEAK eraser')]:
    print(f"\n=== {ke_name} (ke={ke}) ===")
    for mode, punc in [('continuous', False), ('PUNCTATE', True)]:
        mk_mb = np.mean([sim(**MB, punctate=punc, ke=ke, seed=s) for s in range(4)])
        mk_gnp = np.mean([sim(**GNP, punctate=punc, ke=ke, seed=s) for s in range(4)])
        hl_mb = np.mean([halflife(**MB, punctate=punc, ke=ke, seed=s) for s in range(3)])
        print(f"  {mode:11}: GNP M={mk_gnp:.3f}  MB M={mk_mb:.3f}  MB/GNP={mk_mb/mk_gnp:.2f}  | MB half-life={hl_mb:.0f} (arb units)")
