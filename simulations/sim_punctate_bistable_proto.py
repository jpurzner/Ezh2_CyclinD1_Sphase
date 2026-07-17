"""Prototype 2: BISTABLE mark + STRONG punctate bursts. Cooperative read-write (M^n) + suppressed nucleation (small
a0) makes the locus bistable (latched ON=high-mark or OFF=low-mark). STRONG, brief transcriptional bursts can switch
a locus OFF; between bursts the cooperative read-write re-latches it ON. MB (frequent bursts, high transcription)
spends more time OFF -> lower mark; GNP (rare bursts) stays latched ON -> higher mark. Test: does this give (a)
MB<GNP without a strong eraser, and (b) LONG persistence (a latched ON mark is stable for long stretches)?
Report the population-fraction-in-ON and the mean mark, continuous vs punctate, across cooperativity n and burst strength."""
import numpy as np
KW, DEL, MU = 0.025, 0.0008, 0.20
DT, T = 0.5, 400000.0

def sim(EZH2, Gli, lam, punctate, ke, a0, a_rw, n, g, seed, write_on=True, M0=0.6):
    rng = np.random.default_rng(seed); M = M0; B = 0
    mean_evict = g * lam / (lam + MU)
    nsteps = int(T / DT); burn = nsteps // 2; acc = []
    for i in range(nsteps):
        if punctate:
            if B == 0:
                if rng.random() < lam * DT: B = 1
            elif rng.random() < MU * DT: B = 0
            evict = g * B
        else:
            evict = mean_evict
        w = (KW * EZH2 * (a0 + a_rw * M**n) * (1 - M)) if write_on else 0.0
        M += (w - (DEL + ke * Gli) * M - evict * M) * DT
        M = min(1.0, max(0.0, M))
        if i >= burn: acc.append(M)
    acc = np.array(acc)
    return float(np.mean(acc)), float(np.mean(acc > 0.4))   # mean mark, fraction-of-time ON

MB  = dict(EZH2=3.6, Gli=0.19, lam=0.06)
GNP = dict(EZH2=2.0, Gli=0.027, lam=0.012)

print("BISTABLE + strong bursts. ke=0.03 (WEAK eraser). scan cooperativity n and burst strength g.")
for n in (2, 4):
    for g in (0.15, 0.35):
        a0, a_rw = 0.004, 0.9
        print(f"\n n={n} g={g} a0={a0} a_rw={a_rw}:")
        for mode, punc in [('continuous', False), ('PUNCTATE', True)]:
            mb = [sim(**MB, punctate=punc, ke=0.03, a0=a0, a_rw=a_rw, n=n, g=g, seed=s) for s in range(4)]
            gnp = [sim(**GNP, punctate=punc, ke=0.03, a0=a0, a_rw=a_rw, n=n, g=g, seed=s) for s in range(4)]
            mb_m, mb_on = np.mean([x[0] for x in mb]), np.mean([x[1] for x in mb])
            gn_m, gn_on = np.mean([x[0] for x in gnp]), np.mean([x[1] for x in gnp])
            ratio = mb_m / gn_m if gn_m > 1e-6 else np.nan
            print(f"   {mode:11}: GNP M={gn_m:.2f}(ON {gn_on:.0%})  MB M={mb_m:.2f}(ON {mb_on:.0%})  MB/GNP={ratio:.2f}")
