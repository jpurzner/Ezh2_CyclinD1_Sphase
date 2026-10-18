# H3K27me3 as a division speed-limiter — JP's framework (2026-07-21)

*JP's synthesis of the EZH2/H3K27me3→CyclinD1 axis. Demonstrated in the reduced mark equation:
`simulations/sim_speed_limiter.py` + `fig_speed_limiter.py`. Conceptual — the working v44 default is unchanged.*

## The idea
H3K27me3 is **not** a bistable latch — it is a **graded writing-vs-dilution rate-sensor**. The bistable element
is the **division decision** (go / no-go, the R-point), and CyclinD1-being-a-mark-target turns the sensor into a
**negative-feedback governor** on division speed.

## The self-consistent model (all from the equation)
```
Cd = Cd_drive(mitogen) · R(P)          R = mark repression of CyclinD1  (the governor)
r  = rate(Cd)                          division rate 1/Tc, commitment-thresholded
w  = W0·(WBAS + (1-WBAS)·r/(KR+r))      writer EZH2 = baseline + PROLIFERATION-COUPLED, SATURATING
P  = w / (δ + ln2·r)                    mark = writing / (turnover + replicative dilution)
```

## The four results (fig panels A–D)
- **A — speed limit.** P peaks at a moderate rate; the writer **saturates**, so dividing too fast lets `ln2/Tc`
  outrun `w` and P collapses below **θ_diff** → the H3K27me3-marked **differentiation genes de-repress**. There is a
  **max sustainable division rate**. A fast GNP still accumulates the mark *because* the writer is proliferation-coupled
  — up to that ceiling.
- **B — governor.** With the mark⊣CyclinD1 brake, division is held in the sustainable window and P is maintained;
  **without** it, mitogen drives division straight to max → P collapses → differentiation. The governor is the brake.
- **C — bistable division.** At intermediate mitogen the same input gives two outcomes by mark history
  (governed-slow vs runaway-fast). The bistability is in **divide-or-not**, not the mark (which is graded). Low mitogen
  → **locked out** (mark holds CyclinD1 below commitment).
- **D — regime map** (mitogen × oncogenic CyclinD1 drive): **locked out / sustainable proliferation / too-fast →
  differentiation**. **GNP** sits at the sustainable ceiling — pushed hardest (half the brain's neurons), right at the
  edge. **MB** overrides the governor via constitutive Gli/MYCN into the **too-fast / low-mark** regime — consistent
  with **ChIP MB<GNP** H3K27me3; MB escapes differentiation by oncogenic means, not the governor.

## How it reframes the model
The v44 pieces already exist — a monostable mark, the Skp2-p27-Rb R-point (the bistable commitment switch), the
mark⊣CyclinD1 that raises the commitment threshold, and the graded G1-lengthening (the rate-governing). This framework
**names what they collectively do: a self-limiting proliferation governor**. The one explicit piece still missing is the
H3K27me3 mark on **differentiation genes** (only CyclinD1 is marked today); a passive "differentiation-sentinel" locus
with the same writing/dilution balance would make the "too fast → differentiation" readout explicit.
