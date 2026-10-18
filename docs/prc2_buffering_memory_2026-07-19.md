# PRC2 / H3K27me3 as a BUFFERING MEMORY for proliferation — conditions (2026-07-19)

*JP overnight: "explore under what conditions the PRC2 inhibition on CyclinD1 would produce a buffering memory
so cell division persists through fluctuations and mitogen withdrawal. Note: NO real transient G0s in GNP, so
H3K27me3 ACCUMULATES DURING MAXIMAL DIVISION (ChIP-seq maximal at max division, declining with differentiation)."*

Harness: `simulations/sim_prc2_buffering_memory.py` (reduced, scale-free, writer proliferation-coupled) +
`fig_prc2_buffering_memory.py`. Grounded in the full v44 model via the (regenerated) `sim_ezh2_mitogen_buffering`
and `sim_mitogen_withdrawal`.

## Bottom line
There are **two different jobs** a PRC2/H3K27me3 mark on the CyclinD1 axis can do, and they need **opposite circuit
topologies**. The one JP literally describes (mark REPRESSES CyclinD1) gives buffering of the **dose-response**, but
it **cannot** produce the persistence memory. The persistence memory (division continues through fluctuations +
withdrawal) requires the mark to **SUSTAIN the proliferative state** — i.e. PRC2 silencing anti-proliferative /
differentiation genes — as a bistable **latch**. Both are consistent with "mark maximal during division"; only the
latch gives a memory, and its decline-with-differentiation is an **active erasure**, distinct from mitogen loss.

## Topology A — mark ⊣ CyclinD1 (negative feedback; what the current v44 model has)
- **Buffers the steady-state DOSE-RESPONSE.** As mitogen rises, the proliferation-written mark rises and represses
  CyclinD1, flattening CyclinD1(mitogen): reduced-model gain 0.85→0.53; **full v44: 0.85→0.47** (sim_ezh2_mitogen_
  buffering). The mark is maximal at high mitogen/division ✓ (JP's constraint).
- **Does NOT buffer fluctuations.** Division-rate CV *rises* with feedback at every frequency (0.06→0.6) — the
  negative feedback + slow mark adds phase lag near the steep proliferation threshold.
- **Does NOT create withdrawal memory — mathematically can't.** A repressor only *scales* the mitogen-dependent
  drive, so at M=0 CyclinD1 → basal regardless of the mark. Reduced model: feedback *shortens* the proliferation
  half-life after withdrawal (2.1→0.4–1.1). **Full v44 confirms:** feedback ON gives **0.6–1.2 fewer** divisions
  after withdrawal than OFF (it sharpens arrest, it does not sustain division).
- So topology A = *setpoint robustness* (a real, useful buffering of where CyclinD1 sits vs mitogen), NOT a memory.
- **Short-dip persistence is a pure cell-cycle (Rb-E2F R-point) memory that the mark does NOT help — and slightly
  hurts.** The bistable restriction point lets a committed cell ride through a sub-cycle dip and still divide, and it
  does so *with or without* the mark (identical at feedback-off). Adding the topology-A mark makes CyclinD1 lower
  (repressed), so the committed cell sits *closer* to de-committing: in a fair comparison (feedback-ON vs
  feedback-OFF, not vs the un-repressed raw drive) the mark gives *fewer* banked divisions through a dip, and reduces
  the R-point rescue. The apparent "buffering benefit" in sim_ezh2_mitogen_buffering is relative to the *un-repressed
  drive* (f0_prc2=1), i.e. it shows CyclinD1 is flatter than the raw Gli/MYCN drive — not that division is more
  persistent than with no mark. So all sub-cycle-dip persistence is the **R-point (cell-cycle) memory**, independent
  of chromatin; the chromatin mark contributes nothing positive to it.

**Adversarial verification (independent agent, 14 tests): CLAIM CONFIRMED — and it generalizes.** No topology-A
regime beat feedback-off on withdrawal persistence across 4 Hill exponents × 4 mark timescales × 4 basal levels × 3
Krep × 3 read-write × 4 gains (incl. n=8 cooperative, basal 0.5, half-life to 125): best banked-division gain of any
g>0 over its own g=0 = **−0.0001** (a numerical tie in the degenerate limit; never positive). The bound is not just
the multiplicative form: for ANY net-repressive coupling with rep(Mk=0)=1, feedback-off holds the drive at its
maximum (rep=1) at *every instant*, so **D_on(t) ≤ D_off(t) pointwise** after withdrawal → banked divisions and
proliferation half-life are always ≤ off; no transient overshoot (worst +4×10⁻⁵), delay/oscillation doesn't help,
and the R-point can't be rescued (E_on ≤ E_off everywhere). A mark that boosts CyclinD1 *degradation* is even worse
(more repressive). The ONLY couplings that create withdrawal memory are those whose *net* effect on CyclinD1 is
positive — i.e. topology C by construction (the agent's `sign(dD_ss/dMk)` classifier flips to +). Deep reason:
topology A closes a **negative** loop (setpoint regulation, never a latch); memory needs the **positive** loop only a
net-pro-proliferative sustain/double-negative mark provides.

## Topology C — mark SUSTAINS proliferation (double-negative: PRC2 silences differentiation genes)
- **This is the buffering MEMORY.** The mark adds a mitogen-INDEPENDENT pro-proliferative drive, so the loop
  proliferation → EZH2/writer → mark → (silences the brake) → proliferation is **positive feedback**. Above a
  threshold sustain gain it **latches**: once proliferating (mark ON), division **persists through full mitogen
  withdrawal and through transient dips** (Fig C, E) — a genuine memory.
- **It is BISTABLE, not graded.** The transition is a sharp bifurcation: below sₙ no memory, above sₙ a *permanent*
  latch (Fig D). This holds even with a non-cooperative (monostable) mark, because the *sustain loop itself* is the
  positive feedback — read-write cooperativity is not required. So the memory is all-or-none per allele/cell; graded
  population behaviour would come from cell-to-cell variation in the loop gain (population-of-models).
- **Reconciles the ChIP-seq via a two-input logic** (Fig F): the latched mark stays ON through mitogen withdrawal
  (→ division persists = the memory) but is released by an **active differentiation erasure** (a demethylase burst /
  EZH2 shutdown that is part of the differentiation program, not triggered by mitogen loss). Hence H3K27me3 is
  maximal during division, persists through mitogen fluctuations, and **declines with differentiation** — exactly the
  ChIP-seq — while still buffering division against mitogen withdrawal.

## Conditions for the buffering memory (topology C latch)
1. **Sustain-loop gain supercritical** — the mark's pro-proliferative drive × proliferation→writer coupling × mark
   response must exceed 1 (bistability). Below it: no memory; far above: a rigid latch.
2. **Slow mark** — the memory duration through a dip/withdrawal is set by the mark's persistence; the mark must be
   slower than the fluctuations it buffers (chromatin timescale, hours–days).
3. **Proliferation-coupled writer** — the mark is written in proportion to division (E2F/EZH2), so it is maximal at
   max division and self-reinforces the proliferative state (JP's new constraint is *required*, not just allowed).
4. **Release is a separate input** — differentiation actively erases the mark; mitogen withdrawal alone does not.

## Implications / what to test
- The current v44 model implements **only topology A** (EZH2 ⊣ CyclinD1). It therefore *should* show dose-response
  buffering (it does) but *should not* rescue division through withdrawal by the mark (it doesn't — feedback reduces
  post-withdrawal divisions). This is a clean, already-grounded prediction, not a deficiency to "fix" — unless the
  biology calls for the persistence memory.
- If GNP division genuinely **persists through mitogen withdrawal via H3K27me3**, the responsible circuit is
  topology C (PRC2 maintaining the proliferative state by silencing differentiation genes), and the model would need
  a topology-C module to capture it. Experimental discriminators: (a) ChIP-seq H3K27me3 on *differentiation-gene*
  loci (not just CyclinD1) should be high in dividing GNPs and fall with differentiation; (b) **EZH2i during mitogen
  withdrawal should COLLAPSE the persistence** if topology C is operating (the latch is lost), whereas under topology
  A EZH2i would if anything *raise* CyclinD1 (de-repression) — an opposite-sign prediction that cleanly separates them.

Files: `simulations/sim_prc2_buffering_memory.py`, `simulations/fig_prc2_buffering_memory.py`
(`fig_prc2_buffering_memory.pdf/png`).

## Addendum — the writer(EZH2)-vs-dilution(division) balance and the two STATES (JP, same night)
JP's refinement: in fast-dividing GNP (no transient G0), H3K27me3 is set by EZH2 **writing** vs replication
**dilution**, and BOTH scale with proliferation — high EZH2 (E2F-driven) + heavy per-cycle dilution. When division
stops, dilution falls away but EZH2 collapses, so the mark is not replaced → differentiation genes de-repress. This
is the concrete **writer side** of the topology-C sustain latch. Model + fig:
`sim_prc2_writer_dilution_states.py`, `fig_prc2_writer_dilution_states.py`
(`dMk/dt = kw·EZH2·(1−Mk) − (del + μ·p)·Mk`, EZH2 tracks proliferation, mark sustains proliferation).

Findings:
- **The two states are real and bistable** (Fig A/B/D). ON = fast-dividing, EZH2 high, mark maintained,
  differentiation genes silenced; OFF = arrested, EZH2→0, mark low (0.19), differentiation genes de-repressed. A
  mitogen up/down ramp shows **hysteresis** — the memory — and withdrawal from fast division **latches** proliferation
  (Fig C). This is the same topology-C latch, now with the EZH2/dilution mechanism filled in.
- **Why the mark is high DURING division (not diluted away):** writing must outpace dilution, `kw·kE > μ`. When it
  does, the mark is maintained (even maximal) through fast division; if dilution won, fast dividers would be mark-LOW
  (contradicting the ChIP-seq). So "lots of EZH2 + lots of dilution → net high mark" is a specific, satisfiable
  balance condition, and it is what makes H3K27me3 maximal at max division.
- **One correction to the dilution intuition:** stopping division does NOT preserve the mark via reduced dilution.
  When division stops EZH2 collapses (no writer), and the low-EZH2 effect dominates the low-dilution effect → the mark
  is lost anyway → differentiation genes fire (exactly JP's stated conclusion). In fact the division-coupled dilution
  during fast division slightly **weakens** the ON state (it competes with writing), making the latch a bit *less*
  sticky than pure turnover would (μ=0.9 relaxes toward p≈0.30 after withdrawal vs μ=0 latching at p≈0.51). So the
  dilution is a **cost** to the memory, not its source — the memory is the sustain loop (mark → differentiation-off →
  proliferation → EZH2 → mark). The population-average mark peaks when the proliferating FRACTION is highest (P7) and
  declines as cells drop to the OFF branch (P14), matching the ChIP-seq without any single cell's mark being graded.
