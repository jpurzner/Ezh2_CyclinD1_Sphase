# Memory function on the current model — mitogen ramp + fluctuation (2026-07-20)

*JP: re-explore the memory function with mitogen ramp up/down and fluctuation, on the current default (Cdh1-EZH2 +
graded G1-lengthening + Skp2). Files: `sim_memory_ramp_fluctuation.py`, `fig_memory_ramp_fluctuation.py`.
Single-cell deterministic, full v44 model, piecewise Shh(t).*

## Two memories, distinguishable
The model has **two** memory mechanisms and this exploration separates them:

1. **Chromatin low-pass buffer (the H3K27me3 mark).** Under a fluctuating Shh (period 14–80 h), CyclinD1 and Skp2
   track the drive (CV 0.57–0.73 and 0.85–1.6), EZH2 is attenuated ~3× (CV 0.20–0.27), and **the mark is nearly flat
   (CV 0.002–0.07)** — it filters out the mitogen noise almost completely, and more so at fast periods (a genuine
   low-pass filter). So the hierarchy of buffering is Skp2 < CyclinD1 < EZH2 < **mark (deepest, most stable memory)**.
   The mark also creates a small CyclinD1 **dose-response hysteresis** under ramps (loop area +0.09 with the mark,
   **0.00 with the mark's CyclinD1 repression off** — the loop *is* the mark).

2. **Cell-cycle commitment hysteresis (Skp2-p27-Rb-E2F R-point + graded G1-lengthening).** Under a mitogen ramp, the
   **proliferation rate shows a clear hysteresis loop: on the down-ramp division persists to LOWER Shh than it
   started on the up-ramp** (Fig C) — a bistable-commitment memory. This is a *cell-cycle* memory, not chromatin; it
   is the same R-point bistability that carries a committed cell through a sub-cycle dip.

## Reconciles the earlier reduced-model result
Last week's reduced-model finding was that a CyclinD1-*repressing* mark gives dose-response buffering but not
withdrawal-persistence memory. The full model is consistent: the **persistence memory (proliferation hysteresis) is
the cell-cycle commitment**, while the **mark's role is the low-pass buffer + a small dose-response hysteresis**. The
graded G1-lengthening deepens the buffering (subthreshold dips lengthen the cycle, slowing EZH2's S-phase supply).

## Notable, mechanism-level
- **Skp2 is now the FASTEST tracker** (highest CV) because of the new mitogen-dose Skp2 induction — it follows Shh
  even more tightly than CyclinD1, and shows its own ramp hysteresis (Fig F). A concrete, testable consequence of the
  Skp2 change.
- **Rate-dependent sign of the CyclinD1 loop**: fast ramp → CyclinD1 propped up on the down-ramp (mark relief lags);
  slow ramp → CyclinD1 lower on the down-ramp (mark had time to accumulate and represses more). The mark's
  accumulation-vs-relief race is visible in the loop sign.

## Caveats / next
- Single-cell deterministic; the population version (heterogeneous cd_scale / birth-p27) would smear the sharp
  commitment staircase into a graded population hysteresis — worth running for the realistic readout.
- Ramp panels are cycle-averaged (Shh-binned) to see through the cell-cycle oscillation.
- Open question this sets up: an EZH2i-during-ramp run should REMOVE the CyclinD1 dose-response hysteresis (mark
  loop) while LEAVING the proliferation-commitment hysteresis — a clean in-silico separation of the two memories,
  and a wet-lab discriminator.
