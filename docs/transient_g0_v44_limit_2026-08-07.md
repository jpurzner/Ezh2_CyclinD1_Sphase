# #3 Periodic high-p27 transient G0 — v44 deterministic core cannot do it (2026-08-07)

**JP direction:** "we need some form of transient G0 with high p27 periodically in cells" (MB-specific; GNP has none).
**Finding:** the v44 deterministic core produces a **clean commitment switch, not a transient G0** — G0-fraction is
robustly **0%** across every lever I tested. Transient G0 needs a *dedicated* mechanism (a population/toggle layer),
not a parameter tweak. This matches the stochastic-commitment finding, `cdki-baked-g0-g1-behavior`, the two-quiescence
roadmap (Module A = population layer), and the existing v45 model (built for exactly this).

## What was tested (all give MB G0-frac = 0, cell either cycles cleanly or arrests permanently)
1. **Birth-p27 (`P21_div`) in the default (option A):** raising P21_DIV_MB 1.4→2.1 does *nothing* — MB identical
   (6 div, 25.9h, 30/32). p27 is buffer-dominant and inert (CdP21 active), so birth-p27 is sponged up.
2. **Birth-p27 under option B (p27-inhibitory, 31/32 params):** raising 0.6→2.5 *still* does nothing — MB identical
   (6 div, 25.9h, 31/32). Reason: **p27 is cleared far too fast** (trajectory: born-p27 2.5 → 0.007 within ~30h)
   for the pulse to hold a G0 dwell; the cell commits before p27 can arrest it.
3. **Slower p27 clearance (kDeKPC 0.02→0.002) + high birth-p27 + option B:** G0-frac *still* 0 (DMSO phase
   fractions 0/74/19/7 throughout). Slowing clearance doesn't open a dwell.
4. **Single high-p27 birth (probe):** below a threshold → cycles cleanly (no delay); above → **permanent arrest**
   (no re-entry, mass runaway) — a bistable OFF lock, never a transient dwell. Thresholds GNP ~1.0, MB ~2.0.

## Mechanistic reason
- The commitment is a **sharp bistable switch** (two-step Rb + Skp2-p27 double-negative). A newborn cell either
  crosses the R-point quickly (commits — pre-S is brief, pRb rises fast → G0-frac ≈ 0) or falls into the bistable OFF
  state (**permanent** arrest: high p27 → low CDK2 → no Skp2 → p27 stays high; CyclinD1 alone only mono-phosphorylates
  Rb and cannot bootstrap the CyclinE-CDK2 loop, so it never re-enters). There is **no intermediate "reversible dwell
  then re-enter" regime** reachable by birth-p27 or p27 stability.
- A **transient** G0 requires holding cells in a *reversible* pRb-hypophospho state and then **releasing** them —
  i.e., a bistable toggle with a **mitogen/CyclinD1-gated escape** (CyclinD1 accumulates during G0 until it flips the
  toggle → re-entry) plus a **growth cap** (no mass runaway during the dwell). The v44 core has none of these.

## The path (a design fork for JP — this is real mechanism work, not a tweak)
Transient G0 is a **population-layer** phenomenon (roadmap Module A). Three ways to deliver it:
1. **Build a transient-G0 module into v44:** a self-sustaining p27-CDK2 (or Rb-E2F) toggle that HOLDS a reversible G0,
   + a mitogen/CyclinD1-gated re-entry, + a growth cap. This is what makes the dwell transient and MB-specific
   (MB's high CyclinD1/CDK6 drive — now two-armed via #1 — is exactly the escape signal; GNP with low drive would
   lock or never enter). Substantial build (≈ what the v45 thread was).
2. **Use the existing v45 model** (`v45-stochastic-commitment`, 20/20): already has the bistable CDK2-p27 toggle +
   mitogen-gated bootstrap + mass-cap — built for transient commitment. Adopt it as the transient-G0 layer.
3. **Reduced population model** (`simulations/stoch_commit_g0/`): a fast vectorised ensemble with birth-p27
   heterogeneity + a mitogen-gated escape → demonstrate the MB-specific transient-G0 fraction (cycling in/out) and
   GNP ~0%, as the Module-A population layer on top of the v44 core.

**MB-specificity comes for free once there is a re-entry mechanism:** re-entry is mitogen/CyclinD1-driven, and MB's
drive is now two-armed (CyclinD1 + CDK6, #1) — so MB born high-p27 can escape (transient G0), while a GNP with the
same high birth-p27 (lower drive) would lock (no re-entry) or, at its normal low birth-p27, never enter. That is the
clean "MB transient G0, GNP none" story — but it needs the escape route, which the v44 core lacks.

## Status
Nothing baked. #1 (CDK6) and #2 (graded GNP G1) are done at 31/32; #3 is blocked on a design decision about *how* to
add the transient-G0 layer (build into v44 / adopt v45 / reduced population model). Recommend deciding the approach
before building, since it's the largest of the three pieces.

Probes: `simulations/birthp27_probe.py` + the P21_DIV_MB / kDeKPC harness sweeps (this session).
