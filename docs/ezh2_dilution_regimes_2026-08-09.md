# EZH2 generation vs replicative dilution of H3K27me3: failing / sufficient / excess (2026-08-09)

JP: "how does the model behave if the amount of EZH2 being generated is either sufficient to compensate for
dilution, in excess, or failing to keep up?" Knob = `kTlEZ` (EZH2 translation = amount of EZH2 generated),
`ez_scale` × baked (baked default = 1.0). On the baked default. Scripts: `recon_ezh2_dilution_balance.py`,
`fig_v44_ezh2_dilution_regimes.py`.

## The three regimes (steady cycling, mean Mk)

| regime | ez_scale | EZH2 | Mk | rate-sensitivity* | CyclinD1 | character |
|---|---|---|---|---|---|---|
| **FAILING** | ~0.1 | 0.5–0.9 | 0.43–0.61 (eroded) | **0.11–0.13** | de-repressed (12–41) | dilution wins; mark low & strongly division-rate-dependent; weak repressor |
| **SUFFICIENT** | 0.25–0.5 | 1.2–3.0 | 0.7–0.88 (balanced) | 0.05–0.065 | governed (6–22) | writing ≈ dilution; mark tracks division rate |
| **EXCESS** (**baked = 1**) | ≥0.75 | ≥2.3 | 0.94–0.98 (saturated) | 0.01–0.02 | floored (~7.5 MB / ~2.2 GNP) | writing >> dilution; mark pinned near 1, insensitive to rate |

\* rate-sensitivity = Mk(slow, 4× CDKI) − Mk(fast, 1× CDKI): how much the mark responds to the cell slowing
down (the dilution/speed-limiter effect).

## Findings

1. **The baked model is in EXCESS** for both GNP and MB (Mk ≈ 0.94, rate-sensitivity ~0.02). EZH2 generation
   is more than enough to overwrite dilution, so the mark saturates near 1. (The EZH2 MB/GNP fold ≈ 2.0 is
   *preserved across ez_scale* — both scale together — so the fold does not distinguish the regime; the
   absolute Mk level does, and the calibration puts it high.)

2. **The mark is a division-rate sensor ONLY when the writer is not in excess.** rate-sensitivity falls
   0.11 (failing) → 0.06 (sufficient) → 0.02 (baked, excess) → 0.01 (deep excess). This is *why* replicative
   dilution does not feed the permanent-vs-transient fate at the baked setpoint (see
   `permanent_vs_transient_2026-08-09.md`): the mark is saturated, so dilution can only nudge it. **The
   "speed-limiter" / rate-sensor behaviour lives in the SUFFICIENT regime**, not where the model currently
   sits. Panel D: Mk rises with cycle period (slower → less dilution → higher mark) steeply in failing,
   moderately in sufficient, flat in excess.

3. **EZH2 generation sets the governor strength — with diminishing returns in excess.** More EZH2 → more
   mark → lower CyclinD1, but the CyclinD1 curve flattens once in excess (MB Cd 7.8→7.45 as ez 1→3). So the
   baked model sits *past the useful-repression knee*: adding EZH2 barely represses further, and removing a
   lot (toward sufficient) is needed before CyclinD1 climbs. Conversely, in the failing regime CyclinD1 is
   strongly de-repressed (Cd 12–41) — the mark has lost control of the drive.

4. **GNP tips into the regimes at lower EZH2 than MB.** At ez 0.1, GNP is FAILING (Mk 0.43) while MB is still
   SUFFICIENT (Mk 0.61) — MB's higher EZH2 (2× fold) keeps its mark higher at every scale. Both are in
   EXCESS at the baked setpoint.

## Interpretation / open question

The model is deliberately calibrated (EZH2 MB/GNP 2.05, H3K27me3 ChIP MB<GNP, EZH2i/vismo folds) into the
**excess** regime, where H3K27me3 is a saturated, largely static repressor and a weak rate-sensor. That is
internally consistent (it's why the governor is a threshold-setter and why dilution doesn't switch fate).

But JP's mechanistic intuition — H3K27me3 as a **replicative-dilution rate-sensor** (accumulates when slow,
lost when fast) — is a property of the **sufficient** regime (ez ~0.3–0.5, Mk ~0.7–0.85). If the biology
sits there, then: (a) the mark would genuinely track proliferation rate (a speed-limiter), (b) dilution
would gain leverage on the drive and could feed the graded slowdown / fate near the boundary, and (c) the
CyclinD1 governor would be more responsive to EZH2 dose (no diminishing returns). **The discriminator is the
absolute H3K27me3 occupancy at the CyclinD1/CDK6 loci** in cycling cells: near-saturated (broad, ~1) → the
model's excess regime; mid-range and division-rate-dependent → the sufficient/rate-sensor regime. This is
directly measurable (CUT&RUN quantitation + a proliferation-rate perturbation).

Nothing baked/changed — this is an exploration of the current model's regimes.
