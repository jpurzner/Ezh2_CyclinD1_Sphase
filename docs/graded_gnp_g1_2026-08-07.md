# Variable GNP G1 by mitogen dose (#2, 2026-08-07)

**JP direction:** "we need variable G1 duration in GNPs based on mitogen dose." Deliver a graded GNP G1 dose-response
(lower Shh → longer G1), keeping it simple and holding validation. Nothing baked — presented as a param set for review.

## Mechanism (already in the model, was latent)
`mu_eff := mu/(1+k_mu_cki·INK4) · (mu_min_frac + (1-mu_min_frac)·Cd^n_g1len/(K_g1len^n_g1len + Cd^n_g1len))`
As CyclinD1 (mitogen) falls, growth slows → the **size-gated** cell takes longer to reach S-entry → **G1 lengthens
while S/G2 stay fixed** (Fan-Meyer shape: density lengthens G1, S/G2 flat). Latent at baseline because K_g1len=0.6 sat
far below GNP's Cd (~2.9), so the slowdown never engaged at validation doses.

## What it took (two coupled things, not one)
K_g1len alone lengthens the baseline but does **not** widen the graded window — the window is capped by the **arrest
threshold** (baseline model arrests at SHH≈0.42, leaving no dose range to grade over). So graded-G1-by-dose needs:
1. **a wider GNP cycling range** (mild commitment-threshold drop) to create the dose range, AND
2. **K_g1len engaged** to lengthen G1 across it,
3. **while restoring the SHH=0 arrest** (the wider range erodes Shh-dependence — the same tension as option B) via
   the **low basal Cd** lever (`k_Cd_tx_basal`), and restoring the periods via `mu`/`k_mu_cki`.

A single-lever attempt (candidate A: kPhRbCd 0.22 + K_g1len only) gave a beautiful wide curve (G1 14→50h) but broke
the SHH=0 arrest (cycled at SHH=0) and slowed the baseline to 24h → 29/32. The joint re-balance fixes both.

## Result (`g1_optimize.py`, ~140-candidate search with a graded-G1 objective so the re-balance can't disengage it)
**31/32** (only the structural MB-HU-G2 fails; all arrests + Shh-dependence intact). Param set (`g1_optimize_best.json`):
`mu 0.000856, k_mu_cki 0.603, kPhRbCd 0.159, K_CdRb 0.412, k_Cd_tx_basal 0.000464, K_g1len 1.557, mu_min_frac 0.221,
n_g1len 3`. CyclinD1 MB/GNP 4.81.

Graded GNP G1 dose-response (`g1_dose_sweep.py`):
| SHH | period | G1 | S | G2M | ndiv |
|---|---|---|---|---|---|
| 0.50 | 20.7 | 11.3 | 5.9 | 3.5 | 16 |
| 0.42 | 22.5 | 12.7 | 6.3 | 3.5 | 14 |
| 0.35 | 25.3 | 15.3 | 6.4 | 3.5 | 14 |
| 0.30 | 28.4 | 18.1 | 6.9 | 3.4 | 12 |
| 0.25 | 29.0 | 26.6 | (exit transition) | | 3 |
| ≤0.20 | arrest | | | | |

**G1 lengthens 11→18h over the robust cycling range (SHH 0.5→0.3, S/G2 flat), extends to ~27h at the exit
transition, then the cell exits below SHH≈0.22.** No discrete G0 — graded G1 lengthening, exactly the biology.

## Caveat / trade-off
- **Baseline period is ~20.4h at SHH=0.5** (vs the ~16h target) — it passes the validation band but has drifted up,
  because K_g1len is partially engaged already at SHH=0.5 (Cd≈2.4 is not far above K_g1len≈1.56). This is inherent:
  `baseline 16h` and `a wide graded window` can't both hold unless Cd ranges much higher at saturating Shh than the
  grading threshold. If SHH=0.5 is treated as a mid dose (not saturating), the 16h would sit at higher Shh — worth
  deciding how SHH=0.5 maps to "full mitogen."
- The graded window (SHH 0.5→0.3) is narrower than the single-lever candidate A (which reached 50h) because the
  re-balance had to keep the SHH=0 arrest; the exit now happens at SHH≈0.22 rather than ≈0.05. That's the right
  trade-off — a real dose-response that still exits at low mitogen.

## Status / next
Nothing baked (these are re-tuned values of existing baked params). Baking = updating `_ts_bake` with the param set
above — JP's call. Combines naturally with #1 (CDK6 arm) and #3 (periodic high-p27 transient G0). If we promote this
+ CDK6 + option B, the shared prerequisite remains a partial-Shh dose-response validation target (this graded-G1
curve is exactly that target).

Files: `simulations/g1_dose_sweep.py`, `g1_optimize.py`, `g1_optimize_best.json`, `g1_optimize_log.jsonl`.
