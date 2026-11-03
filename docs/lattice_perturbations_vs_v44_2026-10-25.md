# Perturbing Gli1, EZH2 and KDM6B — lattice vs the v44 targets

2026-10-25. SLURM 12386415, 12386468, 12386556. Ccnd1, 50 cells per condition,
400 h, scored as fold change in Pol II initiation rate — which is what the
measured CyclinD1 ratios are.

## Scorecard

| test | target (tol) | v44 model | lattice A | lattice B |
|---|---|---|---|---|
| EZH2i CyclinD1 fold (GNP) | 2.20 (±42%) | **1.63** | 1.41 | **2.03** |
| Ezh2 cKO CyclinD1 | 2.46 (±30%) | — | 1.36 | **2.06** |
| CyclinD1 GNP+HHi / GNP | 0.157 (±40%) | — | **0.15** | **0.14** |
| Ezh2-OE → CyclinD1 down | direction | — | 0.57 (2×), 0.28 (3×) | — |

Baseline A = KDM6B 30, PRC2 24 (0.784 init/h, 71% of the de-repressed
ceiling). Baseline B = KDM6B 10, PRC2 24 (0.528 init/h, 48% of ceiling).

**Baseline B passes all three.** It also clears the EZH2i target that the
existing v44 model misses (2.03 against 1.63, target 2.20).

## The constraint this exposes

The two targets pull in opposite directions and pin the baseline:

- the ceiling (no PRC2) is ~1.1 init/h, so EZH2i = 2.2 requires the baseline
  at 1.1/2.2 = **0.50 init/h**;
- the HHi floor is ~0.08 init/h, so HHi = 0.157 requires the baseline at
  0.08/0.157 = **0.51 init/h**.

Those agree, which is why baseline B works and baseline A (0.78) cannot: at
A the gene is already at 71% of maximum and EZH2i has only 1.4× of headroom
however completely PRC2 is removed. **The measured EZH2i fold is therefore a
direct measurement of how repressed Ccnd1 is at baseline in GNP** — about
half-maximal, not near-maximal.

That is a stronger statement than the v44 shortfall alone supported, and it
reframes that shortfall: 1.63 is not a missing mechanism so much as a baseline
that sits too far up the dose-response.

## What is fit and what is predicted

**Fit.** Baseline B was chosen after baseline A failed, and the GNP+HHi
condition uses a 25% EZH2 reduction — the low end of the measured 25–47%
range for GNP+HHi. Both are selections within measured bounds, but they are
selections, so the three passing scores are a calibration, not a prediction.
Across the measured EZH2-drop range the HHi fold runs 0.14 (25%) → 0.22 (35%)
→ 0.39 (47%), so only the low end fits.

**Predicted.** Nothing was tuned to these:

| | lattice | note |
|---|---|---|
| Ezh2-OE 2× | 0.57 | correct direction (Fig 3E, qualitative) |
| Ezh2-OE 3× | 0.28 | |
| **Kdm6b KO** | **0.04** | no measured target |
| **Kdm6b OE 2×** | **1.09** | no measured target |

## The sharpest prediction

**Losing KDM6B is equivalent to losing Gli1.** Kdm6b KO gives 0.04 and Gli1
withdrawal alone gives 0.04 — the same number. In the model Gli1's entire
effect on Ccnd1 runs through KDM6B recruitment, so a Kdm6b knockout should
phenocopy Hedgehog inhibition at this locus.

**But KDM6B overexpression does almost nothing** (1.09). The eraser is already
saturating at baseline, so the system is sensitive to losing it and insensitive
to adding more. That asymmetry is testable and is the cleanest discriminating
prediction the model currently makes.

## Caveats

- Ccnd1 only, one parameter set for the Pol II arm.
- Initiation rate is compared to measured mRNA ratios; equal mRNA half-life
  across conditions is assumed.
- `Ezh2_OE` at baseline B did not finish, so the direction test is from
  baseline A only.
- EZH2 dose-response is smooth and monotone throughout (0.5 → 72 molecules),
  consistent with the monostability established earlier — there is no switch
  anywhere in this panel.
