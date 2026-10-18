# Cdh1-gated EZH2 degradation — restoring the S-elongation → EZH2 accumulation (2026-07-19)

*JP (from the regenerated figures): "elongation of S does not increase EZH2 abundance — it should, we have
experimentally verified this." This is the fix. JP-approved; baked into the v44 default.*

## The problem
The model showed **no S-elongation → EZH2 accumulation**: the `HU EZH2-in-S boost` (EZH2 in HU-arrested S ÷ EZH2
in control S) sat at **1.009** (target 1.31), and the cycle-phase tracking was flat too (EZH2 protein G2/G0 = 1.04
vs 1.48; EZH2 G0/cycling = 0.44 vs 0.6). All three "passed" only on wide bands. Experimentally, arresting/elongating
S raises EZH2 per cell.

## Why (diagnosed, not guessed)
Three candidate explanations ruled out by direct test:
1. **Not synthesis rate.** Raising `kEZE2f` 6× moved the boost only 1.009→1.031 (and broke other targets).
2. **Not a cyclin collapse.** In control-S vs HU-arrested-S the EZH2 drivers are *identical* (E2f 0.98/0.99,
   CyclinA 1.64/1.71, Cd ~10.4) — the arrested cell reaches the *same* S-phase state.
3. **Root cause: the concentration convention + a single turnover rate.** `with_ezh2_conc=True` preserved EZH2
   across division, so EZH2 was a steady pool already near its S-phase level — arrest could not push it higher.
   For a longer S to accumulate *more* EZH2, division must **dilute** it (so S-*duration* sets how much rebuilds).

## The fix (biologically grounded)
**EZH2 is a bona fide APC/C-Cdh1 substrate** — degraded in G0/G1, stabilized in S/G2. The model already had Cdh1
(the `C1` species, used to degrade Skp2). Three coupled changes:
1. **Dilution convention** (`with_ezh2_conc=False`, default flipped): EZH2 = EZH2/2 at mitosis → dividing cells
   reset it, arrested cells (no division) bank it. This is what makes elongated/arrested S accumulate EZH2.
2. **Low baseline `kDeEZ`** (0.00117 → 0.000281): committed/S/arrest cells hold EZH2 STABLE so it accumulates
   over the long arrest.
3. **Cdh1(C1)-gated degradation** `EZH2_deg: EZH2 => ; Cell*(kDeEZ + kDeEZ_C1*C1)*EZH2` (new `kDeEZ_C1`): clears
   EZH2 in G0 (C1 active) so G0/cycling stays ~0.6 (not driven to ~0). This resolves the tension a single rate
   cannot: **accumulate in committed/HU-arrested S (C1 off) yet clear in G0 (C1 on)**.

A joint re-optimization (`simulations/optimize_cdh1_ezh2.py`, 19 knobs, halving convention) re-balanced the coupled
EZH2/PRC2/transcription knobs so all other targets hold.

## Result (pure default, no overrides)

| Metric | Old default | Cdh1-EZH2 default | Target |
|---|---|---|---|
| HU EZH2-in-S boost | 1.009 | **1.24** | 1.31 |
| EZH2 protein G2/G0 | 1.04 | **1.86** | 1.48 |
| EZH2 G0/cycling | 0.44 | **0.62** | 0.6 |
| EZH2 MB/GNP | 1.64 | 1.95 | 2.05 |
| Vismodegib (MB+HHi/MB) | 0.139 | 0.128 | 0.144 |
| CyclinD1 fold | 4.87 | 4.24 | 5.07 |
| Score | 27/28 | **27/28** | — |
| Lone fail | HU-G2 0.356 | HU-G2 0.352 | (same structural near-miss) |

**Direct evidence it's genuine accumulation** (not a classifier artifact), same drivers, EZH2 protein:
control-S 5.43 → HU-arrested-S **6.50** (ratio 1.20). Before the fix it was 2.37 → 2.39 (flat).

## What was baked
`src/build_model_v44_heldt.py`: builder default `with_ezh2_conc=False`; `EZH2_deg` reaction gains `+ kDeEZ_C1*C1`
(default 0 in the string, set by `_ts_bake`); `_ts_bake` updated with the re-fit config incl. `kDeEZ`, `kDeEZ_C1`,
`K_E2f_EZ`. `simulations/validate_v44.py`: `EZH2_CONC` env default flipped to `0`. Cell-cycle machinery unchanged.
Eraser still OFF; mark half-life ~2.8h preserved. HU-G2 stays the lone structural near-miss (see
`v44-hu-s-g2-structural-limit`). Config: `simulations/optimize_cdh1_ezh2_best.json`. Diagnostic:
`simulations/_diag_ezh2_hu.py [fix]`.
