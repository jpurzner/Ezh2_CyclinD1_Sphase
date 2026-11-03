# Depth-matched re-call of all mouse H3K27me3 libraries

2026-10-24. SLURM 12380557 (downsample), 12380558 (call), 12380606 (titration).
All 8 mouse H3K27me3 libraries downsampled to a common 3.80M reads and
re-called with epic2: `--bin-size 100 --gaps-allowed 3 -egf 0.8881 -fdr 0.01`.
Output in `~/me3_recall/`.

## The libraries

| sample | native ChIP reads | input | input reads |
|---|---|---|---|
| Scott GNP rep1 | 96,428,049 | GNP_input_0513 | 180,337,266 |
| Scott MB rep1 | 76,700,720 | MB_input_0513 | 48,338,218 |
| Scott MB rep2 | 14,672,452 | MB_input_0513 | 48,338,218 |
| Scott GNP rep2 | 8,808,239 | GNP_input_0513 | 82,574,349 |
| EL Tumor4 | 7,042,101 | Tumor4_Input | 6,188,681 |
| EL Tumor1 | 4,564,987 | Tumor1_Input_1 | 4,948,958 |
| EL GNP1 | 4,221,855 | GNP1_Input_1 | 4,088,592 |
| EL Tumor3 | 3,822,030 | Tumor3_Input | 1,704,863 |

A 25-fold depth range. Inputs were left at native depth — downsampling the
input was separately shown to change coverage by under 0.1 percentage points,
and a deeper input is a better background estimate, not a biased one.

## Result: matched depth works, but 3.8M is the wrong depth to match at

All eight hit the target within 0.05%:

| sample | n peaks | median bp | Mb | genome % |
|---|---|---|---|---|
| Scott GNP rep1 | 6,378 | 2,199 | 17.1 | **0.657** |
| Scott GNP rep2 | 13,383 | 1,599 | 34.4 | 1.322 |
| Scott MB rep1 | 15,575 | 2,499 | 59.1 | 2.275 |
| Scott MB rep2 | 26,010 | 1,899 | 68.5 | 2.636 |
| EL GNP1 | 8,840 | 2,499 | 36.6 | 1.408 |
| EL Tumor1 | 14,062 | 2,099 | 40.9 | 1.572 |
| EL Tumor3 | 14,367 | 2,199 | 42.0 | 1.617 |
| EL Tumor4 | 19,862 | 2,399 | 76.8 | **2.954** |
| *Scott GNP rep1, full 96M* | *50,709* | *2,299* | *138.8* | ***5.337*** |

**What matching bought.** Median peak width is now consistent at 1.6-2.5 kb
across every library. Before, the stored calls ranged 99-699 bp and were not
comparable at all. Residual spread in coverage is 0.66-2.95%, a 4.5x range
that is now genuine between-library variation rather than a depth artefact.

**What it cost.** Every library lands below the ~5% expected for a broad mark.
Scott GNP rep1 falls from 5.34% at its native 96M to **0.657%** at 3.8M — an
8-fold loss — and becomes the *worst* library in the set despite being the
best one at full depth. At 3.8M the calls are dominated by sampling noise.

## Why: H3K27me3 detection has not saturated even at 96M reads

Titrating Scott GNP rep1 down from its full depth:

| depth | lambda | k | n peaks | genome % | change |
|---|---|---|---|---|---|
| 3.8M | 0.157 | 1 | 6,378 | 0.657 | |
| 8M | 0.330 | **2** | 9,275 | 0.526 | **-0.131** |
| 16M | 0.660 | 2 | 15,702 | 1.487 | +0.961 |
| 24M | 0.991 | **3** | 22,565 | 1.456 | **-0.031** |
| 32M | 1.321 | 3 | 25,820 | 2.318 | +0.862 |
| 48M | 1.981 | **4** | 34,663 | 2.892 | +0.574 |
| 64M | 2.641 | **5** | 43,368 | 3.331 | +0.439 |
| 96M | 3.962 | **7** | 50,764 | 5.290 | +1.959 |

Two things to take from this.

**There is no plateau.** Coverage is still climbing steeply at 96M. Whatever
depth is "enough" for H3K27me3 in this system is beyond anything we have, so
depth is not a solved variable — it is the dominant one.

**The dips are the integer threshold, exactly as predicted.** Both intervals
where coverage *falls* with more sequencing — 3.8M to 8M and 16M to 24M — are
precisely the intervals where the eligibility threshold k steps up (1 to 2, and
2 to 3). Every interval where k holds constant shows a gain. This is the step
function from `src/statistics.pyx` showing up directly in the data: adding
reads raises lambda, which raises the integer number of reads a bin needs to
count at all, and the stricter bar can outweigh the extra signal.

## Revised recommendation

Downsampling to the shallowest library is the wrong move here, and this run
demonstrates it rather than assuming it.

1. **Call the deep Scott libraries at native depth** with `-egf 0.8881`,
   bin 100, gap 3. Those calls are already on disk
   (`scott_GNP_rep1_FULLDEPTH_b100_g3_egf.bed`, 5.337%).
2. **Do not use the four EL libraries (3.8-7.0M) for peak-based analysis.**
   At that depth epic2 recovers under a third of the domain landscape. Use
   them as coverage/signal tracks only.
3. **If cross-sample peak comparison is required**, match within the deep set
   only — Scott GNP rep1 (96.4M) and MB rep1 (76.7M) share a usable common
   depth of 76.7M.
4. **Avoid depths near a k boundary.** With bin 100 the steps fall at roughly
   6M, 20M, 29M, 43M and 62M reads; landing just above one costs coverage
   relative to just below it.
5. Always pass `-egf 0.8881`; without it epic2 silently uses hg19's 0.9020.
