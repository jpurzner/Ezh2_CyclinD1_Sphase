# Porting the chiptools thresholding to Python, and what it fixes
2026-10-24. Code: `src/lattice/chipnorm.py`. Source of truth:
[github.com/jpurzner/chip_tools](https://github.com/jpurzner/chip_tools)
(`chiptools`) plus the newer `chipseq_norm` R package.

Ported because the lattice dataset is npz/python and the Frontenac node has no R.

---

## 1. Why the shipped binary layers could not be used

- `me3_bin` is **saturated**: core plateau = 1.000 in 11/20 ENCODE mouse
  tissues. The index shows why — `cut_me3` runs as low as **−2.035**.
- `state` is **not comparable across samples**: genome-wide Active fraction
  spans **0.001 to 0.577**, because `cut_ac` ≈ 0.005 and `cut_k4 = inf` in
  several samples. Where ac/k4 failed, everything me3⁺ collapses to "Polycomb".

## 2. The canonical pipeline, and why every constant is what it is

Reading the repo resolved something that had looked arbitrary. The defaults only
make sense **after winsorise-and-scale**, which is the step I had been skipping:

```
wins_norm_histones()            winsorise by POINT COUNT (500 low, 20 high), scale 0-1
  -> trim_upper_tail()          rank/value knee; threshold 3e-4 is defined ON [0,1]
  -> chip_binarize_tail_removed()   binarize_counts() on the remainder,
                                    exclude_extreme = 0.1 < x < 1, mu = c(0.5, max-0.2)
  -> chip_categorize_and_collapse()  unbound / bound / tail
```

So `3e-4`, `exclude_extreme`, and `mu = c(0.3, 0.65)` are all [0,1]-scale
constants. Applied to log-ratio or rlog data without the winsorise step they are
meaningless — which is the likely origin of the bad shipped cutoffs.

## 3. Validation against the real R

Run locally (Rscript + mixtools 2.0.0.1) against the identical input vector:

| function | python | R | verdict |
|---|---|---|---|
| `binarize_counts` cutoff_50 | 0.9563 | 0.9563 | **exact** |
| cutoff_75 | 0.7269 | 0.7269 | **exact** |
| mu_lower / mu_upper | −0.0296 / 2.6713 | −0.0296 / 2.6713 | **exact** |
| `wins_norm` (min/max/mean) | 0.0000/1.0000/0.1532 | same | **max abs diff 6.6e−05** |
| `trim_upper_tail` on [0,1] | 0.3939, 0.3926 | 0.3946, 0.4020 | agree ~2% |
| `trim_upper_tail` raw scale, `rescale=False` | **NA** | **NA** | **reproduces R exactly** |

The ~2% residual on the tail trimmer is loess degree-2 (R) versus degree-1
(port, evaluated on a rank grid so it runs on 10^6+ values).

Two bugs found while porting, both mine, both worth recording because a naive
port hits them:
- the slope must be **per unit rank**; evaluating on a grid without dividing by
  the grid spacing inflates it by the grid step.
- walking down the rank axis the spacing is **negative**; dividing by the signed
  step flips every slope and the crossing fires at rank 1.

## 4. soften_zero is required for the expression channel

`chiptools::soften_zero` thins a dominant zero bin to at most 1.5x the tallest
other bin. log2tpm is **~35% exact zeros**, and without it the mixture puts an
entire component on the spike:

| | cutoff_50 |
|---|---|
| synthetic zero-inflated, k=2 raw | **0.000** |
| same, after `soften_zero` | **2.982** |

Confirmed on the real data — every sample so far:

| cell | shipped cut_me3 | NEW me3 cut | tpm cut RAW | tpm cut SOFTENED |
|---|---|---|---|---|
| GM12878 | **−0.78** | **1.139** | 0.000 | 0.126 |
| H1 | 1.30 | 1.371 | 0.000 | 0.164 |
| AK199 | 3.67 | 3.200 | 0.000 | 0.158 |
| AK213 | 5.41 | 3.104 | 0.000 | 0.159 |

**Every raw expression cutoff collapses to 0.000.** That is the pathology, and
it means my earlier "expressed = log2tpm > 0" split was the same arbitrary
choice by another route. The softened cutoffs (0.13–0.22) call ~37–43% of genes
expressed.

The me3 corrections are large where the shipped cutoff was implausible
(GM12878 −0.78 → 1.139) and modest where it was already sensible.

## 5. What is in the port

| function | source | notes |
|---|---|---|
| `binarize_counts` | chipseq_norm | hand-written 1-D 2-component EM, no sklearn, so identical code runs locally and on the cluster |
| `trim_upper_tail` | chipseq_norm | + `rescale` so the 3e-4 default transfers to arbitrary scales |
| `soften_zero` | chip_tools | |
| `split_counts` | chip_tools | BIC over k = 2..4 |
| `wins_norm` | chip_tools | winsorise by point count + unit scale |
| `chip_segment` | chip_tools | the full broad-mark chain, returning unbound/bound/tail on the original scale |

`binarize_counts` differs between the two R packages: chip_tools uses
`arbvar = FALSE` with hard-coded `mu = c(0.3, 0.65)`; chipseq_norm uses
`arbvar = TRUE` with quantile-based inits and states in its own comments that
the hard-coded values "only worked on the old [0,1] scale". The port follows
chipseq_norm.

## 6. Still to do

- Re-run the 47-cell threshold analysis on the **re-derived** calls. Everything
  in `h3k27me3_vs_transcription_2026-10-23.md` used continuous me3 and an
  arbitrary expression split, so the *shape* results should survive, but the
  expression axis should be re-cut at the softened cutoff.
- Apply `chip_segment` (tail trimmer included) rather than bare
  `binarize_counts` — H3K27me3 is a broad mark and the tail is the point.
