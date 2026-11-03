# H3K27me3 vs transcription, re-cut with the chiptools thresholds

2026-10-24. Re-run of the 47-cell core-clearance analysis on calls re-derived
with `chip_segment` (winsorise → `trim_upper_tail` → two-component GMM) for
H3K27me3 and `soften_zero` → `binarize_counts` for expression, replacing the
dataset's shipped binary layers and the arbitrary `log2tpm > 0` expression
split used in the first pass.

46 of 47 samples (one skipped: too few marked genes). Script `recut3.py`,
submitted as SLURM job 12378891; output `recut3.json`.

## Headline: one conclusion strengthens, one does not survive

| | first pass | re-cut (continuous) | re-cut (binary calls) |
|---|---|---|---|
| threshold, expression percentile | 80.3 | 81.5 | 83.1 |
| …its CV across samples | 0.17 | **0.07** | **0.07** |
| threshold, absolute log2tpm | 1.88 | 2.10 | 2.32 |
| …its CV across samples | 0.62 | 0.44 | 0.45 |
| GBM transition width | 33.3 | 13.0 | 11.2 |
| normal transition width | 10.2 | 12.4 | 11.0 |
| GBM / normal | 3.25× | **1.05×** | **1.01×** |
| Mann–Whitney p | 0.0003 | **0.56** | **0.63** |

**Survives, and tightens.** The clearance threshold is universal as a *rank*
and cell-type specific in absolute units. Re-cutting halves the rank CV
(0.17 → 0.07) while the absolute-unit CV stays six times larger (0.44). The
first pass had this right and the re-cut makes it sharper.

**Does not survive.** The 4× broader GBM transition was an artefact of the
expression split. It is now 1.05×, p = 0.56. I reported that result; it is
withdrawn.

## Why the width result was an artefact

The first pass called a gene expressed if `log2tpm > 0`, i.e. any nonzero
count. The re-derived cutoff sits at a median of **0.191** (range
0.115–0.794) and calls a median of **41.6%** of genes expressed.

Everything between the two — low-but-nonzero genes — was being spread across
the bottom deciles of the first pass. Bulk GBM tumours have more of that than
a sorted or cultured normal population does, so their deciles were padded at
the bottom with genes that are not really on, which stretched the measured
transition. Per-cell the correction is large and one-directional: AK089
49.0 → 7.2, AK139 56.2 → 13.5, AK071 38.6 → 10.6, while the threshold
position for those same samples barely moves (80.3 → 80.7, 89.6 → 86.6,
79.4 → 79.5).

So the expression axis, not the me3 axis, was carrying the error. That also
explains why the first pass showed a *wider* CV on the threshold: the same
padding was displacing the crossing point by a variable amount per sample.

## The binarisation was not the problem

Re-deriving the me3 calls barely changes anything, which is worth saying
plainly because it was the motivation for the whole exercise:

- threshold percentile, continuous vs binary core: r = 0.931, median
  difference 0.9 percentile points
- transition width: r = 0.806, median difference 1.7 percentile points

The continuous core-mean was a faithful proxy for the bound fraction all
along. What the re-cut bought was not a better me3 axis — it was discovering
that the *expression* axis was wrong.

## Raw cutoffs: why `soften_zero` was necessary

Every single sample's raw `log2tpm` cutoff came out at **0.000**. The zero
spike (33–47% of genes, depending on sample) owns the lower mixture
component outright, so the fitted boundary collapses onto zero and the GMM
reports a "cutoff" that is really just the spike. Softening the spike first
moves it to 0.115–0.794 and makes the mixture identifiable.

On the me3 side the shipped index cutoffs were also unusable — GM12878 shipped
−0.78, re-derived 1.14 — consistent with the earlier QC finding that
`me3_bin` was saturated (core plateau 1.000 in 11/20 ENCODE mouse) and the
`state` layer's Active fraction ranged 0.001–0.577.

## Baseline differences between cohorts, unresolved

The binary calls expose something the continuous mean hid. Bound fraction of
core sites at the silent end of the expression range:

| group | n | silent end | top decile | continuous amp |
|---|---|---|---|---|
| GNP | 1 | 0.351 | 0.101 | 1.166 |
| ENCODE mouse | 19 | 0.785 | 0.699 | 0.227 |
| ENCODE human | 6 | 0.709 | 0.558 | 0.189 |
| GBM | 20 | 0.316 | 0.205 | 0.845 |

ENCODE samples call ~78% of core sites bound even at silent genes and only
drop to ~70% at the top decile — a dynamic range of 9 percentage points.
GNP and GBM start near 0.32 and clear to 0.10–0.21. Either ENCODE's
H3K27me3 is genuinely much flatter, or its signal-to-noise compresses the
range and the GMM cutoff lands low. This is the same ENCODE-amplitude
question raised in the first pass and the re-cut has not settled it; the
threshold *position* is nonetheless consistent across all four groups, which
argues the compression is in amplitude rather than in the location of the
transition.

ENCODE human also sits lower on threshold percentile (72.2) than mouse
(81.7), on n = 6.

## What to carry forward

1. The ~80th-percentile clearance threshold is the robust result. Quote it as
   a rank, not as log2tpm ≈ 3 — the absolute value varies 6× more.
2. Do not quote the GBM broadening. It is gone.
3. Use `soften_zero` before any mixture fit on log2tpm in this dataset.
4. The ENCODE amplitude compression is still unexplained and is the thing most
   likely to be a data-quality artefact rather than biology.
