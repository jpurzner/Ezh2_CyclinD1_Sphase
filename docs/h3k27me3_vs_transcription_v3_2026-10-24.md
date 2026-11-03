# 47-cell threshold analysis, rewritten-cutoff edition

2026-10-24. Third pass, SLURM job 12379010, 46 samples. me3 cutoffs now from
`chip_segment` with `background_scale` + `background_null_cutoff` — no
winsorisation, no two-component fit. Expression unchanged (`soften_zero` →
`binarize_counts`). Output `recut4.json`.

## Headline: the conclusions are unchanged, and that is the point

| | v1 original | v2 winsorise + GMM | **v3 background + null** |
|---|---|---|---|
| threshold, expression percentile | 80.3 | 81.5 | **81.5** |
| …CV | 0.17 | 0.07 | **0.07** |
| threshold, absolute log2tpm | 1.88 | 2.10 | **2.10** |
| …CV | 0.62 | 0.44 | **0.44** |
| binary threshold, percentile | — | 82.9 | **82.1** |
| …CV (log2tpm) | — | 0.45 | **0.39** |
| GBM / normal transition width | 3.25× | 1.05× | **1.05×** |
| …p | 0.0003 | 0.56 | **0.56** |

The continuous numbers are *identical* to v2, which is the expected result and
a useful check: the continuous analysis never touches a me3 cutoff, so a
change to the cutoff method must leave it alone. It does.

The GBM broadening stays retracted, now at 0.90× / p = 0.87 on the binary
calls — if anything more clearly null than before.

## What the new cutoffs actually changed

The bound cutoff roughly doubled, from a median of 0.843 to **1.685**, and
the per-site positive fraction fell to a median of **13.8%** (range
4.1–51.8%). The old winsorise+GMM cutoff was calling ~36% of all core sites
H3K27me3-positive, which was never plausible for a mark that covers a minority
of promoters.

The tail collapsed in **9 of 46** samples — in those the knee the trimmer
finds is the marked population itself, with no separate shelf beyond it, so
there is nothing to strip. Where a tail does survive it is a median 2.74% of
sites.

## Two anomalies I previously flagged largely resolve

**ENCODE human was an outlier on threshold position.** On the binary calls
the gap to the other three cohorts closes almost completely:

| | ENCODE human | everyone else | gap |
|---|---|---|---|
| v2 binary | 70.3 | 83.2 | 12.9 |
| **v3 binary** | **80.7** | **82.1** | **1.4** |

So that outlier was mostly a cutoff artefact, not biology. Note the
*continuous* gap is unchanged at 9.4 — it cannot move, for the same reason the
headline numbers cannot. Whatever remains in the continuous signal is either
real or a property of the ENCODE continuous scale, and this analysis does not
separate those.

**The ENCODE amplitude compression shrinks.** Bound fraction of core sites at
the silent end of the expression range:

| group | v2 | v3 |
|---|---|---|
| GNP | 0.404 | 0.284 |
| ENCODE mouse | 0.832 | **0.581** |
| ENCODE human | 0.721 | **0.358** |
| GBM | 0.439 | 0.263 |

The "ENCODE calls ~78% of core sites bound even at silent genes" figure I
reported becomes 58% (mouse) and 36% (human). ENCODE still sits above GNP and
GBM, so the compression is reduced rather than explained away, but a large
part of it was the winsorised cutoff.

## Caveat on the validation statistic

Within this run, restricted to `me3_marked` genes, `corr(marked fraction,
per-site bound fraction)` is −0.189. That is not the +0.216 reported from the
rewrite test, because that test used *all* genes while this analysis is
conditioned on the marked subset — conditioning on being marked removes most
of the variation the correlation was measuring. The two numbers are not
comparable; the all-gene figure is the one that speaks to the winsorisation
bug.

Separately, `me3_marked` is itself derived from the dataset's shipped GMM
cutoffs, which are the ones shown to be unreliable, so it is a weak reference
in either direction.

## What to carry forward

1. The ~80th-percentile clearance threshold stands across all three cutoff
   methods. Quote it as a rank, not as log2tpm — the absolute value still
   varies 6× more across samples.
2. The GBM broadening remains withdrawn.
3. Use the v3 cutoffs for anything that needs binary calls: positive
   fractions are plausible, the ENCODE threshold outlier is gone, and the
   cutoffs are exactly scale-equivariant.
4. Still unresolved: the residual ENCODE elevation in continuous core me3,
   and the fact that 9 samples have no separable tail.
