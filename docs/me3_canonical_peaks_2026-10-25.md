# Canonical H3K27me3 peak sets (mouse)

2026-10-25. Decision record. SLURM 12382617 / 12382626 / 12382627.

## The sets to use

All called with `epic2 --bin-size 100 --gaps-allowed 3 -egf 0.8881
--false-discovery-rate-cutoff 0.01` at **native depth**, mm39, in
`~/me3_recall/`:

| file | source | reads | Mb | genome % | intervals |
|---|---|---|---|---|---|
| `GNP_H3K27me3_canonical_mm39.bed` | Scott GNP rep1 | 96.4M | 138.8 | **5.337** | 50,709 |
| `MB_H3K27me3_canonical_mm39.bed` | Scott MB rep1 | 76.7M | 186.7 | **7.181** | 32,999 |
| `MB_H3K27me3_consensus_mm39.bed` | Scott MB rep1 n rep2 | — | 89.8 | 3.454 | 33,612 |

The same rule for both cell types: **rep1 at native depth**. GNP lands at
5.34% and MB at 7.18%, both in the plausible range for a broad mark against
the project's ~5% figure. The MB rep1-n-rep2 consensus is provided as a
conservative alternative but is depth-limited by rep2 (14.7M), not by biology.

## What was dropped, and why

**The four EL libraries are excluded from peak-based analysis**: GNP1 (4.2M),
Tumor1 (4.6M), Tumor3 (3.8M), Tumor4 (7.0M). A depth titration on Scott GNP
rep1 showed H3K27me3 coverage still climbing steeply at 96M reads with no
plateau, so at 4-7M epic2 recovers well under a third of the domain landscape.
These remain usable as coverage/signal tracks; they are not peak sources.

Also superseded:

- `EL_GNP1_K27me3_b100_g6.bed` (epic2_all_test) — 140,501 peaks at median
  1,004 bp. Does not match any BAM we can identify: re-calling
  `EL_ChIPseq/output_mm39/GNP1_H3K27me3_1` at the same nominal settings gives
  5,610 peaks at median 10,016 bp. **Provenance unknown; do not use.** This is
  the set behind the first version of the Ccnd1/Ccnd2/Cdk6 figure.
- `GNP_H3K27me3_rep1_fdr*_epic2peak.bed` (Scott chip_quant_mm39) — 1,195
  peaks, 0.11% genome. Re-running the command recorded in its own
  `peak_calling.sh` gives 39,256 peaks and 9.00%, a 78x difference from
  identical inputs. Almost certainly a failed run.

## Replicate agreement

| | rep1 | rep2 | shared | % of rep1 in rep2 | % of rep2 in rep1 |
|---|---|---|---|---|---|
| GNP | 138.7 Mb | 41.0 Mb | 24.3 Mb | **17.5%** | 59.2% |
| MB | 186.7 Mb | 129.2 Mb | 89.8 Mb | 48.1% | 69.5% |

GNP agreement looks poor but is a depth artefact: rep2 is 8.8M against rep1's
96.4M, so it recovers a subset (59% of rep2 is inside rep1) rather than
disagreeing. MB, where the two reps are closer in depth, agrees better. Neither
is evidence against rep1.

## Ccnd1 / Ccnd2 / Cdk6, canonical set

Fraction of each gene's H3K27me3 within +/- 500 kb that a window captures:

| gene | kb of me3 in 500 kb | +/-20 kb | +/-50 kb | +/-100 kb | +/-250 kb |
|---|---|---|---|---|---|
| Ccnd1 | 143.0 | **4.3%** | 13.4% | 33.5% | 71.6% |
| Ccnd2 | 151.0 | **5.8%** | 19.0% | 27.1% | 56.9% |
| Cdk6 | 44.1 | **27.7%** | 46.2% | 72.3% | 92.1% |

Essentially unchanged from the uncorrected gap-3 call (4.3 / 5.4 / 26.9%), as
expected — the egf fix is a 1.6% effect. The conclusion stands and is if
anything firmer: **the lattice window at TSS +/- 20 kb sees 4-6% of the
H3K27me3 landscape at the two cyclin D genes.**

## What still needs doing

- The lattice dataset is still built at TSS +/- 20 kb. Rebuilding at +/- 100
  or +/- 250 kb is a `build_lattice_dataset.py <binbp> <half>` re-run.
- Human (GBM, ENCODE) peak sets have not been touched; they need their own
  egf (hg19/hg38) and their own depth assessment.
- Anything downstream that consumed either superseded file should be re-checked.
