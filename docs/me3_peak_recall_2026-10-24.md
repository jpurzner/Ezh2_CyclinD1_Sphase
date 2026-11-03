# Re-calling GNP H3K27me3 with epic2: a gap sweep

2026-10-24. SLURM array 12379549, 18 runs. mm39, effective genome 2.6e9.
epic2 `--bin-size 100 --false-discovery-rate-cutoff 0.01`, sweeping
`--gaps-allowed`.

## 1. The stored Scott peak file is wrong, not just badly parameterised

Re-running the **exact command recorded in `GNP_H3K27me3_rep1_peak_calling.sh`**
— same BAM, same input, bin 100, gaps 6, FDR 0.01 — gives:

| | n peaks | mean bp | total Mb | genome % |
|---|---|---|---|---|
| stored `GNP_H3K27me3_rep1_fdr2_epic2peak.bed` | 1,195 | 2,474 | 3.0 | **0.11** |
| **re-run, identical parameters** | **39,256** | **5,960** | **234.0** | **9.00** |

**78× more coverage from the same command.** The stored file is the product of
a failed or truncated run. No amount of parameter tuning explains it, and any
analysis that used it was working from essentially no peaks — consistent with
it containing zero peaks within 250 kb of Ccnd1, Ccnd2 or Cdk6 while the
signal at those loci is obvious.

## 2. The input-depth imbalance was not the cause

Scott's input is 180.3M reads against a 96.4M ChIP (1.87x deeper), which
looked like a plausible way to inflate epic2's background. Downsampling the
input to 96M makes essentially no difference:

| gap | full input, genome % | input downsampled to 96M, genome % |
|---|---|---|
| 1 | 3.364 | 3.340 |
| 3 | 5.620 | 5.584 |
| 6 | 9.000 | 8.927 |
| 10 | 12.952 | 12.908 |
| 20 | 18.652 | 18.590 |

Hypothesis tested and rejected — epic2 is already scaling by library size.

## 3. The gap sweep

**Scott rep1** (96.4M reads):

| gap | gap bp | n peaks | mean bp | Mb | genome % |
|---|---|---|---|---|---|
| 1 | 100 | 70,670 | 1,237 | 87.5 | 3.36 |
| **3** | **300** | **55,117** | **2,651** | **146.1** | **5.62** |
| 6 | 600 | 39,256 | 5,960 | 234.0 | 9.00 |
| 10 | 1 kb | 23,322 | 14,439 | 336.8 | 12.95 |
| 20 | 2 kb | 6,224 | 77,917 | 485.0 | 18.65 |

**EL GNP1** (4.2M reads):

| gap | n peaks | mean bp | Mb | genome % |
|---|---|---|---|---|
| 1 | 13,432 | 1,867 | 25.1 | 0.97 |
| 3 | 8,341 | 4,633 | 38.6 | 1.49 |
| 6 | 5,610 | 10,016 | 56.2 | 2.16 |
| 10 | 3,921 | 18,940 | 74.3 | 2.86 |
| 20 | 1,552 | 63,067 | 97.9 | 3.77 |

Coverage rises monotonically with the gap and there is no plateau in this
range, so the gap cannot be chosen by looking for one — it has to be set
against an expected coverage. Against the project's own ~5% figure for a broad
mark, **Scott gap 3 (5.62%, mean domain 2.65 kb)** is the operating point, with
gap 6 (9.00%, 5.96 kb) as the permissive alternative. Gap 20 is clearly
over-merged: 78 kb mean domains over 18.7% of the genome.

EL never reaches 5% at any gap. At 4.2M reads it is too shallow to support a
broad-domain call; its apparent adequacy earlier came from a stored file that
does not match this BAM (see below).

## 4. The stored EL file does not match its own BAM either

| | n peaks | mean bp | Mb | genome % |
|---|---|---|---|---|
| stored `EL_GNP1_K27me3_b100_g6.bed` | 140,501 | 1,004 | 141.0 | 5.42 |
| re-run of `GNP1_H3K27me3_1` at b100 g6 | 5,610 | 10,016 | 56.2 | 2.16 |

Same nominal parameters, 25x more peaks at a tenth the mean width. So
`EL_GNP1_K27me3_b100_g6.bed` in `epic2_all_test` was not produced from
`EL_ChIPseq/output_mm39/GNP1_H3K27me3_1`, or not with these settings. **This is
the file I used for the wide panels in the Ccnd1/Ccnd2/Cdk6 figure**, and I
cannot currently account for its provenance.

## What to do

1. Discard the stored Scott me3 peak files and re-call. The re-run at gap 3 is
   already on disk at `~/me3_sweep/scott_full_b100_g3.bed`.
2. Establish what `EL_GNP1_K27me3_b100_g6.bed` actually is before using it
   again.
3. Re-check any downstream result that consumed either file.
4. gap=50 runs were still going when this was written; they will only extend
   the over-merging trend.
