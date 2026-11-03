# Which substrate should the H3K27me3 cutoff be derived from?

2026-10-24. Three candidate substrates, 46 samples, SLURM jobs 12378968 (per-gene
vs per-site) and local analysis of the peak-intersect-TSS tables.

## The three substrates

| substrate | values per sample | available for |
|---|---|---|
| per-site binned bigWig | 1.6–4M (genes × 200 sites) | all 47 |
| per-gene summary of the same bigWig | ~58k (1 per gene) | all 47 |
| **peak ∩ TSS read counts** | ~18k (1 per gene) | **6 ENCODE human + 9 West only** |

The project's own QC doc is explicit that the third is the intended one —
*"Gate 2 — quantification substrate (peak ∩ TSS, NOT raw TSS windows)"* — but
it was only ever generated for ENCODE human and the West cerebellum timecourse.
The GBM tumours, the 20 ENCODE mouse tissues and GNP have no peak ∩ TSS table,
so the lattice dataset is built on raw TSS windows throughout.

## Per-gene vs per-site, all 46 samples

| variant | n/sample | median separation | median bound fraction |
|---|---|---|---|
| per-site, me3_marked genes only | 1.75M | 1.99 | 36.1% |
| per-site, all genes | 4M | 2.64 | 18.2% |
| **per-gene mean, all genes** | 58k | 2.23 | 22.1% |
| per-gene mean, marked only | 8.8k | 2.17 | 27.8% |
| per-gene core mean, all genes | 58k | 2.15 | 19.5% |
| per-gene max, all genes | 58k | 2.27 | 19.2% |

The substrate matters a lot: ENCODE mouse bound fraction goes from 46.1%
(per-site) to 19.6% (per-gene). That reduction is in the right direction — the
78%-bound-at-silent-genes figure I reported for ENCODE was partly a per-site
artefact.

## Correction: the bimodality statistic was misleading

I first scored each fit by the depth of the dip between the two component
means, and reported that per-gene recovers a valley in 60% of samples against
19% for per-site. **The histograms do not support that reading.** Liver scores
a valley of 0.339 on a distribution that is visibly a single broad bell; the
statistic is picking up a shallow numerical dip between two wide overlapping
components. Treat `substrate_comparison.png` as the evidence, not the number.

What the histograms actually show:

- Per-gene distributions are better *shaped* for a mixture — sharply peaked
  near the low end with a long right tail, so the cutoff lands on a real
  shoulder (GNP, AK100, AK231, H1).
- Per-site distributions of the same samples are broader and flatter.
- Neither is convincingly bimodal in most samples.
- ENCODE mouse cerebellum and ENCODE human liver / GM12878 are near-symmetric
  broad bells on the winsorised scale, per-gene and per-site alike. That shape
  is what low signal-to-noise looks like, not what a domain mark looks like.

## The peak ∩ TSS substrate

For the 6 ENCODE human lines (`peaktss_encode.png`), on log2(x+1):

| | GM12878 | H1 | HepG2 | K562 | liver | SK-N-SH |
|---|---|---|---|---|---|---|
| % genes with no overlapping peak | 10.3 | 10.7 | 1.9 | 2.3 | 9.3 | 2.6 |
| bound cutoff | 8.99 | 11.84 | 10.35 | 8.88 | 8.40 | 7.72 |
| separation | 1.99 | 4.52 | 3.67 | 3.02 | 1.96 | 0.52 |
| bound fraction | 28.6% | 13.5% | 16.0% | 23.5% | 31.7% | 51.1% |

**H1 is the one genuinely bimodal sample in this whole exercise** — two clear
modes with the cutoff sitting in the valley between them. That is biologically
coherent: ES cells have a sharp bivalent-promoter class, differentiated cells a
continuum. The other five remain unimodal even on the peak-aware substrate.

Note that only 2–11% of genes lack an overlapping peak, so "peak ∩ TSS" is far
less selective than the name suggests — H3K27me3 peak sets are broad enough to
touch almost every promoter. SK-N-SH's separation of 0.52 means its fit is
effectively degenerate.

## The problem that outranks substrate choice

Across all 46 samples the per-gene bound fraction **anti-correlates** with the
dataset's own `me3_marked` flag (r = −0.415; per-site is worse at −0.593; GBM
alone is −0.84). The cohort ratios swing from 0.37 (GNP) to 2.69 (ENCODE mouse).

A gene-level cutoff and `me3_marked` are answering the same question, so they
should agree. The likely cause is the per-sample winsorise-to-[0,1] step inside
`chip_segment`: it normalises away exactly the cross-sample differences in
H3K27me3 abundance, leaving a shape statistic rather than an abundance one. A
sample where many genes are marked has its marked population spread across the
rescaled range, and the GMM places the cutoff higher relative to it.

**So these per-sample cutoffs should not be compared across samples as they
stand.** That caveat applies to the per-site version I used in the re-cut too.

## Recommendation

1. For deciding *which genes are me3-marked*, use the per-gene summary, not
   per-site — it is the chiptools idiom and the distributions suit it.
2. For the lattice metagene, keep a *per-site* call: a gene-level cutoff
   applied site-by-site is a category error, since the two answer different
   questions ("is this gene marked" vs "is this site occupied").
3. Do not compare bound fractions across samples until the winsorisation is
   replaced by a shared scale.
4. The headline results are unaffected either way: the ~80th-percentile
   clearance threshold and the absence of GBM broadening both come from
   *continuous* core me3 against an expression ranking, and never touch these
   cutoffs.
5. Generating peak ∩ TSS tables for GBM, ENCODE mouse and GNP would mean peak
   calling on 41 samples. Worth doing if gene-level calls become load-bearing;
   not needed for the continuous analysis.
