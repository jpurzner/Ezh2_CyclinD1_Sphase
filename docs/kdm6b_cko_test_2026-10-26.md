# Testing the lattice against Kdm6b cKO — the prediction fails

2026-10-26. SLURM 12395073, 12395146, 12395180, 12395191.
Data: `promoter_segmentation/tables/kdm6b_eraser.tsv` (built by
`scripts/kdm6b_eraser_test.py`, 2026-10-21) — Kdm6b cKO vs WT, 77,455 genes,
with promoter H3K27me3 (RPGC) and expression (log2CPM) in both genotypes.

I asserted that the lattice's Kdm6b-KO prediction could not be tested because
no data existed. That was wrong: this table is exactly the test.

## The prediction and the result

The lattice predicted Kdm6b KO → Ccnd1 **0.04** (a 25-fold repression), on the
grounds that losing the Gli1-recruited eraser lets me3 encroach on the
proximal promoter.

| gene | promoter me3 WT → cKO | me3 fold | expression fold |
|---|---|---|---|
| Ccnd1 | 6.31 → 13.70 | **2.17×** | **1.06** |
| Ccnd2 | 3.43 → 7.30 | **2.13×** | **1.95** |
| Cdk6 | 6.50 → 10.77 | **1.66×** | **1.72** |

**The chromatin arm is right.** Losing KDM6B raises promoter H3K27me3 roughly
two-fold at all three genes, so KDM6B *is* demethylating these promoters —
which, for these loci, runs against the De Santa generalisation I cited
yesterday.

**The transcriptional arm is refuted.** Predicted 0.04, measured 1.06 at
Ccnd1 — wrong by a factor of 25 — and wrong in *direction* at Ccnd2 and Cdk6,
which go UP while gaining me3.

## Genome-wide: the coupling is real but ~30x weaker than modelled

Across 14,089 expressed genes:

- `corr(Δ promoter me3, Δ expression)` = **−0.114** (Pearson, p = 9.3e−42);
  Spearman −0.093.
- The octile dose-response is monotone in the repressive direction, but the
  entire range spans 0.08 log2 — **6% of expression** end to end.
- Genes that more than doubled promoter me3 (n = 141): median expression
  change **−0.063 log2 = 0.96 fold**, i.e. a 4% decrease. Among them 6.4% fall
  >1.5× against 0.7% rising — a 9:1 enrichment, so the direction is genuine.

So doubling promoter H3K27me3 costs about 4% of expression. The lattice's Hill
coupling (h = 3, Km = 0.10) takes the promoter from fully on to essentially
off across a comparable range.

## What this says about the model

**The lattice made the mark the repressor.** Pol II loading is a Hill function
of promoter me3, so me3 level *is* repression. The data says promoter me3 is
only weakly repressive — the direction holds, the magnitude does not.

**The v44 two-compartment model already had this right.** In
`two_compartment_prc2_2026-08-09.md` the repressor is `P_prox`, PRC2 *complex*
occupancy at the TSS, while the mark is the distal reservoir that loads it.
That structure predicts exactly what is seen here: remove the eraser, the mark
rises, but transcription barely moves, because the mark is not what represses.

The fix is structural, not a parameter change: **repression should be a
function of PRC2 occupancy at the promoter, with the mark acting to recruit and
retain PRC2 rather than to block Pol II directly.** That also preserves the
reservoir result, since the reservoir's job becomes loading P_prox faster.

## What this does and does not overturn

| result | status |
|---|---|
| Kdm6b KO = Gli1 loss (0.04) | **refuted** — measured 1.06 |
| Bursting vs encroachment | **suspect** — same Hill coupling drives it |
| v44 perturbation scorecard | **suspect** — EZH2i/HHi folds come from the same coupling |
| Reservoir gives faster off-kinetics | **probably survives** — it is about mark kinetics, not the repression law, but needs re-deriving |
| me3 domain geometry, nucleation at neighbours | **unaffected** — chromatin-side results |
| The monostability proof | **unaffected** — it is about the writing/loss balance |

## A second, independent test available

`West_ChIPseq/cultured_CGN_CUTRUN/` holds **GSK-J4 vs DMSO H3K27me3 CUT&RUN,
3 replicates each** (SRR21377411-16) — a pharmacological KDM6 inhibition in
cultured granule neurons. Still raw fastq, so it would need processing, and
CGN are post-mitotic rather than GNP. It would test the chromatin arm a second
way but not the transcriptional one.
