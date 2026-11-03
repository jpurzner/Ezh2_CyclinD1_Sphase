# GSK-J4 CUT&RUN: acute KDM6 inhibition does not raise promoter H3K27me3

2026-10-26. SLURM 12395476 (alignment), 12395494 (quantification), 12395559,
12395590. Data: `West_ChIPseq/cultured_CGN_CUTRUN/`, processed from raw fastq.

## What was processed

12 libraries aligned to mm39 (bowtie2, `--end-to-end --very-sensitive
--no-mixed --no-discordant -I 10 -X 700`, MAPQ ≥ 20, duplicates removed):

- **Kdm6bi arm** — DMSO ×3 vs **GSK-J4** ×3 (KDM6A/B inhibitor)
- **Ezh2i arm** — DMSO ×3 vs **GSK126** ×3 (EZH2 inhibitor)

Alignment 95.6–97.4%, 14.2–23.7M deduplicated pairs per library. Output in
`~/lattice/cutrun/`.

## GSK126 is a working positive control

| gene | DMSO → GSK126 (promoter CPM) | fold |
|---|---|---|
| Ccnd1 | 133.6 → 65.4 | **0.49** |
| Ccnd2 | 91.7 → 37.1 | **0.41** |
| Cdk6 | 98.5 → 17.9 | **0.19** |

Global: fraction of reads in the canonical GNP me3 domains falls from 0.307
to 0.268, a ~13% reduction.

So the assay, the pipeline and the promoter windows can all detect a 2–5×
change in promoter H3K27me3. That matters for interpreting the other arm.

## GSK-J4 does nothing

| gene | DMSO → GSK-J4 (promoter CPM) | fold |
|---|---|---|
| Ccnd1 | 111.6 → 112.8 | 1.01 |
| Ccnd2 | 96.1 → 80.0 | 0.83 |
| Cdk6 | 84.2 → 81.5 | 0.97 |

Genome-wide over 64,072 promoters: median log2 **+0.014**, with 14.4% up
>1.5× against 14.8% down >1.5× — balanced, i.e. noise. Global domain
occupancy is unchanged (0.302 → 0.303).

## The two KDM6B perturbations disagree

| | Ccnd1 | Ccnd2 | Cdk6 |
|---|---|---|---|
| Kdm6b cKO (genetic) | **2.17×** | **2.13×** | **1.66×** |
| GSK-J4 (acute) | 1.01× | 0.83× | 0.97× |

Candidate reconciliations, in the order I would rank them:

1. **Timescale.** H3K27me3 turns over slowly (t½ ~158 h). An acute drug
   treatment may simply be too short for the mark to accumulate, whereas a
   genetic cKO has had the whole developmental window. Checkable from the
   treatment duration in the source paper.
2. **Cell-cycle context.** These are cultured, post-mitotic granule neurons;
   the cKO table is from a proliferative context. The project's own audit
   argued that replication is what carves the promoter core, so without
   division there may be no route for me3 to redistribute on any timescale.
3. **Target engagement.** GSK-J4 is a prodrug (GSK-J1 is active) and hits
   KDM6A as well as KDM6B; engagement was not verified here.

## A normalisation caveat that matters for one arm but not the other

CPM is the wrong normalisation when a perturbation changes the mark globally.
GSK126's genome-wide median of **+0.46** is an artefact: reads lost from true
me3 domains are redistributed to background, so most promoters appear to rise
while the real targets fall. **The GSK126 per-gene folds are reliable in
direction but their magnitudes are compressed**, and the right global readout
is the in-domain fraction above.

The GSK-J4 null is **not** affected, because nothing changed globally either —
the in-domain fraction is flat, so the relative and absolute measures agree.
That null is solid.

## Where this leaves the KDM6B arm of the model

| arm of the prediction | cKO | GSK-J4 |
|---|---|---|
| me3 rises at the promoter | supported (1.7–2.2×) | **not supported** |
| transcription falls | **refuted** (1.06, 1.95, 1.72) | not testable here |

So each arm of the lattice's KDM6B mechanism is contradicted by at least one
dataset. Combined with the abundance problem — KDM6B undetectable in the
Itzhak proteome, and the model needing 10–30 molecules per domain against a
bound of 0.02–0.32 — the arm should be demoted exactly as
`locked_parameters_2026-10-23.md` recommended.

## Still unused

The same directory holds a **DIV1 / DIV3 / DIV5 / DIV7 differentiation
timecourse**, 3 replicates each, already downloaded. That is a direct
measurement of how H3K27me3 redistributes as granule neurons exit the cycle
and differentiate — which is the process the lattice is actually trying to
describe, and a better target than either drug arm.
