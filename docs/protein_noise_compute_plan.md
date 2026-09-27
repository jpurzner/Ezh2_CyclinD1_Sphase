# Where to run the burst-imputation analysis: measured, not guessed

Question: server, AWS, or this machine? Benchmarked on 2026-09-27.

**Answer: do it locally. You do not need AWS, and you only need Frontenac for your own raw reads.**

---

## The thing that decides it: we never touch raw reads

The heavy step people assume — FASTQ → alignment → BAM — **is already done**. Both ENCODE and
ChIP-Atlas serve *processed* signal (bigWig) and *peak calls* (BED). The job is aggregating those
to gene level. No alignment, no read processing, no genome indices.

The end product is a **genes × features matrix**: ~20,000 × ~50 = a few MB. Model fitting on that
is seconds.

## Measured numbers

| | measured |
|---|---|
| this machine | **10 cores, 64 GB RAM, 1.2 TB free**; `pyBigWig` 0.3.26 already installed |
| bandwidth here | **4.3 MB/s** (100 MB in 23.3 s) — the one weak link |
| ChIP-Atlas bigWig, one experiment | 718 MB → **~2.8 min** to download |
| ChIP-Atlas peak BED (Q<1e-5) | **10.2 MB** (3.7 MB at Q<1e-10, 0.7 MB at Q<1e-20) |
| remote bigWig streaming | works (`Accept-Ranges: bytes`), but **336 ms/query** → ~112 min per file for 20k genes |
| ENCODE portal | **unreachable from here** (HTTP 000, 45 s timeout) |
| ENCODE S3 (`s3://encode-public`) | **publicly readable over plain HTTPS, no credentials, no AWS CLI** |

## Two phases, sized

**Phase 1 — peak-based features. ~400 MB, minutes.**
40 antigens × 10.2 MB BED = 400 MB, about 1.5 minutes at measured bandwidth. Build gene × feature
matrices from peak overlap (presence, summed score, fraction of the promoter/body covered). This
is essentially what ChIP-Atlas's own *Target Genes* tables do, and that already gave a real signal
(Suz12 vs burst frequency, marginal r = −0.23). Enough for a first cross-cell-type R².

**Phase 2 — continuous RPM, only if Phase 1 justifies it. ~28 GB, ~2 hours.**
40 bigWigs × 718 MB. Matches Ochiai's measurement (RPM over −1000/+100 and gene body). Still
trivially local: 28 GB against 1.2 TB free.

**Do not stream bigWigs remotely for this.** It works, but at 336 ms/query it is ~112 min per file
against ~2.8 min to download, and firing 20,000 range requests per file at a public academic
server 40 times over is antisocial. Download once, query locally.

## Why not AWS

ENCODE *is* on AWS Open Data and I confirmed `s3://encode-public` is publicly listable and
readable with no credentials. Computing in us-west-2 would give zero egress and near-zero latency.

It still is not worth it here:

* the one-time transfer is 0.4 GB (Phase 1) or 28 GB (Phase 2), not terabytes;
* the compute is a regression on a few-MB matrix, not a cluster job;
* setup, IAM and cost-tracking overhead exceeds the ~2 hours of download it would save.

AWS would become the right answer if the scope changed to many cell types, or to anything needing
the raw reads.

## Why not Frontenac — for *this* part

Frontenac is the right place for **your own GNP/MB ChIP-seq**: that is genuine FASTQ → alignment →
signal work, needs the cluster, and should live in a reproducible pipeline next to your other
sequencing. But the public training data is already processed, so pushing it through a cluster
queue adds scheduling latency and data-transfer steps for no gain.

**One caveat that could flip this:** local bandwidth measured **4.3 MB/s**, which is slow. If that
is this machine's connection rather than ChIP-Atlas throttling, Frontenac likely has far better
throughput and is worth using purely as a download host for Phase 2 — fetch there, aggregate
there, and bring back the few-MB matrix. Worth testing the same download on Frontenac before
committing to Phase 2 locally.

## Recommended split

| step | where | size |
|---|---|---|
| Phase 1 peak features, public cell types | **local** | 400 MB |
| Train + cross-cell-type validation (mESC ↔ GM12878) | **local** | seconds |
| Phase 2 continuous RPM, if justified | local, or Frontenac if its bandwidth is better | 28 GB |
| **Your own GNP/MB ChIP from raw reads** | **Frontenac** | real alignment work |
| Apply the fitted model to GNP/MB | **local** | trivial |

## Practical notes

* `pyBigWig` is installed; `deeptools`, `pybedtools` are not. For Phase 1, plain interval
  intersection in Python or `bedtools` is enough — no need for the full deepTools stack.
* Cache raw downloads outside the repo (`data/priors/raw/` is already gitignored) and commit only
  the derived gene × feature matrices, which are small and are the reproducible artefact.
* Filter experiments to **wild-type only** before aggregating. This mattered: 37% of ChIP-Atlas's
  333 mESC Suz12 experiments are KO-flagged, and removing them changed Ezh2's frequency effect
  from 0.96× to 0.51× (`protein_noise_bursting_chromatin.md` §2c).
* ENCODE's portal was down during this benchmark. Resolve accessions via the S3 bucket or retry
  later; do not assume the REST API is available.
