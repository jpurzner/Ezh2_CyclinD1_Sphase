# Data requests

What I need to take the noise work past prediction. Ordered by how much it unblocks, not by
how hard it is to get.

**Short version:** only **two papers** are genuinely blocked (Eden 2011, McShane 2016) — everything
else in the literature I can fetch myself, see §C. The real bottleneck is **§A, your own data**.
Items A1–A3 are blocking: without them the model can produce noise but cannot be tested or even
put on an absolute scale.

---

## A. From the Purzner lab — nobody else has these

### A1. The HCI single-cell tables ⬅ **the single biggest gap**

Everything in `protein_noise_data_exploration.md` is a *prediction* of per-protein CV. There is
currently nothing to test it against. This is also the brief's own T8, where the loader is a
placeholder ("data format to be supplied").

**Format:** one row per cell, tidy CSV/parquet. Any column naming is fine — I'll map it.

| column | why it is needed |
|---|---|
| **mean intensity** of the stained protein | this is the **concentration** — the quantity the ODE species represent |
| **integrated intensity** of the same | the **amount**; the two differ by cell size and the difference is itself informative |
| **DNA content** (integrated DAPI/Hoechst) | cycle staging — CV must be computed *within* a phase, not across |
| **EdU / BrdU** if stained | resolves S-phase properly; DNA content alone cannot separate early-S from late-G1 |
| **nuclear area** (and cell area if segmented) | size regression — the brief's §2.3 point, and the Rb-dilution term |
| **condition** | GNP / MB, ± vismodegib, ± EZH2i, ± palbociclib, ± HU |
| **plate / well / replicate / batch** | batch is a variance component and has to be separable from biology |
| **protein stained** | one file per antibody is fine |

**Proteins, in priority order:** CyclinD1, EZH2, H3K27me3, p21/p27, Gli1, MYCN, phospho-Rb,
Ki67. The first three are the model's core; CyclinD1 and EZH2 are where the predictions are
sharpest (predicted CV ≈ 0.87 and ≈ 0.89 respectively).

**Raw, not summarised.** Per-cell rows, not means ± SD. The whole question is the shape of the
distribution, so means destroy the signal.

### A2. Something to set the absolute scale (Ω)

`docs/DATA_PROVENANCE.md` says the model's species are in **arbitrary units**, calibrated to
ratios. Discrete bursts (T6) and Gillespie (T2) both need molecule **counts**. Any *one* of these
is enough to anchor the whole model:

- GNP nuclear volume, or cell volume (µm³) — even a rough mean from your imaging;
- an absolute copy number for any single species, by any method;
- a quantitative western against a recombinant standard for one protein.

Cheapest is probably nuclear volume from the HCI segmentation you already have. Without this, the
Poisson term `1/⟨p⟩` and the burst discretisation have no defined magnitude and T6 cannot be
parameterised at all — only made self-consistent.

### A3. Measured cell-cycle length for GNP and MB

µ = ln2/T_cc scales the entire dilution analysis, and the extrinsic correlation time is quoted in
**generations**, so T_cc converts it to hours. The model produces ~22–23 h; I'd rather this were
an input than an output. Live count, cumulative EdU, or dilution of a label — whatever you have.

### A4. mRNA half-lives in GNP — `Ccnd1` and `Ezh2`

Only if such an experiment already exists (actinomycin-D chase, SLAM-seq, TT-seq). These are the
two specific numbers that are currently doing the most damage:

- ***Ccnd1*:** Schwanhäusser says 14.0 h, the actinomycin-D literature says <1 h, and after the
  `/60` fix the model implies 52 min. This 12× spread moves CyclinD1's noise correlation time from
  1.7 h to 20 h — the difference between "fluctuates several times per G1" and "frozen for longer
  than a cycle". Opposite predictions for G0 escape.
- ***Ezh2*:** the model needs ~1.2 h to reproduce the EZH2 phase gradients it is calibrated on;
  Schwanhäusser measures 9.85 h. Those cannot both be true of the same system (`protein_noise_fixes.md`
  §F3). It also moves EZH2's CV between 0.42 and 0.89.

### A5. The scRNA-seq you already have (Fig 4A,B)

`DATA_PROVENANCE.md` lists scRNA-seq (Sepp, Vladoiu, Ocasio, Luo) used for EZH2 transcript by
phase. Even with the zero-inflation caveat that motivated this whole exercise, the **per-cell count
distributions** for our genes would give:

- ⟨m⟩ in the actual system, rather than borrowed from NIH3T3 or fibroblasts;
- a Fano factor / CV² per gene, which constrains `(1+B)/⟨m⟩` directly;
- coverage for `Gli1`, `Ptch1`, `Smo`, `Mycn` — **absent from all three external datasets**.

A gene × cell count matrix (raw counts, not normalised) plus the phase assignment per cell is what
I'd want. This is not a substitute for burst inference, but it beats borrowing from fibroblasts.

### A6. A technical-noise floor

To compare a measured CV to a predicted one I have to subtract the measurement noise. Whatever you
have: secondary-only/isotype control wells, technical replicate stains of the same sample, or
beads. Without this, any excess CV is uninterpretable — it could be biology or the microscope.

---

## B. Paywalled — all now supplied by the user; kept for provenance

### B1. Eden et al. 2011, *Science* 331:764 — **Supporting Online Material, Table S1**
doi:10.1126/science.1199784 · "Proteome half-life dynamics in living human cells"

Bleach-chase protein half-lives, ~100 proteins, human H1299. science.org returns 403 on the
article, the `doi/suppl/` listing and every `downloadSupplement` pattern; not open access, no PMC
record, and the 2012 *Nat Protoc* follow-up is closed too.

**Why:** an *arbiter*, not a per-gene resource. It measures half-lives in **living** cells by a
completely different route from SILAC pulse-chase, so it is the cleanest independent check on
whether Schwanhäusser's half-lives run systematically long — the question behind the *Ccnd1* 14 h
vs <1 h conflict. It is small, human and a cancer line, so it will not fill the Hedgehog gaps.

### B2. McShane et al. 2016, *Cell* 167:803 — **supplementary tables**
PMID 27720452 · "Kinetic analysis of protein stability reveals age-dependent degradation"

No PMC record; cell.com returns 403.

**Why:** ~10% of proteins are degraded **non-exponentially** — newly made or unassembled copies are
cleared fast, assembled copies are stable. v44 carries a lot of complexes (`CeP21`, `CaP21`,
`RbE2f`, `RbmE2f`, `E1C1`, `iPcna`, `iRc`) and currently degrades bound and free pools at the
**same rate**. If any of our species are in their non-exponential class, that is a structural
error in the model, not a parameter tweak. Ideally the table listing which proteins show
age-dependent degradation, with the fast/slow rate pair.

*If either is easier as a PDF than as a table, that's fine — I can pull numbers out of a PDF.*

### ~~B3. Battich et al. 2020~~ — **supplied 2026-09-26**, see `protein_noise_cellcycle_rates.md`

### B3 (for the record)
doi:10.1126/science.aax3072 · "Sequencing metabolically labeled transcripts in single cells
reveals mRNA turnover strategies" (scEU-seq)

No PMC record; science.org 403s.

**Why:** the only dataset that gives **synthesis and degradation rates simultaneously, resolved
across the cell cycle** — and it was done in **RPE1**, the line Heldt 2018 was fitted in. Every
other stability dataset we hold is a single steady-state number per gene; this one says how
turnover itself changes through the cycle, in a model whose entire subject is cycle-dependent
regulation.

*Status of B1/B2: both supplied 2026-09-24 and ingested — see `protein_noise_data_exploration.md` §7.*

---

## C. Don't bother — I can get these myself

Checked and reachable; I'll fetch them unless you say otherwise.

| paper | route | what it gives |
|---|---|---|
| Zatulovskiy 2020 *Science* | PMC7489475 | Rb dilution by growth; size-independent synthesis |
| Padovan-Merhar 2015 *Mol Cell* | PMC4402149 | burst frequency vs cell size; S-phase dosage compensation |
| Lanz 2022 *Mol Cell* | PMC9444988 | proteome-wide concentration scaling with cell size |
| Li, Bickel & Biggin 2014 *PeerJ* | PMC3940484 (fully OA) | the Schwanhäusser reanalysis — abundances and translation rates several-fold higher |
| Beentjes 2020 *Phys Rev E* | arXiv:1912.04337 | exact solutions with bursting + replication + division; the analytic test cases for T6 |

---

## D. Beyond the ask — one experiment that would settle the central question

Flagging this because it is cheap and it targets the thing nothing in the data can currently give.

The brief's core problem is that a snapshot cannot distinguish frozen heterogeneity from fast
fluctuation, and that correlation times therefore have to come from priors. The priors now say
**99% of proteins have τ > 8 h and 0% have τ < 2 h** — but that is borrowed from NIH3T3 and H1299,
not measured in GNPs.

You can measure it in **fixed** cells, without live imaging, using **sister-cell pairs**: in
recently divided doublets, compare the CV of the sister-to-sister ratio against the CV of the
whole population. The ratio of those two variances is the fraction of the variance that is
heritable, i.e. slow relative to the time since division. Repeating at two post-mitotic times
gives the decay directly.

For CyclinD1 this would discriminate τ ≈ 1.7 h from τ ≈ 20 h in one experiment — which is exactly
the 12× uncertainty that A4 is otherwise trying to resolve indirectly. Sister pairs can be
identified from recently divided doublets by morphology, or with a short EdU/lineage window.

Entirely your call whether that is worth the bench time — mentioning it because it converts the
project's biggest assumption into a measurement.
