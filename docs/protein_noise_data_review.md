# Data review: what we hold, what it gives us, what is missing

Stock-take before starting the analysis. Regenerate the per-gene table with
`python3 simulations/audit_gene_coverage.py` (matrix also at
`data/priors/gene_coverage_matrix.csv`).

---

## 1. The ten datasets

24 raw files, ~80 MB. Raw supplements are gitignored; `data/priors/SOURCES.md` carries URL, date,
sha256 and **dilution convention** for each.

| dataset | system | n | what it uniquely gives |
|---|---|---|---|
| **Schwanhäusser 2011** | mouse NIH3T3 | 5029 | the **only** source of absolute **protein copy number**, plus mRNA copies, both half-lives, transcription and translation rates |
| **Herzog 2017** (SLAM-seq) | mouse ESC | 5969 | **best mRNA half-lives** — metabolic labelling, no transcription block |
| **Larsson 2019** (S1/S3/S6) | mouse fibroblast, ESC, human fibro | 7870 / 7309 / 1609 | **burst size and frequency**, allele-resolved. S6 pairs mouse↔human, which gave the transfer asymmetry |
| **Larsson S2** | — | 3 regressions | promoter sequence → burst kinetics (TATA) |
| **McShane 2016** | NIH3T3 **and RPE** | 3293 / 3134 | **two-state** degradation: fast unassembled `kA` → stable assembled `kB`, plus ED/NED class. RPE is Heldt's line |
| **Mathieson 2018** | 5 primary cell types | 8805 | post-mitotic → **dilution-free** protein half-lives; the class-(b) anchors |
| **Battich 2020** (scEU-seq) | RPE1-FUCCI | 528 | **cycle-resolved** synthesis and degradation. The only source of *how turnover changes through the cycle* |
| **Dörrbaum 2018** | rat neuron / glia | 3778 / 4126 | neural-lineage cell-type pair, one method |
| **Eden 2011** | human H1299 | 100 | bleach-chase in **living** cells; an independent method, used as arbiter |
| **Sigal 2006** | human H1299 | 40 clones | **extrinsic** noise: mixing time and CV |
| **Ochiai 2020** | mouse ESC | count matrices | bursting × ChIP-seq. **Fitted kinetics are not deposited** — only raw counts |

### Structure of the derived tables

`data/priors/turnover_burst_priors.csv` — **924 rows, long/tidy**: one row per
`gene × quantity × source × cell_type`, with `value`, `units`, `table`, `growth_corrected` and
`notes`. **313 numeric slots filled, 512 explicit `NA`** — nothing imputed. Long format because
the same quantity has several conflicting sources and they must stay separable.

Also `sigal2006_mixing_times.csv` (40 rows), `eden2011_bleachchase_halflives.csv` (100 rows,
parsed from the SOM PDF), `gene_coverage_matrix.csv` (this audit).

---

## 2. Per-gene coverage

33 genes mapped to v44 species or parameters. Coverage by quantity:

| quantity | coverage |
|---|---|
| **mRNA half-life** (Herzog) | **26/33 — 79%** |
| burst size (Larsson, either cell type) | 20/33 — 61% |
| Battich cycle strategy | 15/33 — 45% |
| protein copies / protein t½ (Schwanhäusser) | 13/33 — 39% |
| protein t½ (McShane) | 12/33 — 36% |
| translation rate | 11/33 — 33% |
| protein t½, post-mitotic (Mathieson) | 7/33 — 21% |

### Can each gene actually be parameterised?

The noise budget needs five things: `⟨p⟩`, `⟨m⟩`, `B`, `t½(protein)`, `t½(mRNA)`.

**Complete (5/5) — 9 genes:** `Ccnd1`, `Ccna2`, `Ezh2`, `Skp2`, `Cdk2`, `Cdk4`, `Chek1`, `Eed`,
`Suz12`.

Of those, only **four are explicit v44 species** — `Ccnd1`, `Ccna2`, `Ezh2`, `Skp2`. The rest map
to lumped parameters (`Cdk4` into `Cd`, `Cdk2` into `Ce`/`Ca`) or are implicit (`Eed`, `Suz12` are
PRC2 partners, not species).

Fortunately those four include **both burst candidates and the epigenetic brake** — the nodes the
noise work most depends on.

**Partial (3–4/5) — 17 genes.** And the pattern is systematic: nine of them are exactly
`.mB.T` — mRNA copies, burst size and mRNA half-life present, **both protein columns missing**
(`Ccne1`, `Cdkn1a`, `E2f1`, `Trp53`, `Mycn`, `Gli2`, `Dtl`, `Cdc25a`, `Wee1`).

**Essentially absent (≤2/5) — 7 genes:** `Ccnb1`, `Cdkn2d`, `Fbxo5`, `Gli1`, and — with **nothing
at all** — `Gli3`, `Ptch1`, `Smo`.

---

## 3. The systematic gap is the protein side

The mRNA and bursting side is in reasonable shape (79% and 61%). **Protein abundance and protein
half-life are the bottleneck at 21–39%**, and that one gap is what demotes most genes from
complete to partial.

The cause is structural: **Schwanhäusser is the only source of absolute protein copy number**, and
it covers 5029 abundant proteins — which systematically misses low-abundance regulators. Exactly
our genes. `Cdkn1a`, `E2f1`, `Trp53`, `Ccne1` and `Mycn` are all absent from it.

That matters unevenly, though. In `CV² = 1/⟨p⟩ + (p_cap+B)/⟨m⟩ · filter`, the `1/⟨p⟩` Poisson term
is **negligible above ~10³ copies**, which nearly all these proteins exceed. So a missing `⟨p⟩`
costs us little. **A missing protein half-life is the real loss** — it sets the filter factor and
the correlation time, and without it τ collapses to the dilution-only assumption `1/µ` (flagged
`[!]` in `analyze_noise_budget.py`).

### The Hedgehog arm

| gene | what we have |
|---|---|
| `Gli1` | mRNA t½ only (3.05 h, Herzog) |
| `Gli2` | mRNA t½, burst size, ⟨m⟩ |
| `Gli3` | **nothing** (a Battich strategy call only) |
| `Ptch1` | **nothing** |
| `Smo` | **nothing** |
| `Mycn` | mRNA t½, burst (ESC only), ⟨m⟩ — **no protein half-life anywhere** |

And zero BioModels entries annotate EZH2, MYCN, GLI2 or GLI3, so there is no published model to
borrow structure from either. **The Hedgehog arm will run on assumed values.** Worth stating in
the methods rather than discovering per-parameter.

---

## 4. What is missing, by consequence

**Blocked — nothing in the literature fixes these:**

1. **HCI per-cell tables.** Every predicted CV is untested. Nothing to compare against.
2. **An absolute scale (Ω).** Model species are in arbitrary units (`DATA_PROVENANCE.md`); discrete
   bursts and Gillespie both need counts. One number would anchor it: GNP nuclear volume, or one
   absolute copy number.
3. **Measured T_cc** for GNP and MB. Scales µ and converts Sigal's extrinsic time from generations
   to hours.

**Would materially improve the parameterisation:**

4. **Protein half-lives for the nine `.mB.T` genes** — the single highest-leverage gap, and it is
   what stands between "partial" and "complete" for half the gene set.
5. **Anything at all for `Ptch1`, `Smo`, `Gli3`.**
6. **A proliferating neural progenitor dataset.** Every neural dataset held is post-mitotic, so the
   mismatch with GNP is cycle state, not lineage.

**Known-unusable rather than missing:**

7. **Ochiai's fitted burst kinetics are not deposited** — only the count matrices. Getting per-gene
   values means re-inferring from Tables S1/S2 (the `txburst` route), or reading them off figures.
8. **Battich is relative only** — log2 of each gene's own cycle median. It constrains *modulation*,
   not level, so it cannot supply a mean rate.

---

## 5. Bottom line

* The **four species the noise work most depends on** — `Ccnd1`, `Ccna2`, `Ezh2`, `Skp2` — are
  fully parameterisable today.
* **mRNA and bursting are in good shape**; the **protein side is the bottleneck**, and it is
  concentrated in exactly the low-abundance regulators.
* The **Hedgehog arm has essentially nothing** and will be assumption-driven regardless of what
  else we obtain.
* Everything still blocking is **lab data**, not literature. There is nothing left to download.
