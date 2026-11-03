# Absolute abundance and subcomplex stoichiometry of PRC2 and related chromatin enzymes

Literature census, compiled 2026-10-23.

**Standard applied:** every number below was read in the Methods, Results, figure legend, or
a machine-readable supplementary table of the primary paper, except where the row is explicitly
marked `SECONDARY` (an introduction/discussion statement citing someone else) or `ABSTRACT-ONLY`.
Nothing here comes from a review's summary figure.

**Confidence grades**

| Grade | Meaning |
|---|---|
| **A** | Direct absolute measurement with a physical standard (SRM/PRM with isotope-labelled peptides, calibrated FCS/single-molecule counting), read from the primary text. |
| **B** | Label-free absolute estimate (iBAQ / proteomic ruler / LFQ stoichiometry) read from the primary text or the paper's own supplementary data table. Typical accuracy 1.5–3x. |
| **C** | Primary paper, but the number is relative, derived, or a statement in the Intro/Discussion citing another paper (`SECONDARY`). |
| **X** | Abstract-only, or value could not be read in the primary text. |

---

## 0. Headline findings

1. **There is no SRM/PRM, isotope-dilution, or calibrated-FCS absolute copy number for any PRC2
   subunit anywhere in the indexed literature.** I searched EuropePMC for
   `"selected reaction monitoring" AND EZH2`, `"parallel reaction monitoring" AND (EZH2 OR SUZ12)
   AND absolute`, `"EZH2 molecules per cell"`, `"PRC2 molecules per cell"`, `recombinant AND
   "quantitative western" AND (EZH2 OR SUZ12) AND "nuclear extract"`, and
   `"absolute quantification" AND chromatin AND "copies per cell"` — zero hits returning a PRC2
   copy number. **PRC2 has no Holzmann-equivalent paper.** Every absolute number available is
   label-free (iBAQ / proteomic ruler), i.e. Grade B at best.
2. The best single independent cross-check is **Itzhak et al. 2016 eLife** (HeLa, proteomic ruler,
   six replicate maps), whose own supplementary table gives **EZH2 = 142,635 copies/cell**
   — **2.8x higher than your PaxDb-derived ~50,700**. See §1 and §6.
3. The SUZ12 / EED excess over EZH2 is **real in the data but modest**: EZH2 : SUZ12 : EED =
   1 : 1.37 : 2.47 in the same HeLa dataset — not the 2x for both that a casual read of PaxDb
   suggests. Affinity-purification stoichiometry says the assembled core is strictly 1:1:1
   (Smits 2013), so the excess is free subunit, iBAQ protein-specific bias, or both (§2).
4. **Only ~20% of PRC2 is chromatin-bound at any instant** (Youmans 2018, single-particle
   tracking, Grade A for the fraction). Composing this with Itzhak gives
   **~26,000–29,000 chromatin-engaged EZH2 molecules per cell ≈ 1 per ~1,600 nucleosomes**.
5. Your check on PaxDb vs RAD21 does not reproduce in the Itzhak data: the same histone-normalised
   approach there gives RAD21 **1.43x above**, and SMC3 **2.96x above**, the Holzmann SRM values —
   i.e. over-reporting, not 5.2x under-reporting. The 5.2x is most likely an artefact of PaxDb's
   cross-dataset integration rather than a property of histone normalisation (§6).

---

## 1. Absolute copy numbers (Q1, Q4)

### 1a. Itzhak, Tyanova, Cox & Borner 2016, *eLife* 5:e16950 (PMC4959882)

Read directly from the paper's **Supplementary file 1, sheet "Compact HeLa Spatial Proteome"**
(column *"Estimated Copy number per cell"*). Method stated in the Notes sheet: *"Estimates were
calculated using the Proteomic Ruler approach (Wisniewski et al., MCP 13:3497–3506; 2014), and
represent medians from six replicate measurements."* Cell type **HeLa**, asynchronous (no cell-cycle
sorting). These are **TOTAL CELLULAR** copies. n = 6 replicate organellar maps.

The *"Cytosolic Pool"* column is defined in the Notes sheet as *"the proportion of a protein that
remains in the cytosol after lysis and centrifugation at 78,400 x g for 30 min"* — so
`1 − cytosolic` is **nuclear + organellar pellet**, which for these nuclear proteins is a
reasonable nuclear proxy but is **not** chromatin-bound.

| Protein | **TOTAL** copies/cell | Conc. (nM) | Cytosolic fraction | ⇒ non-cytosolic (≈nuclear) | Abundance percentile | Grade |
|---|---:|---:|---:|---:|---:|:--:|
| **EZH2** | **142,635** | 88.2 | 0.083 | ~130,800 | 58.7 | B |
| EZH1 | 9,699 | 6.0 | 0.137 | ~8,400 | 18.5 | B |
| **SUZ12** | **194,816** | 120.5 | 0.084 | ~178,400 | 63.7 | B |
| **EED** | **351,948** | 217.7 | 0.107 | ~314,300 | 72.9 | B |
| RBBP4 | 1,285,129 | 795.0 | 0.382 | ~794,300 | 88.0 | B (see note) |
| RBBP7 (Q16576 canonical) | 702,123 | 434.3 | 0.629 | ~260,300 | 82.1 | B (see note) |
| RBBP7 (isoform Q16576-2, separate row) | 6,459,071 | 3,995.7 | 0.531 | — | 96.8 | **X** — implausibly large, almost certainly an isoform/protein-group assignment artefact; do not use |
| AEBP2 | 41,548 | 25.7 | 0.699 | ~12,500 | 38.5 | B |
| JARID2 | 8,861 | 5.5 | 0.000 | ~8,900 | 17.3 | B |
| MTF2 (PCL2) | 138,218 | 85.5 | 0.000 | ~138,200 | 58.3 | B |
| PHF1 (PCL1) | 13,853 | 8.6 | 0.000 | ~13,900 | 22.6 | B |
| PHF19 (PCL3) | 20,227 | 12.5 | 0.020 | ~19,800 | 27.6 | B |
| EPOP (C17orf96) | 6,482 | 4.0 | 0.000 | ~6,500 | — | B |
| PALI1 (LCOR) | 9,364 | 5.8 | 0.000 | ~9,400 | — | B |
| **EP300** | **27,961** | 17.3 | 0.046 | ~26,700 | 32.5 | B |
| **CREBBP** | **54,919** | 34.0 | 0.054 | ~52,000 | 42.9 | B |
| **KDM6A** | **4,183** | 2.6 | 0.203 | ~3,300 | 10.1 | B |
| **KDM6B** | not detected | — | — | — | — | **X** |
| **HDAC1** | **461,170** | 285.3 | 0.126 | ~402,900 | 76.7 | B |
| **HDAC2** | **1,542,148** | 954.0 | 0.150 | ~1,310,900 | 89.5 | B |
| HDAC3 | 103,977 | 64.4 | 0.576 | ~44,100 | — | B |
| RNF2 (RING1B) | 133,236 | 82.4 | 0.033 | ~128,800 | — | B |
| RING1 (RING1A) | 188,100 | 116.4 | 0.132 | ~163,300 | — | B |
| BMI1 (PCGF4) | 102,982 | 63.7 | 0.000 | ~103,000 | — | B |
| CBX7 | not detected in this dataset | — | — | — | — | X |
| CBX2 / CBX4 / CBX8 | 46,050 / 153,900 / 348,400 | 28.5 / 95.2 / 215.6 | 0.026 / 0.135 / 0.000 | — | — | B |
| KDM2B | 21,300 | 13.2 | 0.181 | — | — | B |
| EHMT2 (G9a) | 127,100 | 78.7 | 0.034 | — | — | B |
| SETDB1 | 40,000 | 24.8 | 0.643 | — | — | B |
| DNMT1 | 830,900 | 514.0 | 0.032 | — | — | B |
| CTCF | 253,400 | 156.8 | 0.001 | — | — | B |
| POLR2A | 342,500 | 211.9 | 0.228 | — | — | B |
| **RAD21** (benchmark) | **501,337** | 310.1 | 0.009 | ~496,900 | 77.7 | B |
| **SMC3** (benchmark) | **1,233,877** | 763.3 | 0.014 | ~1,216,600 | 87.7 | B |
| SMC1A (benchmark) | 1,392,000 | 861.3 | 0.012 | — | — | B |
| STAG1 / STAG2 | 99,420 / 366,800 | 61.5 / 226.9 | 0.002 / 0.020 | — | — | B |
| NIPBL | 152,400 | 94.3 | 0.040 | — | — | B |

*RBBP4/7 note:* MaxQuant protein-grouping of the two near-identical paralogues plus the rogue
Q16576-2 row makes the absolute RBBP4/RBBP7 figures the least trustworthy in the table. Treat
RBBP4+RBBP7 as "~1–2 x 10⁶ total, shared across PRC2, NuRD, CAF-1, Sin3 and others", not as a
PRC2-specific pool.

**Derived quantities (my arithmetic on the above; flagged as derived, not measured):**

- EZH2 : SUZ12 : EED = **1 : 1.37 : 2.47**; (EZH1+EZH2) : EED = 1 : 2.31.
- Nucleosomes in the same dataset (see §5): ~3.9 x 10⁷ (sum of core histones / 8) to
  ~4.6 x 10⁷ (mean of H4, H2B, H2A / 2).
- **EZH2 ≈ 1 per 275–323 nucleosomes** (total cellular), vs your PaxDb figure of 1 per 594
  → **2.2x discrepancy**.
- EZH2 x 0.20 chromatin-bound fraction (§3) = **~28,500 chromatin-engaged EZH2 ≈ 1 per ~1,600
  nucleosomes**. This is the number most relevant to a lattice/spreading model.

### 1b. Everything else I could find

| Source | What it gives | Why it is not a copy number |
|---|---|---|
| Wiśniewski, Hein, Cox & Mann 2014 *MCP* 13:3497 (PMC4256500) | The method itself. Validation: proteomic-ruler copy numbers recapitulated PrEST-SILAC spike-in values **"within an average deviation of 1.5-fold"**, using on average 19.4 peptides/protein. Histone MS signal was **2.07–4.03% of total MS signal** across four human lines; total protein per cell recovered within **1.24 ± 0.29x** of cell-counting. | Benchmark for how much to trust §1a; no PRC2 values. Grade A for the benchmark statement. |
| Beck et al. 2011 *MSB* 7:549 (PMC3261713); Nagaraj et al. 2011 *MSB* 7:548 (PMC3261714) | U2OS / HeLa deep proteomes, iBAQ. | PRC2 subunits are only in the supplementary tables; these datasets are already inside PaxDb, so they are **not independent** of your estimate. |
| Ly et al. 2014 *eLife* 3:e01630 (PMC3936288) | Cell-cycle-resolved proteome of NB4 cells (G1/S/G2M, elutriation), >9,000 genes. | Reports **iBAQ-scaled intensities only, no copies per cell** — I checked the Results text; the only copy-number statement is a reference to Nagaraj for dynamic range. So **no cell-cycle-resolved absolute PRC2 numbers exist**. Grade X for this purpose. |
| Hein et al. 2015 *Cell* 163:712 (PMID 26496610); Schwanhäusser et al. 2011 *Nature* 473:337 (PMID 21593866); Sharma et al. 2015 *Nat Neurosci* 18:1819 (PMC7116867) | HeLa / NIH3T3 / mouse-brain abundances | Not in PMC OA; could not read the primary tables in this pass. Flagged in §8. |

---

## 2. What fraction of EZH2 is in a complete PRC2? Free pools (Q2)

| Finding | Value | Compartment | Method | Cell type | n | Source | Grade |
|---|---|---|---|---|---|---|:--:|
| Assembled PRC2 core is a strict monomer of each subunit | **EED : EZH1/2 : SUZ12 = 1 : 1 : 1** | nuclear extract (soluble) | GFP affinity enrichment of BAC-expressed EED-GFP at near-endogenous level + label-free iBAQ stoichiometry; corrected iBAQ of bait set to 1 | HeLa | 3 GFP pulldowns vs 6 controls (9 total) | Smits, Jansen, Poser, Hyman, Vermeulen 2013 *NAR* 41:e28 (PMC3592467), Results "PRC2" + Fig. 4B,C | B |
| Core stoichiometry is unchanged on differentiation | Ezh1/2 : Eed : Suz12 stable ESC→NPC | nuclear extract | Eed-GFP BAC pulldown, iBAQ relative to Ezh1+2 = 1 | mouse ESC and NPC | 3 pulldowns | Kloet et al. 2016 *NSMB* 23:682 (PMC4939079), Fig. 1d | B |
| Cellular EZH2 : SUZ12 : EED ratio | **1 : 1.37 : 2.47** | TOTAL cellular | proteomic ruler | HeLa | 6 | Itzhak 2016 supp. table (§1a) | B |

**Interpretation.** The complex is 1:1:1 but the cell contains ~2.5x more EED than EZH1+EZH2
protein. Two explanations are not separated by any paper I found:

- **a genuine free EED (and modest free SUZ12) pool.** EED is a WD40 protein that also sits in
  EED–EZH1 and in PRC2-independent contexts; a 1.4–2.5x excess is entirely plausible.
- **iBAQ protein-specific bias.** The proteomic-ruler benchmark is 1.5x average deviation
  (Wiśniewski 2014), so a 1.37x SUZ12:EZH2 ratio is *within noise* and only the 2.5x EED figure
  exceeds it. EED is a small (50 kDa), peptide-rich, well-digested protein; EZH2 (85 kDa, SET
  domain) is not — the bias runs in exactly the direction observed.

**Conclusion for your model:** treat the number of catalytically competent PRC2 complexes as set by
**EZH1+EZH2, not by EED or SUZ12** — i.e. ~1.5 x 10⁵ per HeLa cell, with a free EED pool of
~2 x 10⁵ that should be regarded as unproven. No paper measures the free pool directly (e.g. by
size-exclusion of nuclear extract with absolute quantification of each fraction) — see §8.

---

## 3. Chromatin-bound fraction (the TOTAL vs CHROMATIN distinction)

This is the PRC2 literature's best-measured quantity and the one your model most needs.

| Quantity | Value | Compartment | Method | Cell type | n | Source | Grade |
|---|---|---|---|---|---|---|:--:|
| EZH2 chromatin-bound fraction | **20% ± 3%** | fraction of **NUCLEAR** molecules | live-cell single-particle tracking of endogenous CRISPR knock-in 3xFlag-HaloTag-EZH2, JF646, 97.5 Hz, Spot-On model fitting | human U2OS | 3 biological replicates | Youmans, Schmidt & Cech 2018 *Genes Dev* 32:794 (PMC6049511), Results + Fig. 2D | **A** |
| SUZ12 chromatin-bound fraction | **25% ± 3%** | NUCLEAR | as above, HaloTag-SUZ12 knock-in | U2OS | 3 | same, Fig. 2D | **A** |
| EZH2 D_free | 2.09 ± 0.12 µm²/s | nucleoplasm | SPT | U2OS | 3 | same, Fig. 2C | A |
| SUZ12 D_free | 1.99 ± 0.13 µm²/s | nucleoplasm | SPT | U2OS | 3 | same | A |
| EZH2 bound fraction **at an H3K27me3-rich lacO array** | **~40%** | locus-local | FRAP immobile fraction | U2OS F42B8 | ≥3 | same, Discussion + Fig. 3 | A |
| Residual bound fraction when SUZ12 cannot bind PCL/AEBP2 | **0.08** (≈ HaloTag-alone background ~0.10) | NUCLEAR | SPT, SUZ12(VEFS) / SUZ12(PHF1+ABHmut) | U2OS | 3 | same, Fig. 5D,E | A |
| Overexpressed eGFP-EZH2 D_free | ~7 µm²/s | nucleoplasm | FCS | — | — | Hetey et al. 2017, cited in Youmans 2018 Discussion — `SECONDARY`; the Youmans authors attribute the discrepancy to a free, non-complexed EZH2 population created by overexpression | C |

**Note the internal tension, exactly analogous to Holzmann's cohesin case:** SPT says 20% of PRC2
is chromatin-bound *genome-wide*, but ChIP-based and biochemical-fractionation pictures of PRC2
imply a much more concentrated, CpG-island-focused pool. Youmans 2018 also shows this bound
fraction is **almost entirely dependent on the accessory subunits (PCL + AEBP2/JARID2), not on
EED–H3K27me3 read-out**: 30 min of A-395 (EED–H3K27me3 inhibitor) changed neither the bound
fraction nor the residence time.

---

## 4. PRC2.1 vs PRC2.2 stoichiometry (Q3)

### 4a. Human, HeLa — relative to the PRC2 core = 1.0

Smits et al. 2013 (PMC3592467), Results section "PRC2", read verbatim from the primary text.
Method: EED-GFP BAC pulldown from HeLa nuclear extract, iBAQ stoichiometry, n = 3.

| Accessory subunit | Stoichiometry vs core | Subcomplex |
|---|---:|---|
| RBBP4/RBBP7 (RbAp48/46) | **0.60** | both (core-adjacent) |
| **EPOP (C17orf96)** | **0.47** | PRC2.1 |
| **PCL1/2/3 (PHF1/MTF2/PHF19)** | **0.35** | PRC2.1 |
| **AEBP2** | **0.20** | PRC2.2 |
| **JARID2** | **0.12** | PRC2.2 |
| PALI1 (C10orf12) | **0.06** | PRC2.1 |

⇒ **PRC2.1 : PRC2.2 ≈ (0.35 + 0.06) : (0.20 + 0.12) ≈ 0.41 : 0.32**, i.e. roughly **55:45** in
HeLa nuclear extract, with ~27% of the core carrying neither PCL nor AEBP2/JARID2.
**Caveat stated by the authors:** they cannot distinguish genuinely substoichiometric binding from
partial wash-off during affinity purification, so all accessory values are **lower bounds** and the
"naked core" fraction is an upper bound.

### 4b. Mouse, ESC vs NPC — the cell-type switch (also answers part of Q5)

Kloet et al. 2016 *NSMB* 23:682 (PMC4939079). Eed-GFP BAC in mESC, differentiated to NPC;
iBAQ of each interactor divided by iBAQ of Ezh1+Ezh2 (set to 1); n = 3 pulldowns; Fig. 1d,e.
Numbers below are the fold-changes stated in the Results text (the absolute ratios are in Fig. 1d,
a heatmap I could not read numerically).

| Subunit | ESC vs NPC | Subcomplex |
|---|---|---|
| **Mtf2** | ~10x or more enriched in **ESC** | PRC2.1 |
| **Epop (C17orf96)** | ~10x or more enriched in **ESC** | PRC2.1 |
| Tceb1/Tceb2 (Elob/Eloc) | ~10x or more enriched in ESC | PRC2.1/EPOP module |
| **Phf19** | ~10x enriched in **NPC** | PRC2.1 |
| **Pali1 (Gm340)** | ~10x enriched in **NPC** | PRC2.1 |
| Phf1, Ehmt2 | NPC-specific (detected only in NPC) | PRC2.1 |
| Core (Ezh1/2, Eed, Suz12) | unchanged ratio | — |

**Therefore the PCL identity, not the PRC2.1/2.2 split per se, is what switches in the neural
lineage: MTF2→PHF19/PHF1.** That is directly relevant to a cerebellar/GNP system.

### 4c. Human quiescent vs cycling fibroblasts

McCole et al. 2025 *Cell Rep* 44:115192 (PMC11931288). TIG3T hTERT fibroblasts, 96 h serum
starvation; IP-MS of endogenous EZH1, EZH2, EED; iBAQ-based stoichiometry relative to the bait
after subtracting mean IgG iBAQ; n = 3 IPs. Read from Results + Fig. 1E,F legend.

- In quiescence **EZH1 replaces EZH2, and PHF1 replaces PHF19** in the PRC2 complex.
- **MTF2, EPOP, PALI1 and JARID2 were undetectable in PRC2 IP-MS from these fibroblasts**,
  consistent with their low expression — i.e. in a quiescent diploid fibroblast PRC2 is
  essentially *all* PRC2.1-PHF1, with no measurable PRC2.2.
- In quiescent cells **only 14% of H3K27me3 lies in "targeted regions" (EED-bound promoters)**;
  the other 86% is in dispersed regions. By contrast **90% of canonical PRC1 is at targeted
  regions.** (ChIP-Rx, Fig. 4H,I.) Grade B.

### 4e. Cross-check: AP-MS stoichiometry vs whole-cell abundance — EPOP does not reconcile

My arithmetic, dividing the Itzhak 2016 HeLa total-cellular copy numbers (§1a) by EZH2 = 142,635,
and comparing with the Smits 2013 HeLa AP-MS stoichiometry (§4a). Both are HeLa, so this is a
like-for-like test.

| Accessory | whole-cell copies / EZH2 copies | Smits AP-MS stoichiometry | Agreement |
|---|---:|---:|---|
| MTF2 | 0.97 | — | — |
| PHF1 | 0.10 | — | — |
| PHF19 | 0.14 | — | — |
| **PCL total (MTF2+PHF1+PHF19)** | **1.21** | 0.35 | **cellular pool is 3.5x larger than the complexed pool** |
| AEBP2 | 0.29 | 0.20 | good |
| JARID2 | 0.06 | 0.12 | good (same order) |
| PALI1 | 0.07 | 0.06 | excellent |
| **EPOP** | **0.045** | **0.47** | **10x contradiction** |

Two conclusions:

- **The PCL proteins are not limiting.** There is slightly *more* PCL protein in the cell than
  EZH2, yet AP-MS recovers only 0.35 PCL per core. Either two-thirds of PCL protein is free of
  PRC2, or PCL is washed off during purification. Either way the "55:45 PRC2.1:PRC2.2" split in
  §4a is **weakly constrained** and should not be used as a hard model parameter.
- **The EPOP stoichiometry of 0.47 cannot be right as a cellular fraction.** EPOP is 22x less
  abundant than EZH2 in whole cells, so at most ~4.5% of PRC2 can carry EPOP. Smits' 0.47 is
  presumably an iBAQ artefact for a small, peptide-poor protein, or reflects enrichment of an
  EPOP-containing subpopulation by the EED-GFP bait. **Do not model EPOP as present on half of
  PRC2.**

### 4d. Mechanistic context (not stoichiometry, but constrains the competition)

Youmans, Gooding, Dowell & Cech 2021 *Mol Cell* 81:488 (PMC7867654). Flag-SUZ12 knock-in hiPSC
(WTC-11), separation-of-function SUZ12 point mutants, tandem-affinity IP-MS with iBAQ normalised
to SUZ12 (n = 2 biological replicates, heatmap Fig. 1F / Table S1 — numeric values not readable
from the text). Key measured quantity in the primary text: **PRC2.1 binds dsDNA with
K_d,apparent ≈ 20–40 nM and higher affinity than PRC2.2** (EMSA, n = 3 independent experiments,
two protein preps). PRC2.1 and PRC2.2 **compete** for a limiting SUZ12 pool; forcing PRC2 into
PRC2.1 raises SUZ12 occupancy and H3K27me3 (~50% increase at HOXA2) at canonical Polycomb targets,
not at new loci. Grade B/C.

---

## 5. Nucleosomes per cell (Q6) — there is no direct count

**Short answer: nobody counts nucleosomes.** Every published "nucleosomes per cell" number is
one of two calculations, and the proteomic-ruler route is *circular* with respect to this question.

**Route 1 — histone mass = DNA mass (this is what the proteomic ruler assumes).**
Wiśniewski et al. 2014 (PMC4256500), Results: *"In eukaryotic cells, DNA is packaged in chromatin
by histones, and the mass of the DNA is about equal to the combined mass of histones"*; they use
**6.5 pg DNA for a diploid human cell** and a genome of 3.2 x 10⁹ bp at 615.9 Da/bp.
6.5 pg of core histone at a mean core-histone mass of ~13.7 kDa = 2.9 x 10⁸ core histones =
**~3.6 x 10⁷ octamers per diploid human cell**, equivalent to one nucleosome per **~178 bp**.
*Because this ratio is an input, any nucleosome count read out of an iBAQ/ruler dataset simply
returns the assumption.* Grade C for the purpose of Q6.

Consistency check on the Itzhak HeLa table (my arithmetic, derived):

| Histone | copies/cell | implied nucleosomes (÷2) |
|---|---:|---:|
| H4 | 1.09 x 10⁸ | 5.4 x 10⁷ |
| H2B (all groups) | 9.7 x 10⁷ | 4.8 x 10⁷ |
| H2A (incl. H2AX, H2AZ, macroH2A) | 7.1 x 10⁷ | 3.6 x 10⁷ |
| H3 (all groups) | 3.8 x 10⁷ | 1.9 x 10⁷ |
| **sum of core histones ÷ 8** | — | **3.9 x 10⁷** |

The 2.8x spread between H4 and H3 is a real and important warning for your normalisation:
**H3 is systematically under-recovered in shotgun MS** (heavy PTM load, poor tryptic behaviour of
the K27/K36 peptide). Your "sum of core histones / 8" happens to be robust to this (3.9 x 10⁷ vs
4.6 x 10⁷ from the H4/H2B/H2A mean, only 17% apart), but a per-histone normalisation would not be.

**Route 2 — genome size / nucleosome repeat length (NRL), measured by MNase digestion.**
This is the only genuinely *measured* handle, and NRL is strongly cell-type-dependent — including
within the brain, which matters for your system.

| System | NRL | Method | Source | Grade |
|---|---|---|---|---|
| Rabbit **cerebral cortex neurons** | **~160 bp** | micrococcal nuclease ladder on isolated nuclei | Thomas & Thompson 1977 *Cell* 10:633 (PMID 862024) — read from the abstract; full text not accessible | **ABSTRACT-ONLY / A-quality method** |
| Rabbit **cerebellar neuronal nuclei** | **~200 bp** | same | same | ABSTRACT-ONLY |
| Rabbit non-astrocytic glia (same cortex), liver | ~200 bp | same | same | ABSTRACT-ONLY |
| Mouse dorsal root ganglion neurons | ~165 bp | MNase-seq | cited in Shtumpf/Teif et al. 2026 *NAR* "Nucleosome spacing across cell types, diseases, and ages" (PMC12961432), §"NRL in different genomic regions" | C (`SECONDARY`) |
| Mouse cortical OPC / astrocytes | ~182 / ~183 bp | MNase-seq | same review | C (`SECONDARY`) |
| Mouse cortical neurons during development | **200 bp before birth → 170 bp at P14** | MNase | same review | C (`SECONDARY`) |
| Mammalian heterochromatin, genome average | ~190 bp | MNase-seq / ATAC re-analysis | same review | C (`SECONDARY`) |
| Human telomeric repeats | ~157 bp | MNase | same review | C (`SECONDARY`) |

**Direct consequence for your parameterisation.** The cerebellar result is the useful one:
Thomas & Thompson found **cerebellar neurons have the long (~200 bp) repeat**, unlike cortical
neurons. For mouse the diploid genome is **5.4 Gb, not 6.4 Gb** (2.7 Gb haploid) — so at 200 bp
the mouse cerebellar nucleosome count is **~2.7 x 10⁷ per cell**, about **16% below** the
6.4 Gb/200 bp = 3.2 x 10⁷ you are currently using, and your figure appears to be the *human*
diploid number. If GNPs behave like proliferating progenitors with the shorter progenitor-type
repeat (~180 bp), the count rises to ~3.0 x 10⁷.

---

## 6. Does anything contradict EZH2 ≈ 50,000 copies/cell? (your explicit flag request)

| Comparison | Value | Ratio to your 50,700 | Verdict |
|---|---:|---:|---|
| Itzhak 2016, HeLa, proteomic ruler, TOTAL | 142,635 | **2.8x** | **CONTRADICTS — flag.** Same normalisation family (histone-anchored), one deep dataset rather than an integration. |
| Itzhak 2016, HeLa, non-cytosolic (≈nuclear) | ~130,800 | 2.6x | contradicts |
| Itzhak 2016 x Youmans 20% chromatin-bound | ~28,500 | 0.56x | **Your number sits between total and chromatin-bound.** If your 50,700 was meant as total cellular it is ~3x low; if it happens to approximate the chromatin-engaged pool it is ~2x high. |
| EZH2 per nucleosome, Itzhak | 1 per 275–323 | vs your 1 per 594 | **2.2x discrepancy** |
| EZH2 per nucleosome, chromatin-bound only | 1 per ~1,600 | — | the model-relevant figure |

### Checking your RAD21 reasoning

You report that PaxDb under-reports RAD21 relative to histone by **5.2x** against the Holzmann SRM
value of 350,000. **The same comparison in the Itzhak dataset goes the other way:**

| Protein | Itzhak HeLa (proteomic ruler) | Holzmann 2019 SRM (G1 HeLa) | Ratio |
|---|---:|---:|---:|
| RAD21/SCC1 | 501,337 | 350,000 | **1.43x over** |
| SMC3 | 1,233,877 | 417,000 | **2.96x over** |
| SMC1A | 1,392,000 | (expected ≈ SMC3) | internally consistent with SMC3 |

Three things follow:

1. **Histone-anchored normalisation is not intrinsically 5x low for RAD21.** A single deep,
   six-replicate dataset normalised the same way comes out 1.4x *high*. Your 5.2x is therefore far
   more likely a property of **PaxDb's integration step** — it averages ppm across dozens of
   heterogeneous studies of wildly differing depth, and shallow studies systematically
   under-sample low-abundance proteins while over-weighting histones. That is a different and
   larger error than the method's intrinsic one.
2. **The method's intrinsic error is ~1.5x, not 5x.** Wiśniewski 2014 benchmarked the proteomic
   ruler against PrEST-SILAC spike-ins and recovered copy numbers "within an average deviation of
   1.5-fold" (Results, Fig. 5A). Anything beyond ~2–3x needs a different explanation.
3. **Anchoring on a single protein propagates that protein's error to everything else.** By
   anchoring on RAD21 = 350,000 and correcting a putative 5.2x, you may have moved EZH2 down by a
   factor that belongs to RAD21's PaxDb entry, not to EZH2's. I would recommend re-anchoring on
   the histone ruler directly (with the per-histone caveat in §5) and reporting a 3x uncertainty
   band rather than a point estimate.

**Recommended working values for the model** (human, proliferating, asynchronous), stated with
honest error bars:

- EZH1+EZH2 total: **1.5 x 10⁵ (range 0.5–3 x 10⁵)**
- Chromatin-engaged EZH1/2 at any instant: **3 x 10⁴ (range 1–6 x 10⁴)**
- Nucleosomes: **4 x 10⁷ human diploid-equivalent HeLa; 2.7 x 10⁷ mouse diploid at NRL 200 bp**
- Chromatin-engaged PRC2 per nucleosome: **~1 : 1,000–2,000**

---

## 7. H3K27 methylation stoichiometry (Q7)

Honest state of the field: the widely quoted percentages all trace back to **three primary MS
papers, none of which is openly accessible**, and almost every modern paper restates them as
intro-sentence citations. I could not read the primary numbers; everything below marked
`SECONDARY` is a statement in another primary paper's text that cites them.

| Statement | Value | Cell type | Method | Primary source | Where I read it | Grade |
|---|---|---|---|---|---|:--:|
| H3K27me2 is the most abundant state | **~50% of all H3** | mouse ESC | histone MS | Peters et al. 2003 *Mol Cell* 12:1577 (PMID 14690609) — **not accessible** | Højfeldt et al. 2018 *NSMB* 25:225 (PMC5842896), Introduction, citing refs 2/16/17 | C `SECONDARY` |
| H3K27me1 and H3K27me3 | **each 10–15% of H3** | mouse ESC | histone MS | same cluster | Højfeldt 2018 Intro (PMC5842896) | C `SECONDARY` |
| H3K27me2 coverage | **>70% of total H3** | mouse ESC | histone MS | Ferrari et al. 2014 *Mol Cell* 53:49 (PMID 24289921) — **not in PMC, cell.com Cloudflare-blocked** | Lavarone et al. 2019 *Nat Commun* 10:1679 (PMC6459869), Introduction + Discussion, citing ref 17 | C `SECONDARY` |
| H3K27me2 45–70% of all H3; H3K27me3 5–10% | range across flies and mammals | fly S2 + mammals | histone MS | Peters 2003; Ebert 2004; Jung 2010; Voigt 2012; Ferrari 2014 | Lee et al. 2015 *Genome Res* 25:1170 (PMC4510001), Introduction | C `SECONDARY` |
| K27me2/K36me2 peptides = **55–62%** of all observed H3.2+H3.3 (aa 27–40) peptides | 55–62% | mouse ESC (control) | SILAC + RP-HPLC + LTQ-Orbitrap, "GP" relative abundance of modified peptide species | **Jung, Pasini, Helin & Jensen 2010 *MCP* 9:838 (PMC2871418)** | **primary Results text, read directly** | **B** |
| H3.3K27me3 | **~10%** in WT mESC | mouse ESC | as above | same | same Results text | B |
| K27me2/K36me2 remains **~12%** after *Suz12* deletion | 12% | *Suz12*-KO mESC | as above | same | same Results text | B |
| K27me1-K36me2 is the dominant combinatorial species (60% of signal) in MMSET-high myeloma; in MMSET-low, K27me1-K36me2 = 18% and K27me2-K36me0 = 29% | — | NTKO / TKO myeloma lines | SILAC + triple-quadrupole MRM of the H3 27–40 peptide, 84 SILAC intermediates monitored, kinetic modelling | **Zheng et al. 2012 *PNAS* 109:13549 (PMC3427122)** | primary Results text, read directly | **B** |
| Only **14% of H3K27me3 is at EED-bound promoters**; 86% is dispersed. 90% of cPRC1 is at targeted regions | 14% / 86% / 90% | quiescent human TIG3T fibroblasts | ChIP-Rx with spike-in normalisation, genome partitioned into targeted vs dispersed | McCole et al. 2025 *Cell Rep* (PMC11931288) | primary Results, Fig. 4H,I | B |
| Absolute occupancy % of individual H3 marks, interphase vs mitosis, measured against **known concentrations of synthetic internal-standard peptides** | values are in Fig. 2B only | HeLa S3 | targeted MS with synthetic standards (Broad "global chromatin profiling") | Javasky et al. 2018 *Genome Res* 28:1455 (PMC6169886) | primary Results — method read, **numbers not extractable from text** | **X** (recoverable from their Fig. 2B / supplement) |

**Best defensible numbers for your model**, stated as a range because the primary sources are
inaccessible: of total H3K27 in a mouse ES/neural-progenitor-like cell,
**me0 ≈ 10–20%, me1 ≈ 10–15%, me2 ≈ 50–70%, me3 ≈ 5–15%**.

**H3K27me3-marked nucleosomes per cell (derived, not measured).** With 2 H3 per nucleosome and
H3K27me3 largely symmetric on marked nucleosomes (Voigt et al. 2012, PMC3498816 — not read):

| me3 fraction of H3 | marked H3 molecules | marked nucleosomes (symmetric) | marked nucleosomes (fully asymmetric upper bound) |
|---|---:|---:|---:|
| 5% | 4.6 x 10⁶ | 2.3 x 10⁶ | 4.6 x 10⁶ |
| 10% | 9.2 x 10⁶ | 4.6 x 10⁶ | 9.2 x 10⁶ |
| 15% | 1.4 x 10⁷ | 6.9 x 10⁶ | 1.4 x 10⁷ |

(using 4.6 x 10⁷ nucleosomes/cell). So **2–7 x 10⁶ H3K27me3 nucleosomes per cell**, against
~3 x 10⁴ chromatin-engaged PRC2 complexes — **roughly 100–200 marked nucleosomes per bound
enzyme**, which is the single most important ratio this census produces for a maintenance model.

---

## 8. What is still missing

**Hard gaps — no paper exists, as far as I can determine**

1. **Any Grade-A absolute copy number for any PRC2 subunit.** No SRM/PRM with isotope-labelled
   reference peptides, no quantitative western against a recombinant standard with stated
   recovery, no calibrated FCS or single-molecule photobleaching count. PRC2 simply has no
   Holzmann-2019 equivalent. Every number in §1 is Grade B.
2. **Cell-cycle-resolved absolute PRC2 abundance.** Ly 2014 (the obvious candidate) reports only
   iBAQ intensities. No G1 / S / G2M copy numbers for EZH2 exist. Given that your project is about
   S-phase dilution of EZH2, this is the most important missing measurement — and it may be worth
   generating rather than finding.
3. **Direct nucleosome counting.** Both available routes are calculations (§5). The highest-value
   missing experiment is a quantitative amino-acid-analysis or isotope-dilution measurement of
   total core histone per nucleus in a defined cell type, paired with its own MNase NRL.
4. **Absolute free-pool measurement.** Nobody has size-fractionated nuclear extract and absolutely
   quantified EZH2/SUZ12/EED in monomeric vs complexed fractions. The 2.5x EED excess remains
   unexplained.
5. **Cerebellar granule neuron progenitor or medulloblastoma proteome with copy numbers.** Searched
   `"cerebellar granule neuron progenitor" AND proteome` and `medulloblastoma AND proteomic AND
   "copy number" AND "mass spectrometry"` — nothing with per-cell copy numbers. The medulloblastoma
   proteogenomics papers are TMT-relative.
6. **KDM6B (JMJD3) absolute abundance anywhere.** Not detected in the Itzhak HeLa proteome.

**Soft gaps — the number exists but I could not read it in this pass**

7. **Kloet 2016 Fig. 1a bottom table** explicitly *"lists the absolute abundance of Eed, Ezh2 and
   Suz12 in nuclear extracts from each cell type"* (mESC and NPC), taken from **Spruijt et al. 2013
   *Cell* 152:1146**. **This is the single highest-value thing left to retrieve**: it is mouse,
   it is ESC *and* neural progenitor, and it is nuclear-extract absolute abundance. The values are
   in a figure panel; get them from the Spruijt 2013 supplementary table or by reading the Kloet
   Fig. 1a image.
8. **Hauri et al. 2016 *Cell Rep* 17:583** "A High-Density Map for Navigating the Human Polycomb
   Complexome" — AP-MS + SEC across the whole Polycomb complexome, the natural place for PRC2.1 vs
   PRC2.2 partitioning. No PMCID; a bioRxiv version exists (EuropePMC PPR7022) whose full text the
   REST API would not serve. Worth retrieving by hand.
9. **Youmans 2021 Table S1** — iBAQ values normalised to SUZ12 for every PRC2 accessory subunit in
   wild-type hiPSC. Would give a hiPSC PRC2.1:PRC2.2 ratio directly comparable to Smits' HeLa
   numbers.
10. **Javasky 2018 Fig. 2B / supplement** — absolute occupancy percentages of H3K27me1/2/3 in HeLa,
    interphase *and* mitosis, measured against synthetic peptide standards. This is the cleanest
    available answer to Q7 and the only cell-cycle-resolved one.
11. **Ferrari et al. 2014 *Mol Cell*** and **Peters et al. 2003 *Mol Cell*** — the actual primary
    H3K27me0/1/2/3 percentages. Neither is in PMC; cell.com is Cloudflare-protected. Try an
    institutional proxy.
12. **Hein 2015, Schwanhäusser 2011 (NIH3T3 — the only mouse absolute proteome), Sharma 2015
    (mouse brain regions incl. cerebellum)** — all would give independent copy numbers, all behind
    non-OA supplementary tables.
13. **Leroy et al. 2013 *Epigenetics Chromatin* 6:20 (PMC3710262), Additional file 1** — actual
    relative abundances of every H3K27/K36 combinatorial state across 24 human cancer cell lines.
    OA and downloadable; would give cell-type variance for Q7.

**Methodological caution to carry forward**

- Always state which of TOTAL / NUCLEAR / CHROMATIN-BOUND a number is. For PRC2 the total→
  chromatin-bound factor is **5x** (Youmans 20%), which is larger than most of the measurement
  disagreements being argued about.
- iBAQ/proteomic-ruler numbers carry ~1.5x intrinsic error (Wiśniewski 2014) but 2–3x real-world
  error against SRM (RAD21 1.4x, SMC3 3.0x here). Quote PRC2 copy numbers to one significant
  figure with a 3x band.
