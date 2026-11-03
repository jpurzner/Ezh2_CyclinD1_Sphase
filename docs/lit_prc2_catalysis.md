# PRC2 catalytic enzymology — verified kinetic constants

Compiled 2026-10-23. Literature census of primary measurements of PRC2 catalysis.

## Standard applied

Every number below was read in the **Results, Tables, Figure legends or Methods of the primary
paper**. Nothing is taken from an abstract, a review, or a search snippet. Where a tabulated value
exists only inside a rendered figure panel, that is stated explicitly in the Notes column and the
grade is capped at B.

**Confidence grades**

| Grade | Meaning |
|---|---|
| A | Direct measurement stated verbatim in the primary text / table |
| B | Primary text but model-derived (Michaelis–Menten or global fit), **or** a value read off a figure-panel image in the primary paper |
| C | Secondary / second-hand (review, citing paper, model input) |
| X | Searched for, not found |

**Reading note.** Essentially every `kcat` in this field is a Michaelis–Menten nonlinear-regression
output, so most kinetic rows are grade B by construction. That is not a criticism of the data; it is
a statement about what kind of number it is. Grade A is reserved for things actually measured at a
point — fold-stimulation from densitometry, ITC K_d, single-point apparent turnover.

**Substrate is flagged in every row.** This is the single largest source of error in this
literature.

---

## 1. Master table — PRC2 kinetic constants

### 1a. kcat and K_M on NUCLEOSOMAL substrates

This is the section that was previously missing. Five independent primary datasets exist.

| # | kcat (min⁻¹) | K_M | kcat/K_M | Assay | Substrate | Enzyme construct | Species | Buffer / T | [SAM] | n | Direct/fit | Source | Grade | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N1 | **5.71 ± 0.46 ×10⁻³** (H3K36me0) | 89 ± 20 nM (per 601 site) | 6.4 ± 1.4 ×10⁻⁵ min⁻¹ nM⁻¹ ≈ 1.1×10³ M⁻¹s⁻¹ | ³H-SAM, P81 filter, scintillation | **12-mer 601 nucleosome array**, H3K36 methyl-lysine analogue | PRC2 **5-mer**: FLAG-EZH2 + EED + SUZ12 + RbAp48 + **AEBP2** (Sf9/MultiBac) | human | 50 mM HEPES pH 7.5, 0.5 mM MgCl₂, 1 mM DTT, 30 °C, 45 min, 25 nM PRC2 | 1 µM (**≈1.6× SAM K_M — NOT saturating**) | 3 (SEM) | fit | Jani 2019 PNAS 116:8295, **Table 1** | B | Table 1 header is "kcat, min⁻¹ × 10⁻⁵". No exogenous H3K27me3 peptide for array assays |
| N2 | 1.36 ± 0.19 ×10⁻³ (H3K36me2) | 44 ± 22 nM | 3.1 ± 1.5 ×10⁻⁵ | " | " | " | " | " | " | 3 | fit | Jani 2019, Table 1 | B | K_M error is ~50% |
| N3 | 1.15 ± 0.09 ×10⁻³ (H3K36me3) | 63 ± 15 nM | 1.8 ± 0.4 ×10⁻⁵ | " | " | " | " | " | " | 3 | fit | Jani 2019, Table 1 | B | me0/me3 kcat ratio = 4.97; authors state activation coefficient "γ = 5" |
| N4 | **9.6 ×10⁻⁴** (1.6 ± 0.2 ×10⁻⁵ s⁻¹), unmodified | 120 ± 60 nM | 1.4 ± 0.7 ×10² M⁻¹s⁻¹ | ³H-SAM scintillation | **12-mer 601 array, 30 bp linker** (Ar-3) | PRC2 **core 4-mer**: EZH2/EED/SUZ12/RBBP4 (MultiBac, insect) | **not stated in main text** | 30 °C, 40 min, 15 nM PRC2, 10 µL; buffer in SI only | 0.66 µM, **not stated saturating** | 3–4 | fit | Ge, Jani, Diehl, Müller, Muir, JACS 2019;141:15029, **Table 1** | B | Paper writes **kcat,app** throughout, never kcat |
| N5 | 5.58 ×10⁻³ (9.3 ± 0.5 ×10⁻⁵ s⁻¹), **H3K27me1 array** | 230 ± 30 nM | 4.0 ± 0.6 ×10² M⁻¹s⁻¹ | " | 12-mer H3K27me1 array (Ar-9) | PRC2 core 4-mer | " | " | " | 3–4 | fit | JACS 2019, Table 1 | B | me1→me2 step is ~5.8× faster than me0→me1 on arrays |
| N6 | **N.D.** (me2 array) | N.D. | 3.6 ± 1.1 ×10¹ M⁻¹s⁻¹ | " | 12-mer H3K27me2 array (Ar-10) | PRC2 core 4-mer | " | " | " | 3–4 | fit | JACS 2019, Table 1 | B | "we were unable to achieve substrate saturation for the H3K27me2 substrates" — **the me2→me3 step could not be saturated on nucleosomes** |
| N7 | 7.2 ×10⁻³ (1.2 ± 0.1 ×10⁻⁴ s⁻¹) | 62 ± 23 nM | 2.0 ± 0.8 ×10³ M⁻¹s⁻¹ | " | 12-mer unmodified array | **PRC2 core 4-mer + full-length JARID2** (MS-confirmed K116me2/3) | " | " | " | 3–4 | fit | JACS 2019, Table 1 | B | JARID2: kcat,app ×7.5, K_M ÷1.9, kcat/K_M ×14 vs core |
| N8 | 2.4 ×10⁻³ (4.0 ± 0.1 ×10⁻⁵ s⁻¹) | 30 ± 6 nM | 1.3 ± 0.2 ×10³ M⁻¹s⁻¹ | " | 12-mer **H3K4me3** array (Ar-11) | PRC2 core + JARID2 | " | " | " | 3–4 | fit | JACS 2019, Table 1 | B | Core alone on Ar-11: kcat/K_M = 1.6 ×10¹, i.e. **8.5-fold inhibited**; with JARID2 only 1.5-fold |
| N9 | **7.10 ×10⁻³** (0.426 h⁻¹) | 141 ± 31 nM | 30.2 ± 6.9 ×10⁻⁴ h⁻¹ nM⁻¹ | ³H-SAM, glass-fibre filter plate | **chicken erythrocyte oligonucleosomes** (native) | PRC2 **5-mer**: EZH2 + EED + SUZ12 + RbAp48 + **AEBP2** | human enzyme / **avian substrate** | 20 mM Bicine pH 7.6, 1 mM DTT, 100 mM KCl, **temperature not stated** | 2.0 µM (SAM K_M 208–304 nM → **saturating**) | 6 reactions/timepoint | fit | Sneeringer 2010 PNAS 107:20980, **Table 1** | B | Native substrate = heterogeneous starting methylation state |
| N10 | **0.11 ± 0.01** | 0.08 ± 0.00 µM (80 nM) | 1.30 ± 0.08 µM⁻¹min⁻¹ | ³H-SAM 96-well filter plate | **HeLa mono- and dinucleosomes** (native) | PRC2 **5-mer**: EZH2 + EED1 + SUZ12 + RBBP4 + **AEBP2** | human enzyme / human native substrate | 50 mM Tris pH 8.0, 2 mM MgCl₂, 4 mM DTT, 0.002% Tween, **ambient temperature** | fitted axis; SAM K_M 0.13 ± 0.03 µM | "replicate experiments", n not stated | fit (two-substrate sequential global fit) | McCabe 2012 PNAS 109:2989, **Table S2** | B | Highest nucleosomal kcat in the WT set. Nucleosome K_M is **68× tighter** than the peptide K_M in the same paper |
| N11 | **0.0417** (2.5 ± 0.06 h⁻¹) | not determined | — | ³H-SAM filter assay | "human nucleosome" (5 µM), **mono- vs oligo- not specified** | PRC2 **trimer**: EZH2/EED/SUZ12 | human | 20 mM Tris pH 8.0, 5 mM DTT, 0.01% Triton X-100, **23 °C**, 30 min, 200 nM PRC2 | 2 µM | 3 (SD) | **single fixed [S] apparent turnover, not a MM fit** | Huang 2019 Sci Rep 9:6524, **Table 1** (DMSO control row) | B | Paper self-contradicts on the paired peptide value (text 23 h⁻¹ vs Table 24 h⁻¹) |
| N12 | **0.0336** (5.6 ± 0.24 ×10⁻⁴ s⁻¹) | 111 ± 15.8 nM | 5.1 ± 0.7 ×10⁻⁶ s⁻¹nM⁻¹ ≈ 5.1×10³ M⁻¹s⁻¹ | **MTase-Glo** luminescent SAH | **chromatinised 3.6 kbp ATOH1 genomic DNA** (not a 601 array), **conc. per octamer** | PRC2 **core 4-mer**: EZH2/SUZ12/EED/RBBP4 (all full length) | human | 50 mM Tris pH 8.0, 0.5 mM MgCl₂, 0.1% Tween-20, 5 mM DTT, **35 mM KCl**, 30 °C, 120 min, 50 nM PRC2 | 25 µM, **not stated saturating, no SAM K_M** | 3 (SD; fit SE on constants) | fit | Zhang et al. (Davidovich), Nat Commun 2021;12:4592, **Fig. 7c inset table** | B | **Read off the figure image**; not in text, caption or Source Data. Cross-checks against text-stated K_M values and the "7-fold" statement |
| N13 | **N.D.** | N.D. | N.D. | " | same ATOH1 chromatin, **100 mM KCl** | PRC2 core 4-mer | human | as N12 but 100 mM KCl | 25 µM | — | — | Nat Commun 2021, Fig. 7c | X | **The physiologically salted core-complex baseline could not be measured.** Real hole in the dataset |
| N14 | **0.240** (4.0 ± 1.3 ×10⁻³ s⁻¹) | 224 ± 20.6 nM | 17.9 ± 1.7 ×10⁻⁶ s⁻¹nM⁻¹ | " | same ATOH1 chromatin, 35 mM KCl | **PRC2 core + PALI1 "PIR" (aa 1058–1250)** = 5-mer | human | as N12, 15 nM enzyme | 25 µM | 3 | fit | Nat Commun 2021, Fig. 7c | B | kcat values read from figure image; K_M values confirmed in running text |
| N15 | **0.474** (7.9 ± 4.7 ×10⁻³ s⁻¹) | 392 ± 57.8 nM | 20.1 ± 3.2 ×10⁻⁶ s⁻¹nM⁻¹ | " | same, **100 mM KCl** | PRC2 core + PALI1 PIR | human | " | 25 µM | 3 | fit | Nat Commun 2021, Fig. 7c | B | **Highest verified PRC2 kcat on a nucleosomal substrate: 0.47 min⁻¹ ≈ 28 methyl transfers/enzyme/hour** |
| N16 | **0.0708** (1.18 ± 0.03 ×10⁻³ s⁻¹) | 37.7 ± 3.4 nM | 3.1 ± 0.3 ×10⁻⁵ nM⁻¹s⁻¹ ≈ 3.1×10⁴ M⁻¹s⁻¹ | MTase-Glo | same ATOH1 chromatin, **conc. per NCP equivalent** | **6-mer**: EZH2 + SUZ12(70–685) + EED(80–411) + RBBP4 + **AEBP2ᴸ (WT)** + **JARID2** | human | 50 mM Tris pH 8.15, 100 mM KCl, 0.5 mM MgCl₂, 0.1% Tween-20, 5 mM DTT, 30 °C, 60 min, 20 nM enzyme | 25 µM, not stated saturating | 3 (SD; fit SE) | fit | *Auto-inhibition of PRC2 by the long isoform of AEBP2*, EMBO J 2025, **Fig. 6I** | B | Read off figure image; cross-checked against the deposited Source Data rates. **Note SUZ12 and EED are truncated here** |
| N17 | 0.143 / 0.117 / 0.151 (AEBP2ᴸ acidic-tract mutants mt1K2K / mt1A2A / mt1Δ2Δ) | 34.8 / 24.9 / 30.6 nM | 6.9 / 7.8 / 8.2 ×10⁻⁵ nM⁻¹s⁻¹ | " | " | as N16 with mutant AEBP2ᴸ | human | " | 25 µM | 3 | fit | EMBO J 2025, Fig. 6I | B | Text: "All the mutants increased the catalytic efficiency of PRC2 by 2.2- to 2.6-fold" — effect is **almost entirely kcat**, K_M barely moves |

### 1b. kcat and K_M on FREE PEPTIDE substrates (for contrast — do not mix with 1a)

| # | Step | kcat (min⁻¹) | K_M or K₁/₂ | Assay | Substrate | Enzyme | Species | Buffer / T | [SAM] | n | Source | Grade |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P1 | me0→me1 | **0.0140** (0.840 h⁻¹) | **K₁/₂** 184 ± 10 nM | ³H-SAM streptavidin Flashplate | **H3 21–44 peptide**, C-term biotin | PRC2 5-mer (+AEBP2) | human | 20 mM Bicine pH 7.6, 1 mM DTT, **T not stated** | 2.0 µM (saturating) | 6 rxn/pt | Sneeringer 2010, Table 1 | B |
| P2 | me1→me2 | 0.0109 (0.654 h⁻¹) | K₁/₂ 436 ± 42 nM | " | H3K27me1 21–44 peptide | " | " | " | " | " | Sneeringer 2010, Table 1 | B |
| P3 | **me2→me3** | **0.00100** (0.060 h⁻¹) | K₁/₂ 178 ± 16 nM | " | H3K27me2 21–44 peptide | " | " | " | " | " | Sneeringer 2010, Table 1 | B |
| P4 | me0→me1 | **1.01 ± 0.39** | 5.40 ± 1.17 µM | ³H-SAM filter plate | H3 21–44 peptide | PRC2 5-mer (+AEBP2) | human | 50 mM Tris pH 8.0, 2 mM MgCl₂, **ambient T** | fitted; SAM K_M 0.81 ± 0.47 µM | not stated | McCabe 2012, Table S2 | B |
| P5 | me1→me2 | 0.57 ± 0.16 | 6.07 ± 1.37 µM | " | " | " | " | " | SAM K_M 0.79 ± 0.28 | " | McCabe 2012, Table S2 | B |
| P6 | **me2→me3** | **0.03 ± 0.01** | 1.99 ± 0.44 µM | " | " | " | " | " | SAM K_M 0.34 ± 0.13 | " | McCabe 2012, Table S2 | B |
| P7 | me0→me1 | 0.0800 (4.80 ± 0.20 h⁻¹) | K₁/₂ 157 ± 12 nM | ³H-SAM SPA | H3 21–44 peptide | PRC2 **4-mer, no AEBP2** | human | by reference to Sneeringer; **T not stated** | **1.5 µM ≈ SAM K_M → half-saturating** | 3 | Wigle 2011 FEBS Lett 585:3011, Table 1 | B |
| P8 | me1→me2 | 0.0555 (3.33 ± 0.21 h⁻¹) | K₁/₂ 337 ± 26 nM | " | " | " | " | " | " | 3 | Wigle 2011, Table 1 | B |
| P9 | **me2→me3** | 0.0180 (1.08 ± 0.04 h⁻¹) | K₁/₂ 144 ± 11 nM | " | " | " | " | " | " | 3 | Wigle 2011, Table 1 | B |
| P10 | overall (unmod H3) | **2.53 ± 0.21** | 0.84 ± 0.21 µM | **LC-MS/MS of SAH** with d4-SAH standard | **H3 1–45-biotin peptide** | PRC2 4-mer, **Drosophila** (inferred from context, not stated in the Fig. 3A legend) | *D. melanogaster* | 20 mM Tris pH 8.5, 0.5 mM DTT, 0.1% BSA, 7 nM PRC2, **T not stated** | **2 µM vs SAM K_M 5.42 µM → 0.37× K_M, sub-saturating** | **not stated** | Schmitges 2011 Mol Cell 42:330, Results text + Fig. 3A | B |
| P11 | overall (**H3K4me3** peptide) | **0.32 ± 0.08** | 0.36 ± 0.07 µM | " | H3K4me3 1–45 peptide | " | " | " | SAM K_M 10.04 ± 1.56 µM | not stated | Schmitges 2011 | B |
| P12 | apparent, single [S] | **0.400** (24 ± 1 h⁻¹) | — | ³H-SAM SPA | H3 21–44 peptide, 1 µM | PRC2 **trimer** EZH2/EED/SUZ12 | human | 20 mM Tris pH 8.0, 5 mM DTT, **23 °C**, 20 nM PRC2, 1 h | 2 µM | 3 (SD) | Huang 2019 Sci Rep 9:6524, Table 1 | B |

### 1c. SAM K_M

| Value | Enzyme | Substrate context | Source | Grade |
|---|---|---|---|---|
| **0.62 ± 0.2 µM** (unmod array), 0.67 ± 0.2 µM (H3K36me3 array) | PRC2 core + JARID2 | 12-mer 601 nucleosome array, 400 nM | Jani 2019 PNAS, Table 1 footnote †, Fig. 2F | B |
| 0.208 ± 0.050 – 0.304 ± 0.064 µM (range across WT + four Y641 mutants) | PRC2 5-mer | peptide | Sneeringer 2010, Results text | B |
| WT me0 peptide **0.81 ± 0.47 µM**; WT **nucleosome 0.13 ± 0.03 µM**; whole-table range 0.13–3.08 µM | PRC2 5-mer | peptide and HeLa nucleosomes | McCabe 2012, Table S2 | B |
| ~1.5 µM (WT), ~4 µM (Y641C) | PRC2 **4-mer** | peptide | Wigle 2011, Discussion | B |
| 5.42 ± 0.65 µM (unmod H3), **10.04 ± 1.56 µM (H3K4me3)** | *Drosophila* PRC2 4-mer | H3 1–45 peptide | Schmitges 2011, Results | B |
| ~0.5 µM (WT) → ~2.9 µM (SUZ12-S583A) | PRC2 5-mer | peptide | Zhang 2022 Nat Commun 13:6671 (CK2/SUZ12), Results text | B |

**Conclusion on SAM K_M:** no consensus value. Verified range spans ~20-fold (0.13 – 10 µM). The two
nucleosome-context values are the lowest (0.13 and 0.62 µM). The two MTase-Glo nucleosome-array
papers (Davidovich lab) use 25 µM SAM and report **no SAM K_M at all** — so their kcat values are
formally apparent.

**kcat/K_M for SAM.** Only McCabe 2012 Table S2 reports it directly: WT 5-mer PRC2, **0.87 ± 0.12
µM⁻¹min⁻¹ on HeLa nucleosomes**, 1.46 ± 0.36 µM⁻¹min⁻¹ for the me0→me1 peptide step. Everything else
in the field reports kcat/K_M against the histone substrate, not SAM. Grade B.

---

## 2. Allosteric stimulation by H3K27me3 — PINNED

Four primary measurements, all on nucleosomal substrates. The "up to 3-fold" figure you were
carrying is **too low**; the primary literature says 7–12-fold.

| Fold | kcat vs K_M | Substrate | Activator and how delivered | Enzyme | Assay | Source | Grade |
|---|---|---|---|---|---|---|---|
| **~7-fold** | **kcat.** "in both cases the substrate concentration required to achieve the half-maximal reaction rate (K_M) is not significantly affected by the incubation with peptides" (Fig. 4B SAM titration, Fig. 4C nucleosome titration) | **recombinant oligonucleosomes** | free H3K27me3 peptide in *trans*, 10 and 40 µM; half-max stimulation **30–40 µM**, plateau ~100 µM | recombinant PRC2-Ezh2 | ³H-SAM HKMT, gel/autoradiography | **Margueron 2009 Nature 461:762, Fig. 4A + Results** | **A** |
| **~12-fold (PRC2-EZH2)**; **~2-fold (PRC2-EZH1)** | not resolved | **core nucleosomes**, 300 nM | H3K27me3 peptide in *trans* | human PRC2-EZH1 or -EZH2 4-mer (EZH1/2+EED+SUZ12+RBAP48), 15–60 nM | ³H-SAM, 50 mM Tris pH 8.5, 5 mM MgCl₂, 4 mM DTT, **500 nM SAM**, 60 min, **30 °C**, n=3 | **Lee CH, Holder M, … Reinberg, Armache. Mol Cell 2018;70:435, Fig. 2A** | **A** |
| ~8.8-fold proximal, ~3.2-fold distal (vs H3K27R control) | not resolved | **12-mer 601 array** with a pre-installed H3K27me3 tetranucleosome block **in cis** | in *cis*, same array | PRC2 core 4-mer | ³H-SAM, 0.5 mM Mg²⁺, 30 °C | Ge/Muir JACS 2019, **Fig. 3d** | **B** (figure-image read) |
| qualitative only | — | mononucleosomes | — | PRC2-AEBP2-JARID2 | western | Cookis/Nogales NSMB 2025;32:393 | C for magnitude |

**Three things to carry forward:**

1. **The allosteric effect is on kcat, not K_M.** Margueron 2009 states this explicitly for
   both SAM and nucleosome titrations. Ge/Muir reach the same conclusion for the JARID2 analogue:
   activation acts "primarily through an increase in kcat,app".
2. **Stimulation is in *cis* and short-range.** Ge/Muir: "spreading is restricted to between 1–4
   nucleosomes in our experimental system"; mixing separately-modified arrays in *trans* gave no
   stimulation at 0 or 1.5 mM Mg²⁺. They frame Ar-6 vs Ar-7 as "short-range (0.1–0.7 kbp) versus
   medium-range (≥1 kbp)".
3. **Occupancy half-max requires 30–40 µM free H3K27me3 peptide** (Margueron 2009), consistent with
   the weak EED K_d (below). In a cell this is only achievable by *cis* presentation on a
   neighbouring nucleosome — which is exactly the geometric result.

**EED affinities for allosteric ligands** (ITC, ΔEED 77–441, *direct* measurements, grade A):

| Ligand | K_d | Source |
|---|---|---|
| JARID2-K116me3 | **3.4 µM** | Sanulli 2015 Mol Cell 57:769, Fig. 5C |
| JARID2-K116me2 | 10.4 µM | " |
| H3K27me3 | **36.4 µM** | " |
| H3K27me2 | 68 µM | " |
| JARID2 unmodified, JARID2-K116me3-F117A | not detectable | " |
| H3K27me3 (fluorescence anisotropy, EED 40–441) | 41.7 ± 2.8 µM | Nat Commun 2021 (PALI1), Fig. 3a |
| H3K9me3 / H4K20me3 / H1K26me3 / H3K27me3 | 10–45 µM; **~4-fold weaker per methyl group removed** | Margueron 2009, Supp Table S1 |
| H3K4me3, H3K36me3, H3K79me3 | **no appreciable binding** | Margueron 2009 |

Text: "Jarid2-K116me3 peptide exhibited ten times more affinity for ΔEed than H3K27me3 peptide"
(36.4 / 3.4 = 10.7).

---

## 3. Allosteric stimulation by JARID2-K116me3

| Fold | kcat vs K_M | Substrate | Enzyme | Source | Grade |
|---|---|---|---|---|---|
| **"5- to 6-fold"** vs K116A / K116R fragments | not resolved | **recombinant chromatin (oligonucleosomes)** — NOT peptide | PRC2-Ezh2 + equimolar JARID2 109–450 (K116me3 by native chemical ligation) | Sanulli 2015, Results text + **Fig. 5F** | **A** |
| **lost entirely** on native chromatin | — | native chromatin | same | Sanulli 2015, Fig. 5F | A |
| **~10-fold** vs PRC2 core alone (K116me2 ~8-fold; unmodified tJARID2 ~50% of the me2/me3 level) | — | **unmodified 12-mer 601 array** | PRC2 core 4-mer + equimolar semisynthetic tJARID2 108/137–450 | Ge/Muir JACS 2019, Fig. 5d, text | **A** for the stated folds |
| **kcat,app ×7.5, K_M ÷1.9, kcat/K_M ×14** | **kcat** | 12-mer unmodified array | PRC2 core vs PRC2 + full-length JARID2 | Ge/Muir JACS 2019, Table 1 (rows N4 vs N7 above) — ratio is my arithmetic from their tabulated values | B |
| ~30-fold max in a JARID2 titration, plateau at **1:1 stoichiometry** | — | 12-mer array Ar-3 | PRC2 core + JARID2 titration | Ge/Muir JACS 2019, Fig. 3f | B (figure-image read) |

**JARID2 overrides the geometric preference.** With JARID2 present "the stimulating effect of having
preinstalled blocks of H3K27me3 is largely lost"; the authors' conclusion is that "PRC2-JARID2 is
optimized for H3K27me3 **establishment**, while the PRC2 core complex efficiently **propagates**".

**Serine phosphorylation gates it.** JARID2 pS120 in the K116me3 context "reduced HMT activity back
nearly to the level associated with … tJARID2-K116R"; pS124 and pS126 did not (Ge/Muir Fig. 5d, grade
A for the statement, B for the bar magnitudes).

**Structural basis** (Cookis/Nogales NSMB 2025;32:393, qualitative, grade A as mechanism): JARID2-K116me3
and H3K4me3 **compete for the same EED regulatory site**; H3K4me3 binds the aromatic cage but fails
to fold the SRM helix, so it is an *allosteric antagonist* rather than a weak agonist.

---

## 4. Accessory subunits — rate vs binding

This is the cleanest-answerable question in the set, because three papers separate the two.

| Subunit | Effect on catalytic rate | Effect on binding | Substrate | Source | Grade |
|---|---|---|---|---|---|
| **JARID2** (full length, K116me2/3) | **kcat,app ×7.5** (9.6×10⁻⁴ → 7.2×10⁻³ min⁻¹) | K_M ÷1.9 (120 → 62 nM) | 12-mer 601 array | Ge/Muir JACS 2019, Table 1 | B |
| **PALI1 (PIR, aa 1058–1250)** | **kcat ×7.1** (0.0336 → 0.240 min⁻¹ at 35 mM KCl) | **K_M gets 2× WORSE** (111 → 224 nM). Net kcat/K_M only ×3.5 | ATOH1 chromatin array | Nat Commun 2021, Fig. 7c | B |
| **PALI1 (PIR)** — equilibrium binding, separate assay | — | **K_d ÷53** for NCP₁₈₂ (1409 → 17.1 nM); ÷53 for NCP₁₄₇ (2560 → 48.5 nM); DNA affinity >25-fold | EMSA, mononucleosome probes | Nat Commun 2021, Fig. 7b + abstract | A |
| **PALI1 separation of function** | K1241A mutant loses allosteric activation **yet still binds nucleosomes within 2-fold of WT K_d**, and is still ~5× more active than core PRC2 | — | — | Nat Commun 2021, Fig. 6b | A |
| **AEBP2 (short isoform)** | stimulates PRC2-EZH1 and PRC2-EZH2 by **similar fold** (magnitude **not stated numerically**); **additive with, and independent of, H3K27me3 allosteric activation** (unaffected by EZH2-P132S, which abolishes allostery) | increases nucleosome binding by EMSA; KR-motif (aa 165–188) required for both | chromatinised plasmids, core nucleosomes | Lee CH 2018 Mol Cell 70:435, Figs. 4A–4F, 5 | A for the mechanism, **X for the fold** |
| **AEBP2 long isoform (AEBP2ᴸ)** acidic tracts | **auto-inhibitory**: removing/charge-swapping them raises **kcat 1.65–2.13×** and kcat/K_M 2.2–2.6× | K_M barely moves (37.7 → 24.9–34.8 nM); CpG-DNA K_d improves >15-fold; nucleosome K_d improves only ≤2-fold | ATOH1 chromatin array | EMBO J 2025, Fig. 6I + text | B (values), A (folds stated in text) |
| **AEBP2ᴸ — is the active site itself inhibited?** | **No.** TEV cleavage of the AEBP2ᴸ N-terminal region gave "the same HMTase activity irrespective of whether the N-terminal region was cleaved or not" on an **H3 tail peptide** substrate. "This data indicates that the catalytic site of PRC2 is not inhibited" | — | H3 peptide | EMBO J 2025, Fig. 6D | A |
| **AEBP2 long vs short isoforms** | long isoform's DE-motifs inhibit both H3K27 methylation **and EZH2 automethylation**; effect assigned to catalysis, not binding ("may modestly attenuate, but does not abolish, nucleosome association") | both isoforms enhance nucleosome binding; short slightly stronger | 601 di-/mono-nucleosomes | Li, … Lee C-H. Nucleic Acids Res 2026 (PMC13347268) | A for the statements; **no kcat/K_M anywhere in that paper** |
| **PHF1** | not separated: "increases the binding affinity and residence time of PHF1-PRC2 on nucleosomes about two- to three-fold, resulting in more efficient H3K27 methylation" | 2–3× | mononucleosomes | Finogenova 2020 eLife, citing Choi 2017 | C (restated from Choi 2017; not re-measured) |
| **MTF2 / PHF19 / EPOP** | — | — | — | — | **X — no primary kcat or K_M found** |

**Summary of the rate-vs-binding question:** for JARID2, PALI1 and AEBP2ᴸ, the measured effect on
Michaelis–Menten parameters is **predominantly kcat, not K_M**. PALI1 is the sharpest case — it
improves nucleosome K_d by >50-fold in EMSA while making K_M *worse* on an array substrate, and still
gains 7-fold in kcat. Equilibrium K_d and kinetic K_M are decoupled in this system, and papers that
quote only one of them are not interchangeable. The EMBO J authors note this explicitly, citing
"a substantial defect in HMTase but only a moderate change in the apparent affinity for nucleosome
probes during EMSA".

---

## 5. Inhibition by H3K4me3 and H3K36me2/3

### H3K4me3

| Effect | kcat vs K_M | Substrate | Enzyme | Source | Grade |
|---|---|---|---|---|---|
| **kcat ÷7.9 (reported as 8-fold): 2.53 → 0.32 min⁻¹** | **kcat.** Peptide K_M actually *tightens* (0.84 → 0.36 µM); **SAM K_M roughly doubles** (5.42 → 10.04 µM) | **free H3 1–45 peptide** | *Drosophila* PRC2 4-mer | Schmitges 2011, Results + Fig. 3A | B |
| kcat/K_M ÷8.5 (1.4×10² → 1.6×10¹ M⁻¹s⁻¹) | — | **12-mer 601 nucleosome array** | PRC2 core 4-mer | Ge/Muir JACS 2019, Table 1 | B |
| only ÷1.5 (2.0×10³ → 1.3×10³ M⁻¹s⁻¹) | kcat ÷3.0, **K_M ÷2.1 (tighter)** | 12-mer array | **PRC2 + JARID2** | Ge/Muir JACS 2019, Table 1 | B |
| mechanism: binds EED aromatic cage but **fails to fold the SRM** → allosteric antagonist competing with JARID2-K116me3 | — | mononucleosomes | PRC2-AEBP2-JARID2 | Cookis/Nogales NSMB 2025 | A (mechanism), X (fold) |

**Caveat on Schmitges' quoted specificity constant.** The printed kcat/K_M values (7.8×10³ →
0.53×10³ M⁻¹s⁻¹, "~15-fold") reproduce exactly only if computed against the **SAM** K_M, not the
peptide K_M. Against the peptide K_M the drop is ~3.4-fold. The paper does not say which K_M it used.
Use the kcat ratio (8-fold) as the citable number.

### H3K36me2 / H3K36me3

| Effect | kcat vs K_M | Substrate | Enzyme | Source | Grade |
|---|---|---|---|---|---|
| **kcat ÷4.97 (me3), ÷4.2 (me2); K_M essentially unchanged** (89 → 63 nM) | **pure kcat.** Verbatim: "the difference in PRC2 activity … is primarily a kcat effect". Authors state activation coefficient **γ = 5** | 12-mer 601 array, H3K36 methyl-lysine analogue | human PRC2 5-mer (+AEBP2) | Jani 2019 PNAS, Table 1 + Results | B |
| **">10-fold inhibited"** for both H3K27me1 and H3K27me3 formation; authors' own summary: **"about 10- to 20-fold less effective"** | not resolved (no kinetic fit) | **H3Kc36me3 mononucleosomes, 446 nM** (*Drosophila* nucleosomes) | full-length **human PRC2**, and separately full-length **PHF1-PRC2** | Finogenova 2020 eLife 9:e61964, Fig. 3C + author response | **A** for the ">10-fold"; the 10–20-fold refinement is in the peer-review response |
| ~5-fold inhibition | not resolved | **H3K36R and H3K36A mononucleosomes** | human PRC2 | Finogenova 2020, Fig. 3C + author response | A |
| **no inhibition on free peptide** — "no significant difference" in H3K27me1 formation on H3 18–42 vs H3 18–42 K36me3 (Welch's t-test, technical triplicate, ESI-MS) | — | **free H3 18–42 peptide** | human PRC2 | Finogenova 2020, Fig. 3D | **A** |
| **no inhibition on free peptide** — "no decrease in enzymatic turnover is observed" with H3 1–45 K36me3 | — | free H3 1–45 peptide | *Drosophila* PRC2 | Schmitges 2011, Fig. S6D legend | **A** |
| **no kcat or K_M for H3K36me2/3 exists in Schmitges 2011 at all** | — | — | — | — | X (genuine absence, not an access failure) |

**This is the cleanest substrate-dependence result in the whole census:** H3K36me3 inhibition is
**entirely nucleosome-dependent** and vanishes on free peptide, in two independent papers from two
labs on two species. Finogenova's framing: "The allosteric inhibition of PRC2 by H3K36me3 therefore
only occurs in the context of the geometric constraints of the nucleosome." Binding is unaffected —
PRC2 binds H3Kc36me3 and unmodified mononucleosomes with indistinguishable affinity (Fig. 3A,B).

Finogenova HMTase conditions: 446 nM mononucleosomes or 223 nM dinucleosomes, 20 mM HEPES pH 7.8,
50 mM NaCl, 2.5 mM MgCl₂, 5% glycerol, 0.25 mM EDTA, 0.5 mM DTT, **80 µM SAM**, 90 min at RT,
western-blot densitometry, n = 3, SEM.

---

## 6. Methylation events per PRC2 per unit time; processivity

**Processivity / distributivity: grade X.** No primary quantitative measurement found in any paper
read. Searched in Jani 2019, Ge/Muir JACS 2019, Sanulli 2015, Oksuz 2018, Schmitges 2011,
Finogenova 2020, Lee 2018, Margueron 2009, and both Davidovich papers. The strings "processiv" /
"distributiv" occur in only one of them — Schmitges 2011, and only as Discussion speculation
("limiting processivity of the enzyme"), with no coefficient and no assay. Ge/Muir explicitly ran
their assays under conditions where "only single methylation state turnover … is observed", which is
an assay design choice made *to avoid* the processivity question, not an answer to it.

**Methylation events per enzyme per unit time** must therefore be derived from kcat. Converting the
nucleosomal rows of §1a:

| Source | Enzyme / substrate | kcat (min⁻¹) | **methyl transfers per PRC2 per hour** |
|---|---|---|---|
| Ge/Muir JACS 2019 | core 4-mer, unmodified 12-mer array | 9.6×10⁻⁴ | **0.058** |
| Jani 2019 | 5-mer +AEBP2, unmodified 12-mer array | 5.7×10⁻³ | 0.34 |
| Ge/Muir JACS 2019 | +JARID2, unmodified array | 7.2×10⁻³ | 0.43 |
| Sneeringer 2010 | 5-mer, chicken oligonucleosome | 7.1×10⁻³ | 0.43 |
| Nat Commun 2021 | core 4-mer, ATOH1 chromatin, 35 mM KCl | 0.0336 | 2.0 |
| Huang 2019 | trimer, "human nucleosome" | 0.0417 | 2.5 |
| EMBO J 2025 | 6-mer AEBP2ᴸ-JARID2, ATOH1 chromatin, 100 mM KCl | 0.0708 | 4.2 |
| McCabe 2012 | 5-mer, HeLa mono/dinucleosome | 0.11 | 6.6 |
| Nat Commun 2021 | +PALI1 PIR, ATOH1 chromatin, 100 mM KCl | 0.474 | **28** |

**Range: 0.06 to 28 methyl transfers per PRC2 per hour — a ~480-fold spread.** This is the honest
answer to "how fast is PRC2 on chromatin". If you need a single number for a model, the defensible
move is to pick one paper whose enzyme composition and salt match your biology and carry its error,
not to average across the table.

A third-party model input worth knowing about and **not** citing as a measurement: a 2023 PNAS paper
on PRC2/G-quadruplex RNA sets `kcat = 10⁻¹ s⁻¹` (= 6 min⁻¹) as a fixed simulation parameter, stating
it "was taken from prior PRC2 literature". That is **~13× above the fastest verified nucleosomal
kcat** in this table and ~6000× above the slowest. Grade C; do not propagate.

### Closest thing to an in vivo rate

Oksuz, Narendra, Lee, … Reinberg. Mol Cell 2018;70:1149 (PMC7700016), mouse ESCs, EED-knockout cells
with 4-OHT-inducible WT EED, TetR-EZH2 tethered to a 5×TetO cassette:

- "Recruitment of PRC2 and expression of WT EED at 24 hr resulted in the formation of H3K27me3
  domains of **~5 kb** in size at both artificial recruitment sites" (Figs. S6C–D). Grade A as a
  domain size.
- **The paper never states a rate.** ~5 kb / 24 h ≈ 0.2 kb/h ≈ 1 nucleosome/h is *my* arithmetic from
  a single timepoint and assumes linear growth. Grade C, upper-bound order of magnitude only.
- Timeline, grade A: H3K27me2 precedes H3K27me3; at 12 h the mark is at nucleation sites only, by
  36 h it has reached spreading sites and me2→me3 conversion has occurred at nucleation sites.
- Global stoichiometry by quantitative histone MS: H3K27me3 is **20% of H3 in WT mESCs**, ~1% in EED
  cage mutants, 0% in EED KO (Fig. 1B). Grade A.

---

## 7. EZH1 vs EZH2

| Observation | Value | Substrate | Source | Grade |
|---|---|---|---|---|
| Basal activity | **similar** for PRC2-EZH1 and PRC2-EZH2 in the absence of activator | core nucleosomes | Lee 2018 Mol Cell 70:435, Fig. 2A lanes 1 and 5 | A |
| **Allosteric activation by H3K27me3** | **EZH2 ~12-fold; EZH1 only ~2-fold** | core nucleosomes (also reproduced on chromatinised plasmids and dinucleosomes) | Lee 2018, Fig. 2A | **A** |
| Linker-DNA preference | PRC2-EZH2 on **40 bp linker** dinucleosomes is **5–6× more active** than on 30 bp linker; corroborated on 12N oligonucleosomes | dinucleosomes / 12N oligonucleosomes, 300 nM | Lee 2018, Fig. 1D–E | A |
| Response to free DNA, EZH1 vs EZH2 | differential (PRC2-EZH1 shows marginal stimulation on 40 bp linker dinucleosomes) | dinucleosomes | Lee 2018, Figs. 1D, S2B | A |
| AEBP2 stimulation | **similar fold for both** paralogues | chromatinised plasmids | Lee 2018, Fig. 4A | A (direction), X (magnitude) |
| Nucleosome binding | AEBP2 enhances binding of both; the enhancement is **smaller for EZH1 because PRC2-EZH1 already has higher intrinsic nucleosome binding** | core nucleosomes, 200 nM, EMSA | Lee 2018, Fig. 4F | A |
| SET domain | "EZH1 and EZH2 share a nearly identical SET domain" — structural work on EZH1 is argued to transfer to EZH2 | — | Zhang 2022 Nat Commun 13:6671 (CK2/SUZ12) | A |
| **kcat or K_M for PRC2-EZH1** | — | — | — | **X — no Michaelis–Menten analysis of an EZH1-containing complex found anywhere** |

Lee 2018 assay conditions (all the fold numbers above): 15 µL, 50 mM Tris-HCl pH 8.5, 5 mM MgCl₂,
4 mM DTT, **500 nM ³H-SAM**, 300 nM nucleosomes, 15–60 nM PRC2, **60 min at 30 °C**, autoradiography
with Coomassie loading control, n = 3, mean ± SD. Note the SAM concentration is **500 nM**, at or
below every reported SAM K_M — these are relative-activity measurements, not kcat measurements, and
the enzyme titration was explicitly chosen to stay in the linear range.

---

## 8. me2 → me3 is rate-limiting: the quantitative backing

Three independent datasets, **all on free peptide**, all agreeing in direction:

| Source | Enzyme | kcat/K ratio me0 : me1 : me2 | kcat fall me0 → me2 |
|---|---|---|---|
| Sneeringer 2010 | human 5-mer | **13.4 : 4.4 : 1** | 14-fold |
| McCabe 2012 | human 5-mer | **9 : 6 : 1** (stated in Results) | 34-fold |
| Wigle 2011 | human **4-mer** | 4.1 : 1.3 : 1 | 4.4-fold |

The effect sits in **kcat**, not in substrate binding: peptide K₁/₂ varies less than 3.5-fold across
all three methylation states in Sneeringer. Sneeringer's footnote matters — on peptides the constant
is **K₁/₂, not K_M**, because "the peptidic substrates displayed sigmoidal binding behavior".

**On nucleosomes the step resolution has never been measured.** Every nucleosomal kcat in §1a is
either a single lumped value on a heterogeneously modified native preparation, or (Ge/Muir) a
me0→me1 and me1→me2 pair where **the me2→me3 array could not be saturated** so kcat and K_M were
undeterminable. Ge/Muir's kcat/K_M for the me2 array (3.6 ± 1.1 ×10¹ M⁻¹s⁻¹) is 11-fold below the
me1 array (4.0×10²) and 3.9-fold below the unmodified array — consistent in direction with the
peptide data, and the only nucleosome-substrate evidence that exists for the claim.

---

## 9. Cross-paper comparability — read before quoting any single number

1. **Peptide kcat spans ~70-fold across three papers using the same H3 21–44 substrate**: Sneeringer
   0.014, Wigle 0.080, McCabe 1.01 min⁻¹ for the me0→me1 step. The peptide binding constant spans
   ~30-fold (184 nM to 5.40 µM). There is no single "PRC2 peptide kcat".
2. **Nucleosomal kcat spans ~480-fold** (§6 table). Driven by: enzyme composition (trimer / 4-mer /
   5-mer / 6-mer; which accessory subunit), substrate architecture (601 12-mer array vs native
   oligonucleosome vs chromatinised genomic DNA vs HeLa mono/dinucleosome), salt (35 vs 100 mM KCl
   changes the PALI1 complex's kcat by 2-fold and makes the core complex unmeasurable), pH (7.5 to
   8.5), and readout (³H filter binding vs MTase-Glo luminescence vs LC-MS).
3. **Within a single paper, the peptide-to-nucleosome kcat gap is only 2–10-fold**, not orders of
   magnitude: Sneeringer 2×, McCabe 9×, Huang 9.6×. The orders-of-magnitude discrepancies in this
   field are **between papers, not between substrates**. This is the opposite of the usual
   assumption and is worth stating in any write-up.
4. **K_M is reported on inconsistent bases.** Jani's K_M is "per 601 site"; the Davidovich papers'
   are per octamer / per NCP-equivalent (explicitly stated, and confirmed against their Source Data);
   **Ge/Muir's JACS Table 1 K_M is ambiguous** — the paper says "varying array substrate
   concentrations", so 1.2×10⁻⁷ M is most likely per *array*, a 12-fold difference from the others.
   Do not place Ge/Muir's K_M on the same axis as Jani's without resolving this.
5. **Sub-saturating SAM is endemic.** Jani ran at 1 µM against a SAM K_M of 0.62 µM (~62% Vmax);
   Schmitges at 2 µM against 5.42 µM (~27%); Wigle at exactly the SAM K_M (~50%); Ge/Muir at 0.66 µM
   with no SAM K_M measured (and they correctly write kcat,**app** throughout); the two MTase-Glo
   papers at 25 µM with no SAM K_M reported. Only Sneeringer and McCabe are clearly saturating.
   Most published PRC2 kcat values are therefore systematic underestimates of true Vmax-derived kcat.
6. **Enzyme concentration confounds the PALI1 7-fold claim**: 50 nM core vs 15 nM PRC2-PALI1_PIR.
   kcat normalises for this, but any error in active-fraction determination maps straight onto the
   7-fold.
7. **Species mixing is common and under-reported.** Sneeringer: human enzyme on chicken oligonucleosomes.
   Finogenova: human PRC2 on *Drosophila* nucleosomes. Cookis: human PRC2 on *Xenopus* nucleosomes.
   Davidovich: human PRC2 on *Xenopus*-type recombinant octamers. Ge/Muir: **species not stated in
   the main text at all**.

---

## 10. Corrections to identifiers encountered during this census

Recorded so they are not propagated:

- McCabe 2012 PNAS 109:2989 is **PMC3287005**, not PMC3286974 (which is a different PNAS paper on
  antiangiogenic agents and breast cancer stem cells).
- Oksuz 2018 Mol Cell 70:1149 is **PMID 29932905 / PMC7700016**.
- Ge/Muir JACS 2019 is **141(38):15029–15039**.
- EMBO J 2025 AEBP2ᴸ paper, Fig. 6 panel I caption misreferences its own panels: it cites "(E)" and
  "(G)" where the correct source panels are **(F)** and **(H)**.
- Sanulli 2015 Fig. 5F: the plotted bars imply a K116me3/K116A ratio of ~13×, larger than the
  authors' stated "5- to 6-fold". Cite the authors' number; note the discrepancy.
- Huang 2019 Sci Rep: Results text (23 and 6 h⁻¹) disagrees with Table 1 DMSO rows (24 and 9 h⁻¹).

---

## What is still missing

**Gaps that are genuine absences in the literature, not access failures:**

1. **Step-resolved (me0/me1/me2) kinetics on a nucleosomal substrate.** Does not exist. The
   me2→me3-is-rate-limiting claim rests entirely on free-peptide data from one lab lineage
   (Copeland/Epizyme, three papers). The only nucleosomal hint is Ge/Muir's unsaturable H3K27me2
   array, which tells you the step is slow but cannot give kcat or K_M for it. **If the me2→me3 rate
   on chromatin matters to your model, it is unmeasured.**
2. **Any kcat or K_M for an EZH1-containing complex.** None found. All EZH1 vs EZH2 comparisons are
   relative-activity gels.
3. **Any kcat or K_M for PRC2 with MTF2, PHF19 or EPOP.** None found. PHF1 has a binding/residence-time
   number (2–3× from Choi 2017) but no catalytic constant.
4. **A PRC2 core-complex kcat at physiological salt on a chromatin array.** Reported as **N.D.** in
   Nat Commun 2021 Fig. 7c — the core complex's activity at 100 mM KCl was too low to fit. The most
   physiologically relevant baseline is the one that could not be measured.
5. **Processivity or distributivity, quantitatively.** No coefficient, no single-molecule
   measurement, no assay designed to answer it. Several papers deliberately avoid the regime.
6. **An in vivo spreading rate in kb/h or kb/cell-cycle.** Oksuz 2018 gives one distance-and-time
   pair (~5 kb at 24 h) and reports it as a domain size, not a rate.
7. **A SAM K_M for the two MTase-Glo nucleosome-array papers** (Nat Commun 2021, EMBO J 2025). Both
   run at 25 µM SAM and neither measures or asserts saturation, so both kcat sets are formally
   apparent.
8. **A fold-stimulation by H3K27me3 measured on a nucleosomal substrate with kcat and K_M separated.**
   Margueron 2009 states K_M is unaffected (so the effect is kcat) but gives no fitted constants;
   Lee 2018's 12-fold is a 500 nM-SAM relative-activity measurement. Nobody has run a full
   Michaelis–Menten series ± H3K27me3 on nucleosomes.
9. **AEBP2 (short isoform) fold-stimulation as a number.** Lee 2018 shows it in three figures and
   never states a fold.

**Gaps that are access failures and could still be closed:**

10. **Swalm BM et al., "Convergent evolution of chromatin modification by structurally distinct
    enzymes: comparative enzymology of histone H3 Lys27 methylation by vSET and EZH2." Biochem J
    2013;453:241.** No PMCID, subscription required, portlandpress.com served a Cloudflare
    interstitial. From the title this is a direct comparative-enzymology paper on PRC2 kcat and SAM
    K_M. **Highest-value remaining target — worth pulling through institutional access.**
11. **Swalm BM et al., "Reaction coupling between wild-type and disease-associated mutant EZH2."
    ACS Chem Biol 2014.** No PMCID, subscription required. Likely contains WT/mutant coupling kinetics.
12. **Jani 2019 PNAS SI Appendix.** PNAS returns 403; the PMC supplementary link sits behind a
    JavaScript download gate. Fig. S2 in particular explains *why* peptide, histone and mononucleosome
    substrates failed Michaelis–Menten analysis and 12-mer arrays were required — directly relevant
    to anyone choosing a substrate.
13. **Ge/Muir JACS 2019 supplemental .docx (171 MB).** Contains buffer, pH, salt, species and the
    enzyme-concentration series. **Without it, the species of the enzyme in Table 1 is unknown.**
14. **Nat Commun 2021 Fig. 7c replicate count.** Caption says three independent experiments; the
    deposited Source Data has four repeat columns for three of the four conditions.
15. **NAR 2026 (PMC13347268) supplementary PDF.** Bot-blocked. The assay format (fixed substrate,
    enzyme titration, gel endpoint) cannot yield Michaelis–Menten constants, so "no kinetics" is
    near-certain, but the PDF itself is unverified.

---

## Papers read in full for this census

Primary, with kinetic content:
Margueron et al. Nature 2009;461:762 (PMC3772642) ·
Sneeringer et al. PNAS 2010;107:20980 (PMC3000297) ·
Schmitges et al. Mol Cell 2011;42:330 ·
Wigle et al. FEBS Lett 2011;585:3011 ·
McCabe et al. PNAS 2012;109:2989 (PMC3287005) ·
Sanulli et al. Mol Cell 2015;57:769 (PMC4352895) ·
Lee CH, Holder M, et al. Mol Cell 2018;70:435 (PMC5949877) ·
Lee CH, Yu JR, et al. Mol Cell 2018;70:422 (PMC5935545) ·
Oksuz et al. Mol Cell 2018;70:1149 (PMC7700016) ·
Jani et al. PNAS 2019;116:8295 (PMC6486736) ·
Ge, Jani, Diehl, Müller, Muir. JACS 2019;141:15029 (PMC6941475) ·
Huang et al. Sci Rep 2019;9:6524 (PMC6484020) ·
Finogenova et al. eLife 2020;9:e61964 (PMC7725500) ·
Zhang et al. Nat Commun 2021;12:4592 (PMC8319299) ·
Zhang et al. Nat Commun 2022;13:6671 (PMC9645763) ·
Cookis, Lydecker, Sauer, Kasinath, Nogales. Nat Struct Mol Biol 2025;32:393 (PMC11832421) ·
EMBO J 2025, AEBP2 long isoform (PMC12669776) ·
Li et al. Nucleic Acids Res 2026 (PMC13347268)

Checked and found to contain **no** kinetic constants:
Martin, Cao, Zhang. JBC 2006;281:8365 (zero occurrences of kcat / K_M / kinetic) ·
Poepsel, Kasinath, Nogales. NSMB 2018;25:154 (PMC5805599) ·
Mol Cell 2024, automethylated PRC2 dimerization (PMC11980035)

Checked and rejected as *primary* sources:
PNAS 2023 PRC2/G-quadruplex (PMC10266057) — uses kcat = 10⁻¹ s⁻¹ as a fixed literature-derived model
input · "Forecasting histone methylation by Polycomb complexes with minute-scale precision" 2023
(PMC10745708) — modelling; cites Sneeringer's 13:4:1, measures nothing new.
