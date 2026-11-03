# Census of PRC2-antagonising enzymes: KDM6A/KDM6B, p300/CBP, HDACs

Compiled 2026-10-23. Scope: kinetics, abundance and binding dynamics needed to parameterise the
nucleosome-lattice model of the H3K27me3-depleted core at active promoters.

**Standard applied.** Every number below was read in the Methods, Results, or a numbered Table of the
*primary* paper unless explicitly flagged. Values that exist only in a figure, only in an abstract, or
only as a citation inside someone else's paper are flagged as such and graded down.

## Confidence grades

| Grade | Meaning |
|---|---|
| **A** | Number read in the Methods/Results/Table of the primary paper, with assay, substrate, construct and conditions all stated |
| **B** | Primary paper read, but the number itself lives in a figure, or conditions are only partly stated |
| **C** | Primary source *identified* but full text inaccessible; number known only second-hand |
| **X** | **No primary measurement exists.** Not "I could not find it" — searched and the measurement has not been made |

---

## HEADLINE

**KDM6B kinetics do exist — one primary measurement, on a peptide, with a catalytic-domain construct.**
Walport et al. 2016 *Nat Commun* 7:11974 Table 1 reports kcat = 1.17 × 10⁻² s⁻¹ (0.70 min⁻¹) and
K_M = 6.7 µM for KDM6B(1141–1590) on an H3(14–34)K27me3 peptide. This appears to be the only
published steady-state kcat for KDM6B. It was previously hard to find because the paper is titled and
indexed as a study of *arginine* demethylation; the H3K27me3 measurement is the lysine control.

**But: no KDM6 enzyme has ever been assayed kinetically on a nucleosome.** Not KDM6A, not KDM6B, not
full-length, not catalytic domain. There is also no KDM6–nucleosome structure. Every KDM6 number in
the literature is a free-peptide number from a truncated construct.

**And both KDM6 enzymes sit at ~10³ copies per HeLa cell.** In the deepest available HeLa proteome
(Bekker-Jensen 2017, 14,200 protein groups, iBAQ copy numbers calibrated against established HeLa
standards, log-space Pearson r = 0.95), **KDM6A = 974 and KDM6B = 859 copies/cell** — against HDAC1
at 851,000, HDAC2 at 515,000, EP300 at 37,000 and EZH2 at 31,000. KDM6B is therefore **~1,000× rarer
than HDAC1 and ~36× rarer than EZH2**, and is essentially equimolar with KDM6A. Any model in which
KDM6B is the *dominant* antagonist setting the promoter core is asserting a mechanism carried by an
enzyme that is (i) present at ~10³ copies/cell, (ii) never measured on its physiological substrate,
and (iii) argued in the best in-vivo dataset to work largely *independently* of its demethylase
activity.

---

## Q1. kcat / K_M for KDM6A and KDM6B on H3K27me3

| Enzyme | Construct | Substrate | Parameter | Value | Assay & conditions | n | Source | Grade |
|---|---|---|---|---|---|---|---|---|
| **KDM6B** | 1141–1590 (JmjC catalytic domain), His6, *E. coli* | H3(14–34)K27me3 **peptide** | **kcat** | **(11.7 ± 0.1) × 10⁻³ s⁻¹** = **0.70 min⁻¹** | Formaldehyde-dehydrogenase/NAD⁺-coupled fluorescence; Fe(II)(NH₄)₂(SO₄)₂ 10 µM, Na-ascorbate 100 µM, 2OG 200 µM, NAD⁺ 500 µM, enzyme 0.2–4 µM, FDH 0.025 U, 30 µL, **37 °C**, 25 min | 3 | Walport et al. 2016 *Nat Commun* 7:11974, **Table 1**; PMC4931022 | **A** |
| **KDM6B** | as above | as above | **K_M** | **6.7 ± 0.6 µM** | as above | 3 | same | **A** |
| **KDM6B** | as above | as above | **kcat/K_M** | **≈ 1.75 × 10³ M⁻¹ s⁻¹** | computed by authors (table column mislabelled "× 10⁻⁶"; the printed 1,750.0 reconciles with kcat/K_M in M⁻¹ s⁻¹) | 3 | same | **A** |
| KDM6B | as above | H3(14–34)K27**R**me2a peptide (unnatural) | kcat / K_M | 28.8 ± 3.6 × 10⁻³ s⁻¹ / 676 ± 118 µM | as above | 3 | same | A |
| **KDM6B** | truncated JMJD3 (JmjC + C-terminal Zn²⁺ GATA-like domain) | H3K27me3 peptide | **K_M** | **22 µM** | MS quantification of H3K27me2 product, 6 min reaction, peptide 0.6–250 µM | n.s. | Kruidenier et al. 2012 *Nature* 488:404, Results; PMC4691848 (author ms.) | **B** |
| **KDM6A** | **unknown** | **unknown** | "overall catalysed rate constant" | **0.29 min⁻¹** (= 4.8 × 10⁻³ s⁻¹) | **unknown** | — | **Primary source identified: Kristensen JBL et al. 2011 *FEBS Lett* 585:1951–1956, PMID 21575637**, "Enzyme kinetic studies of histone demethylases KDM4C and KDM6A". Cited as ref [23] by Ramanan/Christov 2023 *Chemistry* 29:e202300655 (PMC10526731), which converts it to a 22 kcal/mol barrier | **C** |
| KDM6A | 940–1401 (catalytic domain) | H3K27me3 / me2 / me1 peptides | *qualitative* | active on me3 and me2, **not** me1 | MALDI-TOF MS endpoint | — | Hopkinson et al. 2014 *Epigenetics* 9:306, Table 1; PMC4623018 | A |
| KDM6B | 1141–1590 | H3K27me3 / me2 / me1 peptides | *qualitative* | active on **all three** states | as above | — | same | A |
| KDM6B | catalytic domain | H3K27me3 peptide | **apparent K_D** | **0.8 µM** | AlphaScreen, 500 nM His-enzyme, biotin-H3K27me3 peptide, 15 min RT then 1 h, 22 °C, 50 mM HEPES pH 7.5 / 0.01% Tween-20 / 0.1% BSA | n.s. | Walport et al. 2014 *JBC* 289:18302, Results; PMC4140284 | **A** |
| KDM6A | catalytic domain | H3K27me3 peptide | apparent K_D | **2 µM** | as above | n.s. | same | A |
| UTY/KDM6C | 818–1347 | H3K27me3 peptide | apparent K_D | 78 µM (P1214I variant: 3.5 µM) | as above | n.s. | same | A |
| KDM6A / KDM6B / UTY | catalytic domains | 2OG (with 23-mer H3K27me3 peptide) | 2OG K_M and kcat | reported "similar for all five enzymes/variants" — **numerical values appear only in the Fig. 4B bar chart, not in the text or any table** | MALDI-MS coupled assay; enzyme 1–10 µM, Fe(II) 10 µM, 2OG 20–500 µM, ascorbate 100 µM, peptide 10–50 µM, 50 mM HEPES pH 7.5 / 50 mM NaCl / 5% glycerol, **37 °C** | n.s. | same | **B** |

### Reading of Q1

- **Only one kcat exists for the whole KDM6 family in a form you can cite with conditions: KDM6B, 0.70 min⁻¹, peptide, catalytic domain, 37 °C.**
- The KDM6A 0.29 min⁻¹ you were carrying is **real and traceable to Kristensen 2011 FEBS Lett**, not invented by the QM/MM authors. But FEBS Letters 2011 is paywalled (Wiley/Elsevier; blocked from ScienceDirect and Wiley; not OA per Unpaywall; no PMCID). Substrate, construct, cofactor concentrations, temperature and n are **unverified**. Treat as grade C.
- Taken at face value the two numbers put **KDM6A ≈ 2.4× slower than KDM6B**, but they come from different labs, different assays and (probably) different constructs, so the ratio is not defensible.
- **Full-length vs catalytic domain: no comparison exists for either enzyme.** This matters: Hopkinson 2014 and the KDM6A literature both note the TPR domain is required for efficient H3K27me1 demethylation, i.e. the non-catalytic domains demonstrably change the substrate spectrum — so the truncated-construct kcat is not safely transferable. Walport 2016 did assay *full-length* Flag-KDM6B(1–1636) from HEK293T, but only in qualitative on-bead endpoint assays, never kinetically.

---

## Q2. Is KDM6B catalytically active on nucleosomes at all?

**There is no kinetic measurement of KDM6A or KDM6B on any nucleosomal substrate. Grade X.**
Searches for KDM6 + nucleosome/nucleosomal/designer nucleosome/recombinant nucleosome return no
biochemical kinetics and no cryo-EM or crystal structure of a KDM6 enzyme on a nucleosome. Contrast
LSD1/CoREST, which has had a nucleosome-bound crystal structure since 2020.

### What the in-vivo evidence actually says

| Observation | Value / statement | Source | Grade |
|---|---|---|---|
| Jmjd3 genomic distribution in LPS+IFNγ macrophages | 54% of Jmjd3 peaks within ±0.5 kb of a TSS; >70% within 2.5 kb | De Santa et al. 2009 *EMBO J* 28:3341, Results; PMC2752025 | A |
| **Jmjd3 targets mostly lack H3K27me3 to begin with** | **2963/3339 (88%) of Jmjd3 target genes had no H3K27me3 peak within ±1 kb before LPS** | same | **A** |
| **Jmjd3 binding rarely coincides with H3K27me3** | only **511/4398 (11.6%)** of Jmjd3 peaks lay within ±500 bp of an H3K27me3 peak | same | **A** |
| **Where they do coincide, H3K27me3 does not preferentially fall** | **83/511 (16.3%)** of those peaks lost >2-fold H3K27me3 at 4 h LPS — **versus a genome-wide background of 10,385/59,173 = 17.6%** | same | **A** |
| Bulk marks in Jmjd3-null macrophages | bulk H3K27me3/me2/me1 and H3K4me3/me2/me1 "not grossly affected" by absence of Jmjd3 (western) | same | B |
| Authors' own conclusion | Jmjd3's role in tuning the LPS transcriptional output "is largely independent of H3K27me3 demethylation" | same, Abstract/Results | A |
| Jmjd3/UTX demethylase-independent chromatin remodelling | Jmjd3 required for Brg1 co-precipitation with T-bet and for *Ifng* remodelling in Th1 cells; "independent from its H3K27-demethylase potential" | Miller, Mohn & Weinmann 2010 *Mol Cell* 40:594; PMC3032266 | A (qualitative) |
| UTY is demethylase-dead yet redundant with UTX in vivo | compound X^Utx− Y^Uty− males phenocopy Utx-null females; "UTX H3K27 demethylation is non-essential for embryonic viability" | Shpargel et al. 2012 *PLoS Genet* 8:e1002964; PMC3459986 | A (qualitative) |
| Enzyme-dead UTX knock-in rescues mesoderm | UTX controls *Brachyury*/mesoderm differentiation "independent of its enzymatic activity"; UTX KO has no effect on global H3K27me3 | Wang C et al. 2012 *PNAS* 109:15324; PMC3458330 | A (qualitative) |

### Counter-evidence that the catalytic activity is real

- KDM6B(1141–1590) demonstrably turns over H3K27me3 peptide at 0.70 min⁻¹ with a 6.7 µM K_M
  (Walport 2016) — the enzyme is not catalytically dead.
- KDM6B acts on all three methylation states on peptides (Hopkinson 2014).
- GSK-J1 inhibits JMJD3 with IC₅₀ 60 nM (AlphaScreen) and GSK-J4 raises cellular H3K27me3
  (Kruidenier 2012) — so cellular H3K27me3 does depend on KDM6 catalysis to some degree.

### Verdict

**KDM6B is a genuine, catalytically competent H3K27me3 demethylase on peptides. It is not established
as the enzyme that creates the H3K27me3-depleted core at active promoters, and the strongest in-vivo
dataset argues against that role**: in LPS macrophages KDM6B lands overwhelmingly on TSSs that are
*already* H3K27me3-free, and where it does overlap H3K27me3, the mark is lost at exactly the
genome-background rate. Making KDM6B the dominant antagonist in the lattice model is currently a
modelling choice, not a measurement.

---

## Q3. Binding kinetics / residence times for KDM6A and KDM6B on chromatin

**No residence time, no bound fraction, no diffusion coefficient has been published for either
enzyme. Grade X.** No single-molecule tracking study of KDM6A or KDM6B exists.

The closest data:

| Measurement | What was actually reported | Source | Grade |
|---|---|---|---|
| FRAP of endogenous mEGFP-UTX knock-in, mouse ESCs | normalised recovery **curve** only (Fig. 2b), n = 13 cells. **No t½, no mobile fraction, no residence time given numerically.** | Shi B et al. 2021 *Nature* 597:726; PMC9008583 | **B** (curve exists) / **X** (number does not) |
| FRAP of transfected UTX-mEGFP vs His→Tyr IDR mutant, 293T | recovery curves, n = 7 cells each; HY mutant "recovered slower and less completely" — **no numbers** | same | X |
| UTX condensation | UTX forms cIDR(549–848)-dependent nuclear foci / in-vitro droplets at 30 µM EGFP-cIDR | same | A (qualitative) |

The UTX-condensate result is relevant to the lattice model in a different way: if KDM6A partitions
into condensates, its *effective local* concentration at a target site is decoupled from the
whole-nucleus concentration, and the ~10³ copies/cell bulk number below is a lower bound on local
availability.

---

## Q4. Abundance of KDM6A / KDM6B, and the magnitude of KDM6B induction

### Absolute copy numbers — two independent HeLa datasets

Both are primary, both read directly from the supplementary tables, and **they disagree by up to
4–5× on individual proteins.** Report the pair, not one of them.

**Dataset 1 (deeper; the one to use for KDM6B).** Bekker-Jensen DB et al. 2017 *Cell Syst* 4:587,
PMC5493283, **Table S3 (mmc4.xlsx), column "Copy Number (iBAQ estimated)"**; 14,237 protein groups /
12,200 coding genes in growing HeLa. Copy numbers are iBAQ values calibrated against previously
established HeLa copy-number standards, log-space Pearson r = 0.95 (Results, Fig. 4F).

**Dataset 2 (shallower, but gives concentrations and localisation).** Itzhak DN, Tyanova S, Cox J,
Borner GHH. 2016 *eLife* 5:e16950, PMC4959882, **Supplementary file 1, sheet "Compact HeLa Spatial
Proteome"**; 8,710 proteins; columns "Estimated Copy number per cell" and "Median cellular
con-centration [nM]". (File verified byte-identical, md5 81aac26f…, between the eLife CDN copy and
the copy inside the PMC supplementary archive.)

| Protein | **Bekker-Jensen 2017** copies/cell | **Itzhak 2016** copies/cell | Itzhak conc. | Itzhak percentile / MS-MS | Grade |
|---|---|---|---|---|---|
| **KDM6A** | **974** | 4,183 | 2.59 nM | 10.1 / 28 counts | **A** (both) |
| **KDM6B** | **859** | **not detected** (absent from all 8,710) | — | — | **A** (both) |
| EP300 | 37,256 | 27,961 | 17.3 nM | 32.5 / 374 | A |
| CREBBP | 13,901 | 54,919 | 34.0 nM | 42.9 / 800 | A |
| HDAC1 | 850,903 | 461,170 | 285 nM | 76.7 / 399 | A |
| HDAC2 | 514,720 | 1,542,148 | 954 nM | 89.5 / 998 | A |
| EZH2 | 30,794 | 142,635 | 88.2 nM | 58.7 / 277 | A |
| SUZ12 | 58,010 | 194,816 | 120.5 nM | 63.7 / 296 | A |
| EED | 88,382 | 351,948 | 217.7 nM | 72.9 / 340 | A |

**KDM6B is not absent from the cell — it is absent from the shallower dataset.** At 859 copies/cell
it sits just below KDM6A (974) and right at the depth limit of an 8,710-protein experiment, which is
exactly why Itzhak missed it. The two enzymes are, to within the error of these methods, **equimolar
at ~10³ copies/cell in HeLa**. Treating KDM6B as the dominant antagonist and KDM6A as a bit player is
not supported by abundance.

**Ratios for the model (Bekker-Jensen, whole-cell):**
HDAC1 : HDAC2 : EP300 : EZH2 : CREBBP : KDM6A : KDM6B ≈ **874 : 529 : 38 : 32 : 14 : 1.0 : 0.9**.

Caveats to carry into the model:

- These are **whole-cell** numbers. KDM6A, EP300, CREBBP, HDAC1 and EZH2 are annotated "Mostly
  Nuclear" (Itzhak), so nuclear concentrations are higher by roughly the cell:nucleus volume ratio
  (~3–10× in HeLa). Apply that factor uniformly — the model is sensitive to the ratios, not the
  absolute scale.
- **The cross-dataset spread is the real error bar**: 4.3× for KDM6A, 4.6× for EZH2, 3.9× for CREBBP,
  3.0× for HDAC2, 1.8× for HDAC1. Any parameter taken from a single MS dataset should carry a ±3–5×
  uncertainty, not the quoted precision.
- Itzhak's EED > SUZ12 > EZH2 ordering violates PRC2 stoichiometry; Bekker-Jensen gives the same
  qualitative ordering. Both are wrong about a complex of known 1:1:1 composition, which is a useful
  calibration of how far to trust any single number here.
- HeLa is a poor reporter cell for KDM6B specifically, since it is an immediate-early gene (below).

### KDM6B induction magnitude

**No quantitative baseline-vs-induced measurement of KDM6B protein or mRNA in absolute units has been
published. Grade X for a number.** What exists:

| Observation | Source | Grade |
|---|---|---|
| Jmjd3 induction in primary mouse BMDM peaks at **2 h** after LPS+IFNγ; shown as a western blot plus an RT-qPCR **panel (Fig. 1A)** — fold change not stated in the text | De Santa et al. 2009 *EMBO J*; PMC2752025 | **B** (timing A; magnitude X) |
| Jmjd3 induction by LPS depends on multiple conserved NF-κB sites | De Santa et al. 2007 *Cell* 130:1083 (paywalled; statement as relayed in De Santa 2009) | C |
| PMA (80 nM, 16 h) induces JMJD3 in HL-60 — western blot only | Kruidenier et al. 2012 *Nature*; PMC4691848 | B |

This is a real hole. Baseline KDM6B is now known (**859 copies/cell in unstimulated HeLa**,
Bekker-Jensen 2017), but because KDM6B is an immediate-early gene the model's KDM6B concentration is
not a constant — it is a step function whose **step size nobody has measured in protein units**. Any
steady-state lattice parameterisation using a single KDM6B concentration should be treated as a
scenario, not a measurement; and note that at 859 copies/cell even a 10-fold induction leaves KDM6B
~100× below HDAC1.

---

## Q5. p300/CBP: kcat for H3K27 acetylation, substrate preference, nucleosome kinetics

### The key negative result

Source: Henry RA, Kuo Y-M, Andrews AJ. 2013 *Biochemistry* 52:5746, PMC3756530, Results.

> "We observed low levels of p300 acetylation of H3K27 and CBP acetylation of H3K9; however, these
> sites did not appear in our steady-state experiments because they are not acetylated before 10% of
> the total substrate is acetylated. Specifically, this means that although they are detectable at
> long time points (24 h), **kinetic parameters could not be determined for these sites**."

**There is no kcat or kcat/K_M for p300 acetylation of H3K27 — on peptide, free histone, or
nucleosome. Grade X.** On free recombinant H3, H3K27 is a *minor, kinetically unmeasurable* p300 site.

### What p300/CBP kinetics do exist (free recombinant histone H3, not nucleosome)

Henry 2013 Table 2. Conditions: **37 °C**, 50 mM HEPES pH 7.8 + 100 mM ammonium bicarbonate,
saturating acetyl-CoA 200 µM, H3 0.25–15 µM, p300 1–150 nM or CBP 0.5–10 nM, SRM-MS quantification of
each site, analysis restricted to <10% total substrate acetylated.

| Enzyme | Site | kcat (× 10⁻³ s⁻¹) | K₁/₂ (µM) | kcat/K₁/₂ (× 10³ M⁻¹ s⁻¹) | n_H | Grade |
|---|---|---|---|---|---|---|
| p300 | H3K9 | 5.24 ± 0.39 | 3.52 ± 0.54 | 1.49 ± 0.25 | 1.65 ± 0.25 | A |
| p300 | H3K14 | 12.27 ± 0.47 | 3.33 ± 0.16 | 3.68 ± 0.23 | 3.10 ± 0.40 | A |
| **p300** | **H3K18** | **51.89 ± 1.14** | **2.07 ± 0.10** | **25.10 ± 1.34** | 2.82 ± 0.25 | **A** |
| p300 | H3K23 | 9.14 ± 0.83 | 2.23 ± 0.39 | 4.10 ± 0.81 | n.a. | A |
| **p300** | **H3K27** | **not determinable** (detected only after 24 h) | — | — | — | **A** (negative) |
| CBP | H3K14 | 31.58 ± 1.35 | 0.03 ± 0.07 | 95.35 ± 20.00 | n.a. | A |
| CBP | H3K18 | 106.50 ± 4.05 | 1.12 ± 0.15 | 94.84 ± 13.34 | n.a. | A |
| CBP | H4K16 | 19.38 ± 2.21 | 1.53 ± 0.40 | 12.67 ± 3.65 | n.a. | A |

Additional CBP peptide kinetics: Le Bihan et al. 2025 *JBC* 301:110465, PMC12305233 — CBP catalytic
core (res. 1095–1773), FAM-H3K18 peptide, **kcat/K_M for acetylation = 7.86 ± 1.41 × 10³ M⁻¹ s⁻¹**,
K_M(acetyl-CoA) = 10.13 ± 1.66 µM, n = 3 biological replicates. They confirmed by UFLC-MS that CBP
acetylates a FAM-H3K27 peptide (FAM-TKAARKSAPAT) but **did not determine kinetics for it**. Grade A.

### The in-vivo picture is the opposite of the in-vitro picture

Source: Weinert BT et al. 2018 *Cell* 174:231, PMC6078418 (SILAC acetylome in conditional
*Crebbp*fl/fl;*ep300*fl/fl MEFs; interventions = Cre KO, Compound-R, A-485, CBP112).

| Observation | Value | Grade |
|---|---|---|
| **H3K18, H3K27, H3K36 and N-terminal H2B acetylation are 95–99% reduced by CBP/p300 KO, Cmpd-R or A-485** — "almost exclusively acetylated by CBP/p300" | 95–99% | **A** |
| H3K14 and H3K23 acetylation are **unaffected** by loss of CBP/p300 | — | A |
| CBP/p300 "showed remarkable selectivity for sites on H3, targeting K18 and K27" | — | A |

So: on free recombinant H3 in vitro, p300 prefers K18 ≫ K23 ≈ K14 > K9 and barely touches K27; in
cells, K18 and K27 are the two sites that are *exclusively* CBP/p300-dependent while K14/K23 are not
CBP/p300 sites at all. **The in-vitro substrate preference does not predict the in-vivo one.** The most
likely explanation is nucleosomal/chromatin context plus recruitment, neither of which is captured by
the free-histone assay — which is exactly the gap the lattice model needs filled and which nobody has
filled.

**No p300 or CBP kcat on a nucleosome substrate with an H3K27 readout exists. Grade X.**

---

## Q6. p300 autoacetylation and bromodomain-mediated read-write

| Question | What exists | Source | Grade |
|---|---|---|---|
| Autoacetylation rate constants | **Primary source identified but inaccessible**: Karanam B et al. 2006 *JBC* 281:40292, "Kinetic and mass spectrometric analysis of p300 histone acetyltransferase domain autoacetylation" (PMID 17065153; no PMCID; Elsevier/ASBMB, 403 from every route tried). Related: Thompson PR et al. 2004 *Nat Struct Mol Biol* 11:308 (activation loop). | — | **C** |
| Trans-autoacetylation is activator-triggered | p300s autoacetylates slowly alone; TBK1-phosphorylated **IRF3 dimers** (not monomers, not IRF3ΔC) trigger a rapid burst of autoacetylation of the lysine-rich autoinhibitory loop (AIL). Same for pTyr701-STAT1 dimers. Crystal structure shows the AIL projecting into a neighbouring p300's active site. | Ortega E et al. 2018 *Nature* 562:538; PMC6914384 | A (qualitative) |
| Quantitative effect of activation on histone acetylation rate | H4 acetylation **V_max = 43.8 ± 5.3 cpm/min with pIRF3 vs 22.5 ± 2.8 cpm/min without** (≈2-fold), [³H]acetyl-CoA titration, equimolar pIRF3 | same | **B** (relative units; no kcat) |
| Bromodomain ligand modulates **nucleosome** acetylation | I-CBP112 stimulates p300/CBP nucleosome acetylation **up to 3-fold**; CBP30 does not; **no stimulation on isolated histone H3**, and none with the isolated HAT domain — i.e. the effect requires both the bromodomain and a nucleosome. The stimulated site is **H3K18**. | Wang Z-A et al. 2016 *Biochemistry* 55:4859; PMC5007619 | **A** (fold-change) |
| Bromodomain requirement for H3K27ac at enhancers in cells | Raisner R et al. 2018 *Cell Rep* 24:1722, "Enhancer activity requires CBP/P300 bromodomain-dependent histone H3K27 acetylation" (PMID 30110629; paywalled, not read) | — | C |
| **Quantitative read-write: rate enhancement as a function of pre-existing acetylation; intranucleosomal vs internucleosomal spreading** | **Nothing. No measurement, no rate enhancement factor, no processivity number, no nucleosome-array spreading kinetics for p300.** | — | **X** |

The I-CBP112 result is the closest thing to a quantitative handle on bromodomain-mediated allostery
(3-fold, nucleosome-dependent), but it is driven by a **synthetic ligand**, not by an acetyl-lysine
mark, so it does not license a read-write feedback term. If the lattice model has a p300 read-write
loop, its gain is currently unconstrained by data.

---

## Q7. HDAC kinetics on nucleosomal H3K27ac

### The core dataset — grade A, exactly the right substrate

Source: Wang ZA, Millard CJ, Lin C-L, Gurnett JE, Wu M, Lee K, Fairall L, Schwabe JWR, Cole PA. 2020
*eLife* 9:e57663, PMC7316510, **Table 2**.
Conditions: semisynthetic 146 bp Widom-601 nucleosomes carrying a single site-specific acetyl-lysine
on full-length H3, **100 nM nucleosome**; 30 nM CoREST, 5 nM MiDAC, 80 nM NuRD, 45 nM Sin3B, 60 nM
SMRT, 60 nM free HDAC1; 50 mM HEPES pH 7.5, 100 mM KCl, **100 µM IP6**, 0.2 mg/mL BSA, **37 °C**,
timepoints 0–120 min, site-specific western-blot readout, **n = 2–5**.

| Substrate | HDAC1 (free) | CoREST (RCOR1) | MiDAC (MIDEAS) | NuRD (MTA1) | Sin3B | SMRT (NCoR2, HDAC3) |
|---|---|---|---|---|---|---|
| **H3K27ac nucleosome, V/[E] (min⁻¹)** | N.D. | **0.093 ± 0.0088** | **0.18 ± 0.036** | **0.0077 ± 0.0028** | **0.0089 ± 0.00053** | **0.039 ± 0.0027** |
| H3K9ac nucleosome | N.D. | 0.10 ± 0.012 | 1.2 ± 0.19 | 0.012 ± 0.0032 | 0.0050 ± 0.0043 | 0.046 ± 0.0045 |
| H3K14ac nucleosome | N.D. | <0.005 | 0.52 ± 0.034 | 0.0091 ± 0.00080 | <0.005 | 0.019 ± 0.0022 |
| H3K18ac nucleosome | N.D. | 0.062 ± 0.0019 | 0.17 ± 0.0051 | 0.0046 ± 0.00088 | 0.0082 ± 0.0030 | 0.021 ± 0.0039 |
| H3K23ac nucleosome | N.D. | 0.067 ± 0.010 | 0.048 ± 0.0071 | 0.0071 ± 0.00083 | 0.0084 ± 0.0033 | 0.0081 ± 0.0028 |
| **H3K27ac on isolated (non-nucleosomal) H3**, V/[E] (min⁻¹) — Table 1 | **1.1 ± 0.25** | 2.0 ± 0.35 | 0.43 ± 0.038 | 1.5 ± 0.48 | 1.1 ± 0.31 | 1.8 ± 0.80 |

Two facts the model should absorb:

1. **Nucleosome packaging costs HDAC complexes 10–200×.** For H3K27ac, CoREST drops from 2.0 min⁻¹ on
   free H3 to 0.093 min⁻¹ on nucleosomes (21×); NuRD drops from 1.5 to 0.0077 min⁻¹ (195×).
   Any parameterisation using free-histone or peptide HDAC rates will overestimate erasure by
   1–2 orders of magnitude.
2. **H3K27ac is not a strongly discriminated site.** Across the five complexes, H3K27ac nucleosome
   rates span 0.0077–0.18 min⁻¹ and sit at or near the top of the per-complex range. The dramatic
   site selectivity in this dataset is the *rejection of H3K14ac* (Gly13-dependent), not anything
   special about K27.

### Transcription-factor stimulation of nucleosomal H3K27ac deacetylation

Source: 2026 *Nucleic Acids Res* (Cole lab), "Multivalent interactions mediate SNAIL transcription
factor stimulation of the nucleosome deacetylase activity of the CoREST complex", PMC13316430.
Conditions: 100 nM 185 bp semisynthetic nucleosome, 100 nM SNAIL, LSD1–HDAC1–CoREST (LHC) complex,
50 mM HEPES pH 7.5, 100 mM KCl, 100 µM IP6, 0.2 mg/mL BSA, 37 °C, **n = 4**.

| Substrate | V/[E] without SNAIL (min⁻¹) | V/[E] with SNAIL (min⁻¹) | Grade |
|---|---|---|---|
| **185 bp E-box-entry H3K27ac nucleosome (N9)** | **0.017 ± 0.0040** | **0.35 ± 0.035** (≈20×) | **A** |
| 185 bp E-box-entry H3K14ac nucleosome (N7) | 0.0080 ± 0.0028 | 0.15 ± 0.0040 (≈20×) | A |
| 185 bp E-box-entry H3K9ac nucleosome (N4) | 0.021 ± 0.0027 | 0.81 ± 0.22 (≈40×) | A |
| 185 bp H3K9ac, no E-box (N3) | 0.018 ± 0.0033 | 0.13 ± 0.051 (≈7×) | A |

**This is the most load-bearing HDAC result for a promoter-boundary model**: a sequence-specific
transcription factor bound at the nucleosome entry site raises HDAC1/CoREST deacetylation of
nucleosomal **H3K27ac by ~20-fold**, from 0.017 to 0.35 min⁻¹. Erasure rate at a promoter is not a
constant — it is set by what is bound there.

### HDAC1/2 residence times

**None found. Grade X.** No single-molecule tracking or FRAP-derived residence time for HDAC1 or
HDAC2 was located. (Dissociation kinetics of class-I HDAC *inhibitors* in cells have been measured,
but that is drug–target residence, not enzyme–chromatin residence, and is not the quantity needed.)

---

## Q8. Measured H3K27ac half-life in vivo

### Primary measurement — grade A

Source: Zheng Y, Thomas PM, Kelleher NL. 2013 *Nat Commun* 4:2203, PMC3929392, **Table 1**.
Method: HeLa S3; cells labelled to >90% with ¹³C₆-Arg for 7 d, then switched to light Arg + **¹³C₆-glucose**
so that new acetyl groups carry +2 Da; propionic-anhydride derivatisation, tryptic digest, SRM on a
triple quadrupole; 30 h time course; **mean of two independent fittings of two biological replicates**;
labelling efficiency ~80%, corrected by modelling.

| Site | Turnover rate (h⁻¹) | **t½ (h)** | Labelled population at plateau | Ac stoichiometry | Response to Na-butyrate, 0 vs 6 h |
|---|---|---|---|---|---|
| **H3.1/2-K27ac** | **0.32 ± 0.08** | **2.3 ± 0.6** | **0.308 ± 0.009** | 1.5% | 4.5× |
| **H3.3-K27ac** | **0.61 ± 0.02** | **1.1 ± 0.1** | **0.187 ± 0.008** | 2.7% | 3.9× |
| H3K18ac (reference) | 0.65 ± 0.07 | 1.1 ± 0.1 | 0.698 ± 0.005 | 4.1% | 5.6× |
| H3K14ac (reference) | 0.67 ± 0.08 | 1.0 ± 0.1 | 0.591 ± 0.006 | 27.8% | 2.49× |
| H3K9ac (reference) | 0.4 ± 0.1 | 2.0 ± 0.5 | 0.4927 ± 0.008 | 2.2% | 4.5× |

**The second column is the headline, but the fourth column is the one the lattice model needs.** The
labelled plateau for H3.1/2-K27ac is only 0.31 and for H3.3-K27ac only 0.19, against 0.59–0.70 for
H3K14ac/H3K18ac. The authors state this directly: "H3K27ac from both canonical H3.1/2 and the H3.3
variant showed low labeling levels (30% and 20%, respectively), **suggesting more than half of K27ac
could be stable**." So H3K27ac is kinetically **two populations**: a fast pool with t½ ≈ 1–2 h and a
large (≥ ~70% for H3.1/2) pool that does not turn over within 30 h. A single-exponential H3K27ac
lifetime in the model will misrepresent this.

Also from the same table: H3K27 **methylation** stoichiometry is 70% (H3.1/2) and 58% (H3.3), against
H3K27 **acetylation** stoichiometry of 1.5% and 2.7% — i.e. bulk H3K27ac is a ~2% mark competing with
a ~60–70% methyl mark.

### Secondary measurement — grade B

Source: Weinert BT et al. 2018 *Cell* 174:231, PMC6078418.
MEFs treated with the p300/CBP catalytic inhibitor **A-485**, six timepoints, SILAC acetylome;
deacetylation half-lives computed for 619 sites quantified at all six timepoints.

| Quantity | Value | Grade |
|---|---|---|
| Median deacetylation t½ across 619 CBP/p300-regulated sites | **93.6 min** | A |
| Histone sites: half-lives computed for 14 of 20 CBP/p300-regulated sites; **12 of them < 60 min** | <60 min | A (as a bound) |
| **H3K27ac-specific half-life** | plotted in **Fig. 5B** and listed in **Table S4**; **not stated in the running text**, and the PMC supplementary files are behind a proof-of-work download gate that could not be passed | **B** |
| Deacetylation, not protein turnover, dominates | acetylation t½ ≪ protein t½ (comparison to McShane 2016) | A |

**These two numbers measure different things and should not be averaged.** Zheng 2013 measures
*appearance* of new acetyl groups with writers and erasers both running (a writer-limited exchange
time, 1–2 h). Weinert 2018 measures *decay* after the writer is switched off (an eraser-limited decay
time, <1 h). For a lattice model the Weinert number is the right one for k_deacetylation and the Zheng
number is the right one for steady-state exchange — and the Zheng stable fraction is a constraint on
both.

### Bonus: the only in-vivo H3K27 *demethylation* rate constants

Zheng Y, Sweet SMM, et al., Kelleher NL. 2012 *PNAS* 109:13549 (M4K), PMC3427122 — dual-SILAC
labelling in KMS11 myeloma cells, 44 effective rate constants for H3K27/H3K36 methylation **and
demethylation**, units of **d⁻¹**, conditioned on the neighbouring-site methylation state. The
demethylation constants are reported **only in the Fig. 3A/3B matrices**, not in the text. The authors
note these are compound "effective" constants related to kcat by the scalar [E]/K_M and lump all
demethylases together. **Grade B.** This is the only published in-vivo H3K27 demethylation rate in any
units and is worth extracting from the figure if the model needs an order of magnitude — but it cannot
be attributed to KDM6B.

---

# What is still missing

Ordered by how badly it hurts the lattice model.

1. **Any kinetics for KDM6A or KDM6B on a nucleosome.** Not slow, not fast — *absent*. There is no
   kcat, no K_M, no V/[E], no structure. Every KDM6 number in existence is a free-peptide number. By
   analogy with HDAC complexes (10–200× slower on nucleosomes than on free H3, Wang 2020) and with
   p300 (free-histone site preference that inverts in cells, Henry 2013 vs Weinert 2018), the peptide
   kcat of 0.70 min⁻¹ is probably an **upper bound by one to two orders of magnitude**, but that is an
   inference, not a measurement. *This is the single largest unconstrained parameter in the model.*

2. **KDM6B induction magnitude.** Baseline abundance is now pinned (859 copies/cell in HeLa,
   Bekker-Jensen 2017), but there is still **no measured induction fold in protein units** in any
   system, and no absolute quantification in a cell type where KDM6B is actually doing something.
   Since KDM6B is an immediate-early gene, the induced concentration — not the 859 baseline — is what
   the model needs, and it remains a free parameter spanning at least one to two orders of magnitude.

3. **Residence time / bound fraction for either KDM6 enzyme.** No SMT exists; the two FRAP experiments
   that exist (Shi 2021) report curves without fitted parameters. The p300 numbers you already have
   (Ferrie 2023: f_bound 0.55 ± 0.02 and 0.52 ± 0.03, FRAP residence ~26 s) have no KDM6 counterpart,
   so the two arms of the model are parameterised on incommensurate evidence.

4. **A p300/CBP kcat for H3K27 specifically.** Henry 2013 explicitly could not fit it on free H3;
   Le Bihan 2025 confirmed the activity on an H3K27 peptide but characterised only H3K18. Nobody has
   measured p300 on a nucleosome with an H3K27 readout. The H3K18 numbers (p300 kcat 0.052 s⁻¹,
   kcat/K₁/₂ 2.5 × 10⁴ M⁻¹ s⁻¹) are the nearest available proxy and should be labelled as a proxy.

5. **Quantitative p300 read-write.** No rate-enhancement factor, no processivity, no intranucleosomal
   vs internucleosomal comparison, no nucleosome-array spreading kinetics. The only quantitative
   bromodomain number is a 3-fold *ligand*-driven allosteric stimulation on nucleosomes (Wang 2016),
   which is not a read-write measurement.

6. **Kristensen 2011 FEBS Lett full text** — the primary source of the KDM6A 0.29 min⁻¹. Paywalled
   (Wiley and ScienceDirect both refuse; no PMCID; not OA per Unpaywall). Needs an institutional
   PDF pull to recover substrate, construct, cofactor concentrations, temperature and n. Until then
   the number is grade C and should not be used alongside the grade-A KDM6B value to compute a ratio.

7. **The numeric 2OG K_M and kcat values in Walport 2014 Fig. 4B** for KDM6A, KDM6B and UTY. These
   would give a same-lab, same-assay KDM6A-vs-KDM6B comparison — the only one that would exist — but
   they were never written into the text or a table. Recoverable only by digitising the bar chart or
   contacting the authors.

8. **H3K27ac deacetylation half-life, site-specific, from Weinert 2018 Table S4 / Fig. 5B.** Blocked
   by the PMC supplementary download gate. Obtainable from the Cell website with institutional access,
   or from the authors' web resource.

9. **HDAC1/2 chromatin residence times.** No live-cell measurement located. The model currently has
   enzymatic rates for HDAC complexes on nucleosomal H3K27ac (grade A) but no occupancy term.

10. **A catalytically-dead *Kdm6b* knock-in animal or cell line.** Both *Utx* (Wang 2012, enzyme-dead
    knock-in ES cells) and *Uty* (Shpargel 2012, naturally demethylase-dead) have been used to
    separate catalytic from non-catalytic function. The equivalent experiment for KDM6B was not found.
    Without it, the question "how much of KDM6B's effect on H3K27me3 at promoters is catalytic?"
    cannot be answered in vivo.

---

## Sources used (all read directly)

| Paper | Access route | Used for |
|---|---|---|
| Walport LJ et al. 2016 *Nat Commun* 7:11974, PMC4931022 | EuropePMC fullTextXML | **KDM6B kcat/K_M**, assay conditions, construct |
| Walport LJ et al. 2014 *JBC* 289:18302, PMC4140284 | EuropePMC fullTextXML | KDM6A/B/UTY peptide K_D; MALDI assay conditions |
| Hopkinson RJ et al. 2014 *Epigenetics* 9:306, PMC4623018 | EuropePMC fullTextXML | KDM6A/B methylation-state selectivity, constructs |
| Kruidenier L et al. 2012 *Nature* 488:404, PMC4691848 | PMC HTML (author ms.) | JMJD3 peptide K_M 22 µM; GSK-J1 IC₅₀ 60 nM |
| Ramanan R, Christov CZ et al. 2023 *Chemistry* 29, PMC10526731 | PMC HTML | traced 0.29 min⁻¹ → Kristensen 2011; KDM6B 0.0117 s⁻¹ → Walport 2016 |
| De Santa F et al. 2009 *EMBO J* 28:3341, PMC2752025 | EuropePMC fullTextXML | Jmjd3 TSS occupancy vs H3K27me3 in vivo |
| Miller SA et al. 2010 *Mol Cell* 40:594, PMC3032266 | PMC HTML | demethylase-independent Jmjd3 function |
| Shpargel KB et al. 2012 *PLoS Genet* 8:e1002964, PMC3459986 | PMC HTML | UTX/UTY demethylase-independent redundancy |
| Wang C et al. 2012 *PNAS* 109:15324, PMC3458330 | PMC HTML | enzyme-dead UTX knock-in |
| Shi B et al. 2021 *Nature* 597:726, PMC9008583 | PMC HTML | UTX FRAP (curves only), cIDR condensation |
| **Bekker-Jensen DB et al. 2017 *Cell Syst* 4:587, PMC5493283** | EuropePMC supplementaryFiles archive, Table S3 (mmc4.xlsx), parsed with openpyxl; text via EuropePMC fullTextXML | **KDM6B and KDM6A absolute copy numbers**; iBAQ calibration (r = 0.95) |
| Itzhak DN et al. 2016 *eLife* 5:e16950, PMC4959882 | elifesciences CDN + PMC supplementary archive (md5-verified identical), Supp. file 1 (xlsx), parsed with openpyxl | absolute copy numbers / nM / localisation |
| Henry RA, Kuo Y-M, Andrews AJ 2013 *Biochemistry* 52:5746, PMC3756530 | EuropePMC fullTextXML | p300/CBP site kinetics; **H3K27 not determinable** |
| Le Bihan Y-V et al. 2025 *JBC* 301:110465, PMC12305233 | EuropePMC fullTextXML | CBP H3K18 kcat/K_M; H3K27 peptide activity |
| Weinert BT et al. 2018 *Cell* 174:231, PMC6078418 | PMC HTML | H3K18/K27ac 95–99% CBP/p300-dependent; median t½ 94 min |
| Ortega E et al. 2018 *Nature* 562:538, PMC6914384 | EuropePMC fullTextXML | p300 trans-autoacetylation, ~2× V_max stimulation |
| Wang Z-A et al. 2016 *Biochemistry* 55:4859, PMC5007619 | PMC HTML | bromodomain ligand, 3× nucleosome-specific stimulation |
| Wang ZA et al. 2020 *eLife* 9:e57663, PMC7316510 | EuropePMC fullTextXML | **HDAC complex V/[E] on H3K27ac nucleosomes** |
| 2026 *Nucleic Acids Res*, PMC13316430 | EuropePMC fullTextXML | SNAIL → 20× stimulation of LHC on H3K27ac nucleosomes |
| Whedon SD et al. 2024 *JACS*, PMC11638967 | EuropePMC fullTextXML | sortase nucleosome toolkit; pointer to deacetylation rate papers |
| Zheng Y, Thomas PM, Kelleher NL 2013 *Nat Commun* 4:2203, PMC3929392 | PMC HTML | **H3K27ac t½ in vivo**, stable fraction, stoichiometry |
| Zheng Y et al. 2012 *PNAS* 109:13549, PMC3427122 | PMC HTML | in-vivo H3K27 (de)methylation effective rate constants (figure only) |

### Not obtained

| Paper | Why it matters | Blocker |
|---|---|---|
| Kristensen JBL et al. 2011 *FEBS Lett* 585:1951, PMID 21575637 | **primary source of KDM6A 0.29 min⁻¹** | Wiley 403, ScienceDirect 403, no PMCID, not OA |
| Karanam B et al. 2006 *JBC* 281:40292, PMID 17065153 | p300 autoacetylation rate constants | Elsevier/ASBMB 403, no PMCID |
| Raisner R et al. 2018 *Cell Rep* 24:1722, PMID 30110629 | bromodomain requirement for enhancer H3K27ac | no PMCID |
| Lee MG et al. 2007 *Science* 318:447 | original UTX nucleosome claim | no PMCID |
| Weinert 2018 Table S4 | site-specific H3K27ac deacetylation t½ | PMC proof-of-work download gate |
