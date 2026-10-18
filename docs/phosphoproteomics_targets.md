# Phosphoproteomics / TMT target list (P1, P7, P14 GNP + SHH-MB)

Purpose: pull specific **total-protein** and **phosphosite** quantities from the TMT proteomic/phosphoproteomic
data to resolve open questions in the v44 model that RNA-seq (transcript) cannot answer. Tier 1 = resolves a current
blocker; Tier 2 = high value; Tier 3 = if depth allows.

**Numbering:** residues below are **human UniProt canonical** (most-cited). The data are mouse (P1/P7/P14 GNP + Ptch1
MB) — sites are conserved but mouse numbering may differ by a few residues; the search IDs by phosphopeptide, so map
by sequence context. Key: **read caveats at the bottom first** — several of the highest-value items are low-abundance
and may not be captured without targeted enrichment.

---

## Group A — CIP/KIP inhibitors (dynamic-CDKI module + the p27 dual-function question)

**Total protein (Tier 1):** CDKN1B (p27), CDKN1A (p21), CDKN1C (p57)
- *Resolves:* re-ground the CDKI module's CKI abundances (currently transcript-derived). **The single most important
  number: p21 protein in MB vs P7 GNP** — does p21 protein stay low despite Cdkn1a transcript rising ~2.8×? That tests
  the miR-17/92 "transcript overstates protein" hypothesis and the "can p21 be at abundance" question blocking Module A.
  Also: is p57 protein genuinely low (justifying carrying it off)? Is p27 the dominant CIP/KIP at protein level?

**Phosphosites:**
| Protein | Site (human) | Kinase | What it resolves | Tier |
|---|---|---|---|---|
| p27 (CDKN1B) | **pTyr88** (also pTyr74, pTyr89) | Src/Abl/BRK/Lyn | **THE dual-function readout** — pY88-p27 *activates* CDK4 (Guiley 2019), palbo-insensitive → "p27-high but cycling." Enables Reparam 1. pY88/total-p27 = occupancy. | **1** |
| p27 | pThr187 | CyclinE-CDK2 | SCF-Skp2 degron (p27 turnover term) | 1 |
| p27 | pSer10 | hKIS/DYRK1B | CRM1 nuclear export → cytoplasmic p27 (localization term; Bhatia GNP axis) | 2 |
| p27 | pThr198 | AKT/RSK | 14-3-3 binding, cytoplasmic/stability | 2 |
| p21 (CDKN1A) | pThr145 | AKT | cytoplasmic p21, blocks PCNA binding (ties to the p21-PCNA un-map) | 2 |
| p21 | pSer114 | (GSK3?) | phosphomimetic enhances CRL4^Cdt2 poly-Ub (degron rate) | 3 |
| p57 (CDKN1C) | pThr310 | CDK2 | Skp2 degron (p57 turnover), analogous to p27-T187 | 3 |

---

## Group B — Rb + pocket proteins (two-step Rb, G0 marker, CDK4-vs-CDK2)

**Total protein (Tier 1):** RB1 (Rb); **(Tier 2)** RBL1 (p107), RBL2 (p130) — p130/DREAM = deep-quiescence (Module B / depth)
**Phosphosites (RB1):**
| Site (human) | Kinase | What it resolves | Tier |
|---|---|---|---|
| **pSer807 / pSer811** | CDK4/6 **and** CDK2 | the D20B12 site you stain; hyperphospho/committed marker — validates the two-step Rb G0 marker | **1** |
| **pSer780** | CDK4/6-**specific** | Ser780 vs Ser807/811 ratio separates CDK4 from CDK2 drive — the CDK4/6i-resistance-via-CDK2 axis | **1** |
| pSer795, pThr826, pThr373, pSer612 | CDK | coverage of the mono→hyper phospho-Rb ladder | 2 |

---

## Group C — CDK-activity reporters (bulk mean CDK2/CDK1 activity — the "sensor" analogue)

**Phosphosites (Tier 2)** — canonical substrate sites whose phospho ≈ kinase activity:
| Protein | Site | Kinase reported | Note |
|---|---|---|---|
| NPM1 | pThr199 | CyclinE/A-CDK2 | classic **CDK2-activity** reporter |
| CDC6 | pSer54 | CDK2 | S-entry CDK2 activity |
| RB1 | pSer807/811 | CDK2+CDK4 | (see B) |
| CDK1 | pTyr15 | Wee1/Myt1 (inhibitory) | **MPF switch** — low pY15 = active CDK1 (mitosis) |
| CDC25C | pSer216 | CHK1/PLK (14-3-3) | mitotic-entry gate |
- *Resolves:* mean CDK2/CDK1 activity MB vs GNP and across P1→P14 — the bulk analogue of the live CDK2 sensor
  (constrains the commitment/quiescence setpoint; **does NOT** give the single-cell bimodality — see caveats).

---

## Group D — Cyclins / CDKs / INK4 (the CyclinD1 governor + phase state)

**Total protein (Tier 1):** CCND1, CCND2 (D-cyclin pool — the two-cyclin model); CDKN2C (p18), CDKN2D (p19), CDKN2A
(p16), CDKN2B (p15) (INK4 CDK4/6 brakes); **(Tier 2)** CCNE1/2, CCNA2, CCNB1; CDK4, CDK6, CDK2, CDK1
- *Resolves:* the D1/D2 protein-pool balance (we model a D1 transcript reservoir + a buffered D2); the INK4 protein
  ratios that set MB's longer G1; cyclin protein levels for phase.

---

## Group E — EZH2 / PRC2 (the governor; cell-cycle coupling)

**Total protein (Tier 1):** EZH2; **(Tier 2)** EED, SUZ12, KDM6A, KDM6B (erasers)
**Phosphosites (EZH2, Tier 2):**
| Site (human) | Kinase | What it resolves |
|---|---|---|
| pThr350 | CDK1/CDK2 | promotes EZH2 recruitment/activity — tests "EZH2 activity is cell-cycle-coupled" (our core claim) |
| pThr487 | CDK1 | disrupts EZH2–EED, inhibits methyltransferase |
| pSer21 | AKT | inhibits EZH2 (H3K27me3 ↓) |
- *Resolves:* whether EZH2 *activity* (not just level) tracks the cycle — the EZH2-tracks-mitogen/cycle mechanism.
- **⚠ H3K27me3 itself is NOT measured by phospho-TMT** — the mark central to the whole model needs a dedicated
  **histone-PTM MS** run (or ChIP). Flagging because it's the biggest thing TMT phospho can't give.

---

## Group F — Hedgehog / MYCN / identity (mitogen input)

**Total protein (Tier 1):** MYCN; **(Tier 2)** GLI1, GLI2, GLI3, SUFU, MYC, ATOH1, PTCH1, SMO
**Phosphosites:**
| Protein | Site (human) | Kinase | What it resolves | Tier |
|---|---|---|---|---|
| MYCN | pThr58 | GSK3β | degron (destabilizing) — the MYCN-stability mechanism behind "elevated expression" | 1 |
| MYCN | pSer62 | CDK1/ERK | priming/stabilizing | 1 |
| GLI2/GLI3 | PKA/GSK3/CK1 cluster | processing kinases | rep→act processing switch (if detected) | 3 |

---

## Group G — Commitment / degradation machinery + the irreversibility gate (Module A)

**Total protein (Tier 2):** SKP2, CKS1B; FZR1/Cdh1, CDC20; **FBXO5/Emi1** (the Cdh1↔Emi1 hysteresis = Module A
irreversibility gate); DTL/Cdt2, CUL4A, DDB1 (CRL4^Cdt2 = the p21 degron / un-map)
**Phospho (Tier 3):** FZR1 CDK-inactivating sites; Emi1 degron sites
- *Resolves:* Skp2/Cdh1/Cdt2 protein levels for the p27/p21 turnover terms; Emi1 for the commitment hysteresis.

---

## Group H — Stem / quiescence compartment (Module B — the p57/OLIG2 axis)

**Total protein (Tier 2):** OLIG2, SOX2, NES (Nestin), SOX10, ASCL1
**Phospho (Tier 3):** OLIG2 pSer10/pSer13/pSer14 (N-terminal triple-serine; the CT-179 target; quiescence↔activation)
- *Resolves:* whether the OLIG2/SOX2/p57 stem-quiescence program is detectable at protein in bulk MB (grounds Module B).
  **⚠ bulk cannot resolve the <5 % SOX2⁺ compartment** — expect dilution; a strong signal is informative, absence is not.

---

## Group I — mTOR / localization axis (cytoplasmic-p27 confounder; Bhatia 2009)

**Total (Tier 3):** TSC1, TSC2
**Phospho (Tier 2):** RPS6 pS235/236 & pS240/244; EIF4EBP1 pT37/46; AKT1 pS473/pT308; GSK3B pS9
- *Resolves:* mTOR/AKT activity that drives p27 cytoplasmic mislocalization (why a tumor can be p27-high yet cycling) —
  pairs with p27-pS10/pT198.

---

## Caveats of the technique (what will and won't be usable)

1. **pTyr detection is the big one.** The single most valuable site, **p27-pY88**, is a *phospho-tyrosine*. pTyr is
   <1 % of the phosphoproteome and standard IMAC/TiO2 TMT phospho-enrichment is biased to pSer/pThr — **pY88 is likely
   MISSED unless a dedicated pTyr (anti-pY antibody) enrichment was run.** → *Please tell me whether pTyr enrichment was
   done.* If not, we can't get the p27-activator readout from this dataset (would need a targeted pY IP-MS or pY88 IHC).
2. **Low-abundance proteins may not be quantified.** CKIs (especially **p21, p57**), GLI proteins, OLIG2, and
   transcription factors are low-copy; TMT depth often misses them or gives them at high CV. The *phosphosites* on them
   are harder still. Treat non-detection as "below depth," not "absent."
3. **Bulk = population mean.** TMT gives the *average* — it **cannot** resolve the single-cell **bimodality** (the
   CDK2-low *fraction*) that *is* the transient-G0 bifurcation (Module A). It gives mean levels + phospho-states to
   *parameterize* the modules, not the stochastic distribution (that needs single-cell / a CDK2 sensor).
4. **Phospho occupancy needs paired total protein.** pSite/total ratios (e.g., pY88/p27, pS780/Rb) require the total
   proteome from the same samples. If only phospho was run, we get *relative* phosphopeptide change, not stoichiometry.
5. **Species mapping.** Human numbering above → map to the mouse ortholog by peptide sequence.
6. **P7 GNP is the model's "GNP" reference; MB is the tumor endpoint.** P1→P7→P14 is a proliferation→differentiation
   trajectory (bonus: grounds the graded-exit / K_g1len axis).

## What the results would let us do (interpretation plan)
- **p21 protein low in MB** → confirms carrying p21 low-activity is *correct* (not a fudge); **p21 protein high** →
  Module A must engage p21 as a real inhibitor and the p21-abundance/population question reopens.
- **p27-pY88 present in MB** → build Reparam 1 (p27-as-CDK4-activator) as a *calibrated* term, not a hypothesis.
- **Rb pS780 vs pS807/811 split** → apportion the model's CDK4 vs CDK2 Rb-drive; test the CDK4/6i-resistance story.
- **CKI protein ratios (p18/p19/p27/p21/p57)** → replace the transcript-derived abundances in the CDKI `_ts_bake`.
- **EZH2 pT350/pT487 tracking the cycle / cell type** → direct support (or not) for EZH2-activity-tracks-the-cycle.
- **MYCN pT58/pS62** → whether "elevated expression" is stability-driven.
