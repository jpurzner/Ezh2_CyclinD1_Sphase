# v44 model — gene / transcript map for bulk RNA-seq integration

Purpose: a complete list of every gene / transcript represented in the current v44 model
(`src/build_model_v44_heldt.py`, built on the Heldt 2018 core + the EZH2 / Hedgehog / H3K27me3
biology), so bulk RNA-seq can be pulled per gene to quantify transcript differences across time
(cell-cycle phase and/or condition). Companion CSV: `model_gene_transcript_map.csv`.

Compiled from the model source + `DATA_PROVENANCE.md`. Symbols given for **mouse** (primary, GNP/MB
biology) and **human** (SHH-MB). "Explicit mRNA species" = the model carries a transcript variable
whose trajectory you can compare directly to a time-resolved RNA-seq value; where it's "no", only a
protein/activity/state variable exists (you can still overlay RNA-seq, but there is no in-silico
transcript to match against).

## Legend
- **mRNA species** — the model variable that is the transcript, if any.
- **RNA-seq target** — already used as a calibration target in `DATA_PROVENANCE.md` (with the fold used).
- **Module** — HH = Hedgehog→CyclinD1; EZH2 = epigenetic writer; H3K27 = mark/PRC2 layer;
  CC = Heldt cell-cycle engine; MIT = mitotic switch; CKI = CDK-inhibitor panel.

---

## A. Genes with an explicit transcript variable in the model (best for a time-course figure)

| Gene (mouse / human) | mRNA species | Protein/other species | Role in model | RNA-seq target? | Module |
|---|---|---|---|---|---|
| **Ccnd1** / **CCND1** (Cyclin D1) | `Cd_mRNA` | `Cd` (protein) | Convergence node: mitogen (Gli+MYCN) vs EZH2/H3K27me3 brake; drives Rb phosphorylation | **Yes** — MB/GNP 5.07×; +HHi folds; EZH2i ≈2.2× de-repression | HH |
| **Ezh2** / **EZH2** | `EZH2m` | `EZH2` (protein) | PRC2 writer; E2f×(CycE+CycA)-gated, mitogen-dose-scaled; integrates S-phase duration | **Yes** — cycling/G0 ≈2.0×; G2/G0 1.48×; MB/GNP 2.05× | EZH2 |
| **Gli1** / **GLI1** | `Gli1_mRNA` | `Gli1`, `Gli1_epi` (chromatin memory) | HH pathway amplifier/readout; feeds CyclinD1 + the Kdm6b eraser | **Yes** — MB/GNP 6.9×; +HHi >98%↓ | HH |
| **Ptch1** / **PTCH1** | `Ptch1_mRNA` | `Ptch1_free`, `SHH_Ptch` | Gli target; SHH receptor; broken Gli→Ptch1→Smo feedback in MB → high Ptch1 mRNA | partial (copy-number pinned by Gli1 fold) | HH |

> These four are the only genes with a model transcript. For a "transcript vs time" figure they are the
> clean comparisons. Note the **nascent vs mature** caveat for `Cd_mRNA` in Problems (#6).

---

## B. Genes represented as protein / activity only (RNA-seq can overlay, but no model transcript)

| Gene (mouse / human) | Model species | Role in model | RNA-seq target? | Module |
|---|---|---|---|---|
| **Mycn** / **MYCN** | `MYCN` (protein flux) | Mitogen-independent CyclinD1 floor; amplified in MB | **Yes** — MB/GNP 2.86× — *but no mRNA species* (#2) | HH |
| **Shh** / **SHH** | `$SHH` (boundary input) | Pathway ligand input (0.5 default) | no (set as input) | HH |
| **Smo** / **SMO** | `Smo_active` | HH transducer; vismodegib (HHi) target | no | HH |
| **Gli2/Gli3** / **GLI2/GLI3** | `Gli_act`, `Gli_rep` (lumped) | Activator/repressor Gli forms upstream of Gli1 | no — lumped, no transcript (#3) | HH |
| **Rb1** / **RB1** | `Rb`, `Rbm` (mono), `pRb` (hyper), `RbE2f`, `RbmE2f` | Restriction-point switch; pRb(hyper)=Ser807/811 marker | no | CC |
| **E2f1** / **E2F1** (E2F1–3 activators) | `E2f`, in `RbE2f`/`RbmE2f` | Restriction-point TF; drives CycE/A, Skp2, EZH2, Emi1 | palbo qPCR: E2f1 ↓56.7% | CC |
| **Ccne1** / **CCNE1** (Cyclin E) + **Cdk2** | `Ce`, `CeP21` | S-entry kinase (with CDK2) | no | CC |
| **Ccna2** / **CCNA2** (Cyclin A) + **Cdk2** | `Ca`, `CaP21` | S/G2 kinase; APC/Cdc20 substrate | no | CC |
| **Ccnb1** / **CCNB1** (Cyclin B) + **Cdk1** | `MPF`, `preMPF` (`Cb` readout) | Mitotic (MPF) switch | no | MIT |
| **Cdk4 / Cdk6** / **CDK4/CDK6** | implicit in Cd→Rb kinase (`kPhRbCd`) | CyclinD1 partner kinase; palbociclib target | no (implicit) | CC |
| **Skp2** / **SKP2** | `Skp2` | E2F-target SCF F-box; degrades p27; Cdh1-regulated | no | CC |
| **Fzr1 (Cdh1)** / **FZR1** | `C1`, `pC1`, `E1C1` | APC/C-Cdh1; degrades Skp2/CycA | no | CC |
| **Fbxo5 (Emi1)** / **FBXO5** | `E1`, `E1C1` | APC/C-Cdh1 inhibitor; E2F target | no | CC |
| **Cdc20** / **CDC20** | `Cdc20` | Mitotic APC/C activator; degrades CycA/B | no | MIT |
| **Chek1 (Chk1)** / **CHEK1** | `Chk1` (assignment of `aRc`) | Intra-S replication checkpoint | no | MIT |
| **Cdc25a** / **CDC25A** (Cdc25 family) | `Cdc25a` (assignment) | MPF-activating phosphatase | no | MIT |
| **Wee1** / **WEE1** | `Wee1a` (assignment) | MPF-inhibiting kinase | no | MIT |
| **Pcna** / **PCNA** | `aPcna`, `iPcna` | Replication processivity; p21 sink | no | CC |
| **Dtl (Cdt2)** / **DTL** | `Cdt2` (constant) | CRL4-Cdt2 p21 degradation | no (constant) | CC |
| **Trp53 (p53)** / **TP53** | `P53` | DNA-damage arm (p21 synthesis) | no | CC |
| **Kdm6b (Jmjd3)** / **KDM6B** | eraser term `k_jmjd3_gli*Gli1` | Gli-driven H3K27me3 demethylase (MB<GNP mark) | direction only | H3K27 |
| **Eed** / **EED** | inside `PRC2` read-write (`a_rw`) | PRC2 allosteric reader (me3→activate EZH2) | no | H3K27 |
| **Suz12** / **SUZ12** | inside `PRC2` accessory (`a0`) | PRC2 core / recruitment | no | H3K27 |

---

## C. CDK-inhibitor panel — constants, not dynamic species (#9)

The CKIs are fixed parameters (GNP vs MB set-points), so there is **no time course** to compare against,
even though the panel is RNA-seq/DESeq2 measured. Listed for completeness of the gene set.

| Gene (mouse / human) | Model param | GNP → MB | Notes |
|---|---|---|---|
| **Cdkn2a (p16)** / **CDKN2A** | `p16` | 0 → >0 | H3K27me3-silenced in GNP; ~139× induced MB |
| **Cdkn2c (p18)** / **CDKN2C** | `p18` | 0.464 → higher | Constitutive INK4 |
| **Cdkn1b (p27)** / **CDKN1B** | `P21` pool (see #1) | high at birth (`P21_div`) | Dominant CIP/KIP per docs |
| **Cdkn1a (p21)** / **CDKN1A** | `P21` species (Heldt) | — | Heldt's literal p21; conflated with p27 (#1) |

---

## D. Marks / states / drugs (not transcripts — for reference only)

- **H3K27me1 / me2 / me3** at the *Ccnd1* locus → `m1_me`, `m2_me`, `Mk` (chromatin marks, ChIP not RNA-seq).
- **Mki67 / MKI67** — proliferation readout referenced in comments; not a model species (validation only).
- Drugs (inputs, not genes): `HHi` = vismodegib, `EZH2i` = tazemetostat, `HU` = hydroxyurea,
  CDK4/6i = palbociclib (`kPhRbCd`→0).

---

## Problems / issues found (for the modeling agent — flagged, not fixed)

1. **`P21` species is doubly assigned (p21 vs p27).** Heldt's `P21` is Cdkn1a (p21), but the v44 comments
   repeatedly reinterpret the *same* pool as p27/Cdkn1b ("the dominant CIP/KIP in GNP/MB is p27"; `P21_div`,
   `w_p27*P21`). One pool cannot map cleanly to a single transcript — decide whether the RNA-seq comparator
   is *Cdkn1a* or *Cdkn1b* (they differ in the data), or split the species.
2. **MYCN has no mRNA species.** `MYCN` is a direct protein synthesis flux, yet *Mycn* is a measured RNA-seq
   target (2.86× MB/GNP) that the paper will want to show over time. There is no in-silico *Mycn* transcript
   to compare a time course against.
3. **Gli2/Gli3 are lumped into `Gli_act`/`Gli_rep`** with no transcripts. If the RNA-seq includes *Gli2*/*Gli3*
   there is no model variable to map them onto; only *Gli1* has a transcript.
4. **EZH2 transcript trajectory depends on the just-changed concentration convention.** The Phase 0b commit
   (dca584f) dropped the division halving and re-fit `kTlEZ`; `with_ezh2_conc=True/False` give different
   `EZH2m`/`EZH2` levels, and there are **two different baked `kTlEZ` values** (0.011771 in-block vs the
   `_ts_bake` two-step override). Pin one build config before generating an EZH2 transcript figure.
5. **`Cd_mRNA` dynamics are build-flag-dependent.** The H3K27me3 repression term on `Cd_mRNA` transcription is
   swapped across `with_h3k27_memory` / `with_h3k27_dilution` / `with_prc2` / `with_h3k27_chain`, so the
   *Ccnd1* transcript time-course differs per build. The figure build must be pinned (default is dilution+PRC2+chain).
6. **`Cd_mRNA` is used as both nascent and mature transcript.** It is the PRC2-eviction proxy ("nascent output")
   *and* the pool translated to Cd protein. Bulk RNA-seq measures steady-state mRNA; the eviction arm wants
   nascent. Interpret the *Ccnd1* comparison accordingly.
7. **`Ptch1` transcript vs functional copy number.** In the model `Ptch1_mRNA` stays *high* in MB (broken Gli
   feedback), while `Ptch1_copy_number` (a functional fraction) is what's lowered. Compare RNA-seq *Ptch1* to
   `Ptch1_mRNA`, **not** to the copy-number parameter — otherwise the direction looks inverted.
8. **CKIs (p16/p18/p27/p21) are constants (#9), no time course.** The DESeq2 CKI panel is measured, but the
   model has no CKI dynamics to plot against time — only fixed GNP/MB set-points.
9. **Time axis must be defined.** The model runs in minutes with an *emergent* cycle; "across time" in bulk
   RNA-seq is most likely cell-cycle *phase* (needs sorting/scRNA to resolve) or *condition* (GNP / Ptch⁺ᐟ⁻ /
   MB ± vismo). Bulk RNA-seq cannot resolve phase without sorting — clarify which axis the figure uses so the
   model readout (phase-binned vs condition-binned) matches.
</content>
</invoke>
