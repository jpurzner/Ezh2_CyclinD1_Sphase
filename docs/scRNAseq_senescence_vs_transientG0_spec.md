# scRNA-seq analysis spec — is the model's "transient G0" actually senescence?

**For:** the scRNA-seq agent (knows the GNP / SHH-MB datasets).
**From:** the EZH2–CyclinD1 cell-cycle model (Chahin et al. v44).
**One-line ask:** test whether the model's *reversible transient-G0 arrest* in MB is being conflated
with *senescence*, by scoring a senescence program **with all CDK-inhibitor genes removed** and asking
whether it tracks the model's arrest/brake features or is a separate state.

---

## 1. Why we're doing this (the claim under test)

The model reproduces the MB drug responses (vismodegib arrests, EZH2i rescues) by raising the
proliferation **commitment threshold** in MB with the CDK-inhibitor tone — INK4 (p16/p15/p18 → CDK4/6)
+ CIP/KIP (p21/p27 → CDK2). It then calls the resulting arrest a **reversible "transient G0"** (a
p27-high pause that re-enters when CyclinD1 recovers), *not* senescence.

The fair criticism: those CKIs (especially p16/Cdkn2a) are the canonical **senescence** effectors, so
a reviewer may argue we have relabeled senescence as transient quiescence. The model makes a
falsifiable prediction that settles it: in MB the high-CKI cells are **still proliferating** (CyclinD1
overpowers the brake; pauses are transient), so they should be **senescence-LOW**. Senescence, if
present, should be a **separate, non-cycling, SASP⁺** minority that the model does not describe.

---

## 2. ⚠️ CRITICAL methodological rule — pull the CKIs OUT of every signature

The CDK-inhibitor genes are the *model's brake*. They must **not** also be inside the senescence (or
quiescence) signatures, or the correlation is circular (senescence would "correlate" with our brake
only because they share p16/p21).

- **Remove from the senescence signature (and any quiescence signature):**
  `Cdkn2a (p16), Cdkn2b (p15), Cdkn2c (p18), Cdkn2d (p19), Cdkn1a (p21), Cdkn1b (p27), Cdkn1c (p57)`.
- **Score the CKIs only as their own separate modules** (INK4 and CIP/KIP, below).
- The real question becomes: *with the shared CKIs removed, does the rest of the senescence program
  (SASP, DDR, Lmnb1-loss, etc.) light up in the arrested / high-CKI cells — or not?*
  - **Not** → high-CKI MB cells carry CKIs **without** the senescence program → reversible quiescence,
    model supported.
  - **Yes** → a genuine senescent state exists → scope the model to the proliferative compartment.

---

## 3. Datasets & cell subsetting

- Use the GNP-lineage / SHH-MB scRNA-seq you know best (Ptch⁺/⁻ MB and developing cerebellar GNPs;
  include ± vismodegib and ± EZH2i conditions if available).
- **Subset to the malignant SHH-MB / GNP-lineage cells FIRST.** SASP/senescence gene sets are dominated
  by immune/stromal expression — microglia, macrophages, and T cells will score "senescent" if left in.
  Exclude immune/stroma/endothelial before scoring. (Flag the fraction removed.)
- Mouse symbols throughout (orthologize senescence sets to mouse).

---

## 4. Signatures (score per cell; CKIs already removed where noted)

**A. Senescence — CKIs removed (the key readout)**
- **SenMayo** (Saul et al. 2022, *Nat Commun*; mouse list) — minus any CDKN genes.
- **Core SASP** (subset, since this is intra-tumoral): `Il6, Il1a, Il1b, Cxcl1, Cxcl2, Cxcl12, Ccl2,
  Ccl20, Serpine1 (Pai1), Igfbp3, Igfbp4, Igfbp7, Mmp3, Mmp12, Timp1, Tnf, Nfkb1`.
- **Senescence structural/effector:** `Lmnb1` (expect **DOWN** in senescence — score separately or
  invert), `Glb1, Hmga1, Hmga2, Trp53` and p53-target score.
- Optionally also: CellAge / Fridman-Tavazoie / GO:0090398 (cellular senescence) — CDKN-stripped.

**B. CDK inhibitors — scored SEPARATELY (the model brake)**
- **INK4 (→ CDK4/6):** `Cdkn2a (p16), Cdkn2b (p15), Cdkn2c (p18), Cdkn2d (p19)`.
- **CIP/KIP (→ CDK2):** `Cdkn1a (p21), Cdkn1b (p27)` (report `Cdkn1c/p57` but keep aside — separate project).
- Report per-gene **detection rate**; `Cdkn2a` is dropout-prone, so trust the **module score**, and
  note whether p16 specifically is broad vs restricted.

**C. Proliferation / cell cycle**
- `Mki67, Pcna, Top2a, Ccnb1, Ccna2, Cdk1, Birc5, Mcm2-7`.
- `Seurat::CellCycleScoring` (Tirosh S/G2M sets) → `S.Score, G2M.Score, Phase`.
- **E2F-target score** (Hallmark `E2F_TARGETS`) — the model's Rb-E2F commitment readout.

**D. Quiescence / transient-G0 — the model's reversible arrest (CKIs removed)**
- Operational, *model-defined* cell call: **G0/G1 phase AND low E2F-target/low Mki67 BUT mitogen-rich**
  (high `Ccnd1` + high HH/Gli drive, section F) — i.e. *paused but poised to re-enter*.
- Plus a published quiescence signature (e.g. Cheung & Rando, or Coller G0) with CDKN genes removed.
- The discriminator vs senescence: transient-G0 cells **retain proliferative drive** (Ccnd1/Gli high)
  and are **not** SASP⁺.

**E. EZH2 / PRC2 activity**
- `Ezh2, Eed, Suz12, Ezh1` + a PRC2/H3K27me3-target de-repression score (classic senescence = EZH2 loss
  → p16 up, so expect **EZH2-low where senescence-high** — a clean orthogonality test vs MB which is
  EZH2-high).

**F. Hedgehog / CyclinD1 drive (the model's "drive" axis)**
- `Gli1, Gli2, Ptch1, Ptch2, Hhip, Boc, Sfrp1, Mycn, Ccnd1, Ccnd2, Atoh1` (+ a HH-target module score).

**G. Differentiation — EXPLORATORY, DISREGARD FOR THIS PAPER**
- Granule-neuron program: `Neurod1, Rbfox3 (NeuN), Tubb3, Map2, Cntn2, Gabra6, Grin2b, Grin2c, Reln,
  Cntn1, Sema6a`. Score it and note any anti-correlation with proliferation/EZH2, **but flag clearly
  as out-of-scope** for the current paper (the differentiation switch is not in the model).

**Scoring method:** prefer **UCell** or **AUCell** (robust to dropout/depth) over AddModuleScore; report
which. Regress out / report depth and condition as covariates.

---

## 5. Analyses & the model's prediction for each

1. **Module-score correlation matrix** across {senescence(no CKI), SASP, INK4, CIP/KIP, proliferation,
   E2F, quiescence/transient-G0(no CKI), EZH2/PRC2, HH-drive, differentiation}.
   - *Predicted:* senescence(no CKI) **anti-correlates** with proliferation/E2F; INK4 & CIP/KIP
     correlate with quiescence/transient-G0 but **weakly/not** with senescence(no CKI); EZH2 **high
     where senescence low**.
2. **Is there a discrete senescence-high cluster?** Its size, cycling status (Mki67/phase), SASP, EZH2,
   HH-drive. *Predicted:* if present, it's **small, non-cycling, EZH2-low, mitogen-low** — distinct from
   the bulk transient-G0 cyclers.
3. **The high-CKI cells** (top INK4 / CIP-KIP quantile): are they **cycling (Mki67⁺/S-G2M)** and
   **senescence(no CKI)-LOW**? *This is the core test.* Predicted: yes.
4. **Tonic vs subpopulation:** within the *proliferating* compartment, is the INK4 tone **broadly**
   elevated (supports the model's tonic brake) or confined to few cells? Compare **p18 (expected broad)
   vs p16 (possibly restricted/senescence-tail)** — if p16 is restricted, the working brake on cyclers
   is p18/p15, which the model already separates.
5. **Vismodegib (if available):** does HHi shift cells **cycling → transient-G0** (quiescence↑,
   proliferation↓, drive↓) **without** inducing senescence(no CKI)? *Predicted:* yes = "vismo arrest is
   reversible quiescence, not senescence." If **EZH2i** is available, does it push back toward cycling?
6. **GNP developmental axis (if in data):** as GNPs exit (P7→P14-like), do cells move
   cycling → quiescent/differentiating with **EZH2 falling alongside the proliferation drop** (model:
   EZH2 is S-phase-coupled, so it should decline *with* cycle exit, not clearly before) — and crucially
   **without** a senescence program? Useful contrast to MB.

---

## 6. Caveats to honor
- `Cdkn2a`/p16 dropout — rely on module scores; report detection.
- Snapshot ≠ reversibility: argue transient vs terminal via the **poised** proxy (mitogen/Ccnd1/Gli
  retained) and via the trajectory (cycling↔quiescent vs a one-way senescent branch).
- **Immune confound on SASP is the biggest trap** — subset to tumor/GNP lineage before scoring.
- Mouse↔human ortholog mapping for the senescence sets.

---

## 7. Deliverables
1. **Correlation heatmap** of the module scores (the headline; senescence-no-CKI vs everything).
2. **UMAP** colored by: senescence(no CKI), INK4, CIP/KIP, proliferation, quiescence/transient-G0, EZH2,
   HH-drive.
3. **The key panel:** senescence(no CKI) **vs** proliferation (and vs transient-G0) — scatter / 2-D
   density — to show they're orthogonal (the model's claim) or not.
4. **Table:** among high-INK4 / high-CIP-KIP cells, fraction cycling vs fraction senescence-high.
5. (If ± drug) the vismo → quiescence (not senescence) shift, ± EZH2i reversal.
6. One-paragraph verdict: **does the MB proliferative compartment carry the CKI brake without the
   senescence program (model supported), or is there a distinct senescent state (scope the model)?**

---

*Note for the modeler (not the agent): whatever the result, the paper text gets a sentence + one supp
panel — either "the modeled arrest is reversible quiescence, not senescence (senescence program is
orthogonal/absent in the cycling compartment)" or "senescence is a separate minority state; the model
describes the proliferative compartment." The CKI-removal rule (§2) is what makes either statement
non-circular.*
