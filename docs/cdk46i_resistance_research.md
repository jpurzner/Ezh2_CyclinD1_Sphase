# CDK4/6-inhibitor resistance & reversibility — literature research (for the v44 palbo model)

*2026-07-15. Commissioned to reconcile JP's palbociclib data (in culture: pRb⁺ cycling drops to 16–18% under palbo,
rebounds to 50% on washout after 72h) with a resistance mechanism for the model. Full sourced summary; the
modeling-relevant conclusions are in **bold**. Sources are peer-reviewed (PMIDs/DOIs inline).*

## Bottom line (the correction)
The residual pRb-**positive** cycling under palbociclib is **CDK2 (cyclin E/A–CDK2) re-phosphorylating Rb
independently of CDK4/6** — NOT the drug failing to hit CDK4/6, and NOT alternative pocket proteins. It is present
first as a non-genetic adaptive subpopulation and later selectable into a CCNE1/CCNE2-amplified clone.

1. **Alternative pocket proteins (p107/RBL1, p130/RBL2) are NOT a resistance route** — the opposite. All three
   pocket proteins are CDK4/6 substrates, co-dephosphorylated by palbo; p107/p130 build the **DREAM** repressor
   (Litovchick; Mages *eLife* 2017;6:e26876) which CDK4/6i *stabilizes* → **deepens** arrest. That RB1 loss itself
   causes resistance proves p107/p130 do not substitute for RB1. Where they enter resistance: CDK2 *inactivates*
   p130/p107 (Farkas *JBC* 2000; PMID 10906146) — they are victims of the CDK2 node, not a bypass. **Do not build a
   p107/p130 term.**

2. **Resistance landscape, by pRb status** (prevalence from Wander *Cancer Discov* 2020, PMID 32404308; 41 resistant
   HR+ MBC biopsies): cyclin E–CDK2 activation (CCNE2 amp 14.6% + adaptive; the dominant convergent node) →
   **pRb⁺**; non-genetic/adaptive persister (34% had NO genomic driver) → **pRb⁺**; RB1 loss (~5–10%) → **pRb⁻**
   (only common pRb-negative route); AURKA (26.8%), AKT/PI3K, CDK6 amp, RTK/MAPK → all **pRb⁺**. **Nearly every
   common mechanism leaves cells pRb-POSITIVE; a pRb⁺ residual is diagnostic of CDK2/cyclin-E bypass, not RB1 loss.**
   (Asghar *Clin Cancer Res* 2017 PMID 28606920 — single-cell: resistant cells exit mitosis into high-CDK2,
   bypassing the R-point; Herrera-Abreu *Cancer Res* 2016 PMID 26977878 — non-canonical cyclinD1–CDK2/cyclinE1
   *recover* pRb under palbo.)

3. **Incomplete inhibition — in vitro NO, in the brain YES.** Palbo fully engages CDK4/6 at experimental doses
   (pRb IC50 ~60 nM vs ~0.5–1 µM used; Fry *Mol Cancer Ther* 2004 PMID 15542782) → in **culture**, residual pRb⁺ is
   NOT incomplete target inhibition (it's CDK2 bypass). **BUT in the CNS, palbo is a P-gp/BCRP substrate → genuinely
   subtherapeutic brain exposure** (de Gooijer *Invest New Drugs* 2015 PMID 26123925 — ~115× in transporter-KO;
   Parrish 2015 PMC4613960 — subtherapeutic at the invasive edge). So "make the drug worse" is physically real **for
   MB in vivo (BBB efflux)**, a separate axis from the in-vitro CDK2 bypass.

4. **Reversibility (matches JP's washout data).** Short/moderate palbo (1–2 d) → G1 arrest **fully reversible within
   24 h of washout**, near-synchronous re-entry (Crozier/Foy *EMBO J* 2022 e108599; Kim 2020 PMC7653349). Prolonged
   (>3–4 d) progressively converts a **p53/MDM2/mTOR/cell-size-dependent** fraction to irreversible senescence
   (Foy *Mol Cell* 2023; Klein/Koff *Cancer Cell* 2018 PMID 30017615). Reversible-vs-senescent split is
   context-dependent, not fixed.

5. **SHH-medulloblastoma — direct validation of the model's reversible-reserve story.** Cook Sangar
   *Clin Cancer Res* 2017 (PMID 28637687), SHH-MB PDX: palbo → **~97% G1 arrest**, pRb inhibited (Ser780/807/811),
   ~63% tumor reduction — **but 15/20 flank tumors recurred within 60 days of withdrawal**, Ki67⁺ residual cells in
   scar. **Cytostatic, reversible arrest with a quiescent reserve that regrows** — exactly the model's reserve +
   reversibility. **CDK6 amplification is the primary recurrent SHH-MB lesion.** CDK4/6i also drives neuronal
   differentiation + collapses transcriptional heterogeneity (Neuro-Oncol 2023 PMID 37127652). Combo with mTORC1i
   (sapanisertib) + nanoparticle palbo (CNS PK) prolongs survival (Sci Adv 2022, 10.1126/sciadv.abl5838).

## pRb readout — antibody definition (what "pRb" means in JP's data)
JP's pRb stain = **Cell Signaling phospho-Rb Ser807/811, clone D20B12** (#41359). Key properties:
- Detects **Rb HYPER-phosphorylation** ("multiply phosphorylated forms, NOT mono-phosphorylated Rb"); Ser807/811
  phospho is a validated **hyperphosphorylation / E2F-active (committed, cycling) marker**. → maps onto the model's
  **`pRb` (hyper species)**, so the reframe's G0 marker (`pRb` < threshold = hypophospho) is the correct analog.
- **Ser807/811 is phosphorylated by BOTH CDK4/6 (priming) and CDK2 (completion)** — only **Ser780** is
  cyclin-D/CDK4/6-*exclusive*. So a **residual Ser807/811⁺ under palbo = Rb kept hyperphosphorylated by CDK2 =
  the CDK2-bypass signature** (not incomplete CDK4/6 inhibition, which would instead drop Ser780).
- **Shared caveat (antibody AND model):** "hypophospho / D20B12-negative" = un-hyperphosphorylated = true G0 **plus**
  early-G1 mono-phosphorylated (Rbm) committed cells. Neither the stain nor the model cleanly isolates pure G0 from
  early-committed G1; the model uses the same convention, so they are consistent (the "G0" fraction is really
  "hypophospho-Rb").
- **Discriminating experiment:** co-stain **Ser780** (CDK4/6-specific). CDK2-bypass cells = Ser780⁻ / Ser807-811⁺;
  incomplete CDK4/6 inhibition = both⁺. This directly separates the two residual-mechanism hypotheses.
Refs: Cell Signaling #41359; CDK4/6-initiates-Rb / CDK2-commitment, *Sci Rep* 2022 (s41598-022-20769-5);
Cdk4/6 phospho-gradient in G1, PMC8290049.

## Implication for the v44 model
- The escape-sim "resistant = RB1-loss / low-Rb (**pRb-negative**)" label is **the wrong sign** for JP's pRb⁺ data —
  the residual should be **CDK2/cyclin-E bypass (CCNE, pRb-positive)**.
- **But a proper pRb⁺-*dividing* CCNE clone is hard to build in v44** (tested 2026-07-15): a constitutive
  (E2f-independent) CDK2 term raises pRb (→ pRb⁺) but does NOT restore division — a constant Rb-phosphorylation term
  fights the oscillator (Rb never resets), so the cell stalls pRb⁺-but-arrested even with p27 clearance gated on
  CDK2. An E2f-*gated* CyclinE term can't bootstrap (E2f is off under palbo arrest). So capturing pRb⁺ escape needs
  real surgery on the R-point/oscillator — the elaborate resistance JP said not to bother with.
- **"Make the drug worse" (partial CDK4/6 block) is all-or-nothing** in v44 (the R-point is a sharp bistable switch:
  ≤60% residual → all arrest, ≥70% → all cycle), so a uniform partial block can't produce a graded 16–18% via
  CyclinD heterogeneity. It IS physically justified for MB via BBB efflux, but doesn't give the graded residual alone.
- **What v44 captures cleanly + what the MB literature validates: reversible cytostatic arrest + a reversible
  quiescent reserve that regrows on withdrawal** (Cook Sangar). That is the defensible, evidence-backed story; the
  pRb⁺ CDK2-bypass residual is best NOTED as a literature-grounded resistance phenomenon rather than dynamically
  modeled.
