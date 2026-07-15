# Transient G0 in GNP and SHH-medulloblastoma — evidence synthesis

*JP, 2026-07-15. Consolidated from the transient-G0 investigation (GNP/MB quiescence RIS library + Ki67 RIS
library), with the CDK4/6-inhibitor data folded in. Kept to what the evidence supports; **estimates vs measured
values are flagged explicitly**. This is the biological anchor for the corresponding model revision (see "Model
implications" below and `src/build_model_v44_heldt.py`).*

## 1. GNP outer EGL: transient G0 is effectively zero *(firmest conclusion)*

For the bulk transit-amplifying population this is the firmest conclusion. The classic p27 layering (Miyazawa)
puts PCNA and cyclin D1 in the outer two-thirds of the EGL with p27 low there, p27 rising only in the middle and
inner EGL as a prelude to **one-way exit**. The outer EGL runs a high growth fraction (Fujita; Mares & Lodin). The
deterministic outer-to-inner model is looser than the textbook (Hanzel) but not in a way that introduces
reversible pausing. The Spencer CDK2-low transient-G0 framework is real but was built in **cell lines and is
p21-driven**, whereas the GNP-relevant Cip/Kip is **p27**. Direct read of the sections confirms very few p27+ cells
in the outer EGL. So bulk GNP behaviour is **obligate cycling with a graded p27 ramp toward permanent exit**, not a
bistable dip into and out of G0. **Baseline transient G0 in normal GNP homeostasis is essentially nil.**

## 2. MB: G0 is real but hard to quantify, and p27 does NOT measure it

This is where p27 collapses as a quantification tool. Tissue papers (Adesina; Bhatia) framed p27 as a
differentiation marker, but the images show widespread staining — expected, because **p27 IHC reports TOTAL
protein**, which sits in cycling cells as cytoplasmic, phosphorylated, or G1-phase p27 that is *not*
CDK-inhibitory. Once p27 is widespread and Ki67 is substantial, the two populations cannot be disjoint → extensive
p27/Ki67 co-expression. A decent fraction of **Atoh1+ (proliferation-competent) cells are p27-positive**, nailing
this from the progenitor side. **Consequence: neither the presence nor the extent of p27 staining can estimate the
G0 fraction in MB tissue.**

## 3. Ki67 numbers and the partition

- Human MB whole-tumour Ki67 ≈ **30%** across cohorts (Zhao 30.0 ± 11.3%; a pediatric MIB-1 mean near 30%) —
  **measured**. Caveat: SHH ran below Group 4 in Zhao; SHH desmoplastic/nodular tumours are **bimodal** (nodules
  very low, internodular up to ~70%; McManamy & Ellison; Onda ~56% in active areas).
- Murine SHH-MB is far higher and more uniform: ~**54%** by scRNA-seq (Zhang) up to ~**86%** Ki67+ in SmoM2
  (Huang) — the mouse models lack the extensive differentiated-nodule compartment.

**Partition of the human tumour** (reasonable reading): ~30% cycling (Ki67, measured), a large differentiated
fraction (~half), and a **reversible G0 remainder ~20%** — held **exactly as an informed guess, NOT a
measurement**. It is hard *structurally*: Ki67− / p27+ describes *both* differentiated and reversibly-quiescent
cells, and reversibility cannot be seen in fixed tissue. The only thing that splits the non-cycling fraction into
permanent exit vs reversible G0 is **differentiation status**.

## 4. What would actually quantify it *(the panel)*

- **Cycling axis:** Ki67, or better **phospho-Rb** (hyperphospho = cycling; hypophospho = genuinely arrested,
  a more direct functional readout than p27, independent of p27 localisation). **The pRb stain = Cell Signaling
  phospho-Rb Ser807/811, clone D20B12 (#41359)** — detects Rb *hyper*phosphorylation (not mono-phospho), a
  committed/cycling marker that maps onto the model's `pRb` (hyper) species. NB Ser807/811 is a CDK4/6 *and* CDK2
  site (only Ser780 is CDK4/6-exclusive), so a **residual Ser807/811⁺ under CDK4/6i = CDK2-bypass** (Rb kept
  hyperphospho by CDK2); a **Ser780 co-stain** would separate CDK2-bypass (Ser780⁻) from incomplete inhibition
  (Ser780⁺). Caveat: "hypophospho / D20B12⁻" = true G0 **plus** early-G1 mono-phospho committed cells — neither the
  stain nor the model isolates pure G0 (see docs/cdk46i_resistance_research.md).
- **Differentiation axis:** NeuN + NeuroD1, Tag1, MAP2, βIII-tubulin.
- **p27:** a *supporting* stain, not a definer.
- **Permanent exit** = NeuN+/NeuroD1+. **Reversible G0** = the operational intersection: Ki67−, pRb-hypophospho,
  differentiation-marker−, p27+. Scored per cell, this converts the ~20% guess into a real number and separates it
  from the differentiated bulk.

## 5. The inducible reversible G0 (CDK4/6-inhibitor data) — the key result

At baseline the reversible-G0 state is near-zero in GNP and modest in MB, but the **latent capacity is real and
pharmacologically accessible**. CDK4/6-inhibitor data (consistent with the Chahin palbociclib results): cells
**accumulate p27 extensively** (the expected CDK4/6i phenotype — p27 redistributes onto CDK2, Rb goes
hypophosphorylated) and critically they **re-enter the cycle** → true reversible quiescence, not senescence or
differentiation. This maps directly onto the Rb / E2F / cyclin D–CDK4/6 axis. So the correct statement is **not
that MB lacks a transient-G0 program, but that the program is minimally occupied under baseline mitogenic drive and
becomes the dominant, reversible state when CDK4/6 is inhibited.** That reversibility is also the **therapy-escape
liability** — the same logic as the SOX2 / OLIG2 quiescent reserves that survive antimitotics and repopulate.

## 6. Model implications

The reversible p27-high G0 arm belongs in the model, but **gated on CDK4/6–cyclin D activity** and tuned so it is
nearly empty at baseline (**GNP ~0, MB ~20%**) and **fills when CDK4/6 is suppressed**. Two representational points:

1. **Split p27** into a widespread **total pool** (rides along with cycling and Atoh1+ cells; does *not* drive
   arrest) and a **CDK-inhibitory nuclear pool** (the functional G0 driver); only the second couples to cell-cycle
   state.
2. Use **pRb phosphorylation state as the functional cycling/arrest switch**, not p27 level — that is what the
   marker logic and the CDK4/6i data actually track.

Biologically there are **three destinations out of cycle — cycling, reversible CDK4/6-gated G0, and permanent
differentiated exit — and only the first two interconvert.** **Scope note (2026-07-15):** the permanent
differentiation exit (the NeuroD1 / H3K27me3-linked one-way path, the *larger* non-cycling fraction) is
**deliberately NOT modeled here** — it is a separate compartment with its own biology and is out of scope for this
paper. This model represents the **proliferation-competent pool** and the two fates that interconvert within it:
**cycling ↔ reversible CDK4/6-gated G0**. The reversible-G0 fraction (~20% in MB) is therefore reported *within*
that pool, not as a fraction of the whole tumour (which also carries the unmodeled differentiated bulk and dilutes
the whole-tumour Ki67 to ~30%).

### Calibration anchors for the model
- **GNP baseline G0 ≈ 0** and the **CDK4/6i-inducible fill** are the FIRM anchors — calibrate to these.
- **MB baseline G0 ~20%** is a SOFT target (informed guess) — do not over-fit.
- Human MB Ki67 ~30% (measured); differentiated fraction ~half; reversible G0 ~20% (guess).

## 7. Two therapy-escape routes require heterogeneity *(2026-07-15, JP cell-culture data)*

Two observations do **not** fall out of a single homogeneous MB population, and forcing them from one deterministic
cell fails on direct test — each needs a distinct, grounded subpopulation:

1. **MB carries ~20% reversible G0 at baseline.** A homogeneous MB does **not** sit in baseline G0: it is
   **Hh-autonomous** (Ptch1-low + MYCN drive cyclinD), so raising INK4 (p18 up to 8×), lowering SHH, or lowering
   MYCN all keep it cycling (verified sweeps — the mitogen route is a dead end). What **does** produce a reversible,
   pRb-hypophospho, p27-high baseline G0 is a **high CDK-inhibitor tone** (kSyP21 ≈ 0.11). So the reserve is a
   **high-CKI subpopulation** — the SOX2 / OLIG2 quiescent-reserve logic (survives antimitotics, re-enters).
2. **~10–20% of MB cells keep DIVIDING under CDK4/6i** (JP, cell culture). A homogeneous MB fully arrests
   (p27 0→1.38, pRb→0). What escapes is a **low total-Rb subpopulation** (Rb pool ≈ 0.5): with little Rb to hold
   E2F, cells cross the restriction point **independent of cyclinD–CDK4/6** and keep dividing under the drug
   (verified). This is the canonical **RB1-loss / low-Rb CDK4/6i-resistance** route.

**GNP has neither** (no high-CKI reserve, full Rb) → it **arrests cleanly** under CDK4/6i (~2% baseline G0 = the
low-CyclinD1-abundance tail crossing the pRb threshold — *not* differentiation, which is out of scope; ~0% escape).
The MB−GNP gap is the **therapy-escape liability**: under CDK4/6i **~35% of MB persist** — a ~20% quiescent reserve
(arrested but reversible: survives the drug and re-enters when the quiescence signal is released) **plus** a ~15%
resistant fraction that keeps dividing. **Only the ~15% actively proliferate during the drug**; the ~20% is a
surviving regrowth reservoir, not active escape.

**What the model does vs does not claim (important — this is a MIXTURE model, not an emergent-fraction prediction).**
`simulations/sim_cdk46i_escape.py` composes MB as 65% normal / 20% reserve / 15% resistant (GNP 100% normal).
- **The subpopulation SIZES (20% / 15%) are hand-set — JP's culture ESTIMATES, a SOFT target (§3, §6), not a model
  output.** The reported baseline-G0 and still-dividing percentages therefore **equal those inputs**; the per-cell
  CyclinD1-abundance heterogeneity is included for realism but does **not** move the classification (the state
  overrides dominate), so it does not make the magnitudes emergent. The GNP-vs-MB contrast is likewise an
  **assumption** (GNP lacks both subpopulations — normal tissue, no RB1 loss, no tumour quiescent reserve), not a
  derivation.
- **What the model VERIFIES (genuine ODE outputs):** each mechanism produces its claimed phenotype — high-CKI →
  reversible baseline G0; low-Rb → keeps dividing under CDK4/6i; normal → arrests under CDK4/6i; and **reversibility**
  is confirmed by a clean arrest-then-release run (`--reversibility`): the reserve re-enters on CKI release and the
  drug-arrested normal cell re-enters on washout, **both 0 divisions during arrest → 12 after release**. So the
  model asserts these mechanisms **are sufficient**, at the assumed sizes, to reproduce the observed contrast.
Figure `simulations/fig_cdk46i_escape.png`.

**Falsifiable predictions** (these would make the fractions measurable rather than assumed): reserve cells =
p27-high / pRb-hypophospho / SOX2⁺ (or OLIG2⁺), Ki67⁻, re-enter on release; resistant cells = Rb-low / RB1-loss,
Ki67⁺ **during** CDK4/6i. Distinct stains → the two escape routes are separable experimentally.
