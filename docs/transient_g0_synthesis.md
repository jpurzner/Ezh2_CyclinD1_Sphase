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

## 7. CDK4/6i is a DOSE-RESPONSE, not a resistant clone *(2026-07-15, JP palbociclib data)*

**The reframe (supersedes an earlier resistant-clone reading).** JP's pRb-Ser807/811⁺ data: the dividing fraction
falls to **16–18% at 1 µM** (72h, washout → rebounds to 50%: reversible) and **~2% at 5 µM** (saturating; slower to
wash out but still reversible; the extra 5 µM effects are off-target toxicity). This is a **DOSE-RESPONSE**, not a
genetically-resistant subpopulation. Modeling palbo as **graded residual CDK4/6 activity** (`kPhRbCd = KP·resid`)
over a population with CyclinD1 (CV 0.70) + CKI (CV 0.42) heterogeneity, three things fall out of v44 as-is —
**no CDK2-bypass "option B", no senescence compartment** (`simulations/sim_cdk46i_dose_response.py`,
`fig_cdk46i_dose_response.png`):

1. **Dose-limited, pRb-positive residual.** As dose rises the dividing fraction slides down a graded curve, and
   **every dividing cell is pRb-Ser807/811⁺** — the high-CyclinD tail overcomes partial inhibition and keeps
   hyperphosphorylating Rb. The 16–18% residual = a *sub-saturating* dose; it vanishes at saturating dose. This is
   why the residual is pRb-POSITIVE without needing CDK2-bypass. (Consistent with subtherapeutic CNS/BBB exposure
   for MB in vivo — palbo is a P-gp/BCRP substrate.)
2. **Reversible quiescent reserve** (drug-INDEPENDENT): a distinct high-CKI baseline-G0 fraction (MB ~20%, GNP ~2%;
   SOX2/OLIG2-like) that sits out of cycle at every dose and re-enters when its quiescence signal is released
   (Cook Sangar 2017: SHH-MB regrows on withdrawal). This is a **distinct bistable state, NOT a distribution tail**
   — a swept CKI/CyclinD1 distribution only reaches ~20% baseline G0 at an extreme, effectively bimodal width
   (Cip/Kip tail 30–119× median), which *is* a separate mode. GNP lacks it (~2% = its small low-CyclinD1 tail).
3. **Reversibility / re-entry kinetics.** Arrest is reversible at every dose; **deeper arrest (higher dose) → slower
   washout re-entry** (more p27 to clear: ~1180 model-units at 38% inhibition → ~1555 at 100%), which reproduces
   JP's "5 µM takes longer to wash out" with **no senescence state**.

**What the model claims.** The dose-response *shape* and the pRb-positive sign of the residual are **emergent ODE
outputs** (the residual IS the high-CyclinD tail; dividing == pRb⁺ at every dose was verified True). The reserve
FRACTION (~20%) is a hand-set SOFT estimate (§3, §6) — a distinct state, not derived.

**Fit to JP's two doses (2026-07-15).** Only two doses measured (1 µM → 16-18% pRb⁺; 5 µM → ~2%, saturating;
nothing below 1 µM). Mapping concentration → residual CDK4/6 activity by `resid = K/(K+C)` and sweeping heterogeneity
width, the model fits BOTH points at an **empirical K ≈ 1.5 µM** (unrelated to the ~11 nM in-vitro IC50 — JP: the
biochemical value is routinely off from the cell-culture dose; no CDK6 amplification). **KEY caveat/prediction:**
the transition is intrinsically **STEEP** (IC50 ~0.7 µM, essentially saturated by ~2 µM) and **widening the
CyclinD1/CKI heterogeneity does NOT soften it** — the bistable R-point makes each cell near-binary and sets a hard
wall, so CyclinD spread can't spread the population dose-response. This is a **falsifiable prediction**: intermediate
doses (2-3 µM) should already be ~saturated (~2%), not still dropping toward 5 µM. If a real 2-3 µM shows
intermediate values, the true curve is shallower than v44 predicts → there are heterogeneity sources beyond
CyclinD1/CKI (CDK4/6 abundance, drug uptake) that v44 doesn't carry. `fig_cdk46i_dose_response.png` is on the µM axis.

**Falsifiable / discriminating experiments.** (a) A **Ser780 co-stain** separates the two residual mechanisms:
dose-limited residual = Ser780⁻/Ser807-811⁺ at *sub-saturating* dose (partial CDK4/6 still primes)… actually both
sites fall together under partial inhibition, so the cleaner test is the **dose-response curve itself** — a graded,
reversible, pRb⁺ residual that disappears at saturating dose is dose-limitation; a residual that **persists at
saturating dose and grows over weeks** is genuine acquired resistance (CCNE1 amplification → "option B", see
`docs/option_b_cdk2_bypass_reading.md`). (b) Reserve cells = p27-high/pRb-hypophospho/**SOX2⁺**, Ki67⁻, re-enter on
withdrawal. **Note (2026-07-15):** the old `sim_cdk46i_escape.py` resistant-clone framing is RETIRED by this section.
