# Transient G0 & Spontaneous Quiescence (GNP vs SHH-MB)

> The paradigm JP wants the model to hold: **GNPs can *extend/graded-lengthen* G1 but do not spontaneously
> G0-arrest** (they cycle or they differentiate-exit), whereas **SHH-MB is *competent* for reversible G0
> arrest** — a quiescent, therapy-resistant compartment. This page reviews (1) what the biology says the two
> quiescences are, (2) **what the current model actually does** (measured), (3) the gap between them, and (4)
> the concrete module changes needed. It incorporates a commissioned literature review (JP 2026-08-01) as the
> evidence base. See also [Mechanisms Explored](04_mechanisms_explored.md) §8, [Model Architecture](02_model_architecture.md)
> §4H, and memory `cdki-baked-g0-g1-behavior`, `mb-celltype-transient-g0-parameterization`, `transient-g0-literature-constraints`.

## TL;DR

- **There are TWO distinct quiescences and the model conflates neither correctly.** (A) A **shallow, stochastic,
  reversible mitotic-exit G0** that is **p21 (CDKN1A)-specific and p21-necessary-and-sufficient** (Spencer 2013,
  Overton 2014, Barr 2017). (B) A **deep, lineage-programmed dormancy** of the SOX2⁺/OLIG2⁺ MB stem compartment,
  best-associated with **p57 (CDKN1C)/OLIG2 — not p21, not p27** (Vanner 2014; Zhang 2019; Desai/Li 2025; a 2025 p57 preprint).
- **The current model builds neither.** Its commitment is a **sharp p27/CdP21 R-point switch**: GNP and MB both run
  with **G0-fraction = 0 %**; MB's longer cycle is a **longer G1 set by high INK4**, not a p27-high pause; and p21/p57
  are carried **inert**. GNP mitogen response is **go/no-go** (fixed G1, then abrupt full arrest), not graded G1.
- **The model's quiescence CKI is the wrong one.** It routes everything through **p27**, but the literature says the
  stochastic transient-G0 bifurcation is **p21-driven and p27 does not reproduce it** (p27 failed to discriminate the
  two daughter fates in Spencer's single-cell panel). The recently-baked **p21-PCNA un-map** now wires p21 correctly
  (the foundation), but p21 is still held low.
- **Required:** two separate modules — **Module A** (p21 mitotic-exit bifurcation + a population/birth-p21 layer) and
  **Module B** (p57/OLIG2 deep dormancy) — plus **reparameterizing p27 as dual-function** (pY88-p27 *activates* CDK4;
  Guiley 2019) and a **quiescence-depth** state variable.

---

## 1. The biological paradigm — two quiescences, different CKIs

### 1a. GNP: obligate cycling with graded exit, essentially no spontaneous G0
Outer-EGL GNPs proliferate for days under sustained Shh, then exit **one-way** into differentiation (Wechsler-Reya
& Scott 1999). Developmental studies report **graded G1/cycle lengthening** approaching exit — but the one
growth-fraction-unbiased measurement (Legué clone S-phase index) points the other way, and much "graded G1
lengthening" is a growth-fraction artifact (memory `gnp-g1-lengthening-buffering`, wiki 07 §7.2). The operational
GNP picture the model should hold: **cycle-length/G1 is graded-tunable, but there is no stochastic reversible G0**
— cells cycle or they commit to exit.

### 1b. SHH-MB: competent for TWO reversible G0 states
- **(A) The p21 mitotic-exit bifurcation** — shallow, stochastic, reversible, cell-autonomous. Asynchronous cells
  bifurcate at mitotic exit into **CDK2-inc** (immediate S-entry, M→S ≈ 5–7 h, CV≈0.26) vs **CDK2-low** (Rb
  hypophosphorylated, mitogen-sensitive transient G0 until a second restriction window). Sister fates are shared
  98 % of the time → **the decision is set before division**, by **p21 inherited from the mother's G2** (produced by
  endogenous replication stress via p53). p21 is **necessary and sufficient**; CDK2-low fraction is a **graded
  function of p21 dose** (0→100 %). (Spencer 2013; Overton 2014; Barr 2017; Arora 2017; ODE: Heldt/Barr/Novák 2018.)
- **(B) Lineage-programmed deep dormancy of the SOX2⁺/OLIG2⁺ stem compartment.** In Ptch1⁺ᐟ⁻ SHH-MB, <5 % SOX2⁺
  cells are the quiescent, therapy-resistant medulloblastoma-propagating cells (Vanner 2014). OLIG2⁺ progenitors are
  transit-amplifying at onset but become **quiescent, slow-cycling, stem-like** in full tumors (Ki67⁻; Zhang 2019).
  **OLIG2 is the quiescence↔activation switch** (Desai 2025, Li 2025; the OLIG2 inhibitor CT-179 forces
  differentiation/quiescence and synergizes with palbociclib). **Crucially the CKI wiring is inverted vs the p21
  paradigm: OLIG2 directly represses CDKN1A/p21 (Ligon 2007), so this compartment is quiescent *despite low p21*.**
  The best MB-specific evidence points to **p57 (CDKN1C)** as the operative quiescence CKI here (nuclear p57 enriched
  in SOX2⁺/Nestin⁺ cells; inducible p57 → 6× G0 enrichment + vincristine resistance; 2025 preprint, *non-peer-reviewed*).

### 1c. The key CKI verdicts (which the model must respect)
- **p21 = the transient-G0 CKI** (necessary & sufficient for the stochastic bifurcation).
- **p27 ≠ the transient-G0 CKI.** p27 is competent for *mitogen-withdrawal/contact-inhibition* arrest and *deep
  dormancy*, but **did not discriminate the daughter-fate bifurcation** (Spencer panel); no experiment shows p27
  substituting for p21 in a CDK2-sensor bimodal system (a genuine, modelable gap). In MB genetics p27 is
  **progression, not initiation** (Hatton 2013: p27 loss doesn't change tumor frequency; haploinsufficient; not a
  reproducible biomarker).
- **p57 = the deep-dormancy CKI** (HSC/NSC dormancy; the MB stem candidate) — slow kinetics, high exit threshold.
- **p27 is dual-function:** unphosphorylated it inhibits CDK4/6+CDK2; **Tyr88/89-phosphorylated it ACTIVATES
  CDK4–cyclinD1** (Guiley 2019) and is palbociclib-insensitive — so "p27-high but proliferating" is real.

---

## 2. What the current model actually does (measured, baked default 2026-08-01)

Direct single-lineage characterization (`scratchpad/g0_characterize.py`, current baked model = un-map + species-recal):

| | cycle | G1 dwell | **G0-frac** | pRb min | p27 (pre-S) | p21a | CdP21 |
|---|---|---|---|---|---|---|---|
| **GNP** (SHH 0.5) | 17.6 h | 10.8 h | **0 %** | 1.66 | 0.007 | **0.000** | 0.57 |
| **MB** (SHH 0.5) | 26.1 h | 21.1 h | **0 %** | 1.68 | 0.004 | **0.000** | 0.61 |

GNP mitogen dose-response: **fixed** cycle (17.6 h) / G1 (10.8 h) from SHH 1.0 down to 0.35, then **abrupt full
arrest** at SHH ≤ 0.30 (G0-frac jumps 0 %→100 %). No graded G1 lengthening; no transient G0.

**What this means:**
1. **No discrete/transient G0 in either cell type.** pRb never drops below ~1.66 pre-S; the sharp **CdP21 R-point**
   (p27 buffered on CyclinD1–CDK4/6, the Shh-gated commitment lever, [Architecture §4H](02_model_architecture.md)) is a
   near-binary go/no-go. This is a *change* from the pre-CDKI model, which produced a discrete birth-p27 transient G0
   (memory `v44-transient-g0-finding`); the reviewer's fast-KPC clearance + CdP21 buffer removed it
   (`cdki-baked-g0-g1-behavior`).
2. **MB's longer cycle is a longer G1 set by INK4 (p18=1.73), NOT a p27/G0 pause** — pre-S p27 is actually *lower* in
   MB (0.004) than GNP (0.007). So the model does **not** implement "MB competent for G0"; it implements "MB has a
   slower engine (longer G1)."
3. **Commitment runs on p27, and p21/p57 are inert** (p21a ≈ 0 even with the un-map wired in; iRc = 0). The un-map
   restored the *correct p21-PCNA machinery* but p21 is held at `kSyp21a≈2e-5` — a foundation with the switch not engaged.
4. **GNP response is a switch, not graded G1.** The graded μ_eff/G1-lengthening lever exists but is **latent**
   (`K_g1len=0.6` sits below the operating Cd; raising it to 1.8–2.5 recovers graded G1 — `cdki-baked-g0-g1-behavior`,
   `withdrawal-graded-g1-lengthening`).

---

## 3. The gap — model vs biology (five specific mismatches)

| # | Biology (§1) | Model (§2) | Consequence |
|---|---|---|---|
| G1 | GNP: graded G1/cycle lengthening | GNP: fixed G1, abrupt arrest (graded lever latent) | Can't show graded GNP slowdown without raising `K_g1len` |
| G2 | MB: competent for reversible **transient G0** | MB: G0-frac 0 %, longer G1 via INK4 | Model has **no MB G0 competence** at all |
| G3 | Transient G0 is **p21**-driven, stochastic, set by **birth-p21** heterogeneity | Commitment on **p27**; single deterministic lineage; p21 inert | Wrong CKI **and** no population layer → stochastic G0 is structurally impossible |
| G4 | MB stem dormancy = **p57/OLIG2** deep, lineage-programmed | No stem compartment, no OLIG2 node, p57 off | The therapy-resistant MB quiescence is entirely unmodeled |
| G5 | p27 is **dual-function** (pY88 → CDK4 *activator*) | p27 = pure inhibitor (CdP21 buffer) | Can't represent "p27-high but cycling"; palbo-resistance route missing |

**The load-bearing insight:** the model achieves its GNP/MB dichotomy and its Shh-gating with a **p27/CdP21 switch**,
but the *quiescence* biology it is trying to capture lives on **p21 (stochastic transient G0)** and **p57/OLIG2 (deep
stem dormancy)** — two CKIs the model currently carries **inert**. The recently-baked p21-PCNA un-map is the first
correct step (p21 is now wired where the biology puts it), but engaging it as a *functional bifurcation variable*
requires the changes below — and, per the earlier un-map work, engaging p21 at abundance in the current deterministic
architecture just G1-arrests the cell (the degron acts in S, the arrest is in G1; memory `p21-pcna-unmap-finding`).
The missing ingredient is exactly what the literature names: **birth-p21 heterogeneity across a population**, not a
single lineage.

---

## 4. Required changes — two modules + two reparameterizations

Staged, mapping the review's recommendations onto the model's actual state.

### Module A — the p21 mitotic-exit bifurcation (shallow, stochastic, reversible)
- **Variable:** engage **p21 (the now-correctly-wired `p21a`)** as the bifurcation determinant via the
  **stoichiometric CyclinD1-vs-p21 competition** setting an ultrasensitive Rb-E2F switch (Yang 2017 is the natural
  functional form; the model already has the two-step Rb + CdP21 pieces).
- **Noise source:** **birth-p21 inherited from mother G2** (replication-stress/p53 term), *not* a fixed value — this is
  the crux. It requires the **population layer** (`v45`-style ensemble or a birth-p21 distribution), because a
  deterministic single lineage cannot produce a stochastic CDK2-low *fraction*. Anchor the dose-response to Spencer
  2013 (CDK2-low fraction graded 0→100 % with p21).
- **Bistability:** the **p21 ⊣ CDK2 ⊣ p21** double-negative (SCF^Skp2^ + CRL4^Cdt2^) — the un-map already installed
  the CRL4^Cdt2^/PCNA half; add the reciprocal.
- **Irreversibility gate:** **APC/C-Cdh1 ↔ Emi1** hysteresis (Cappell 2016/2018) as the commitment point past which G0
  exit is blocked — the model has APC/Cdh1 (C1) but not the Emi1 dual-negative-feedback bistable form.
- *Status:* the model has the R-point machinery and now the p21-PCNA wiring; missing = birth-p21 heterogeneity +
  population layer + the Emi1 hysteresis. This is the route to **MB transient G0** (gap G2/G3).
- **Now protein-grounded (TMT, [§07.13](07_data_and_evidence.md)):** p21 protein is **undetectable in bulk MB**, so
  the mean is minimal — the bifurcation *must* be birth-p21 heterogeneity in a rare subpopulation, confirming Module A
  is a **population-layer** mechanism and that keeping p21 inert in the deterministic mean is correct.
- **Now TESTED (2026-08-02, `docs/stochastic_commitment_g0_findings.md`):** a reduced stochastic-commitment model
  (p21/p18/p27) shows the GNP-no-G0 / MB-G0 dichotomy **cannot** come from the CKI-vs-CyclinD1 stoichiometry — the
  measured MB changes (CyclinD1 ×4, CDK4/6 ×2, p18 ×3; p27 flat, p21 undetectable) make MB commit **harder** (less
  G0). Adversarially verified (72k-eval hunt + 5 wirings; max MB−GNP gap +0.0 pp). Reason: all three measured
  differences sit on the CDK4/6 arm the ×8 drive dominates; the CDK2 arm (p21/p27) that could bifurcate is flat/shared.
  **MB G0 works only from an MB-specific birth-p21 *tail* (variance, not mean; Spencer/Barr replication-stress) — so
  Module A's noise variable is birth-p21 with an MB-specific tail, NOT birth-p27, and NOT the drive-brake abundances.**

### Module B — p57/OLIG2 lineage-programmed deep dormancy (MB stem compartment)
- A **quiescence↔activation node** gated by **OLIG2** (Desai/Li 2025), *upstream* of the CDK machinery — NOT routed
  through the p21 branch. **p57 (and p27)** as **slow, high-exit-threshold** CKI variables; **p21 held low by OLIG2**
  (Ligon 2007) so this state is quiescent *despite* low p21.
- Slower kinetics + higher exit threshold than Module A (deep vs shallow). This is the therapy-resistant, SOX2⁺/OLIG2⁺
  compartment (gap G4). Entirely new; needs a stem/differentiation axis the current model doesn't have.

### Reparameterization 1 — p27 as dual-function
Make p27 a **stoichiometric CyclinD–CDK4/6 assembly factor** that is **inhibitory when unphosphorylated** and
**activating when Tyr88/89-phosphorylated** (Guiley 2019), with **Ser10-P/CRM1 nuclear export** and **Thr187-P/Skp2
degradation** as localization/turnover terms. Lets the model reproduce "p27-high but proliferating / palbo-resistant"
MB cells (gap G5). The current CdP21 buffer is the inhibitory half only.

### Reparameterization 2 — quiescence-depth state variable
A continuous **Rb-E2F activation threshold** that rises with time-in-G0 and with an ERK-attenuation/lysosomal state
(Yao lab: Kwon 2017, Fujimaki 2019; Ma 2024 p21-degron auto-maintenance), so exit probability falls *continuously* —
converting binary G0 into a graded depth and capturing deep-pool therapy resistance.

---

## 5. Literature evidence base (commissioned review, JP 2026-08-01)

**Routes into transient G0:** inherited mother-G2 p21 from replication stress (Spencer 2013; Barr 2017; Arora 2017;
Heldt 2018) · p21-independent low-CDK4/6-at-mitotic-exit (Yang/Meyer eLife 2020) · lineage quiescence of SOX2⁺/OLIG2⁺
MB stem cells (Vanner 2014; Zhang 2019; Desai/Li 2025) · therapy-induced (vismodegib: Ocasio 2019; radiation; antimitotics:
Vanner 2014) · APC/C-Cdh1↔Emi1 irreversibility (Cappell 2016/2018) · PRC2 represses Cdkn2a/1a/1c (chromatin gate; ties
to this model's EZH2 axis) · mTOR/TSC → cytoplasmic p27 (Bhatia 2009).

**p21 competent, p27 not (for the stochastic bifurcation):** Spencer's marker panel — only phospho-Rb and p21
discriminated CDK2-inc vs CDK2-low; p27 did not. Overton 2014 isolated p21 as *the* gene driving cycling/quiescent
heterogeneity via p21⊣CDK2⊣p21 double-negative feedback. **Gap:** no p27-knock-in-at-CDKN1A / matched-stoichiometry
substitution anywhere.

**Stoichiometry controls the fraction:** Spencer 2013 (p21 dose → CDK2-low fraction 0→100 %) · Overton 2014 (a mitogen
window for bimodality; p21-null → near-uniform cycling) · Ma 2024 (level-dependent switch: low-p21 reversible ↔
high-p21 ERK-attenuated auto-maintenance) · Yao lab (quiescence depth = tunable Rb-E2F threshold + size/growth memory).

**MB CKI compartment data (thin but convergent):** **p57/CDKN1C** marks quiescent SOX2⁺/Nestin⁺ MB stem cells (2025
preprint; *non-peer-reviewed — preliminary*) · **OLIG2 represses p21** (Ligon 2007; Treisman 2019) so the MB stem
compartment is quiescent with low p21 · older 14-tumor IHC: **p21 undetectable in all MB, p27 inverse to MIB1**
(single small series) · **p27 = progression not initiation** (Hatton 2013; Ayrault 2009 haploinsufficiency). **No
published scRNA-seq atlas gives a clean CDKN-by-compartment matrix** (Hovestadt/Vladoiu/Ocasio/Riemondy) — re-analysis
of those GEO datasets is the single highest-value empirical step.

**GNP/MB genetics anchoring this model's CKIs:** p18^Ink4c^ + p27^Kip1^ both restrain GNP proliferation, N-Myc
downregulates both (Zindy 2006) — consistent with this model's INK4-set MB cycle length. Cip/Kip triple-null MEFs still
show a recognizable G0 (Susaki 2012; Cerqueira 2014) → **G0 is not solely CKI-imposed** (a caution for Module B).

**Contested / single-source (flagged):** the p57-in-MB quantitation (one 2025 preprint); "p21 undetectable / p27∝1/MIB1"
(one ~14-tumor series); Overton's exact split percentages (unverified digits; qualitative claim robust).

---

## 6. Open experiments / benchmarks that would change the model

- **Re-analyze the public MB scRNA-seq atlases for a CDKN-gene × compartment matrix** (SOX2/OLIG2-stem vs ATOH1/MKI67-TA
  vs NEUROD1-diff) — replaces the analogy-based p27/p57 assignments with measured values. Highest value.
- A **CDK2-activity sensor in MB cells**: does a bimodal CDK2-low fraction appear under uniform mitogen, and does it
  survive **p21 knockout**? If it persists p21-null, Module A must admit a p27/p57 bifurcation (currently unsupported).
- **Acute p21-degron in MB cells** (never done) — tests the low-p21-reversible ↔ high-p21-auto-maintenance transition.
- **Dual p27 / p27-pY88 IHC in MB** — if pY88-p27 marks the cycling bulk, upweight the p27-activator term (Reparam 1).
- Whether **p57 produces a transient vs deep G0** (never tested) — sets Module B kinetics/threshold.

**Model-side benchmark:** a correct Module A should reproduce Spencer's graded CDK2-low fraction (0→100 % with birth-p21)
and the sister-fate 98 % concordance; a correct Module B should give a small (<5 %), deep, therapy-resistant SOX2⁺/OLIG2⁺
quiescent pool that the p21 route cannot re-activate.

---

*Cross-references: [Mechanisms Explored §8](04_mechanisms_explored.md) (transient G0 / G1-lengthening / withdrawal) ·
[Model Architecture §4F–H](02_model_architecture.md) (Skp2-p27 R-point, two-step Rb, dynamic CDKI, CdP21) ·
[Things Tried & Abandoned §3](05_things_tried_and_abandoned.md) (transient-G0 mechanisms that were structurally absent) ·
[Data & Evidence §7,§10](07_data_and_evidence.md) (cycle timing; senescence-vs-transient-G0 discriminator) ·
[Reviewer & Open Questions §3](10_reviewer_and_open_questions.md) (population layer). Memory: `cdki-baked-g0-g1-behavior`,
`mb-celltype-transient-g0-parameterization`, `p21-pcna-unmap-finding`, `v45-stochastic-commitment` (the population successor).*
