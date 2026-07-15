# Option B — CDK2 bypass of CDK4/6 inhibition, the Cip/Kip sequestration paradox, and how to model it

*Reading guide + mechanism dossier + modeling appendix for JP. Compiled 2026-07-15. Every mechanistic
claim is cited to a peer-reviewed primary source with PMID + DOI. The annotated reading list at the end is
the part meant to be printed and worked through; the modeling appendix is specific to our two-step-Rb Heldt
core (`src/build_model_v44_heldt.py`).*

Related internal docs: `docs/cdk46i_resistance_research.md` (the pRb⁺-residual resistance landscape),
`docs/EZH2_palbo_cycle_gating_tension.md` (why a naïve constitutive CDK2→Rb term breaks the oscillator).

---

## Plain-language summary (one paragraph)

Palbociclib arrests cells by shutting off cyclin D–CDK4/6, which normally puts the first ("mono")
phosphate on Rb and licenses the restriction point. Some cells keep dividing anyway because **cyclin
E/A–CDK2 can phosphorylate and hyperphosphorylate Rb (including Ser807/811) without CDK4/6 having gone
first** — a genuine bypass. Whether a given cell escapes hinges on the paradoxical double life of the
Cip/Kip proteins p21 (CDKN1A) and p27 (CDKN1B): at low, substoichiometric levels they are not pure
inhibitors — they **assemble, stabilize and nuclear-localize cyclin D–CDK4/6**, and that assembled
complex then acts as a **reservoir that sequesters p21/p27 away from CDK2**. JP's hypothesis is that a
cell carrying an unusually *high* total Cip/Kip load can fail to arrest under palbociclib because the
catalytically-dead-but-still-assembled cyclin D–CDK4/6 holds p21/p27 *off* CDK2, leaving CDK2 active
enough to phosphorylate Rb and run S-phase. The logic is a titration argument — the sequestration sink
buffers CDK2 from the inhibitor — and it has a hard limit: once total Cip/Kip exceeds what the sink can
hold, the excess spills onto CDK2 and the cell arrests after all. This document assembles the primary
evidence for each link in that chain, flags where the evidence complicates it, and gives a concrete,
oscillator-safe way to add "Option B" to our two-step Rb ODE.

---

## The mechanism

### 1. The Cip/Kip assembly / sequestration paradox

**p21 and p27 are not pure inhibitors; substoichiometrically they are cyclin D–CDK4/6 *assembly
factors*.** LaBaer et al. showed that p21, p27 and p57 each *promote* the association of CDK4 with
D-type cyclins, that the amount of assembled cyclin D–CDK4 rises with increasing inhibitor levels both
in vitro and in vivo, and proposed that the Cip/Kip family are adaptor proteins that "assemble and
program" kinase complexes rather than only retarding G1 (LaBaer et al., *Genes Dev* 1997; PMID 9106657;
DOI 10.1101/gad.11.7.847). Cheng et al. made the loss-of-function case: primary fibroblasts null for
*both* p21 and p27 **fail to assemble detectable cyclin D–CDK complexes, express cyclin D at reduced
levels, and cannot direct cyclin D to the nucleus** — all three defects reverse when Cip/Kip function is
restored (Cheng et al., *EMBO J* 1999; PMID 10075928; DOI 10.1093/emboj/18.6.1571). This is the
foundation of the "essential activators" view.

**The corollary — cyclin D–CDK4/6 is a Cip/Kip reservoir (sink) that titrates p21/p27 away from CDK2.**
Because assembled cyclin D–CDK4/6 stably holds Cip/Kip, it competes with cyclin E/A–CDK2 for the same
limited pool of inhibitor. Sherr & Roberts frame the whole G1 machine around this: cyclin D–CDK
complexes act as reservoirs that sequester p21/p27, so mitogen-driven cyclin D synthesis *indirectly*
activates CDK2 by pulling inhibitor off it — Cip/Kip are "positive and negative regulators of G1"
depending on which complex they occupy (Sherr & Roberts, *Genes Dev* 1999; PMID 10385618; DOI
10.1101/gad.13.12.1501). This review is the canonical statement of the titration/threshold logic.

**The threshold / titration logic, and where it breaks.** The stoichiometric switch is set by how much
Cip/Kip is bound in the cyclin D–CDK4/6 sink versus free to load onto CDK2, and by post-translational
state: p27 tyrosine-phosphorylation converts it between a CDK4-tolerated (even activating) partner and a
CDK4 inhibitor, and modulates the cyclin D–CDK4 vs CDK2 partitioning (James et al., *Mol Cell Biol*
2008; PMID 17908796; DOI 10.1128/MCB.02171-06; see also the mechanistic dissection in Ray et al., *Mol
Cell Biol* 2009, PMID 19075005 — **note: this paper carries a 2022 Expression of Concern over image
duplication, so treat its figures cautiously while the two-mode conclusion is corroborated
elsewhere**). The **titration reading that underlies JP's hypothesis**: under CDK4/6 inhibition, if the
cyclin D–CDK4/6 sink stays intact and large, raising *total* Cip/Kip can leave the *free* pool available
to CDK2 roughly flat — so more total p21/p27 does **not** translate into more CDK2 inhibition, and a
high-Cip/Kip cell can retain active CDK2 ("sequestration sink" protection). **Where it breaks:** the
sink has finite capacity (set by cyclin D–CDK4/6 abundance). Once total Cip/Kip exceeds sink capacity,
the excess spills onto cyclin E/A–CDK2, CDK2 is inhibited, and the cell arrests. So the relationship
between total Cip/Kip and escape is **non-monotonic**: protective while the sink absorbs the load,
arresting once it saturates. This is exactly a bistable/threshold problem and is why it wants a model.

### 2. Under palbociclib, is the complex disassembled or is the sequestration preserved?

This is the crux: if palbociclib redistributes p21/p27 from cyclin D–CDK4/6 onto CDK2, high Cip/Kip
*sensitizes* (more inhibitor floods CDK2); if the trimeric complex is *maintained* and keeps holding
Cip/Kip, high Cip/Kip can *protect*. The primary data point toward "**largely maintained**," with
important caveats.

- **Palbociclib preferentially binds monomeric CDK4/6, not the cyclin D–CDK4–Cip/Kip trimer, and
  phospho-p27–CDK4–cyclin D1 is drug-insensitive.** Guiley et al. solved the crystal structure showing
  that tyrosine-phosphorylated p27 **allosterically activates** CDK4–cyclin D1 to phosphorylate Rb (by
  remodeling the ATP site and releasing the activation segment), and that purified and endogenous
  phospho-p27–CDK4–cyclin D1 complexes are **insensitive to palbociclib**, which instead primarily
  targets monomeric CDK4/6 in breast-tumor cells (Guiley et al., *Science* 2019; PMID 31831640; DOI
  10.1126/science.aaw2106). Implication: a large fraction of Cip/Kip stays **sequestered** in an active,
  drug-resistant trimer under palbociclib rather than being dumped onto CDK2 — supporting the
  "sequestration preserved" side of JP's hypothesis, and adding a second mechanism (allosteric CDK4
  reactivation) by which high p27 blunts the drug.
- **But arrest does not strictly require the Cip/Kip redistribution model.** Pennycook & Barr deleted
  *both* p21 and p27 (RPE1 and MCF7) and found **palbociclib still initiates and maintains G1 arrest**,
  calling "into question the essentiality of the indirect model." They propose parallel direct
  (CDK4/6-catalytic) *and* indirect (Cip/Kip-redistribution) pathways: high Cip/Kip may add sensitivity
  but is not obligatory for arrest (Pennycook & Barr, *Open Biol* 2021; 11(11):210125; PMID 34784791;
  DOI 10.1098/rsob.210125). This is the key complicating paper and must be cited alongside the
  sequestration story — the "sink" is real but is not the *only* route to arrest.
- **Single-cell readout of the same crux:** cells that inherit high CDK2 activity at mitotic exit are
  the ones that ignore palbociclib. Asghar et al. showed by single-cell CDK2-activity imaging that
  palbociclib-resistant cells **exit mitosis directly into a high-CDK2 (CDK2^high) proliferative state,
  bypassing the restriction point and the requirement for CDK4/6**; tumors enriched for CDK2^high cells
  are resistant (Asghar et al., *Clin Cancer Res* 2017; 23(18):5561–5572; PMID 28606920; DOI
  10.1158/1078-0432.CCR-17-0369). Whether a cell lands in CDK2^high vs CDK2^low is set at birth by
  Cip/Kip (p21) dynamics — the link to Section 1.

**Bottom line for Section 2:** the trimeric cyclin D–CDK4/6–Cip/Kip complex is substantially preserved
under palbociclib (drug hits monomer; phospho-p27 trimer is insensitive), so much of the Cip/Kip pool
stays off CDK2 — the substrate for JP's sink hypothesis — **but** arrest can also proceed without any
Cip/Kip redistribution, so the sink is one mechanism among parallel ones, and its *net* effect (protect
vs sensitize) is a quantitative, threshold-dependent question, i.e. a modeling question.

### 3. CDK2-mediated Rb phosphorylation independent of CDK4/6 (the genuine bypass)

**Two-step Rb: mono-phosphorylation (cyclin D–CDK4/6) then hyperphosphorylation (cyclin E/A–CDK2).**
Narasimha et al. established that in early G1 Rb is *exclusively mono-phosphorylated* by cyclin
D–CDK4/6 — one phosphate on one of ~14 sites, and this mono-P Rb is still active (E2F-binding);
inactivating hyperphosphorylation is a separate, later, cyclin E/A–CDK2 event (Narasimha et al., *eLife*
2014; 3:e02872; PMID 24876129; DOI 10.7554/eLife.02872). Sanidas et al. extended this to a
"mono-phosphorylation code": the 14 mono-P isoforms bind distinct partners and drive distinct outputs,
confirming that mono-P (CDK4/6) and hyper-P (CDK2) are mechanistically distinct steps, not points on one
continuum (Sanidas et al., *Mol Cell* 2019; 73(5):985–1000; PMID 30711375; DOI
10.1016/j.molcel.2019.01.004). **This two-step biology is precisely our model's `Rb → Rbm → pRb`
structure.**

**CDK4/6 and CDK2 dock Rb differently — the bypass is structurally real.** Topacio et al. showed cyclin
D–CDK4/6 uses a specific docking interaction with an α-helix in Rb's **C-terminus that cyclins E, A and
B do *not* recognize**; mutating that helix blocks CDK4/6-driven Rb phosphorylation and arrests cells
(Topacio et al., *Mol Cell* 2019; 74(4):758–770; PMID 30982746; DOI 10.1016/j.molcel.2019.03.020). So
CDK2 reaches Rb by a different (RxL/cyclin-groove) route — it is not simply doing CDK4/6's job through
the same door. That independence is what makes a CDK4/6-inhibitor bypass mechanistically possible.

**Site specificity — can Ser807/811 be CDK2-phosphorylated de novo?** Yes. Ser807/811 (the epitope of
JP's phospho-Rb stain, Cell Signaling D20B12 #41359) is phosphorylated by **both** CDK4/6 and CDK2;
**only Ser780 is cyclin-D/CDK4/6-exclusive** (see `docs/cdk46i_resistance_research.md` for the
antibody/site mapping). So a residual Ser807/811⁺ signal under palbociclib is a **CDK2-bypass signature**
(Ser807/811⁺, Ser780⁻), not incomplete CDK4/6 inhibition (which would keep both⁺). The clean
discriminating experiment is co-staining Ser780.

**The important nuance — CDK2 can *maintain* hyper-P from S-phase onset but does not normally *initiate*
it in G1.** Chung et al. found by single-cell imaging that CDK4/6 activity is required to hyperphosphorylate
Rb *throughout G1*, whereas cyclin E/A–CDK2 can only maintain Rb hyperphosphorylation *starting at the
onset of S phase*; the restriction point is a transient-hysteresis (short-memory) property of CDK4/6, not
an irreversible latch (Chung et al., *Mol Cell* 2019; 76(4):562–573; PMID 31543423; DOI
10.1016/j.molcel.2019.08.020). Read together with Asghar (CDK2^high cells born past the R-point) and the
adaptive-resistance data below, the synthesis is: **the CDK2 bypass is not a property of normal G1 — it
switches on when CDK2 is abnormally elevated (high cyclin E, or low *free* Cip/Kip on CDK2), i.e. exactly
the regime JP's high-sink hypothesis produces.**

**Population/adaptive evidence that CDK2 recovers pRb under palbociclib.** Herrera-Abreu et al. showed
ER⁺ breast-cancer cells adapt to CDK4/6 inhibition within days via **non-canonical cyclin D1–CDK2 (and
cyclin E1–CDK2) that recover Rb phosphorylation and S-phase entry**, and that acquired resistance
selects CCNE1 amplification or RB1 loss (Herrera-Abreu et al., *Cancer Res* 2016; 76(8):2301–2313; PMID
27020857; DOI 10.1158/0008-5472.CAN-15-0728). Zikry et al. (Purvis lab) resolved this at single-cell
resolution as **"fractional resistance"**: a subpopulation completes the cycle *in the presence of the
drug*, showing premature accumulation of E2F1, Rb and CDK2 and **selective sensitivity to CDK2
inhibition** — direct support that the escaping fraction is CDK2-driven and that p21/CDK2 balance is what
sorts cells into escape vs arrest (Zikry et al., *PNAS* 2024; 121(8):e2309261121; PMID 38324568; DOI
10.1073/pnas.2309261121).

---

## Annotated reading list

The 15 papers that matter, grouped by theme. Start with the one starred (★).

### A. Cip/Kip assembly & sequestration (the paradox)

1. **LaBaer J, Garrett MD, Stevenson LF, et al. "New functional activities for the p21 family of CDK
   inhibitors." *Genes Dev* 1997;11(7):847–862.** PMID 9106657 · DOI 10.1101/gad.11.7.847.
   *The original demonstration that p21/p27/p57 promote cyclin D–CDK4 assembly and that assembled complex
   scales with inhibitor level. Read for the founding "adaptor, not just inhibitor" claim.*

2. **Cheng M, Olivier P, Diehl JA, et al. "The p21^Cip1 and p27^Kip1 CDK 'inhibitors' are essential
   activators of cyclin D-dependent kinases in murine fibroblasts." *EMBO J* 1999;18(6):1571–1583.** PMID
   10075928 · DOI 10.1093/emboj/18.6.1571.
   *Genetic proof: p21/p27 double-null fibroblasts can't assemble, stabilize, or nuclear-localize cyclin
   D–CDK. The loss-of-function companion to LaBaer — this is the assembly half of the paradox.*

3. **★ Sherr CJ, Roberts JM. "CDK inhibitors: positive and negative regulators of G1-phase progression."
   *Genes Dev* 1999;13(12):1501–1512.** PMID 10385618 · DOI 10.1101/gad.13.12.1501.
   *The canonical review that states the reservoir/sequestration logic outright — cyclin D–CDK sinks
   titrate Cip/Kip off CDK2. Best single entry point; read this first to load the whole conceptual
   frame before the primary papers.*

4. **James MK, Ray A, Leznova D, Blain SW. "Differential modification of p27^Kip1 controls its cyclin
   D–Cdk4 inhibitory activity." *Mol Cell Biol* 2008;28(1):498–510.** PMID 17908796 · DOI
   10.1128/MCB.02171-06.
   *Stoichiometry/threshold detail: p27 tyrosine-phosphorylation switches it between inhibitor and
   tolerated/activating partner and re-partitions it between CDK4 and CDK2 pools. This is the molecular
   knob behind "sink vs spillover." (Companion Ray et al. Mol Cell Biol 2009, PMID 19075005, dissects two
   inhibitory modes but carries a 2022 Expression of Concern — cite cautiously.)*

### B. Rb phospho-site biology (why the two-step Rb is real, and CDK2 can bypass)

5. **Narasimha AM, Kaulich M, Shapiro GS, et al. "Cyclin D activates the Rb tumor suppressor by
   mono-phosphorylation." *eLife* 2014;3:e02872.** PMID 24876129 · DOI 10.7554/eLife.02872.
   *Establishes mono-P (cyclin D–CDK4/6, active Rb) vs hyper-P (cyclin E/A–CDK2, inactive Rb) as two
   distinct steps. This IS our `Rb→Rbm→pRb` structure; read to justify the two-step topology.*

6. **Sanidas I, Morris R, Fella KA, et al. "A code of mono-phosphorylation modulates the function of
   RB." *Mol Cell* 2019;73(5):985–1000.** PMID 30711375 · DOI 10.1016/j.molcel.2019.01.004.
   *14 mono-P isoforms, distinct partners/outputs; Ser811 is one of them. Confirms mono- vs hyper-P are
   mechanistically separate and site-specific — grounds the model's discretization of Rb states.*

7. **Topacio BR, Zatulovskiy E, Cristea S, et al. "Cyclin D–Cdk4,6 drives cell-cycle progression via the
   retinoblastoma protein's C-terminal helix." *Mol Cell* 2019;74(4):758–770.** PMID 30982746 · DOI
   10.1016/j.molcel.2019.03.020.
   *CDK4/6 docks an Rb C-terminal helix that cyclins E/A/B do NOT use. The structural basis for why CDK2
   is a genuinely independent Rb kinase — the load-bearing paper for "the bypass is real, not redundant."*

8. **Chung M, Liu C, Yang HW, et al. "Transient hysteresis in CDK4/6 activity underlies passage of the
   restriction point in G1." *Mol Cell* 2019;76(4):562–573.** PMID 31543423 · DOI
   10.1016/j.molcel.2019.08.020.
   *The crucial nuance: CDK2 maintains Rb hyper-P only from S-phase onset; CDK4/6 is needed through G1.
   Tells you the bypass is conditional on elevated CDK2, and that the R-point is short-memory hysteresis,
   not a permanent latch. Directly shapes how to gate the CDK2 term in the model.*

### C. CDK4/6-inhibitor bypass & the Cip/Kip–palbociclib crux

9. **★(Section 2) Guiley KZ, Stevenson JW, Lou K, et al. "p27 allosterically activates cyclin-dependent
   kinase 4 and antagonizes palbociclib inhibition." *Science* 2019;366(6471):eaaw2106.** PMID 31831640 ·
   DOI 10.1126/science.aaw2106.
   *Structure + biochemistry: phospho-p27–CDK4–cyclin D1 is catalytically active AND palbociclib-insensitive;
   the drug hits monomeric CDK4/6. The strongest evidence that high p27 keeps Cip/Kip sequestered (and can
   even reactivate CDK4) under palbociclib. Read immediately after Sherr & Roberts.*

10. **Pennycook BR, Barr AR. "Palbociclib-mediated cell cycle arrest can occur in the absence of the CDK
    inhibitors p21 and p27." *Open Biol* 2021;11(11):210125.** PMID 34784791 · DOI 10.1098/rsob.210125.
    *The essential counterweight: p21/p27 double-KO cells still arrest under palbociclib → the indirect
    (redistribution) model is not obligatory; direct + indirect pathways run in parallel. Keeps the sink
    hypothesis honest.*

11. **Herrera-Abreu MT, Palafox M, Asghar U, et al. "Early adaptation and acquired resistance to CDK4/6
    inhibition in ER⁺ breast cancer." *Cancer Res* 2016;76(8):2301–2313.** PMID 27020857 · DOI
    10.1158/0008-5472.CAN-15-0728.
    *Non-canonical cyclin D1–CDK2 / cyclin E1–CDK2 recover pRb and S-phase entry under palbociclib
    (adaptive), then CCNE1-amp / RB1-loss on acquired resistance. Population-level proof of CDK2 bypass.*

12. **Asghar US, Barr AR, Cutts R, et al. "Single-cell dynamics determines response to CDK4/6 inhibition
    in triple-negative breast cancer." *Clin Cancer Res* 2017;23(18):5561–5572.** PMID 28606920 · DOI
    10.1158/1078-0432.CCR-17-0369.
    *Single-cell CDK2-activity imaging: resistant cells exit mitosis into a CDK2^high state that bypasses
    the R-point and CDK4/6 requirement. The single-cell face of Option B — birth CDK2/Cip-Kip sets fate.*

13. **Zikry TM, Wolff SC, Ranek JS, et al. "Cell cycle plasticity underlies fractional resistance to
    palbociclib in ER⁺/HER2⁻ breast tumor cells." *PNAS* 2024;121(8):e2309261121.** PMID 38324568 · DOI
    10.1073/pnas.2309261121.
    *Modern single-cell "fractional resistance": a subpopulation finishes the cycle on-drug with premature
    E2F1/Rb/CDK2 and is selectively CDK2-inhibitor-sensitive. The best contemporary demonstration that the
    escaping fraction is CDK2-driven and heterogeneity-based — the phenomenon Option B must reproduce.*

### D. Cell-cycle ODE / threshold modeling (how to build it)

14. **★(modeling) Heldt FS, Barr AR, Cooper S, Bakal C, Novák B. "A comprehensive model for the
    proliferation–quiescence decision in response to endogenous DNA damage in human cells." *PNAS*
    2018;115(10):2532–2537.** PMID 29463760 · DOI 10.1073/pnas.1715345115.
    *Our core (BioModels BIOMD0000000700). Two coupled bistable switches — the R-point (Rb) and the G1/S
    transition (CDK2–p27 double-negative feedback) — with lumped p21. The base to extend for Option B.*

15. **Barr AR, Heldt FS, Zhang T, Bakal C, Novák B. "A dynamical framework for the all-or-none G1/S
    transition." *Cell Syst* 2016;2(1):27–37.** PMID 27136687 · DOI 10.1016/j.cels.2016.01.001.
    *The CDK2:cyclin ⊣ p27 double-negative-feedback bistable switch in isolation — the exact module JP's
    "free-CDK2 vs sequestered-Cip/Kip" balance feeds into. Read for the switch that Option B tips.*

**Also worth having open (supporting):**

- **Barr AR, Cooper S, Heldt FS, et al. "DNA damage during S-phase mediates the proliferation–quiescence
  decision in the subsequent G1 via p21 expression." *Nat Commun* 2017;8:14728.** PMID 28317845 · DOI
  10.1038/ncomms14728. *The CDK2^low vs CDK2^inc bifurcation at mitotic exit is set by inherited p21 — the
  mechanistic basis for modeling escape as a birth-condition (Cip/Kip) heterogeneity, not a per-cell coast.*
- **Yao G, Lee TJ, Mori S, Nevins JR, You L. "A bistable Rb–E2F switch underlies the restriction point."
  *Nat Cell Biol* 2008;10(4):476–482.** PMID 18364697 · DOI 10.1038/ncb1711. *The canonical bistable
  Rb–E2F threshold model; the minimal switch topology our R-point inherits.*
- **Gérard C, Goldbeter A. "Temporal self-organization of the cyclin/Cdk network driving the mammalian
  cell cycle." *PNAS* 2009;106(51):21643–21648.** PMID 20007375 · DOI 10.1073/pnas.0903827106. *The
  relaxation-oscillator engine of our v42/v43 lineage; useful for the oscillator-reset intuition in the
  appendix.*

---

## Modeling appendix — adding Option B to the two-step Rb ODE

### The current structure and why it fully arrests under palbociclib

Our Heldt-core two-step Rb (`docs/v44_MODEL_DESCRIPTION.md` §A) is:

```
Rb  --(mono; cyclin D–CDK4/6, rate kPhRbCd·Cd)-->  Rbm  --(hyper; cyclin E/A–CDK2, kPhRbCe·Ce + kPhRbCa·Ca)-->  pRb
```

with p21/p27 **lumped** as one inhibitor `P21` that binds/inhibits CDK2 (`CeP21`, `CaP21`), and the
R-point implemented as a sharp Skp2–p27 bistable switch. **Palbociclib = `kPhRbCd → 0`.** In this
topology CDK2 acts only on `Rbm`, and `Rbm` is made only by CDK4/6 — so zeroing `kPhRbCd` starves CDK2
of substrate and the whole population arrests. That is faithful to Chung 2019 (CDK4/6 needed through G1)
but it makes **Option B structurally impossible**: there is no path for CDK2 to touch Rb without CDK4/6
going first.

**The failed naïve fix (already tested — see `docs/cdk46i_resistance_research.md`):** adding a
constitutive `+ k·(Ce+Ca)·Rb` term raises `pRb` (cell reads pRb⁺) but **stalls division** — a constant
Rb-phosphorylation term fights the oscillator, Rb never resets to hypophospho at mitotic exit, and the
cell sits pRb⁺-but-arrested. An E2f-*gated* CDK2 term can't bootstrap either, because under palbociclib
arrest E2f is already off. So the escape needs (a) a CDK4/6-independent CDK2 route **that is cycle-gated,
not constitutive**, and (b) a source of residual CDK2 activity under palbociclib to seed it — which is
exactly what the Cip/Kip sink supplies.

### (a) Minimal structural change: a small CDK4/6-independent CDK2 route on the *mono* step + Cip/Kip partitioning

**Two coupled additions.**

**1. Split the lumped inhibitor into two pools with a maintained cyclin D–CDK4/6 sink.** Introduce (or
reuse the existing `Cd`) a cyclin D–CDK4/6 species `D46` that binds free inhibitor to form a sequestered
trimer `D46_I`, and let free inhibitor `I_free` be what CDK2 sees:

```
association:   D46 + I_free  -->  D46_I     (rate k_seq)
dissociation:  D46_I  -->  D46 + I_free     (rate k_off)
palbociclib:   kPhRbCd -> 0   BUT   k_seq, k_off UNCHANGED   (trimer maintained; Guiley 2019, Pennycook & Barr 2021)
CDK2 inhibition uses I_free (not total I):  free CDK2 activity  ∝ (Ce+Ca) − I_free  (via CeP21/CaP21 as now, but fed by I_free)
```

The key modeling commitment, justified by Guiley 2019, is that **palbociclib does not release the sink**
(`k_off` unchanged when `kPhRbCd→0`). Then the sink capacity `≈ D46·(k_seq/k_off)` buffers `I_free`:
raising total inhibitor while `D46` is abundant leaves `I_free` (hence CDK2 inhibition) nearly flat —
the **sequestration-sink / high-Cip-Kip-protects** regime — until total inhibitor exceeds sink capacity,
at which point `I_free` rises steeply and CDK2 is shut off (arrest). That non-monotonic
`I_free(total inhibitor)` is Option B's core, and it emerges from mass action, not a hand-drawn curve.

**2. Give CDK2 a small CDK4/6-independent mono-phosphorylation of Rb, gated by CDK2 activity (not a
constant):**

```
Rb  --(kPhRbCd·Cd  +  kBypass·(Ce+Ca)·f_gate)-->  Rbm
```

where `kBypass ≪ kPhRbCd` and `f_gate` is a self-consistent cycle gate (see (c)). This lets **residual
free CDK2** (large only in high-sink / high-cyclin-E cells) slowly prime Rb → `Rbm` even with
`kPhRbCd=0`; the existing hyper step then finishes the job. Because the seed is `(Ce+Ca)` — which is high
in CDK2^high cells and low in genuinely arrested cells — the bypass is **conditional on elevated CDK2**,
matching Chung 2019 / Asghar 2017. Keep `kBypass` small enough that a normal low-CDK2 arrested cell never
crosses the Rb threshold; only cells whose sink keeps `I_free` low enough to sustain CDK2 escape.

This reproduces JP's hypothesis end to end: high total Cip/Kip → large maintained `D46_I` sink → low
`I_free` on CDK2 → residual CDK2 activity → `kBypass` route primes Rb → hyper-P → S-phase → escape,
**despite** `kPhRbCd=0`.

### (b) How others partition Cip/Kip and model CDK2 bypass

- **Heldt 2018 (PMID 29463760)** already has the two-switch skeleton (R-point Rb switch + CDK2–p27 G1/S
  switch) and a lumped p21; Option B is the surgical addition of a *second* p21 pool (the D46 sink) and a
  small mono-step CDK2 term. Nothing else in the topology needs to move.
- **Barr et al. 2016 (PMID 27136687)** is the CDK2:cyclin ⊣ p27 double-negative bistable switch — the
  precise module that `I_free` feeds. Option B is literally "let the CDK4/6 sink set the p27 input to the
  Barr switch."
- **Barr et al. 2017 (PMID 28317845)** shows the CDK2^low/CDK2^inc fate is set at *birth* by inherited
  p21. The clean way to get a *fractional* escape (not all-or-none) is to make the escape a **birth-condition
  heterogeneity**: sample inherited total Cip/Kip (and/or cyclin E, `D46`) across daughters; the tail with
  large sink + high CDK2 escapes. This matches Zikry 2024's fractional resistance and avoids needing a
  graded partial drug block.
- **Yao 2008 (PMID 18364697)** and **Gérard–Goldbeter 2009 (PMID 20007375)** are the reference bistable-
  switch and oscillator engines respectively — useful for sanity-checking that the added terms preserve
  bistability (Yao) and clean oscillation/reset (GG).

### (c) The oscillator-reset problem and how to avoid it

**The trap:** any Rb-phosphorylation term that does not collapse at mitotic exit will re-phosphorylate the
freshly reset (hypophospho) Rb of the newborn before the R-point context is re-established → the lineage
reads permanently pRb⁺ and **the DNA-replication licensing / division cannot reset**, so the cell stalls.
Our earlier constitutive-term test hit exactly this.

**The fix — make the bypass cycle-gated by a variable that itself resets at division:**

- Gate `kBypass` on **CDK2 activity itself** (`(Ce+Ca)` or free-CDK2), which is driven low at mitotic exit
  because APC/C degrades cyclin A and the division event resets `Ce`, `Ca`, `E2f`, `Skp2` (model §D). A
  newborn therefore starts with the bypass **off** and only re-engages it if its inherited sink/cyclin-E
  put free CDK2 above threshold — self-consistent, no external clock.
- Do **not** put the bypass on the hyper step (`Rbm→pRb`) as a standalone constant; put it on the **mono
  step** (`Rb→Rbm`) so it is upstream of, and rate-limited by, the same E2f/CDK2 feedback that resets.
- Keep `kBypass` in the "small" regime confirmed by a one-parameter sweep: it must be **below** the value
  at which a WT, mitogen-deprived (or fully palbociclib-arrested, low-CDK2) cell crosses the Rb threshold,
  and **above** the value at which a high-sink/high-cyclin-E cell can.

**Tests (acceptance criteria):**

1. **Baseline unchanged:** with `kPhRbCd` normal, period, phase fractions, and G0 fraction match the
   current 25–28-target calibration (no drift from adding the sink + bypass).
2. **Reset integrity:** over ≥5 simulated divisions, `pRb` returns to ~0 (hypophospho) at each mitotic
   exit and `Dna` resets 1→0 — i.e. no monotonic pRb accumulation, no stall (the direct test the naïve
   term failed).
3. **Graded/fractional escape under palbociclib:** with `kPhRbCd=0` and a *distribution* of inherited
   total Cip/Kip (or `D46`, cyclin E) across the population, a **subpopulation stays pRb⁺ and divides**
   while the rest arrest — target the culture number (residual pRb⁺ cycling ≈ 16–18%,
   `docs/cdk46i_resistance_research.md`). A single homogeneous cell must be all-or-none; the fraction must
   come from birth-condition spread (Barr 2017 / Zikry 2024), not a uniform partial block.
4. **Non-monotonic Cip/Kip dependence:** sweep total Cip/Kip at fixed `D46`; escape fraction should
   **rise then fall** (protective while the sink absorbs, arresting once it saturates) — the falsifiable
   signature of the sequestration-sink hypothesis.
5. **CDK2-inhibitor collapse:** adding a CDK2 inhibition term (raise effective `I_free` on CDK2, or lower
   `Ce+Ca`) should **abolish** the escaping fraction — matching Zikry 2024 (fractional-resistant cells are
   CDK2i-sensitive) and Asghar 2017.
6. **Site signature (qualitative):** the escaping cells should be Ser807/811⁺ but would be Ser780⁻ in a
   two-site model — consistent with the CDK2-bypass discriminator (a note for future data, not currently
   a model species).

**Scope note.** This is a real R-point/oscillator modification, not a cosmetic term — JP has previously
judged the elaborate pRb⁺-dividing clone "not worth building" for the EZH2 story, and the honest fallback
(documented in `cdk46i_resistance_research.md`) is to *note* the CDK2 bypass as a literature-grounded
resistance phenomenon rather than dynamically model it. If Option B is pursued, the sink-partitioning +
small cycle-gated mono-step term above is the minimal structurally-honest way in, and the six tests are
the guardrails.

---

*Compiled from the primary literature via PubMed / journal sources; all PMIDs and DOIs verified against
PubMed records on 2026-07-15.*
