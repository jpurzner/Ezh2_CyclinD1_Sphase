# CDK6 vs CDK4 asymmetry in SHH-MB — why, and what it means (2026-08-07)

Data (JP CSV): CDK6 1508→26003 (**17×**, vismo 0.41×); CDK4 10334→21510 (**2×**, vismo 0.53×). CDK6:CDK4 ratio flips
0.15 (GNP, CDK4-dominated) → 0.29 (P28-Ptch) → **1.2 (MB, CDK6-co-dominant)**. Literature sweep (workflow wvjnh2rfw,
4 angles + synthesis) below.

## Why CDK6 specifically (not CDK4)
- **CDK6 is a DIRECT GLI2 target; CDK4 is not** (Zhang, *JCI* 2018, PMID 29202464). GLI2 ChIP selectively occupies a
  *Cdk6* enhancer ("site 4," a GLI element borrowed from Hh limb development) — necessary + sufficient for Hh-
  responsiveness; **no comparable occupancy at *Cdk4***. Their cohort: ~167× CDK6, ~13.9× *less* differential for CDK4
  — the same "CDK4 high-but-flat / CDK6 low-but-inducible" split JP sees.
- **★ Directly ties CDK6 to the EZH2 governor:** site 4 **gains H3K27ac and LOSES H3K27me3** in MB. So CDK6 induction
  IS de-repression of an EZH2/H3K27me3-controlled Hh enhancer — not a parallel add-on. Warrants `with_cdk6_gli`.
- Paralog logic: CDK6 is a **signal/lineage-gated** gene (inducible TFs); CDK4 is **constitutive** (E2F/Sp1/E2-box) →
  rises only ~2×, in proportion to cycling. A morphogen TF (GLI2) naturally lands on CDK6.
- Coherent feed-forward reinforces the same node: GLI2→CDK6 **and** GLI2→MYCN→E2F5→CDK6 (not CDK4).
- **Epistasis:** CDK4/6i on CDK6-null tumors gives no further growth reduction → abundant CDK4 is *not* the driver.

## What the asymmetry MEANS (functional — CDK6 is not a redundant Rb kinase)
A tumor doesn't need 17× more Rb-kinase (the paralogs are interchangeable on Rb). The CDK6 surge buys CDK6's
**CDK4-absent** functions:
1. **G0-escape rheostat.** CDK6 *level* sets how fast a cell leaves quiescence (Laurenti/Trumpp, *Cell Stem Cell*
   2015, PMID 25704240): deep-dormant cells ~no CDK6; primed cells pre-load high CDK6 for fast re-entry. **More CDK6 =
   lower G0→cycle activation barrier.** ← this is the reversible, mitogen-gated G0-exit dial the transient-G0 work has
   been missing.
2. **Kinase-independent differentiation block / stemness keeper.** CDK6 (not CDK4) acts as a promoter-bound TF cofactor
   blocking differentiation (RUNX1/C-EBPα); kinase-dead CDK6 still deepens stemness (Placke *Blood* 2014 PMID 24764564;
   Kollmann *Cancer Cell* 2013).
3. **INK4-buffer titration.** CDK6 is preferentially held in INK4·CDK6 (OFF) complexes; **forcing CDK6 up 17×
   stoichiometrically overwhelms the p18 brake** — the brake fails by substrate titration, not affinity loss. (Ties to
   JP's p18-restricts point: MB's high CDK6 is how it escapes the p18 brake.)

**So the ratio flip is a lineage-state switch, not a dosage bump:** GNP = CDK4-dominant, INK4-brakeable,
differentiation-competent; MB = CDK6-co-dominant → primed to exit G0 fast, differentiation-blocked, p18 titrated out —
a differentiation-blocked, Hh-addicted engine.

## Implications for the model / transient G0
- Reinforces & mechanistically grounds `with_cdk6_gli` (Gli/EZH2 multiplier on the CyclinD-CDK4/6→Rb Vmax; CDK4 flat),
  and hardens "vismo needs the axis" (CDK6 is the Hh transcriptional readout → collapses harder under vismo).
- **★ Reframes transient G0:** CDK6 = the G0-escape rheostat →
  - **Atoh1⁺ / cycling** compartment = CDK6-high (fast G0 exit + differentiation block → stay Atoh1⁺, keep dividing;
    transient G0 = brief, fast re-entry).
  - **Sox2⁺ quiescent reserve** = CDK6-low / p18-braked / deep-dormant.
  - Suggests giving Gli/EZH2-driven CDK6 a **second, kinase-independent action beyond the Rb-Vmax multiplier: lower the
    G0→cycle barrier + oppose differentiation.** That is the reversible mitogen-gated toggle the deterministic core
    lacks — a mechanistically-grounded route to the Sox2(deep-G0) ↔ Atoh1(primed) reversibility as a population layer.
- p18 (confirmed cycle-restricting in the model: KO cycles at SHH 0.2 where WT arrests) is the Atoh1⁺-compartment G0
  brake; MB overcomes it by CDK6 titration.

## Testable predictions
- **Drug ordering (actionable):** CDK6-co-dominant SHH-MB should be relatively **ribociclib-resistant** (ribociclib
  under-covers CDK6, IC50 ~39 vs ~10 nM) but **palbociclib/abemaciclib-sensitive**. Encodable in the CDK6 arm.
- **Vismo de-represses differentiation** (releasing CDK6's block), not just Rb-kinase → differentiation-marker gain
  under vismo beyond what pure Rb-kinase reduction predicts.
- **Compartment:** CDK6 high in Atoh1⁺ cycling, low in Sox2⁺ quiescent; enforced CDK6 in Sox2⁺ shortens G0; CDK6 KD →
  Sox2⁺/differentiated (Li *J Neurooncol* 2012 PMID 22528790).

Sources: PMID 29202464 (GLI2→Cdk6 site 4), 25704240 (CDK6 quiescence-exit timer), 24764564 / Kollmann 2013 (kinase-
independent CDK6 TF/differentiation block), 16479172 (p18^INK4c MB tumor suppressor), Genes Dev 2000 PMC317144
(INK4·CDK6 sequestration). Caveat: quiescence-timer/kinase-independent data are strongest in HSC/leukemia; a
cerebellar GNP/MB CDK6 ChIP + Sox2/Atoh1-resolved CDK6 protein dataset would be direct confirmation.
