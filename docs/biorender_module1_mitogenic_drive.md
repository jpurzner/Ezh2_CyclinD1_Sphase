# BioRender build — Module 1: Mitogenic drive (Hedgehog → Gli, MYCN → *Ccnd1*)

Source of truth: `HH_MYCN_BLOCK` in `src/build_model_v44_heldt.py`. Every element below is one
species or one flux in that block. Rates are per minute. Bold words are the on-figure labels.

## Context inputs (constants, set per condition)

| symbol | meaning | GNP | Ptch1+/− | MB |
|---|---|---|---|---|
| SHH | ligand | 0.5 | 0.5 | 0.5 |
| HHi | vismodegib (0 off, 1 full Smo block) | 0 | 0 | 0 or 1 |
| c_Ptch1 | functional Ptch1 fraction (gates Ptch1's repression of Smo, not its transcription) | 1.0 | 0.5 | 0.1 |
| A_MYCN | MYCN amplification | 1.0 | 1.0 | 2.8 |

Timescale: every species here turns over in ~1 min except the Gli1 locus memory (~14 h), so the
whole module is quasi-static. It collapses to a single algebraic input to the rest of the model:
the *Ccnd1* transcription drive D (element 10).

## Elements, with equations

**1. SHH** (ligand icon, extracellular). Constant input. Binds Ptch1:
```
Ptch1_free + SHH → SHH:Ptch1      rate k_bind·SHH·Ptch1_free      (k_bind = 5.0)
SHH:Ptch1 → Ptch1_free            rate k_rel·SHH:Ptch1            (k_rel = 0.1)
SHH:Ptch1 → ∅                     rate k_deg·SHH:Ptch1            (k_deg = 0.8)
```
Ligand-bound Ptch1 no longer represses Smo and is degraded. Draw the complex as Ptch1 with SHH docked, fading (internalised).

**2. Ptch1** (12-pass receptor, membrane). A Gli target, which closes the negative feedback:
```
d[Ptch1 mRNA]/dt = k_basal + k_Gli · Gli_act² / (K_P² + Gli_act²) − k_mdeg·[Ptch1 mRNA]
                    (k_basal 0.149, k_Gli 1.87, K_P 0.398, k_mdeg 0.8)
d[Ptch1_free]/dt  = k_tl·[Ptch1 mRNA] − k_pdeg·Ptch1_free − k_bind·SHH·Ptch1_free + k_rel·SHH:Ptch1
                    (k_tl 1.0, k_pdeg 0.5)
```
On-figure: **Ptch1** with a curved return arrow from Gli labelled **Gli → Ptch1 negative feedback**.
Add tag **lost in MB: Ptch1 made but non-functional (c_Ptch1 = 0.1)**. This is why MB has both
constitutive Gli and high Ptch1 mRNA (the SHH-MB marker).

**3. Smo** (7-TM receptor, membrane). Repressed by functional free Ptch1, blocked by vismodegib:
```
d[Smo*]/dt = k_act · (1 − HHi) / (1 + Ptch1_free·c_Ptch1 / K_PS) − k_inact·Smo*
             (k_act 1.57, K_PS 0.1025, k_inact 1.2)
```
On-figure: **Ptch1 ⊣ Smo** (flat bar), **vismodegib ⊣ Smo** (red flat bar, pill icon). Tag on the
Ptch1 bar: **brake scaled by functional Ptch1**.

**4. Gli switch** (one TF icon drawn as a two-state toggle, repressor ⇄ activator; total Gli conserved at 0.6):
```
S(Smo) = Smo*² / (K_S² + Smo*²)                                   (K_S 0.888)
Gli_rep → Gli_act   rate k_on · S(Smo) · Gli_rep                  (k_on 3.0)
Gli_act → Gli_rep   rate k_off · (1 − S(Smo)) · Gli_act           (k_off 1.5)
```
On-figure: **Smo → Gli** with the toggle labelled **Gli-R ⇄ Gli-A (Smo-driven, Hill 2)**.

**5. Gli1** (mRNA → protein; the pathway readout). Activated by Gli-A, repressed by Gli-R:
```
H_act = Gli_act² / (K_a² + Gli_act²)         (K_a 1.18)
H_rep = Gli_rep² / (K_r² + Gli_rep²)         (K_r 0.4)
d[Gli1 mRNA]/dt = Vmax · H_act · (1 − H_rep) · [c_Ptch1 + (1 − c_Ptch1)·0.873]
                  + k_auto·Gli1_epi − k_mdeg·[Gli1 mRNA]        (Vmax 1.27, k_auto 0.012, k_mdeg 0.8)
d[Gli1]/dt      = k_tl·[Gli1 mRNA] − k_deg·Gli1                  (k_tl 1.2, k_deg 0.8)
```
On-figure: **Gli1** as a small node off Gli with tag **pathway readout (MB/GNP ≈ 7×)**. Gli1 also
feeds forward with Gli-A into MYCN and *Ccnd1* (G = Gli_act + Gli1 below).

**6. Gli1 locus memory** (optional; small chromatin icon beside Gli1, MB only). A slow 0–1
chromatin state charged by Smo in cells with broken Ptch1 feedback; it gives the Gli1 residual
that persists ~1 day after vismodegib:
```
d[Gli1_epi]/dt = k_on·(1 − c_Ptch1) · Smo*⁴/(K⁴ + Smo*⁴) · (1 − Gli1_epi) − k_off·Gli1_epi
                 (k_on 0.02, K 0.4, k_off 0.0008 → discharge t½ ≈ 14 h)
```
On-figure if included: **Gli1 locus memory (MB only, t½ 14 h)** with dashed Smo → memory → Gli1.
Silent in GNP because the (1 − c_Ptch1) factor is 0.

**7. MYCN** (TF icon). Gli target plus a mitogen-independent basal term scaled by amplification:
```
G = Gli_act + Gli1
d[MYCN]/dt = k_bas·A_MYCN + k_Gli · G / (K_M + G) − k_deg·MYCN     (k_bas 0.3, k_Gli 0.102, K_M 0.5, k_deg 1.0)
```
On-figure: **Gli → MYCN**, and **MYCN** carrying tag **basal floor × amplification (MB 2.8×)**.
This floor is what makes vismodegib only a partial mitogen withdrawal in MB.

**8. Cell-type badges.** Two mini-cells: **GNP: c_Ptch1 = 1, MYCN ×1** and **MB: c_Ptch1 = 0.1,
MYCN ×2.8**. Ptch1+/− (0.5) can be a third if space allows.

**9. Vismodegib** (red pill). Enters only through element 3 as the factor (1 − HHi).

**10. *Ccnd1* transcription drive D** (the hand-off to Module 2; draw as the arrowhead landing on
the *Ccnd1* promoter). Three additive inputs:
```
D = k_0                                                          basal      (k_0 0.483)
  + k_G · G² / (K_G² + G²) · K_R² / (K_R² + Gli_rep²)             Gli drive  (k_G 59.2, K_G 0.457, K_R 0.3)
  + k_N · MYCN^n / (K_N^n + MYCN^n)                                MYCN drive (k_N 21.7, K_N 1.66, n 3.66)

d[Ccnd1 mRNA]/dt = D · R − k_mdeg·[Ccnd1 mRNA]                    (k_mdeg 0.8)
```
R is the chromatin repression factor from Module 2 (0 < R ≤ 1); Module 1 owns only D.
On-figure: two arrows into the promoter, **Gli → *Ccnd1*** (thick) and **MYCN → *Ccnd1***, with
one tag: **Gli supplies ~86 % of MB Cyclin D1; after vismodegib the residual is ~78 % MYCN, ~22 % basal**.

**11. Cyclin D1 protein** (the module's output node, shared with Module 2):
```
d[Cd]/dt = k_tl·[Ccnd1 mRNA] − k_deg·Cd          (k_tl 0.80, k_deg 1.0 → t½ < 1 min)
```
Protein tracks mRNA instantaneously, so **Cyclin D1** can be drawn as a single node fed by the mRNA.

## Edge list for this module

| from | to | type | label |
|---|---|---|---|
| SHH | Ptch1 | binds (arrow) | |
| Ptch1 | Smo | inhibition | brake × functional Ptch1 |
| vismodegib | Smo | drug inhibition | |
| Smo | Gli toggle | activation | Gli-R → Gli-A |
| Gli-A | Ptch1 | activation, curved back | negative feedback, lost in MB |
| Gli-A | Gli1 | activation | pathway readout |
| Gli-R | Gli1 | inhibition (thin) | |
| Gli-A + Gli1 | MYCN | activation | |
| Gli-A + Gli1 | *Ccnd1* promoter | activation, thick | main axis |
| Gli-R | *Ccnd1* promoter | inhibition (thin) | |
| MYCN | *Ccnd1* promoter | activation | floor × amplification |
| Smo | Gli1 memory | activation, dashed | MB only (optional) |
| Gli1 memory | Gli1 | activation, dashed | t½ 14 h (optional) |
| *Ccnd1* promoter | Cyclin D1 | via mRNA, thick | hand-off to Module 2 |

## Layout note

Membrane strip along the top with SHH outside, Ptch1 and Smo in it. Below the membrane, left to
right: Gli toggle, then Gli1 and MYCN side by side, then the *Ccnd1* promoter at the right edge
of the zone. The Ptch1 feedback arc goes up and left, the only arrow that travels backward.
Vismodegib sits above the membrane over Smo.
