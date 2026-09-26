# Battich 2020 (scEU-seq): mRNA degradation is cycle-modulated, and v44 assumes it is not

Supplied 2026-09-26 (`papers/battish/`); was blocked to automated download. This is the last of
the requested datasets and it turns out to carry the most direct structural implication for v44.

**Why it is different from everything else we hold.** Every other stability dataset gives *one*
half-life per gene. Battich measured synthesis and degradation **simultaneously and resolved
across the cell cycle**, in **RPE1-FUCCI — the line Heldt 2018 was fitted in**. For a model whose
entire subject is cycle-dependent regulation, that is the relevant measurement.

---

## 1. Reading the table

`aax3072_table-s1.csv`: 528 genes × 301 cell-cycle-pseudotime bins.

* `norm_gamma`, `norm_kappa` = **log2(rate ÷ that gene's median rate along the cycle)**. Shape and
  fold-range are informative; absolute magnitude is not. (An earlier pass here misread these as
  min–max normalised and produced nonsense fold-ranges — they are log2 ratios.)
* `cosine_similarity` = cos(norm_gamma, norm_kappa). **Verified by recomputing it from the
  trajectories for all 528 rows: median |difference| = 0.0000.** So the convention is unambiguous:
  * **negative** → degradation *falls* as synthesis rises → the two **cooperate** to swing expression
  * **positive** → degradation *rises* with synthesis → degradation **damps** the swing
* `strategy_group` A–G: SOM clusters of the joint (γ, κ) trajectory shape.

---

## 2. The headline: degradation is not a constant

| | value |
|---|---|
| genes whose **degradation rate** varies **> 2-fold** across the cycle | **90%** |
| median degradation fold-range | **3.3×** |
| p75 / p90 | 5.1× / 8.0× |
| cooperative (cos < −0.5) | 35% |
| neutral | 58% |
| destabilizing (cos > +0.5) | 7% |

**v44 uses constant mRNA degradation rates throughout** — `k_Cd_mRNA_deg`, `kDeEZm`,
`k_Gli1_mRNA_deg`, `k_Ptch1_mRNA_deg` are all fixed numbers. Battich says that is wrong for ~90% of
genes by a median factor of 3.3, and that where it is wrong it is usually wrong *cooperatively*:
degradation drops exactly when synthesis rises, amplifying the cycle swing.

A constant-`k_deg` model therefore has to generate the whole observed cycle amplitude from
synthesis modulation alone, and will need a **larger** synthesis swing than reality to match the
same expression swing. That is a systematic, directional bias — not noise.

### Our genes

| gene | cos.sim | strategy | γ fold | γ peak @ | κ fold | κ peak @ |
|---|---|---|---|---|---|---|
| **EZH2** | **−0.554** | **cooperative** | **2.96** | **0.20** | **4.02** | **0.39** |
| CCNA2 | −0.740 | cooperative | 27.48 | 0.29 | 14.51 | 0.87 |
| CCNE2 | −0.797 | cooperative | 6.00 | 0.12 | 4.64 | 0.53 |
| CDK1 | −0.846 | cooperative | 4.56 | 0.33 | 12.25 | 0.89 |
| CDK2 | −0.721 | cooperative | 5.77 | 0.22 | 3.54 | 0.83 |
| CDC20 | −0.714 | cooperative | 3.02 | 0.34 | 24.14 | 0.87 |
| CDKN2C | −0.646 | cooperative | 5.43 | 0.37 | 10.17 | 0.88 |
| CCNB1 | −0.568 | cooperative | 5.33 | 0.23 | 19.02 | 0.88 |
| CCND3 | −0.212 | neutral | 13.05 | 0.15 | 2.39 | 0.63 |
| EED | −0.468 | neutral | 4.09 | 0.17 | 3.73 | 0.74 |
| GLI3 | −0.483 | neutral | 2.94 | 0.29 | 2.93 | 0.64 |
| CHEK1 | −0.443 | neutral | 4.66 | 0.24 | 2.39 | 0.42 |
| PCNA | −0.400 | neutral | 3.36 | 0.21 | 3.77 | 0.55 |
| DTL | −0.231 | neutral | 5.55 | 0.22 | 5.34 | 0.41 |
| CDKN1A | −0.081 | neutral | 1.69 | 0.82 | 5.73 | 0.69 |
| CDC25A | +0.207 | neutral | 3.17 | 0.26 | 3.49 | 0.45 |
| E2F1 | +0.234 | neutral | 3.54 | 0.67 | 7.34 | 0.47 |

17/38 of our genes are covered. Absent: `Ccnd1`, `Skp2`, `Rb1`, `Mycn`, `Gli1`, `Gli2`, `Ptch1`,
`Smo`, `Cdkn1b`, `Ccne1`, `Cdk4`, `Cdk6`, `Cdt1`, `Wee1`, `Fbxo5`, `Fzr1`, `Trp53`, `Suz12`.
Notably the two nodes the noise work cares about most — `Ccnd1` and `Skp2` — are **not** in it.

---

## 3. This may be the missing mechanism behind v44's worst standing failure

`EZH2 transcript S/G0` is the model's largest standing miss: **1.563 against a target of ~2**
(it passes on tolerance, but it is the one the README singles out). And EZH2 is exactly a
cooperative gene:

* degradation varies **2.96-fold**, peaking **early** (pseudotime 0.20);
* synthesis varies **4.02-fold**, peaking **later** (0.39);
* they are anti-phased (cos = −0.554), so the two **reinforce** the S-phase rise.

v44 has a **constant** `kDeEZm` and must therefore produce the whole S/G0 gradient from
transcription alone — which is precisely where it falls short. Adding cooperative cycle-modulated
degradation would amplify the gradient **in the direction the model is currently missing**.

### It may also dissolve the F3 conflict

`protein_noise_fixes.md` §F3 recorded an apparently irreconcilable conflict: the model needs
`kDeEZm` ≈ 1.2 h to reproduce the EZH2 phase gradients, while Schwanhäusser measures 9.85 h and
Herzog's SLAM-seq measures 3.71 h. Forcing the measured value broke 7 targets.

Battich suggests a third reading: **the model needs an artificially short constant half-life
because it is compensating for missing cycle-modulation.** With degradation free to drop 3-fold
during S, a *longer mean* half-life could still generate the observed gradient. If so, F3 is not a
model-vs-data conflict at all — it is a missing mechanism, and the fix is structural rather than a
parameter the bench has to arbitrate.

**This is a hypothesis, not a result.** It is cheap to test and I have not run it.

### Concrete test

Make `kDeEZm` cycle-dependent behind a flag (`with_cyclic_mrna_deg`, default OFF), phased against
the existing cycle variables, and re-run the suite:

```
kDeEZm_eff := kDeEZm * (1 + a_mk*(gS - <gS>))     // a_mk < 0 = cooperative
```
where `gS` is an existing S-phase indicator (`Ce`, `Ca`, or `aRc`). Then ask whether
`EZH2 transcript S/G0` reaches ~2 **with `kDeEZm` set to Herzog's 3.71 h** rather than the current
1.18 h. Amplitude prior: γ fold-range ≈ 3, phase ≈ a quarter-cycle ahead of synthesis.

If it works, it converts a standing validation miss and a data conflict into one structural fix.
I have not done it — say the word.

---

## 4. Caveats

* **RPE1, not GNP.** Right cell line for the Heldt core, wrong organism and lineage for the
  biology. Cycle *phasing* is likely conserved; fold-ranges less so.
* **528 genes only**, and biased toward well-expressed cycle-regulated ones. `Ccnd1` and `Skp2`
  are missing, which is unlucky given they are the two nodes the noise work most depends on.
* **Relative, not absolute.** The table gives log2 ratios to each gene's own cycle median, so it
  constrains *modulation*, not the mean rate. Pair it with Herzog (mESC) or the model's own
  constant for the mean.
* **Inference-dependent.** Rates come from fitting a kinetic model to labelled/unlabelled counts.
  The authors' own error analysis (fig. S12) reports that mis-calling is asymmetric — true
  destabilizing strategies are not mis-called as cooperative, but the converse can happen — so the
  35% cooperative fraction is arguably a floor.

---

## 5. Status of the request list

All four originally requested turnover datasets, plus the three follow-ups, are now in hand:

| dataset | status |
|---|---|
| Schwanhäusser 2011 | downloaded |
| Mathieson 2018 | downloaded |
| Larsson 2019 (S1, S2, S3) | downloaded |
| Sigal 2006 | downloaded |
| Eden 2011 | **user-supplied** |
| McShane 2016 | **user-supplied** |
| Ochiai 2020, Herzog 2017 | downloaded |
| **Battich 2020** | **user-supplied** |

Nothing from the literature is outstanding. The remaining gaps are all **your own data** —
`protein_noise_data_requests.md` §A: the HCI per-cell tables, an absolute scale for Ω, and a
measured T_cc.
