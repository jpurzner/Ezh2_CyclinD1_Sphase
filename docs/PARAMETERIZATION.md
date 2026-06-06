# Data used to parameterize / validate the v44 model

Every quantitative constraint the model was fit to or checked against, with its source. These are
encoded as the target checks in `simulations/validate_v44.py` (run it to score the current model).
The model is calibrated against **ratios and behaviors**, not absolute concentrations (species are in
arbitrary units). Full target compendium with caveats: `docs/Ezh2_CcnD1_model_targets.md`.

## A. Bulk RNA‑seq — the primary dataset (normalized counts, mean per group)

Cerebellar granule‑lineage / MB bulk RNA‑seq (Purzner lab). These are the canonical numbers the
HH/MYCN→CyclinD1 module is fit to.

| Group | n | Gli1 | Ezh2 | Ccnd1 | Mycn | Mki67 | Atoh1 |
|---|---|---|---|---|---|---|---|
| E15_GNP | 2 | 2,283 | 6,843 | 5,372 | 4,008 | 28,674 | 2,361 |
| P1_GNP | 3 | 3,169 | 4,876 | 7,999 | 3,959 | 25,304 | 1,302 |
| **P7_GNP_wt** (model "GNP") | 5 | **2,591** | **4,996** | **3,639** | **3,077** | 32,614 | 661 |
| P14_GNP | 3 | 1,665 | 3,037 | 2,513 | 1,889 | 14,598 | 347 |
| P7_GNP_Ptch (Ptch+/−) | 2 | 3,775 | 8,250 | 6,942 | 4,251 | 65,326 | 2,391 |
| P28_GNP_Ptch (Ptch+/−) | 4 | 10,666 | 6,641 | 16,056 | 4,904 | 25,561 | 2,339 |
| **MB_Ptch_het** (model "MB") | 2 | **17,881** | **10,172** | **27,590** | **8,808** | 47,839 | 7,018 |
| MB_EZH2ko | 3 | 9,757 | 11,128 | 12,581 | 8,554 | 62,600 | 6,453 |
| ~~MB_JQ1~~ (excluded for now) | 4 | 7,932 | 8,338 | 14,033 | 9,399 | 37,117 | 7,282 |
| **MB_GDC0449** (model "MB+HHi") | 2 | **413** | **6,621** | **3,971** | **7,862** | 55,414 | 10,399 |
| P28_GNP_Ptch_GDC0449 | 3 | 229 | 3,510 | 1,726 | 5,308 | 9,544 | 2,680 |
| **P7_GNP_wt_GDC0449** (model "GNP+HHi") | 4 | **10** | **3,689** | **572** | **2,474** | 12,154 | 381 |

**Derived MB/GNP folds (MB_Ptch_het ÷ P7_GNP_wt) — the cross‑context calibration targets:**

| gene | MB/GNP fold | role |
|---|---|---|
| **Gli1** | **6.90×** | Hedgehog activity — must NOT saturate (model previously gave only ~1.1×) |
| **Ccnd1** | **7.58×** | CyclinD1 — predominantly **Gli‑driven** (see below) |
| Ezh2 | 2.04× | cell‑cycle‑coupled |
| Mycn | 2.86× | `MYCN_amplification` |
| Mki67 | 1.47× | proliferation index |
| Atoh1 | 10.6× | SHH‑MB identity |

**CyclinD1 tracks Gli1 ~1:1 along the Ptch escalation series → CyclinD1 is Gli‑driven, not MYCN‑driven:**

| | Gli1 (fold) | Ccnd1 (fold) | Mycn (fold) |
|---|---|---|---|
| P7_GNP_wt | 1.0× | 1.0× | 1.0× |
| P7_GNP_Ptch | 1.5× | 1.9× | 1.4× |
| P28_GNP_Ptch | 4.1× | 4.4× | 1.6× |
| MB_Ptch_het | 6.9× | 7.6× | 2.9× |

(Gli1 and Ccnd1 rise in lockstep; Mycn rises only ~3×. So the MB CyclinD1 elevation is driven by the
Gli arm. Note: the P28_Ptch point is age‑confounded — same genotype as P7_Ptch but older — so for a
clean Ptch1‑dosage anchor use GNP / P7_Ptch / MB.) **Developmental DECLINE P7→P14** (differentiation):
Ezh2 0.61×, Ccnd1 0.69×, Mki67 0.45×, Gli1 0.64× — all fall together as GNPs exit/differentiate
(confirms the EZH2 drop at exit).

## A2. Bulk RNA‑seq — within‑context ratios (drug / cell‑type)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Cyclin D1: GNP+HHi / GNP | 0.157 (84% ↓) | Gli‑driven fraction of CyclinD1; HH pathway gain | RNA‑seq (P7_GNP_wt_GDC0449) |
| **Cyclin D1: MB+HHi / MB** | **0.144 (86% ↓)** | vismo crashes MB CyclinD1 to **~1.1× a cycling GNP** (≈ threshold); MYCN/basal = the small HHi‑resistant residual | RNA‑seq (MB_GDC0449) |
| Gli1: MB+HHi / MB | 0.023 (98% ↓) | vismo collapses Gli1 in MB | RNA‑seq (MB_GDC0449) |
| **Cyclin D1: MB / GNP** | **7.58×** (raw table; 5.07× by the Fig 4I normalization) | overall MB CyclinD1 drive — Gli‑driven | RNA‑seq |
| **Gli1: MB / GNP** | **6.90×** | HH activity must scale with Ptch1 dosage (de‑saturated) | RNA‑seq |
| Mycn: GNP+HHi / GNP | 0.78 (22% ↓) | Gli‑dependent vs autonomous Mycn synthesis | RNA‑seq |
| Mycn: MB+HHi / MB | 0.86 (14% ↓) | autonomous (amplified) Mycn fraction | RNA‑seq |
| Mycn: MB / GNP | 2.86× | `MYCN_amplification` (MB context) | RNA‑seq |
| Gli1: GNP+HHi | >99% ↓ | HHi (vismodegib) block on Smo→Gli | RNA‑seq |
| Ezh2: MB / GNP | 2.04× (p<0.001) | cell‑cycle‑coupled EZH2 (soft — see note 1) | RNA‑seq, Fig. 4J |
| Ezh2: GNP+HHi | 25–47% ↓ | EZH2 falls when cycling stops | RNA‑seq |

**Why the MB_GDC0449 row matters (vismo on MB).** Vismodegib collapses MB CyclinD1 by 86% (to 0.144
of untreated MB), landing it at **~1.1× a cycling P7 GNP's CyclinD1 — i.e. right at the proliferative
commitment threshold**, not the ~3× I had previously inferred from an over‑estimated 0.40 residual.
This is decisive: at threshold, (i) a heterogeneous tumour partly drops below it → vismo *reduces the
cycling fraction*, and (ii) EZH2i de‑repression (~2×) pushes CyclinD1 back above → *restores cycling*.
So the EZH2i rescue of vismo IS a CyclinD1‑threshold effect after all, now with direct data support;
the calibration target for the HH module is therefore `CyclinD1 MB+HHi/MB = 0.144` (Gli must supply
~86% of MB CyclinD1). *Caveat:* MB_EZH2ko CyclinD1 is *lower* than MB (0.46×) but Gli1 also drops
(0.55×) — that tumour is a differentiated, lower‑Hedgehog state, confounded, and is not used as a
CyclinD1‑repression target (the clean EZH2→CyclinD1 strength comes from the GNP cKO/Taz data, §B).

## A3. CDK inhibitors — the commitment threshold (p16/p21/p27, normalized counts)

| Group | p21 (Cdkn1a) | p27 (Cdkn1b) | p16 (Cdkn2a) |
|---|---|---|---|
| P7_GNP_wt | 890 | 9,168 | **3.0** |
| MB_Ptch_het | 2,511 | 10,893 | **417** |
| MB_GDC0449 (vismo) | 2,339 | **14,261** | **842** |
| P7_GNP_wt_GDC0449 | 893 | **11,632** | 2.8 |
| P7_GNP_EZH2ko | 1,463 | 6,823 | **135** |

**The vismo arrest in MB is multi‑pronged, and only part is EZH2i‑rescuable:**
- **p16 (Cdkn2a)** is silent in GNPs (~3) and ~100–200× induced in MB (~400–840) — a classic
  H3K27me3‑silenced locus (GNP EZH2ko derepresses it 3→135). p16 is THE CDK4/6 inhibitor. Modeled as a
  **competitive** CDK4/6 brake: it raises the CyclinD1 half‑max for Rb phosphorylation,
  `kPhRbCd·Cd/(K_CdRb·(1+p16) + Cd)` (new param `p16`: 0 in GNP, =3 in MB → eff K_CdRb 0.5→2.0). This
  raises the CyclinD1 threshold to commit, so vismo's CyclinD1 drop now crosses it (arrest) — but MORE
  CyclinD1 still overcomes it, so **EZH2i (CyclinD1 up) rescues** a p16‑braked arrest. Mechanistically
  DISTINCT from CDK4/6i (palbociclib = `kPhRbCd=0`, a Vmax block) which raising CyclinD1 cannot bypass
  → **not rescuable** (the model's clean contrast).
- **p27 (Cdkn1b)** rises modestly in MB and is **highest in all vismo (GDC0449) samples** (Shh blockade
  → cell‑cycle‑exit signal) — supports the arrest direction; not separately encoded (the p16 brake +
  CyclinD1 drop already produce it; adding vismo→p27 would over‑arrest).
- **p21 (Cdkn1a)** ~2.5–4× higher in MB; minor.

This is why vismo drops proliferation to ~1/4 (a CyclinD1‑only model only reached ~86%), and why the
EZH2i rescue tops out at ~2/3–3/4 (the non‑rescuable p16/Vmax‑like fraction sets the ceiling).
*Note:* p16's H3K27me3 control is GNP‑specific — in MB p16 has already escaped silencing (high
regardless), so EZH2i in MB does not re‑induce it; the CyclinD1 arm dominates → rescue.

## B. EZH2 perturbation → Cyclin D1 (the feedback strength, `K_EZH2_repression`)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Cyclin D1: **Ezh2 cKO** GNP / WT | log2FC 1.301 → **2.46×** | EZH2 ⊣ CyclinD1 repression strength | RNA‑seq, Fig. 3A |
| Cyclin D1: GNP + Tazemetostat / DMSO | 2.2× | same (catalytic inhibition) | qPCR, Fig. 3C |
| Cyclin D1: MB + Ezh2‑OE / control | significant ↓ | direction of the repression | qPCR, Fig. 3E |
| Cyclin D1: MB + Palbociclib / DMSO | +57% (1.57×) | CDK4/6i ↑ CyclinD1 (via ↓EZH2) | qPCR, Fig. 4F |
| Ezh2, E2f1: MB + Palbociclib | −55.6%, −56.7% | CDK4/6→Rb‑E2F→EZH2 arm | qPCR, Fig. 4F |

> **Note 1 — the cKO 2.46× / Tazemetostat 2.2× are LOWER BOUNDS on the true feedback.** The model
> gives ~3× de‑repression (`K_EZH2_repression = 0.5`). The bulk cKO under‑reads the per‑cell effect
> because (i) Math1‑Cre is incomplete, so WT escapers dilute the signal — at 70–80% recombination the
> bulk 2.46× implies a **per‑cell 2.8–3.1×**, matching the model — and (ii) the Ezh2‑flox deletes the
> SET (methyltransferase) domain but leaves EZH2 protein (residual non‑catalytic repression), and
> H3K27me3 dilutes only with division, so even the per‑cell fold is a floor. The model's 3× is the
> recombination‑corrected value and is plausibly conservative (`fig_v44_ckoo_feedback_strength.py`).

## C. scRNA‑seq — Ezh2 transcript by cell‑cycle phase (the EZH2 cell‑cycle gate)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Ezh2 cycling / G0 (transcript) | ~2.0× (S/G0 ≈ 1.8–2.5 across Sepp, Vladoiu, Ocasio, Luo) | EZH2 transcription gated on E2f×(CycE+CycA) | scRNA‑seq, Fig. 4A,B |

## D. Immunofluorescence — Ezh2 protein by phase / under arrest drugs

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Ezh2 protein G2 / G0 (DMSO) | 1.48× | EZH2 protein accumulation across S/G2 | IF, Supp. 6D |
| Ezh2 in S: HU 10 µM / DMSO | 1.31× (alt 1.22×) | EZH2 integrates S‑phase **duration** (HU→fork speed) | IF, Fig. 4G |
| rShh dose → Ezh2 | dose‑dependent ↑ | Hh→CyclinD1→cycle→EZH2 (whole loop) | IF, Fig. 4H |

## E. Cell‑cycle phase proportions (G0/G1/S/G2)

| Dataset | Values | Constrains | Source |
|---|---|---|---|
| MB, DMSO (renormalized) | G0 24.8 / G1 43.4 / S 15.7 / G2 16.1 % | phase‑duration balance; G0 by p27 marker | Fig. 4G / Supp. 6, Sec. F |
| MB, HU 10 µM (fold vs DMSO) | S ×1.36 ↑, G2 ×0.23 ↓ | HU lengthens S (replication‑fork model) | Sec. F |
| GNP + Tazemetostat (fold) | G0 ×0.92, G1 ×2.03, G2 ×1.44 | EZH2i shrinks G0 / drives cycling | Fig. 3B |
| MB + Ezh2‑OE (fold) | G0 ×2.84, G1 ×0.51 | EZH2‑OE traps in G0 | Fig. 3F |

## F. Dynamics

| Measurement | Value | Constrains | Source |
|---|---|---|---|
| Cell‑cycle period, GNP | ~22–23 h (soft) | overall timescale (growth rate `mu`, fork speed) | Published (Nakashima 2015, ref. 70) |

---

## Optimized / fit parameters

A constrained search (`simulations/v44_recalibrate_gli.py`) over the HH/MYCN→CyclinD1 and EZH2
parameters de‑saturates the Hedgehog arm and fits the bulk RNA‑seq folds **including the
MB_GDC0449 (vismodegib‑on‑MB) point**, with a hard guard requiring the figure‑critical behaviors
(GNP cycles only with SHH; GNP−SHH / GNP+HHi arrest; MB cycles; MB+CDK4/6i arrests and is not rescued
by Ezh2i). Key fit values (full list in `src/build_model_v44_heldt.py`):

- **Gli→Ptch1 negative feedback (NEW).** Ptch1 is a Gli target: `Ptch1_transcription = k_Ptch1_basal
  + k_Ptch1_Gli·Gli_act²/(K_Gli_Ptch² + Gli_act²)`, and functional Ptch1 re‑inhibits Smo
  (`Smo = k_Smo_act/(1 + Ptch1_free·Ptch1_copy_number/K_Ptch_Smo)·(1−GDC0449)`), closing the loop
  Gli→Ptch1⊣Smo⊣Gli. `Ptch1_copy_number` is reinterpreted as the **functional Ptch1 fraction**
  (1.0 GNP, 0.5 Ptch1⁺/⁻, 0.1 MB‑with‑LOH) and now gates Smo‑repression (function), not transcription.
  GNP shows an adaptive Gli1 overshoot to a Shh step; MB is a **broken loop** (induced Ptch1
  non‑functional → constitutive Gli, high Ptch1 mRNA). Fit: `k_Ptch1_basal 0.15, k_Ptch1_Gli 1.87,
  K_Gli_Ptch 0.40`. See `fig_v44_ptch1_feedback.py`. **The broken feedback is the direct cause of high
  MB Gli1:** restoring functional Ptch1 (f→1) in the MB context collapses Gli1 to the GNP level (1.0×),
  even with MYCN still amplified (MYCN is downstream of Gli). The functional‑Ptch1 sweep gives Gli1
  MB/GNP = 1.0 (f=1) → 4.7 (f=0.2) → **7.2 (f=0.1, matches data 6.9×)** → 12 (f=0); so the measured fold
  pins f≈0.1 = ~10% residual feedback (Ptch1⁺/⁻ haploinsufficiency), not a complete null.
- Gli arm de‑saturated (`k_Cd_tx_Gli_max 46`, `K_Gli_act_CycD`, `Vmax_Gli1_tx`, `K_Gli_act_Gli1`,
  `K_Ptch_Smo`, `K_Smo_Gli_switch`, `k_Smo_act`) so Gli1 (7.2×) and CyclinD1 (6.8×) rise GNP→MB and
  vismodegib collapses MB CyclinD1 by ~86% (MB+HHi/MB 0.14 model vs 0.144 data) — Gli supplies ~86%.
- MYCN→CyclinD1 cooperative Hill (`n_MYCN_Cd ≈ 3.7`, `k_Cd_tx_MYCN ≈ 35`, `K_MYCN_Cd ≈ 1.7`) — the
  small HHi‑resistant residual; `MYCN_amplification = 2.86` (MB) from the Mycn MB/GNP ratio.
- `k_Cd_translation = 0.8` — CyclinD1 protein scale set so GNP robustly clears the cycling threshold
  (the de‑saturation dropped GNP CyclinD1 toward a bistable knife‑edge; 0.8 lifts GNP/MB+HHi clearly
  above it). With the saturating CyclinD1→Rb drive (`with_cd_sat`) this stays numerically stable.
- `K_EZH2_repression = 0.5` — EZH2 ⊣ CyclinD1 feedback strength (~3× de‑repression; see note 1).
- EZH2 gate + turnover, and `KmHU_fork` / `vmin_fork` (HU→fork‑speed, S‑phase lengthening) — unchanged.

**The data‑matched rescue (population).** With vismodegib landing MB CyclinD1 at ~1× a cycling GNP
(onto the commitment threshold), a heterogeneous MB ensemble (N=120, CyclinD1/p27 spread;
`fig_v44_fig5_population.py`) cycles 100% (MB) → 86% (MB+HHi) → 93% (MB+HHi+Ezh2i); CDK4/6i holds it
at 0% ± Ezh2i. Single cells still cycle under vismodegib (MB+HHi Ki67 is high), so the rescue is a
population/fractional Cyclin D1‑threshold effect, not an all‑or‑nothing single‑cell arrest; its
magnitude scales with the assumed cell‑to‑cell heterogeneity (the direction is robust).

**Soft / not‑enforced targets**: Ezh2 MB/GNP 2.05× (model ~1.1×; cell‑cycle‑coupled); Ptch1 mRNA
MB/GNP (model 2.3×, elevated but no exact bulk value); and the absolute S/G2 phase proportions
(S is biochemically long in the real‑time core).
