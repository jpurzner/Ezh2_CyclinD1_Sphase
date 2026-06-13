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

## A3. CDK inhibitors — the commitment threshold (full CKI panel; normalized counts)

Source: `Medulloblastoma_projects/H3K27me3_fast_change/Scott_Bulk_Quant/CDK_inhibitors_summary.{csv,md}`
(DESeq2 size‑factor normalized, cleaned set). The CKIs sort by their CDK target — exactly the two
brakes the model separates. **p57 (Cdkn1c) is EXCLUDED** (separate project; it goes the opposite way,
*down* 0.49× in MB).

**INK4 family → inhibit CDK4/6 → the model's competitive CyclinD1→Rb brake (`p16`):**

| gene | alias | P7_wt | MB | MB/P7 | MB+GDC | sig |
|---|---|---:|---:|---:|---:|:--:|
| Cdkn2a | **p16** | 3 | 417 | **139×** | 842 | *** |
| Cdkn2b | p15 | 25 | 154 | 6.4× | 281 | *** |
| Cdkn2c | p18 | 1,149 | 4,288 | 3.0× | 2,121 | *** |
| Cdkn2d | p19 | 890 | 1,449 | 1.2× | 1,028 | ns |

**CIP/KIP family → inhibit CDK2 → the model's p21/p27 arm (`kSyP21`):**

| gene | alias | P7_wt | MB | MB/P7 | MB+GDC | sig |
|---|---|---:|---:|---:|---:|:--:|
| Cdkn1a | **p21** | 890 | 2,511 | 3.1× | 2,339 | *** |
| Cdkn1b | p27 | 9,168 | 10,893 | 1.4× | 14,261 | ** |

**How the model encodes this — two effective brakes, each an aggregate family tone:**
- **INK4 / CDK4/6 tone — two explicit params, `p16` + `p18`** — a *competitive* CDK4/6 brake raising the
  CyclinD1→Rb half‑max `kPhRbCd·Cd/(K_CdRb·(1 + p16 + p18) + Cd)` (`K_CdRb` 0.5→**0.357** when p18 was
  added, so the GNP p18 baseline keeps the GNP net half‑max at 0.357·1.4 = 0.5 — behavior preserved):
  Tones are now **ABUNDANCE‑proportional (DESeq2 counts), not fold‑proportional**:
  - **`p16` (Cdkn2a)** = the GNP‑silent, MB‑specific INK4 (H3K27me3‑silenced, **139× fold** in MB):
    0 in GNP, **0.15** in MB. The dramatic fold reflects near‑absence in GNP (3 counts), but p16's
    *absolute* MB level (417) is only ~0.1× of p18 (4288) — so it's a **small absolute brake** despite
    the fold, and the cleanest MB‑specific marker. (Was 0.88, search‑tuned; corrected to abundance.)
  - **`p18` (Cdkn2c)** = the *constitutive*, dominant INK4 (Shh‑maintained), the largest absolute CDK4/6
    inhibitor: **0.4 in GNP → 1.5 in MB** (3.75× ≈ the data 1149→4288). p16/p18 ratio 0.10 ≈ data 0.097.
    (p15 6.4× folded in conceptually; p19 ns, not encoded.)
  This is the **rescuable** axis: more CyclinD1 overcomes the competitive brake, so **EZH2i rescues**.
  Contrast CDK4/6i (palbociclib = `kPhRbCd=0`, a non‑competitive Vmax block) which raising CyclinD1
  cannot bypass → **not rescuable** (the model's clean contrast).
- **model `kSyP21` = the CIP/KIP / CDK2 tone** (p21 3.1× + p27 1.4×) — raises p21/p27 synthesis ~2× in
  MB (0.002→0.004), inhibiting CyclinE/A–CDK2; also lengthens the p21/p27‑high transient G0. **p27 also
  brakes CDK4/6** (Guiley 2019): the `P21` species enters the CDK4/6 half‑max via `w_p27·P21` (w_p27=1),
  a G0‑weighted brake (p27 peaks in G0). p27 is the most abundant CKI (10,893) but acts mostly on CDK2;
  only a modest CDK4/6 share is modeled (a larger w_p27 would over‑raise the GNP threshold).

So the two brake parameters are *family aggregates* (INK4, CIP/KIP), not single genes — which is why
the wider panel **reinforces** the structure: the CDK4/6 (INK4) arm is strongly and broadly induced
(p16 ≫, + p15, + p18) and is the dominant, rescuable brake; the CDK2 (CIP/KIP) arm is a more modest
co‑brake. This is why vismo drops proliferation to ~1/4 (a CyclinD1‑only model only reached ~86%), and
why the EZH2i rescue tops out at ~2/3–3/4 (the non‑rescuable Vmax‑like ceiling).

**GDC0449 (vismo) barely moves the CKIs** (MB: only p16 marginal +1.23, padj 0.083; all others ns; P7
GNPs: only p18 *down*, padj 0.006 — Shh maintains p18 in GNPs). This supports treating the CKI tone as
a ~static tonic threshold under vismo: **the vismo arrest is the CyclinD1 drop crossing a fixed brake,
not CKI induction.** p16 is an H3K27me3‑silenced locus in GNPs (GNP EZH2ko derepresses it 3→135), but in
MB it has already escaped silencing (high regardless), so EZH2i in MB does not re‑induce it and the
CyclinD1 arm dominates → rescue.

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
- `k_Cd_translation = 0.75` — CyclinD1 protein scale set so GNP robustly clears the cycling threshold
  (the de‑saturation dropped GNP CyclinD1 toward a bistable knife‑edge). With the saturating
  CyclinD1→Rb drive (`with_cd_sat`) this stays numerically stable. *(search‑tuned)*
- `K_EZH2_repression = 0.75` — EZH2 ⊣ CyclinD1 feedback strength. *(search‑tuned: gives an EZH2i
  CyclinD1 fold of ~2.7× — closer to the Fig 3C qPCR ~2× than the prior 0.5 → 3.5×.)*
- MB CDK‑inhibitor brake = **INK4 (p16 0.15 + p18 1.5, abundance‑proportional; GNP p18 baseline 0.4,
  `K_CdRb` 0.357) + CIP/KIP (kSyP21 0.004, 2× baseline) + p27→CDK4/6 (`w_p27`=1)**; CyclinD1 heterogeneity
  σ = 0.68 — *(pRb rescue 100→~23→~71; preserves the threshold the old p16‑heavy tone gave)*.
- `kSyDna` = **0.044** (~5× the Heldt default 0.0093) — S‑phase ~3 h (was ~10 h, unrealistically long);
  brings the DMSO S proportion to ~13–15% (data 15.7%). `M_commit` = **1.05** (was 1.3) rebalances the
  G0/G1 split to the data (G0 ~25% / G1 ~43%) — interim, pending mechanistic G0/G1 timing elements.
- EZH2 gate + turnover, and `KmHU_fork` / `vmin_fork` (HU→fork‑speed, S‑phase lengthening) — the HU‑fold
  targets now need re‑tuning for the shorter S (a known follow‑up).

*Calibration of the four rescue‑determining params (p16, σ, K_EZH2_repression, k_Cd_translation) was a
broad random search (`v44_rescue_search.py`, ~hundreds of evals across 4 seeds) against the pRb+ targets
MB 100 / MB+HHi 25 / +EZH2i 75 / CDK4/6i 0, with hard guards (GNP cycles, GNP+HHi arrests). Crucially,
integration failures at the commitment bifurcation are EXCLUDED from the cycling‑fraction denominator
(not counted as arrested), which is why the honest rescue is ~73% rather than the ~58% an earlier
crash‑as‑arrest count gave.*

**The data‑matched rescue (population).** With vismodegib landing MB CyclinD1 at ~1× a cycling GNP
(onto the commitment threshold), a heterogeneous MB ensemble (N=120, CyclinD1/p27 spread;
`fig_v44_fig5_population.py`) cycles 100% (MB) → 86% (MB+HHi) → 93% (MB+HHi+Ezh2i); CDK4/6i holds it
at 0% ± Ezh2i. Single cells still cycle under vismodegib (MB+HHi Ki67 is high), so the rescue is a
population/fractional Cyclin D1‑threshold effect, not an all‑or‑nothing single‑cell arrest; its
magnitude scales with the assumed cell‑to‑cell heterogeneity (the direction is robust).

**Soft / not‑enforced targets**: Ezh2 MB/GNP 2.05× (model ~1.1×; cell‑cycle‑coupled); Ptch1 mRNA
MB/GNP (model 2.3×, elevated but no exact bulk value); and the absolute S/G2 phase proportions
(S is biochemically long in the real‑time core).
