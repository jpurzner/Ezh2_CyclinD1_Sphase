# v44 Validation Targets — formal specification

**Purpose.** One reviewable place that says exactly what the model is fit to, with **cell type**, **provenance**, and **current model value** on every target. Prevents the MB↔GNP conflation and provenance drift that a 3-source cross-audit (harness vs compendium vs provenance doc, 2026-07-25) surfaced.

- **Canonical machine-readable version:** [`data/validation_targets.json`](../data/validation_targets.json) — edit that first; this file is the readable view.
- **Companion docs:** [`docs/Ezh2_CcnD1_model_targets.md`](Ezh2_CcnD1_model_targets.md) (raw manuscript values, Chahin et al.) · [`docs/DATA_PROVENANCE.md`](DATA_PROVENANCE.md) (parameter→data audit) · [`docs/cdk46i_resistance_research.md`](cdk46i_resistance_research.md) (palbo reversibility).
- **Live harness:** `simulations/validate_v44.py` — **29 scored checks** (the "28/29"). This spec lists those 29 **plus one non-formal target** (the vismo rescue fractions). Two of the 29 are flagged **⚠ harness-mismatch** (period, palbo) because the current check is a proxy that should change.
- **Model values** are from a `validate_v44` run on 2026-07-25 (default model, **28/29**; lone miss = MB HU G2 fold).

## Cell types (keep separate)

| Tag | Meaning |
|---|---|
| **P7 GNP** | postnatal-day-7 cerebellar granule-neuron progenitor |
| **MB** | SHH-medulloblastoma (culture) — *longer cycle, higher CKI, transient G0s* |
| **MB↔GNP** | cross-type ratio; constrains the MB-vs-GNP difference |

## Provenance tags

| Tag | Meaning |
|---|---|
| **manuscript** | Chahin et al. figure (primary lab data) |
| **ext-lit** | published external citation (Contestabile, Fujita, Nakashima, Wechsler-Reya, Heldt) |
| **unpub** | JP lab data not pinned to a figure |
| **prediction** | model output presented as a prediction |
| **assum** | modeling choice |

---

## P7 GNP targets

| Target | Quantity | Target | Tol | Model | ✓ | Src |
|---|---|---|---|---|---|---|
| CyclinD1 GNP+HHi / GNP | Ccnd1, vismodegib vs +Shh | 0.157 | 0.40 | — | ✓ | manuscript §A |
| MYCN GNP+HHi / GNP | Mycn, HHi vs +Shh | 0.78 | 0.30 | — | ✓ | manuscript §A · *model ~0.91, soft* |
| Gli1 GNP+HHi reduction | 1 − Gli1(HHi)/Gli1(+Shh) | 0.99 | 0.05 | — | ✓ | manuscript §A |
| EZH2i CyclinD1 fold | Ccnd1, **Tazemetostat** vs +Shh | 2.2 | 0.42 | **1.63** | ✓* | **manuscript Fig 3C** (qPCR) |
| EZH2 G0/cycling | EZH2 protein, serum-starved vs cycling | 0.6 | 0.40 | 0.61 | ✓ | manuscript §D ⚠ *docs say 0.5* |
| EZH2 Palbo mRNA drop | EZH2 mRNA, CDK4/6i vs +Shh | 0.44 | 0.42 | 0.59 | ✓ | **manuscript Fig 4F** ⚠ *measured in MB* |
| GNP+Shh proliferates | ≥1 division | yes | — | yes | ✓ | ext-lit (Wechsler-Reya '99) |
| GNP−Shh arrests | 0 divisions | 0 | — | 0 | ✓ | ext-lit (Wechsler-Reya '99) |
| GNP+HHi arrests | 0 divisions | 0 | — | 0 | ✓ | inferred (RNA-seq collapse) |
| GNP serum-starve arrests | 0 divisions | 0 | — | 0 | ✓ | assum (mitogen=0 by construction) |
| **Period GNP** ⚠ | cell-cycle length | **16.25 h** | 0.30 | 22.83 | ✗ | **ext-lit (Contestabile 2013)** |

- **EZH2i fold (✓\*)**: passes only on the wide 0.42 band; model 1.63 is well short of the 2.2 point value → **genuine shortfall** (JP confirmed). Separate genetic point: Ezh2 cKO 2.46× (Fig 3A) — the genetic value argues this is a real miss, not an old artifact.
- **EZH2 Palbo drop ⚠**: the 0.444× is measured in **MB** (Fig 4F qPCR), but the harness scores it on the **GNP** condition — a cell-type mismatch to resolve.
- **Period GNP ⚠**: JP adopts **Contestabile 16.25 h** (phase-resolved: G1 7.63/S 7.11/G2 1.18/M 0.33). The harness still checks 22.0 h (model 22.83 passes that); under a 16.25 h target the current model *fails* → the harness change is deferred to the coordinated CKI-split recalibration (needs the MB Tc). This is a **two-citation literature conflict** (Nakashima ~22–23 h vs Contestabile 16.25 h), **not** missing provenance as I first wrote.

---

## MB targets

| Target | Quantity | Target | Tol | Model | ✓ | Src |
|---|---|---|---|---|---|---|
| CyclinD1 MB+HHi / MB | Ccnd1, vismodegib vs MB | 0.144 | 0.45 | 0.11 | ✓ | manuscript §A ⚠ *compendium says 0.40* |
| MYCN MB+HHi / MB | Mycn, HHi vs MB | 0.86 | 0.25 | 0.94 | ✓ | manuscript §A |
| MB proliferates | ≥1 division | yes | — | yes | ✓ | manuscript |
| **MB+CDK4/6i (palbo)** ⚠ | pRb⁺ cycling fraction + reversibility | **16–18%**, reversible (washout→~50%) | — | 0-div (proxy) | — | **manuscript** (JP culture) |
| MB+CDK4/6i+EZH2i | EZH2i does *not* restore cycling | 0 (no rescue) | — | 0 | ✓ | **prediction (Fig S9I,J)** |
| MB 2N (G0+G1) count % | 2N population fraction | 68.2 % | 0.30 | 68.7 | ✓ | manuscript §F (**microscopy**) |
| MB S count % | S population fraction | 15.7 % | 0.42 | 21.9 | ✓ | manuscript §F (**microscopy**) |
| MB G2+M duration | direct (pHH3/BrdU) | 2.5 h | 0.55 | 3.39 | ✓ | ext-lit (Fujita) · *soft (~3.6 h)* |
| MB HU S fold | S fraction, HU vs DMSO | 1.36 | 0.40 | 0.95 | ✓ | manuscript §F (HU) |
| MB HU G2 fold | G2 fraction, HU vs DMSO | 0.23 | 0.40 | **0.36** | **✗** | manuscript §F (HU) |
| EZH2 transcript S/G0 | within-cycle EZH2 mRNA | 2.0 | 0.45 | 1.53 | ✓ | manuscript Fig 4A/4B ⚠ *pooled GNP+MB* |
| EZH2 protein G2/G0 | within-cycle EZH2 protein | 1.48 | 0.45 | 1.86 | ✓ | **manuscript Supp 6D** (IF) |
| HU EZH2-in-S boost | EZH2, HU-arrested S vs cycling S | 1.31 | 0.35 | 1.25 | ✓ | manuscript Fig 4G ⚠ *alt-text 1.22* |

- **MB+CDK4/6i (palbo) ⚠ — this is the row you flagged.** Your culture data: pRb⁺ cycling drops to **16–18%** under palbociclib (a *partial* arrest with a real residual population), and it's **reversible** — rebounds to **~50%** pRb⁺ after 72 h + washout. The harness reduces this to a binary "0 div" on a single representative lineage (a coarse proxy); the reversibility is not scored.
- **MB+CDK4/6i+EZH2i "no rescue"** is a **model prediction** (Fig S9I,J: cycling held ~0% with or without EZH2i, because a Vmax kinase block can't be bypassed by more CyclinD1). The manuscript says it *motivated* the experiments. I had mislabeled it as experimental — corrected.
- **MB HU G2 fold (✗)** is the single current miss (0.36 vs band top 0.32) — the perennial structural near-miss.
- **HU is MB, not GNP.** The HU folds, phase fractions, and within-cycle EZH2 gradients are all MB.

---

## MB ↔ GNP cross-ratios

| Target | Quantity | Target | Tol | Model | ✓ | Src |
|---|---|---|---|---|---|---|
| CyclinD1 MB/GNP | Ccnd1 ratio | 5.07 | 0.20 | 4.61 | ✓ | **manuscript Fig 4I** — load-bearing |
| MYCN MB/GNP | Mycn ratio | 2.80 | 0.30 | 2.71 | ✓ | manuscript §A (~2.86) |
| Gli1 MB/GNP | Gli1 ratio | 6.90 | 0.40 | 7.15 | ✓ | manuscript (raw RNA-seq, *no fig#*) |
| EZH2 MB/GNP | EZH2 ratio | 2.05 | 0.35 | 2.03 | ✓ | **manuscript Fig 4J** · *soft* |
| Skp2 MB/GNP | Skp2 ratio | 2.3 | 0.30 | 2.34 | ✓ | **unpub** (scRNA + MB table) |

---

## Non-formal / candidate targets (measured data not currently scored)

The completeness audit found measured phenotypes in the source docs that the spec should carry — some are cheap candidate checks (the condition already runs).

| Feature | Cell | Value | Src | Note |
|---|---|---|---|---|
| **MB vismo→EZH2i rescue** | MB | pRb⁺ 100 → **~25%** (vismo) → **~75%** (EZH2i) | manuscript Fig 5C | your key rescue result; the RESCUABLE contrast to CDK4/6i; model 100/20/74 |
| Ezh2 cKO GNP de-repression | GNP | CyclinD1 **2.46×** | Fig 3A (RNA-seq) | genetic partner to Taze 2.2× |
| EZH2 GNP+HHi drop | GNP | EZH2 mRNA ↓ 25–47% | §A | EZH2 is Hh-responsive; condition runs → candidate check |
| Palbo CyclinD1↑ / E2f1↓ (MB) | MB | CyclinD1 1.57×, E2f1 0.43× | Fig 4F | E2f1 checks E2f-gating; CyclinD1↑ is orthogonal de-repression |
| GNP+Tazemetostat phase folds | GNP | G1 2.03×, G2 1.44× | Fig 3B | EZH2i *promotes* proliferation |
| Ezh2-OE MB phase folds | MB | **G0 2.84×**, G1 0.51× | Fig 3F | EZH2 GOF → arrest (only UP-direction test) |
| Palbo MB phase folds | MB | S 0.30×, G0 1.96× | §F2 | strongest signals; more discriminating than the binary |
| RO3306 (CDK1i) arrest | MB | S 1.58×, G2 1.26× | §F | only arrest *downstream* of the R-point |
| Qualitative directional | mixed | rShh→EZH2↑; E2f1-OE→EZH2↑; Ezh2-OE→CyclinD1↓ | Figs 4H/4D/3E | non-scored but load-bearing |

---

## Distribution features

| Feature | Cell | Value | Status | Src |
|---|---|---|---|---|
| MB phase distribution | MB | G0 24.8 / G1 43.4 / S 15.7 / G2 16.1 % | active (partial) | manuscript §F (**microscopy**, raw ÷0.874) |
| MB HU phase-fold distribution | MB | G0 1.19 / G1 1.16 / S 1.36 / G2 0.23 | active (partial) | manuscript §F (HU) |
| Per-phase SDs (CV ground-truth) | MB | e.g. DMSO G0 21.7 ± 5.51 SD, N=12 | documented | §F1/F2 — ⚠ the "SE" column is **SD** (note #1) |
| CyclinD1 G0/1 CV | GNP | 0.70 | proposed | JP flow — reconcile with model σ=0.68 (same quantity) |
| CyclinD1 CV under EZH2i | GNP | 0.82 | proposed | model *prediction* (EZH2 buffers variance) |
| Transient-G0 fraction | MB | ~12 % | proposed | JP (GNP outer-EGL ≈ 0) |

---

## Cell-cycle timing

Keep GNP and MB **separate** — the model implicitly sets MB/GNP ≈ 0.99, which is the missing feature.

| | Tc | G1 | S | G2 | M | Src | Status |
|---|---|---|---|---|---|---|---|
| **P7 GNP** | 16.25 h | 7.63 | 7.11 | 1.18 | 0.33 | ext-lit (Contestabile 2013) | **adopted** |
| *(prior GNP estimate)* | ~22–23 h | — | — | — | — | ext-lit (Nakashima 2015, soft) | superseded/conflict |
| **MB** | *TBD (longer)* | — | — | — | — | unpub (JP: ↑CKI + transient G0) | **needs target** |
| model (default) | GNP 22.8 ≈ MB 22.8 | — | — | — | — | assum (mu=0.0005) | current |
| model (CKI-split) | GNP 16.1 / MB 22.0 | — | — | — | — | the `k_mu_cki` feature | candidate |

**Mechanism.** `Tc = ln2 / mu_eff`; the (default-off) `k_mu_cki` term makes it CKI-coupled so MB's higher INK4 slows only MB.

**Expression context** (grounds the CyclinD1 reservoir): CCND1 ≈ 11th pct, CCND2 ≈ top 0.3 % (bulk RNA-seq); CCND1 t½ 1.5–2.5 h (RNAdecayCafe, caveated).

---

## Open issues / to document

1. **Period**: reconcile Nakashima 22–23 h vs Contestabile 16.25 h (JP adopts Contestabile); then decide the harness change (coupled to the CKI-split recal + MB Tc).
2. **MB cycle length**: add a formal target (JP investigating).
3. **CyclinD1 MB+HHi**: 0.144 (harness/provenance, load-bearing) vs 0.40 (compendium §A) — reconcile.
4. **EZH2i fold**: 2.2 target vs 1.63 model = real shortfall to close.
5. **`ezh2_palbo_drop`**: measured in MB, scored on GNP — move the check?
6. **CDK4/6i**: encode the 16–18 % residual + reversibility (not a binary); confirm whether experimental +EZH2i data exists or keep the no-rescue as a prediction.
7. **Modality**: correct "flow cytometry" → cell-culture microscopy in the two existing docs.
8. **Pin figure numbers** for Gli1 MB/GNP and Skp2 MB/GNP.
9. **Distribution**: record per-phase SDs to anchor CV targets; carry the SE-vs-SD caveat.

---

## References

- Chahin et al. (manuscript) — Fig 3A/3B/3C/3E/3F, 4A/4B/4D/4F/4G/4H/4I/4J, Supp 6D, Fig 5, §A–F (RNA-seq, qPCR, scRNA, IF, microscopy).
- Contestabile et al. 2013 — P7 GNP phase durations. · Nakashima 2015 — prior GNP period estimate.
- Fujita — G2/M duration. · Wechsler-Reya & Scott 1999 — Shh mitogen. · Cook Sangar 2017 — SHH-MB palbo reversibility.
- Heldt et al. 2018 (PNAS) — cell-cycle engine. · Vock / RNAdecayCafe — mRNA half-lives.

*Related memory: `validation-targets-formal-spec`, `ezh2-paper-model-figure-map`, `distribution-model-thread`, `recalibration-16h-cycle-attempt`, `v44-hu-s-g2-structural-limit`.*
