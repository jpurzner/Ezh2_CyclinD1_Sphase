# Data provenance — what is constrained vs. assumed in the v44 model

An audit trail of every model parameter/feature and where it comes from, so the manuscript can state
cleanly which parts are **data-constrained** and which are **inherited, soft, or free**. Companion to
`PARAMETERIZATION.md` (full numbers + sources) and `Ezh2_CcnD1_model_targets.md` (the target compendium);
scored by `simulations/validate_v44.py`.

**Guiding principle:** species are in **arbitrary units**; the model is calibrated to experimental
**ratios and qualitative behaviors**, never to absolute concentrations.

## Status legend

| Tag | Meaning |
|---|---|
| **MEASURED** | Constrained by the Purzner-lab experimental data (Chahin et al.): RNA-seq, qPCR, scRNA-seq, IF, flow, CKI panel. |
| **INHERITED** | Taken from external published models/data (Heldt 2018 engine; H3K27me3 chromatin literature). Not fit to our data. |
| **SOFT** | A real target, but loosely enforced / not exactly matched (approximate literature value or a known miss). |
| **FREE** | Not directly measured; set by a constrained search or for numerical stability (direction, not value, is what matters). |

---

## 1. The primary datasets (Chahin et al.)

| Data type | What it measures | Source |
|---|---|---|
| Bulk RNA-seq | Gli1, Ezh2, Ccnd1, Mycn across GNP / Ptch⁺ᐟ⁻ / MB ± vismodegib | normalized counts (Purzner lab) |
| qPCR | CyclinD1 after Ezh2 cKO / Tazemetostat / Ezh2-OE / Palbociclib | Fig 3C,E, 4F |
| scRNA-seq | Ezh2 transcript by cell-cycle phase (Sepp, Vladoiu, Ocasio, Luo) | Fig 4A,B |
| Immunofluorescence | Ezh2 protein by phase; under HU / RO3306 / palbociclib / rShh | Fig 4G,H, S6D |
| Flow cytometry | G0/G1/S/G2 proportions in MB ± arrest drugs | Fig 4G / S6, Sec F |
| CDK-inhibitor panel | INK4 (p16/p18) + CIP/KIP (p21/p27) counts, GNP vs MB | DESeq2 |

---

## 2. Data → parameter → status (the audit table)

### 2.1 Hedgehog → Gli cascade

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| `Ptch1_copy_number` (functional Ptch1 fraction) | Gli1 MB/GNP = 6.9× (RNA-seq) — sweep pins f≈0.1 | 1.0 GNP / 0.5 Ptch⁺ᐟ⁻ / 0.1 MB | **MEASURED** |
| Gli→Ptch1→Smo negative feedback (`k_Ptch1_Gli`, `K_Gli_Ptch`) | The 6.9× fold + the GNP Shh-step Gli1 overshoot; "broken loop" in MB | k_Ptch1_basal 0.15, k_Ptch1_Gli 1.87, K_Gli_Ptch 0.40 | **MEASURED** (structure) / FREE (exact rates) |
| Gli arm de-saturation (`k_Cd_tx_Gli_max`, `K_Gli_act_CycD`, `k_Smo_act`, `K_Ptch_Smo`, `K_Smo_Gli_switch`) | Gli1 6.9× and CyclinD1 6.8× must rise GNP→MB without saturating | search-fit | **MEASURED** (targets) / FREE (individual rates) |
| HHi (vismodegib) block on Smo | Gli1 GNP+HHi >99%↓; MB+HHi 98%↓ | HHi 0→1 | **MEASURED** |
| **CyclinD1 is Gli-driven (not MYCN-driven)** | CyclinD1 tracks Gli1 ~1:1 along the Ptch escalation series | structural conclusion | **MEASURED** (a data-derived design choice) |

### 2.2 CyclinD1 convergence node

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| Gli-vs-MYCN split of the CyclinD1 drive | CyclinD1 MB/GNP = 5.07× (Fig 4I); GNP+HHi/GNP = 0.157; MB+HHi/MB = 0.144 | Gli ~86% of MB CyclinD1 | **MEASURED** |
| MB+HHi CyclinD1 sits at the commitment threshold | MB+HHi/MB = 0.144 → ~1× a cycling GNP | — | **MEASURED** (basis of vismo-arrest / EZH2i-rescue) |
| `k_Cd_translation` (CyclinD1 protein scale) | none — set so GNP robustly clears the cycling threshold | 0.75 | **FREE** |

### 2.3 MYCN

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| `MYCN_amplification` | Mycn MB/GNP = 2.86× (RNA-seq) | 1.0 GNP / 2.86 MB | **MEASURED** |
| MYCN→CyclinD1 cooperative Hill (`n_MYCN_Cd`, `k_Cd_tx_MYCN`, `K_MYCN_Cd`) | the small HHi-resistant CyclinD1 residual (MB+HHi/MB) | n≈3.7, k≈35, K≈1.7 | **MEASURED** (residual) / FREE (Hill shape) |
| Gli-dependent vs autonomous MYCN synthesis | Mycn GNP+HHi 0.78; MB+HHi 0.86 | — | **MEASURED** |

### 2.4 EZH2 production & cell-cycle coupling

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| EZH2 transcription gated on E2f × (CycE+CycA) | Ezh2 cycling/G0 ≈ 2.0× (scRNA-seq, 5 datasets) | `kEZE2f`, `K_E2f_EZ`, `wCe` | **MEASURED** |
| EZH2 protein accumulation across S/G2 (long half-life) | Ezh2 G2/G0 = 1.48× (IF) | `kDeEZ` = 0.00015 | **MEASURED** |
| **EZH2 integrates S-phase duration** | Ezh2 in S under HU = 1.31× (IF) | via `kDeEZ` + S/G2 gate + HU→fork | **MEASURED** |
| Mitogen-dose term `Cd/(Kez_cd+Cd)` | rShh dose→Ezh2 ↑ (IF); Ezh2 MB/GNP = 2.05× | `Kez_cd` = 2.94 | **MEASURED** (dose ↑) / **SOFT** (2.05×; model ~1.9×) |
| CDK4/6→Rb-E2F→EZH2 arm | Palbociclib ↓Ezh2 55.6%, ↓E2f1 56.7% (qPCR) | emergent | **MEASURED** |

### 2.5 EZH2 ⊣ CyclinD1 repression strength

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| Repression depth (AUM leaky floor `f0_mk`, `K_mk`) | Ezh2 cKO 2.46× + Tazemetostat 2.2× de-repression (recombination-corrected ~2.8–3× per-cell) → **EZH2i CyclinD1 fold ≈ 2.2** | f0_mk 0.233, K_mk 0.305 | **MEASURED** |
| Direction of repression | Ezh2-OE ↓CyclinD1 (qPCR, Fig 3E) | — | **MEASURED** |

### 2.6 CDK inhibitors / commitment threshold (the rescuable-vs-not contrast)

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| INK4 / CDK4/6 competitive brake `p16`, `p18` | INK4 panel: p16 139× (small absolute), p18 3.0× (dominant); abundance-proportional | p16 0/0.15, p18 0.4/1.5 | **MEASURED** |
| CIP/KIP / CDK2 brake `kSyP21` | p21 3.1× + p27 1.4× (CIP/KIP panel) | 0.002 GNP / 0.004 MB | **MEASURED** |
| p27→CDK4/6 share `w_p27` | Guiley 2019 (p27 brakes CDK4/6) + G0-weighting | 1 | **INHERITED** (mechanism) |
| Arrest = CyclinD1 drop crossing a **static** brake (not CKI induction) | vismo barely moves the CKIs (RNA-seq) | — | **MEASURED** |
| CDK4/6i (palbociclib) = non-competitive Vmax block, **not rescuable** | the clean contrast to the competitive INK4 brake | `kPhRbCd`=0 | **MEASURED** (behavior) |

### 2.7 Cell-cycle phase timing

| Model parameter / feature | Constrained by (data) | Value | Status |
|---|---|---|---|
| `kSyDna` (fork speed → S-phase length) | S ≈ 3 h (cumulative-BrdU); DMSO S ≈ 15.7% | 0.044 (~5× Heldt default) | **MEASURED** |
| Growth / commitment balance (`M_commit`) | DMSO G0+G1 = 68.2% (resolvable 2N pool) | 1.05 | **MEASURED** (2N pool) / SOFT (G0/G1 split) |
| G2+M duration | direct methods → G2+M ≈ 2.5 h | — | **SOFT** (model ~3.6 h — a standing miss) |
| HU → fork-speed coupling (`KmHU_fork`, `vmin_fork`) | HU: S ×1.36, G2 ×0.23 (flow) | — | **SOFT** (needs re-tune for the shorter S — standing miss) |

### 2.8 Cell-cycle engine (the Heldt core)

| Model component | Source | Status |
|---|---|---|
| Rb–E2F bistable restriction point; cyclin/CDK network; Skp2–p27 commitment; explicit DNA replication; CyclinB/CDK1 mitotic switch; APC/C + division reset; growth/size gate | **Heldt, Barr, Cooper, Bakal & Novák 2018 (PNAS)**, BIOMD0000000700 | **INHERITED** (literature parameters; not fit to our data) |
| Cell-cycle period ~22–23 h | Nakashima 2015 (published estimate for GNPs) | **SOFT** (sets timescale only) |

### 2.9 H3K27me3 AUM module (the mechanistic epigenetic layer)

The mark **dynamics** are anchored to the external chromatin literature; only the mark's steady-state
**repression depth** ties back to our CyclinD1 data. Full sourcing in `H3K27me3_CyclinD1_literature_review.md`.

| Model parameter / feature | Source | Status |
|---|---|---|
| Read-write methylation `k_w_mk` + de-novo floor `k0_mk` | Berry–Howard 2017 A/U/M; Margueron 2009 (EED read-write) | **INHERITED** |
| Turnover `del_mk` (t½ ≈ 16 h de-repression; τ_restore 10–18 h) | Reverón-Gómez 2018, Escobar 2023, Knutson 2012 (EZH2i mark t½ ~1 day) | **INHERITED** |
| Replicative ½ dilution at S; T_cc dependence | Jadhav 2020 (replicational dilution, proliferation-rate-dependent) | **INHERITED** |
| Transcription→PRC2 eviction arm `g_mk`, `K_tx_mk` | Kaneko/Davidovich/Beltran (nascent-RNA evicts PRC2); strength kept weak | **INHERITED** (mechanism) / **FREE** (strength, weak regime) |
| Leaky Hill readout depth `f0_mk`, `K_mk`, `n_mk` | our CyclinD1-fold data (see §2.5) | **MEASURED** |
| EZH2-scaled catalysis; Ser5P Pol II present (leaky floor) | Ser5P Pol II ChIP (our data); Blackledge & Klose 2021 | **MEASURED** + **INHERITED** |

---

## 3. Not data-driven (assumed / free)

- **`k_Cd_translation`** (CyclinD1 protein scale) — numerical, keeps GNP off the bistable knife-edge.
- **CyclinD1 cell-to-cell heterogeneity σ = 0.68** — sets the *magnitude* of the population rescue; only the *direction* is data-robust.
- **`M_commit`** (G0/G1 split) — interim, pending a mechanistic G0/G1 timing element.
- **Eviction-arm strength `g_mk`** — constrained only to preserve the 5.07× fold and stay in the weak regime (per the measured fold; not independently measured).
- The **entire Heldt engine** and the **H3K27me3 kinetic constants** — literature, not our data.

---

## 4. Standing tensions (measured targets not fully matched)

The 5 validation misses (both the default and the AUM model, 22/27):

1. **EZH2 transcript S/G0 gradient** — model too steep (worst miss).
2. **MB HU S-fold & G2-fold** — HU fork-coupling needs re-tune for the shortened S-phase.
3. **MYCN GNP+HHi** — cascade target (0.91 vs 0.78).
4. **MB G2+M duration** — 3.6 h vs 2.5 h (accepted G2 sacrifice from the bake).
5. **EZH2 MB/GNP 2.05×** — soft; model ~1.9× (cell-cycle-coupled).

---

## 5. One-line summary for the manuscript

> The signalling→CyclinD1 cascade (Hedgehog/Ptch1/Gli, MYCN, EZH2 induction and EZH2⊣CyclinD1
> repression) and the CDK-inhibitor commitment threshold are **calibrated to Chahin et al.
> measurements** (bulk/scRNA-seq, qPCR, IF, flow, CKI panel — ratios, not absolute levels). The
> real-time **cell-cycle engine is adopted from Heldt et al. 2018**, and the **H3K27me3 replicative-
> dilution kinetics from the chromatin literature** (Berry–Howard, Jadhav, Reverón-Gómez); only the
> mark's steady-state repression depth is tied to our CyclinD1 data. A small number of scale/heterogeneity
> parameters are set for numerical robustness, with the model's conclusions depending on their direction,
> not their exact values.
