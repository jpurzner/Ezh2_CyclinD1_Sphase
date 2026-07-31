# Origins & Foundations

> **Scope of this page.** How the model was *first assembled* — the exploratory era that ran from a blank slate up through v42, almost entirely inside the **`ezh2_cyclind1_sym`** repository (GitHub: **`jpurzner/Ezh2_Hh_Ccnd1_model`**, "the old foundation"). This covers the systematic review of published BioModels cell-cycle models, the minimal-oscillator search, the Gli/CyclinD1-EZH2 module design, and the long v2→v42 version climb (with its many dead-ends). Where the mature model went from here is covered in [Model Architecture](02_model_architecture.md), [Evolution Timeline](03_evolution_timeline.md), and [Things Tried & Abandoned](05_things_tried_and_abandoned.md).

---

## 1. The starting goal

The project began with a concrete, deceptively simple biological target defined in `EZH2_CyclinD1_Model_Integration_Spec.md`:

Build one ODE model in which

1. the **EZH2 ⊣ CyclinD1 negative feedback loop**,
2. **Hedgehog (SHH → Gli → CyclinD1)** mitogenic signaling, and
3. the **standard mammalian cell-cycle machinery** (G0/G1/S/G2/M),

are coupled through the **Rb–E2F1 pathway**, such that the same model reproduces two regimes:

- **WITH Hedgehog:** sustained oscillations (≈10–15 divisions), stable Gli, ~24 h period.
- **WITHOUT Hedgehog:** G0 arrest (0 divisions) with p27 accumulation.

The scientific context (from the integration spec) is cerebellar granule neuron precursors (GNPs) and SHH-medulloblastoma (SHH-MB), where EZH2 and CyclinD1 form a feedback loop coordinating cell-cycle exit with differentiation. The spec was written for James & Teresa Purzner's differentiation-therapy program and already enumerated the experimental anchors that later became validation targets (Ezh2 cKO → 2.47× CyclinD1 log2FC 1.301; EZH2-OE → ~50% CyclinD1 drop, 2.7× G0; CDK4/6i → EZH2 −55.6%, CyclinD1 +57%; Vismodegib+Tazemetostat rescue). See [Data & Evidence](07_data_and_evidence.md) and [Calibration & Validation](06_calibration_and_validation.md).

The integration spec was explicitly a "buy, don't build" plan: **start from published SBML models** and graft the novel modules on top.

---

## 2. Published models evaluated & integrated (and why each)

Two waves of published-model evaluation happened. The **first** (the integration spec, `EZH2_CyclinD1_Model_Integration_Spec.md`) named three models to download and stitch together; the **second** (the systematic `BIOMODELS_REVIEW_REPORT.md`, `analyze_biomodels.py`, `deep_analysis.py`, `model_comparison_table.csv`, Jan 2026) surveyed five mammalian cell-cycle models to diagnose *why the home-grown oscillator kept failing*.

### 2.1 The original three (integration spec)

| Model | BioModels ID (as cited in spec) | File in repo | Why it was chosen |
|---|---|---|---|
| **Gérard & Goldbeter 2009** — mammalian cyclin/Cdk network | BIOMD0000000109 *(spec)* / BIOMD0000000730 *(review)* | `Gerard2009.xml`, `Gerard2010.xml`, `Gerard2009_only_walkthrough.ipynb` | Core cell-cycle engine: CyclinE/A/B–CDK complexes, all G1/S/G2/M transitions. Known gap: no CyclinD1, no mitogen input. |
| **Swat et al. 2004** — Rb-E2F pathway | BIOMD0000000189 | `swat2004_rb_e2f.xml` | Detailed CyclinD/CDK4/6 → pRb → E2F1 dynamics; intended as the **base** model (has exactly the D-cyclin/Rb/E2F species the EZH2 loop needs). |
| **Schoeberl et al. 2002** — EGF/MAPK signaling | BIOMD0000000048 | `schoeberl2002_egf_mapk.xml` | Template for the *signal → TF → gene-expression* logic, to be adapted for SHH→Gli→CyclinD1. |

The spec's stated build order (its "Notes for Claude CLI"): start from **Swat 2004** (it has CyclinD/CDK4/Rb/E2F1), bolt on **Gérard 2009** for S/G2/M cyclins, then add novel SHH and EZH2 species.

### 2.2 The systematic BioModels review (five models)

`BIOMODELS_REVIEW_REPORT.md` analyzed and scored five models against three requirements — sustained oscillations, proper G0 arrest, and growth-factor dependence at the restriction point. Local SBML copies were kept for all of them (`Gerard2009.xml`, `Novak2004.xml`/`.cps`, `Tyson2001.xml`, `Yao2008.xml`, `Schwarz2018.xml`, plus `yao2008_rb_e2f.xml`, `Abroudi2017.xml`).

| Model | ID | What it contributed to the design |
|---|---|---|
| **Gérard & Goldbeter 2009** | BIOMD0000000730 | The comprehensive 39-variable engine; the eventual answer for "how to get sustained GF-dependent oscillations." Multiple positive-feedback loops (E2F autoactivation, CycE–E2F, CycA, MPF autocatalysis) + APC/Cdc20 & APC/Cdh1 negative feedback + p27. **This is the engine that ultimately won** (see §5, v39+). |
| **Novák & Tyson 2004** | MODEL2006080001 | The **p27/Kip1 restriction-point machinery** — Skp2-mediated p27 degradation, CycE–CDK2 sequestration, mass/µ growth-factor-dependent CycD synthesis. This is where the G0-lock mechanism came from (adopted at v26). |
| **Tyson & Novák 2001** | BIOMD0000000195 | Principles of bistability from mutual antagonism, hysteresis, irreversible transitions; APC/Cdh1 reset for oscillations. Conceptual backbone for the switch design. |
| **Yao et al. 2008** | BIOMD0000000318 | The bistable **Rb–E2F switch** underlying the restriction point (all-or-none E2F, positive feedback E2F→CycE→Rb~P→E2F); justification for Hill n≈4 on E2F. |
| **Schwarz et al. 2018** | BIOMD0000000918 | The quantitative **CDK-activity threshold = 0.84** for R-point passage (95% predictive). Used to inform restriction-point threshold behavior. |

The review's "critical missing mechanisms" verdict (which drove all subsequent versions): the home-grown oscillator lacked **dual positive feedback** (E2F–CycE + Rb–E2F), **APC/Cdh1 negative feedback**, a **CDK threshold**, **CycD degradation when GF absent**, and **p27 accumulation in G0**. `SOURCES_AND_CITATIONS.md` holds the full bibliography and BibTeX for these plus supporting quiescence/APC/CDK-inhibitor literature.

> **Provenance note / discrepancy worth flagging for the paper:** the integration spec cites Gérard 2009 as **BIOMD0000000109** while the systematic review cites **BIOMD0000000730** (and the README cites Gérard & Goldbeter 2009 *Mol Syst Biol* 5:324 vs *PNAS* 106:21643 — these are two different Gérard-Goldbeter 2009 papers). The repo also keeps a `Gerard2010.xml`. Anyone reconstructing the exact SBML lineage should verify which accession each downloaded file actually is.

Model attribution as it was intended to appear in the paper is drafted in `MODEL_DOCUMENTATION_FOR_PAPER.md` (Methods text + a parameter-source table mapping each rate constant to Novak2004 / Gerard2009 / "Round 2 optimization" / "This work").

---

## 3. The minimal-oscillator search

Because the full Gérard/Novak engine (39 variables, multiple bistable switches) proved "difficult to tune… parameters affect WITH/WITHOUT HH equally," a parallel line of work (documented in `MINIMAL_OSCILLATOR_OPTIONS.md`) asked whether a *smaller* oscillator could give clean HH-dependent switching. Four candidates were compared:

| Option | Variables | Source | Coupling idea | Verdict |
|---|---|---|---|---|
| **Tyson 1991** minimal | 2 (u, v) | BIOMD0000000006 | Modulate cyclin synthesis κ via CycD/p27 | Simplest; autocatalytic u² switch; no phase detail. Fallback recommendation. |
| **Goldbeter 1991** cascade | 3 (C, M, X) | classic mitotic oscillator | Modulate synthesis vᵢ via CycD/p27 | *Recommended first* (explicit negative feedback, MM kinetics). **Failed on integration** — see v31. |
| **Custom CycD-driven** | 4 (CycB, APC, Skp2, p27) | novel | Direct CycD gating built in | "Too simple for reviewers" risk. |
| **Gérard skeleton** | 5 | Gérard et al. 2012 | — | Still has multiple bistable switches → same tuning problems as full model. |

This search produced v31 (Goldbeter) and v32 (Tyson) builders. The Goldbeter integration **failed with CVODE stiffness** (see §5); the whole minimal-oscillator branch was eventually abandoned in favor of returning to the full Gérard-Goldbeter engine at v39. `LIVE_IMAGING_INSIGHTS.md` (Spencer/Cappell/Schwarz single-cell CDK2 dynamics) fed the reasoning about what a credible oscillator/commitment readout needed to look like.

---

## 4. The novel modules: Gli and EZH2–CyclinD1

Three modules were original to this work (per `MODEL_DOCUMENTATION_FOR_PAPER.md` §3):

**(a) Hedgehog–Gli signaling** as the growth-factor input, replacing the generic "Mass/GF" term from the published models. The Gli module went through a major biological correction (see `GLI_MODULE_WIREFRAME.md`, `GLI_MODULE_VISUAL_DIAGRAM.txt`, `V28_ARCHITECTURE_SUMMARY.md`):

- **Early (v26–v27) — wrong:** Gli3 → GliR only becomes a repressor after a special "processing" event; `Gli_active` always an activator; in v26 Gli1 *directly* drove CyclinD1.
- **Corrected (v28) — right:** full-length **Gli2/Gli3 are repressors by default** (bound to Gli1/CyclinD1/Nmyc promoters when HH is off); the pathway *converts* repressors → activators via Smo; **Gli1 is an obligate activator and a readout/amplifier, not a driver.** Target genes follow `(basal + Gli_act activation) × (Gli_rep repression)`. A basal CyclinD1 term was added so progenitors divide slowly before HH exposure. Gli biology cited to Hui & Angers 2011, Ruiz i Altaba 2007, Sasaki 1997, Rohatgi & Scott 2007.

**(b) EZH2–bivalent-chromatin regulation of CyclinD1** — E2F-dependent EZH2 deposits H3K27me3 (repressive, slow turnover ~35 h half-life at `k_demethyl_H3K27me3 = 0.02`); Gli recruits H3K4me3 (activating, fast ~1.7 h at 0.4); CyclinD1 transcription depends on the H3K4me3/H3K27me3 ratio, giving epigenetic **memory** of the G0 state. This is the earliest form of what later became the CyclinD1-mRNA-reservoir "Hh memory" story in v44 (see [Mechanisms Explored](04_mechanisms_explored.md)).

**(c) E2F-gated Gli → CycE/CycA transcription** — direct mitogenic input to S-phase cyclins, gated by E2F to preserve restriction-point control.

The CyclinD1 mRNA node became the single integration point where all three inputs (Gli, MYCN, EZH2) converge — the form shown in the README:
```
d[Cd_mRNA]/dt = [ k_basal + k_Gli·(Gli)²/(K²+…)·rep + k_MYCN·MYCN³/(K³+…) ]
                 · K_EZH2/(K_EZH2 + EZH2·(1−EZH2i)) − k_deg·Cd_mRNA
```

A tiny standalone `CyclinD_EZH2_module.xml` and Novak-exploration notebooks (`Novak_Ezh2.ipynb`, `novak_*` sim outputs) are the fossils of these early module experiments.

---

## 5. The v2 → v42 version climb

Every builder from `build_model_v2.py` to `build_model_v42_mycn.py` survives in `ezh2_cyclind1_sym/src/`. This is the single most complete record of "everything tried." The arc:

| Version(s) | Builder file(s) | What it did / what happened |
|---|---|---|
| **v2–v8** | `build_model_v2..v8.py` | Earliest scaffolds. |
| **v9–v10** | `v9`, `v10` | Added Gli1/Ptch1 feedback; stabilized Gli during proliferation. |
| **v11** | `v11` | Fast HH kinetics (6 h Gli response). |
| **v12–v13** | `v12`, `v13` | Cell-cycle speed tuning. |
| **v14** | `v14_bivalent` | **Bivalent chromatin (H3K27me3/H3K4me3) introduced.** |
| **v15** | `v15_stable_gli1` | Gli1 autoregulation = 0.1 → Gli1 explodes (5.5). |
| **v16** | `v16_no_gli1_auto` | Gli1 autoregulation = 0 → signal too weak (1 div). |
| **v17** ⭐ | `v17_balanced_gli1` | **First real milestone.** Gli1 coefficient **0.03** = "perfect balance": 0 div without HH, 23.7 h period with HH, 50× Gli fold-change, bivalent chromatin memory working. Documented in `v17_SUCCESS_SUMMARY.md`. |
| **v18–v20** | `v18_sustained`, `v19_robust`, `v20_positive_feedbacks` | Pushing for sustained multi-division oscillations + more positive feedback. |
| **v21–v25** | `v21_strong_gli1` (0.18, autonomous osc), `v22_multiplicative_gli1`, `v23_gated_gli1`, `v24_no_direct_gli_cyclin`, `v25_gated_gli_cyclin` | Series of Gli1-coupling topologies, all wrestling the same tension: strong enough for signal vs. not self-sustaining. |
| **v26** | `v26_p27_system` | **Added Novák-Tyson p27/Kip1 G0-lock.** But Gli1 directly drove CyclinD1 → **autonomous oscillations, HH had no effect** (Gli1 +1283% in both conditions). Documented in `MODEL_DEVELOPMENT_SUMMARY.md`. |
| **v27** | `v27_gli_architecture`, `v27_1_balanced` | CyclinD1 now driven by `Gli_active` (not Gli1) + Gli3→GliR repressor. **Failed:** GliR dominated in *both* conditions; Ptch1 over-accumulated and suppressed Smo. Triggered a 16-parameter automated search (`parameter_search_v27.py`) that couldn't find a working point. |
| **v28** | `v28_correct_gli`, `v28_1_tuned` | **Biological correction: Gli2/Gli3 are default repressors.** Architecturally right, but a bug (SHH-Ptch1 complex not degraded) meant Ptch1 didn't differentiate (0.99×). `V28_ARCHITECTURE_SUMMARY.md`. |
| **v29** ⭐ | `v29_fixed_ptch1` | **Fixed SHH-Ptch1 degradation + made SHH a boundary species + 18-param differential-evolution search.** Hedgehog pathway now *textbook-perfect*: Ptch1 0.32×, Smo 3.68×, Gli_act 3.20×, Gli_rep 0.38×, Gli1 5.92×, CycD 7.93×. **But cells stuck in S-phase** — CycB couldn't accumulate (Cdh1 too active). `V29_FINAL_STATUS.md`, `V29_PROGRESS_SUMMARY.md`. |
| **v30** | `v30_cycD_dependent_p27` | Added CycD-dependent p27 control (Skp2 = f(E2F, CycD-CDK4²); dual p27 degradation). Achieved **differential *initiation*: 1 division WITH HH vs 0 WITHOUT.** A dedicated 12-param oscillation search then **failed** — it regained cycling in *both* conditions, revealing that "cell-cycle oscillator params are globally coupled" and can't simultaneously give easy cycling WITH HH and hard arrest WITHOUT. `V30_FINAL_RESULTS.md`, `V30_OSCILLATION_SEARCH_RESULTS.md`. |
| **v31** | `v31_goldbeter_oscillator` | Grafted the Goldbeter 3-variable minimal oscillator. **Standalone worked (9 peaks/200 h); integrated it died with `CVODE CV_TOO_MUCH_WORK` at t≈0.74 h** — Michaelis-Menten stiffness (K≈0.01) clashing with slow HH dynamics. `V31_GOLDBETER_INTEGRATION_ISSUES.md`. |
| **v32–v38** | `v32_tyson_oscillator`, `v33_G1_timer`, `v34_v30_plus_timer`, `v35_corrected_chromatin`/`v35_minimal_fix`, `v36_E2F_oscillator`, `v37_modular`/`v37_tuned`, `v38_with_cdc20`/`v38_cdc20_v2` | The recovery/experimentation band: Tyson 2-var oscillator, explicit G1 timers, chromatin fixes, an E2F-based oscillator, a modular rewrite, Cdc20/APC work. |
| **v39** ⭐ | `v39_gerard_hybrid` | **The decisive pivot: abandon the minimal oscillators and adopt the full Gérard-Goldbeter 2009 engine** as the cell-cycle core. This is the architecture that carries through to v42/v44. |
| **v40** | `v40_with_ezh2`, `v40_ezh2_from_v39` | Re-attached the EZH2 module onto the Gérard hybrid. |
| **v41** | `v41_data_constrained` | First systematic constraint to the RNA-seq/functional data. |
| **v42** ⭐ | `v42_mycn` | **Added the MYCN module** (basal synthesis, amplified 2.8× in MB, small Gli term, Hill n=3). **The model submitted with Chahin et al.** — 16/18 validation targets, EZH2i-rescues-HHi-but-not-CDK4/6i prediction. |

The through-line: an **early home-grown chromatin+Gli+p27 oscillator (v9–v38)** that repeatedly hit the same wall — you cannot get sustained HH-dependent oscillations *and* clean G0 arrest from a hand-tuned or minimal oscillator — was **replaced at v39 by the published Gérard-Goldbeter engine**, onto which the novel Gli, EZH2, and MYCN modules were grafted. v17, v29, and v42 are the three "it works" checkpoints along the way.

---

## 6. How this fed into v42 → v44 (the handoff)

The `ezh2_cyclind1_sym` repo's own `README.md` now describes itself as "the v42 model accompanying Chahin et al." — i.e. the old repo *is* where v2–v42 were built, and its README was rewritten to point at v42 as the deliverable.

The forward plan is captured in two places:

- **`checkpoint_model_plan.md`** (in `ezh2_cyclind1_scheckpoint`, May 2026) — the "AMOS" plan to fork `Ezh2_Hh_Ccnd1_model`, freeze the published v42, and build an **intra-S checkpoint** extension (ATR–CHK1–CDC25/Pe axis, Iwamoto/Tyson-Novak family) to test the hypothesis *longer S-phase → more EZH2 → stronger CyclinD1 repression → longer G0/G1*, calibrated to the hydroxyurea EZH2-in-S data (1.31× boost). It specifies `build_model_v43_checkpoint.py` as the next builder and lists concrete phases and predicted experiments (aphidicolin, etoposide).
- The active repo **`Ezh2_CyclinD1_Sphase`** is where v43/v44 and the dynamic-CDKI work actually happened — `src/build_model_v44_heldt.py`, the CyclinD1-mRNA-reservoir memory, the two-cyclin (D1+D2) model, MYCN-elevated-expression, etc.

So the lineage is: **published BioModels (Gérard-Goldbeter 2009, Novák-Tyson 2004, Yao 2008, Tyson 2001, Schwarz 2018, Swat 2004, Schoeberl 2002) → v2–v38 exploratory oscillator → v39 Gérard hybrid → v40 EZH2 → v41 data-constrained → v42 MYCN (paper) → v43 checkpoint plan → v44 (active).** The old `sym` repo is the archive of the assembly; everything from v42 onward lives in the [Sphase repo](09_repos_directories_reproduction.md).

See also: [Evolution Timeline](03_evolution_timeline.md) for the v42→v44 continuation, [Mechanisms Explored](04_mechanisms_explored.md) for the EZH2/chromatin-memory lineage, [Things Tried & Abandoned](05_things_tried_and_abandoned.md) for the minimal-oscillator and Gli3→GliR dead-ends in detail, and [Repos, Directories & Reproduction](09_repos_directories_reproduction.md) for exact file locations.

---

### Key source files for this page (all under `ezh2_cyclind1_sym/` unless noted)

- Design/spec: `EZH2_CyclinD1_Model_Integration_Spec.md`, `MINIMAL_OSCILLATOR_OPTIONS.md`, `GLI_MODULE_WIREFRAME.md`, `GLI_MODULE_VISUAL_DIAGRAM.txt`, `LIVE_IMAGING_INSIGHTS.md`, `IMPLEMENTATION_GUIDE.md`
- Review/attribution: `BIOMODELS_REVIEW_REPORT.md`, `SOURCES_AND_CITATIONS.md`, `MODEL_DOCUMENTATION_FOR_PAPER.md`, `MODEL_DEVELOPMENT_SUMMARY.md`, `analyze_biomodels.py`, `deep_analysis.py`, `model_comparison_table.csv`
- Version status: `v17_SUCCESS_SUMMARY.md`, `V28_ARCHITECTURE_SUMMARY.md`/`.txt`, `V29_FINAL_STATUS.md`, `V29_PROGRESS_SUMMARY.md`, `V30_FINAL_RESULTS.md`, `V30_OSCILLATION_SEARCH_RESULTS.md`, `V31_GOLDBETER_INTEGRATION_ISSUES.md`
- Builders: `src/build_model_v2.py` … `src/build_model_v42_mycn.py` (complete)
- Published SBML: `Gerard2009.xml`/`Gerard2010.xml`, `Novak2004.xml`/`.cps`, `Tyson2001.xml`, `Yao2008.xml`/`yao2008_rb_e2f.xml`, `swat2004_rb_e2f.xml`, `Schwarz2018.xml`, `schoeberl2002_egf_mapk.xml`, `Abroudi2017.xml`, `CyclinD_EZH2_module.xml`
- Forward plan: `ezh2_cyclind1_scheckpoint/checkpoint_model_plan.md`
