# Model figures — which code generates each

## Current status (2026-10-16)

**Baked default** (`build_model_v44()`, `src/build_model_v44_heldt.py`): two-step Rb · EZH2 **dilution + APC/Cdh1 S-phase boost** (the concentration-halving convention was REVERSED, see `docs/cdh1_ezh2_sphase_boost_2026-07-19.md`) · serial me1/me2/me3 H3K27 chain · **two-compartment PRC2** (distal reservoir Mk + proximal occupancy P_prox) · dynamic CDKI (INK4 p18/p19; CIP/KIP p27=P21/p21/p57) with option-B p27 · CDK6–Gli second drive arm · CDK4/6-alone hyper-P **escape** (reversible MB transient G0) · Skp2 mitogen-dose · p21-PIP-degron · GLI1-autoregulation memory · **measured CyclinD1 mRNA half-life (3.4 h; `measured_ccnd1=True`, 2026-10-16)**. **Validation 30/32** (two structural misses, both pre-existing: CyclinD1 GNP+HHi/GNP 0.078 vs 0.157; MB HU-G2 0.456 vs 0.23).

- **CyclinD1 memory is NOT a slow transcript** (2026-10-16): the measured CCND1 mRNA t½ is ~3.4 h (stability-catalog, 135 measurements; 0.5 h cytoplasmic), not the old 11 h reservoir. The ~1-cycle Hh-withdrawal coast comes from the **large CyclinD1 reservoir draining past the commitment threshold** (GNP ~5×, MB ~21× above `K_CdRb`; coast ≈ τ·ln(ratio) ≈ GNP 14 h / MB 17 h). The old 11 h reservoir OVERSHOT to ~3 cycles. See `sim_ccnd1_reservoir_test.py` + `fig_ccnd1_halflife.py`. ⇒ the CyclinD1-memory narrative in `fig_cyclind1_memory` / `fig_pool_vs_halflife` is now **stale (needs caption rework)** even though they re-run.
- **EZH2 protein t½ baked ~42 h** (`kDeEZ`=2.77e-4) + Cdh1-gated G0 clearance; stability-catalog grounds EZH2 at ~20–30 h (dividing NP) / ~57 h (pure degradation). **Older entries citing "t½ ≈ 77 h" / `kDeEZ`=1.5e-4 (e.g. `fig_v44_g0_mark_*`, `sim_ezh2_kdeez_sweep`, `sim_mark_kdeez_ramp_sweep`) are SUPERSEDED.**

> The dated changelog blocks below are **historical** (model-evolution provenance, kept intentionally); the status above is current. Superseded scores: 26/28 → 28/28 → 31/32 → **30/32 (current)**.

## Current-work figures (2026-10) — molecular noise, bursting, stability grounding

Noise/decomposition figures write to the session scratchpad (ephemeral) — rerun the `sim_*` first. Grounded on the AMOS **stability** arm (`stability-catalog` repo; `import stability as S`) and EZH2 IF distributions.

| Output | Generating script | What it shows |
|---|---|---|
| `fig_ezh2_single_trace.png` | `sim_ezh2_single_trace.py` | One cell: EZH2 **mRNA (fast/bursty, τ≈1.7 h) vs protein (stable/slow-drift, τ≈60 h)** — the long protein lifetime low-pass-filters the telegraph (temporal CV 0.45→0.26 = the "more protein than RNA" buffering). |
| `fig_ezh2_bursting.png` | `sim_ezh2_bursting.py` | Transcriptional-bursting EZH2 noise (telegraph on E2f-driven tx): DMSO abundance CV **0.41** vs IF ~0.45; drug directions preserved (CDK4/6i ×0.63, HU ×1.56). Calibrated via `sim_ezh2_calib.py` (duty 0.62). |
| `fig_ezh2_tau_sweep.png` | `sim_ezh2_tau_sweep.py` | Burst **correlation time τ** = the over-time-variance lever: τ 3→20→40 h raises within-cell temporal CV 0.33→0.77 and deep-low (G0-prone) dwell 8%→43%, at the same marginal. |
| `fig_ezh2_variance_decomp.png` | `sim_ezh2_variance_decomp.py` | Variance split: baked ≈ 18% cell-to-cell / 35% cell-cycle / 47% bursting; a no-bursting control isolates the cell cycle. **THE GUESS**: baked & "mostly-frozen" reach the SAME snapshot CV with opposite temporal fractions — the IF snapshot cannot resolve the split. |
| `fig_ccnd1_reservoir_test.png` | `sim_ccnd1_reservoir_test.py` | Measured CcnD1 t½ (3.4 h) + the large reservoir ⇒ ~1-cycle Hh-withdrawal coast (GNP 14 h / MB 17 h); slow 11 h overshoots to 44 h; instant collapses in 3–4 h. Grounds `measured_ccnd1`. |
| `fig_ccnd1_halflife.pdf` | `fig_ccnd1_halflife.py` | CCND1 mRNA half-life vs the model — **UPDATED**: the model now sits AT the measured 3.4 h (thesis inverted from "reservoir too slow"). RNAdecayCafe + stability-catalog. |

> **Known gaps:** `fig_param_heatmap/_settings/_biplot` depend on precomputed `allparams.json`/`sensitivity.json`/`param_meta.json` (no in-repo generator found) — stale, not regenerable here. `fig_ezh2_bursting`/`_variance_decomp` need their scratchpad JSONs regenerated before plotting.

---

> **Model-evolution history** (single→two-step Rb re-opt, population-G0 metric fixes, 2026-08 transient-G0 mechanisms, EZH2/mark recalibrations, figure curation) lives in the dated `docs/*.md` notes and git history — removed here to keep the registry current.

All scripts live in `simulations/` and run against the calibrated **v44 default** (two-step Rb + EZH2
dilution/Cdh1 — concentration convention reversed; measured CyclinD1 mRNA half-life; see Current status;
`src/build_model_v44_heldt.py`). Each writes a `.png` (preview) and a `.pdf` (vector, for the paper)
into `simulations/`. Regenerate any single figure with, e.g.:

```bash
./venv/bin/python simulations/fig_v44_fig5_rescue.py
```

**After any model change, regenerate the whole set with the manifest** (auto-globs every
`fig_v44_*.py` and explicitly runs the slow ensemble figures that don't match that glob, e.g. the
cyclin D1 / birth-p27 bifurcation):

```bash
bash simulations/regenerate_figures.sh
```

The unified validation harness `simulations/validate_v44.py` defines the conditions, the division
counter, the phase (p27) classifier, and the experimental-target checks that the figure scripts import.

---

## Paper figures (model)

| Manuscript panel(s) | Output PDF | Generating script | Notes |
|---|---|---|---|
| **Fig. 5A–C** — Ptch+/− MB / +HHi (arrest) / +HHi+Ezh2i (rescue), Cyclin B traces | `fig_v44_fig5_rescue.pdf` | `fig_v44_fig5_rescue.py` | |
| **Supp. Fig. 7A** — model architecture / wiring diagram | `fig_v44_wiring.pdf` | `fig_v44_wiring.py` | pure schematic (no simulation) |
| **push/pull summary schematic** — what drives EZH2 & CyclinD1 up/down + the H3K27me3 feedback | `fig_v44_pushpull_schematic.svg` | hand-authored SVG (no simulation) | Balance-sheet capstone for the EZH2↔CyclinD1 thread: green = drives up, red = pushes down; EZH2→H3K27me3⊣CyclinD1 leaky-rheostat feedback; three emergent-behavior callouts (progressive repression/buffering; proliferation-rate sets the brake via dilution ∝1/period; **G0 reversible/monostable, not a bistable lock**). Vector; converts to PDF for the manuscript (rsvg-convert/Inkscape). |
| **Supp. Fig. 7B–M** — validation grid (Cyclin B traces; CyclinD1/Mycn/Gli1 model‑vs‑data; EZH2; division counts; time courses) | `fig_v44_validation.pdf` | `fig_v44_validation.py` | |
| **Supp. Fig. 8A–H** — EZH2→CyclinD1 feedback (±feedback: Cyclin B, CyclinD1 mRNA, EZH2 protein, phase relationship) | `fig_v44_feedback.pdf` | `fig_v44_feedback.py` | feedback OFF = `f0_mk = 1.0` (AUM mark cannot repress) |
| **Supp. Fig. 8I,J** — MB+CDK4/6i and MB+CDK4/6i+Ezh2i (no rescue) | `fig_v44_cdk46i_norescue.pdf` | `fig_v44_fig5_rescue.py` | emitted alongside Fig. 5 |
| **Supp. Fig. 7x** — Gli→Ptch1 negative feedback & its disruption in MB (Gli1 step overshoot GNP vs MB; Ptch1 mRNA induction; Smo/Gli broken‑loop) | `fig_v44_ptch1_feedback.pdf` | `fig_v44_ptch1_feedback.py` | Ptch1 = Gli target; `Ptch1_copy_number` = functional fraction (1.0 GNP / 0.1 MB) |
| **Fig. 5C** — population cycling fraction (vismo↓, Ezh2i restores; CDK4/6i no rescue) | `fig_v44_fig5_population.pdf` | `fig_v44_fig5_population.py` | N=120 ensemble, CyclinD1/p27 heterogeneity |
| **Conditions → input parameters + parameter space** (A: input-knob matrix per condition; B: CyclinD1×p16 plane with cycle/arrest boundary + each condition placed) | `fig_v44_condition_paramspace.pdf` | `fig_v44_condition_paramspace.py` | shows MB+HHi arresting at the *same* CyclinD1 a GNP cycles at, because p16 raises the threshold; EZH2i crosses it back |

Drop‑in replacement Results/Methods/captions for the manuscript: `simulations/PAPER_MODEL_SECTIONS_v44.md`.
Data sources used to parameterize the model: `docs/PARAMETERIZATION.md`.

---

## Current model figures

Concise one line each; all consistent with the baked default (see **Current status**). Regenerate with `./venv/bin/python simulations/<script>.py`, or the whole set with `bash simulations/regenerate_figures.sh`.

### Primary — core model behaviour (38)
| Output | Script | What it shows |
|---|---|---|
| `fig_cdk46i_dose_response.png` | `fig_cdk46i_dose_response` | CDK4/6i dose-response on real concentration axis, fit to two palbociclib pRb+ points. |
| `fig_cdk46i_washout.png` | `fig_cdk46i_washout` | CDK4/6i washout in MB: sub-saturating holds cycling residual, saturating arrests; both fully reversible. |
| `fig_celltype_validation.png` | `fig_celltype_validation` | GNP vs MB split validation at the baked cell-type operating point (k_mu_cki), two-step Rb, 30/32. |
| `fig_ezh2_prc2_kinetics.png` | `fig_ezh2_prc2_kinetics` | EZH2i vs EED-KO CyclinD1 de-repression kinetics; serial me-chain as the current model |
| `fig_ezh2_repression_mitogen.png` | `fig_ezh2_repression_mitogen` | EZH2-/CyclinD1 repression shifts the mitogen/commitment threshold up; GNP+MB, with ChIP-correct-mark search. |
| `fig_gnp_hh_mb_vismo_dose.png` | `fig_gnp_hh_mb_vismo_dose` | GNP Hh-dose and MB vismodegib-dose response strips from baked model. |
| `fig_population_mitogen.png` | `fig_population_mitogen` | Many-cell population mitogen dose-response, EZH2->CyclinD1 repression ON vs OFF, GNP and MB. |
| `fig_slide_ezh2_governor.png` | `fig_slide_ezh2_governor` | Talk slide: Hh-dependent cycling, MB transient G0, EZH2 governor threshold; baked model |
| `fig_v44_cd_ezh2_commitment.pdf` | `fig_v44_cd_ezh2_commitment` | Commitment boundary in CyclinD1-drive x EZH2-repression plane; withdrawal and EZH2i rescue unified |
| `fig_v44_cdk46_escape.png` | `fig_v44_cdk46_escape` | CDK4/6-alone Rb-hyper-P escape giving reversible MB transient G0; baked mechanism. |
| `fig_v44_cellcycle_entry.png` | `fig_v44_cellcycle_entry` | EZH2 raises CyclinD1 Hh threshold, delaying cell-cycle entry; EZH2i enters early. |
| `fig_v44_ckoo_feedback_strength.png` | `fig_v44_ckoo_feedback_strength` | Is model ~3x CyclinD1 de-repression supported by Ezh2 cKO/Taze data (per-cell lower bound). |
| `fig_v44_cycle_anatomy.png` | `fig_v44_cycle_anatomy` | Within-cycle anatomy of transient G0 and how CyclinD1 sets its length (p27/pRb markers) |
| `fig_v44_divisions_to_arrest.png` | `fig_v44_divisions_to_arrest` | MB reduced vismo sensitivity via Gli1_epi residual plus resistant proliferating fraction |
| `fig_v44_ezh2_dilution_regimes.png` | `fig_v44_ezh2_dilution_regimes` | EZH2 writer-vs-dilution regimes; baked model sits in EXCESS, mark is rate-sensor. |
| `fig_v44_ezh2i_doseresponse.png` | `fig_v44_ezh2i_doseresponse` | EZH2i graded dose-response across GNP/MB/vismo/palbo contexts; de-repression rescue threshold, p27 G0. |
| `fig_v44_ezh2i_rescue_kinetics.png` | `fig_v44_ezh2i_rescue_kinetics` | EZH2i clears mark, de-represses MYCN-floored CyclinD1, rescues vismo-arrested MB at any dosing time. |
| `fig_v44_feedback_role.png` | `fig_v44_feedback_role` | What the EZH2-/CyclinD1 feedback arm does: buffers CyclinD1/proliferation across mitogen. |
| `fig_v44_g0_arrest_examples.png` | `fig_v44_g0_arrest_examples` | MB lineages dividing then arresting via temporal CDKN random walk; current baked model |
| `fig_v44_g0_mark_dynamics.png` | `fig_v44_g0_mark_dynamics` | G0 mark plateaus, writer-vs-turnover race, monostable reversible (not bistable latch) |
| `fig_v44_g0_mb_specific.png` | `fig_v44_g0_mb_specific` | MB-specific transient G0 as high-CDKI tail; current integrated baked governor. |
| `fig_v44_g0_mechanism.png` | `fig_v44_g0_mechanism` | CyclinD1 to transient-G0 mechanism, p27 marker; EZH2 de-represses CyclinD1 collapsing G0. |
| `fig_v44_g0_temporal.png` | `fig_v44_g0_temporal` | Transient G0 as reversible AR(1) CDKI drift with CDK4/6-escape; labelled current baked model. |
| `fig_v44_gli_kdm6b_eraser.png` | `fig_v44_gli_kdm6b_eraser` | Schematic plus equations of the Gli1->Kdm6b/Jmjd3 H3K27me3 eraser arm. |
| `fig_v44_mitogen_ramp_onoff_population.png` | `fig_v44_mitogen_ramp_onoff_population` | GNP mitogen ramp-on/off population ensemble with vs without mark; entry/exit threshold distributions. |
| `fig_v44_mitogen_sensitivity.png` | `fig_v44_mitogen_sensitivity` | CyclinD1 and proliferation vs Gli1 with EZH2->CyclinD1 feedback on vs off |
| `fig_v44_permanent_transient.png` | `fig_v44_permanent_transient` | Permanent exit vs transient G0 set by CyclinD1/CDK6 drive; baked integrated + CDK4/6 escape, dilution = rate sensor. |
| `fig_v44_progressive_repression.png` | `fig_v44_progressive_repression` | Gas-and-brake: Hh raises CyclinD1 drive and EZH2 brake; threshold Hh to proliferate. |
| `fig_v44_proximal_distal.png` | `fig_v44_proximal_distal` | Two-compartment PRC2 reservoir sustaining reversible transient G0 (30/32) |
| `fig_v44_replicative_dilution.png` | `fig_v44_replicative_dilution` | Replicative dilution ON vs OFF across mitogen sweep; faster cycling ratchets mark down, relieves CyclinD1. |
| `fig_v44_transient_g0_definition.png` | `fig_v44_transient_g0_definition` | Defines transient G0 via CyclinD1/birth-p27 ratio; pRb a bad marker |
| `fig_v44_two_step_rpoint.pdf` | `fig_v44_two_step_rpoint` | Two-step Rb default restores hyper-P Rb as faithful G0/R-point marker |
| `sim_distribution_rescue.png` | `sim_distribution_rescue` | Population-of-models fractional EZH2i rescue phase diagram in CyclinD1 x birth-p27. |
| `fig_v44_mitogen_buffering.png` | `sim_ezh2_mitogen_buffering` | EZH2 tracks mitogen and buffers CyclinD1 (gain 0.85->0.47); R-point commitment memory protects GNP divisions. |
| `sim_ezh2_phase_validation.png` | `sim_ezh2_phase_validation` | Model vs MB IF: EZH2 by phase across DMSO/HU/Palbo/RO3306; stable E2f-gated protein halved at division. |
| `fig_v44_ezh2i_population_dose.png` | `sim_ezh2i_population_dose` | Population EZH2i dose-response: fraction of arrested tumour pushed back into cycle per dose. |
| `fig_v44_mitogen_withdrawal.png` | `sim_mitogen_withdrawal` | GNP mitogen-withdrawal ensemble: depth-gated EZH2->CyclinD1 feedback, sudden vs gradual. |
| `fig_v44_vismo_withdrawal.png` | `sim_vismo_withdrawal` | MB vismodegib ensemble: MYCN-floor CyclinD1, feedback arrests vs no-feedback (EZH2i) coasts, dose-resolved. |

### Exploratory & diagnostic (72)
| Output | Script | What it shows |
|---|---|---|
| `fig_axis_map.pdf` | `fig_axis_map` | Editorial schematic placing v44 parameters on ~4 emergent behavioral axes and couplings |
| `fig_buffering_parameters.png` | `fig_buffering_parameters` | One-at-a-time sensitivity of EZH2/H3K27me3 buffering of CyclinD1; cache present. |
| `fig_diff_gene.png` | `fig_diff_gene` | Poised bivalent diff gene: fast division collapses me3 -> premature activation -> G0 exit. |
| `fig_g1_lengthening_buffer.pdf` | `fig_g1_lengthening_buffer` | Graded G1-lengthening buffer vs abrupt-arrest current model; cache present |
| `fig_governor_noise_filter.png` | `fig_governor_noise_filter` | 6-panel diagnostic: EZH2->CyclinD1 governor as ~15x Shh threshold + dose-response buffer; EZH2_CONC off. |
| `fig_h3k27me3_persistence.png` | `fig_h3k27me3_persistence` | Eraser-vs-ChIP mark-persistence trade-off plus a bistable read-write possibility-space resolution. |
| `fig_latch_consequences.png` | `fig_latch_consequences` | Compares baked (monostable, hi=lo init) vs excluded bistable-latch regime; demonstrates baked G0 is not latched. |
| `fig_latch_dichotomy.png` | `fig_latch_dichotomy` | Didactic monostable (v44 baked) vs bistable-latch mark equation; shows baked case is monostable. |
| `fig_mark_cellcycle_recovery.png` | `fig_mark_cellcycle_recovery` | Cell-cycle duration sets H3K27me3 recovery from S-phase dilution; dilution term consistent with baked model. |
| `fig_matched_withdrawal.png` | `fig_matched_withdrawal` | Transcript-matched ON/OFF/EZH2i withdrawal at varied decline rates; rate-dependent memory vs threshold |
| `fig_memory_parameter_pairs.png` | `fig_memory_parameter_pairs` | Phase diagrams of H3K27me3 buffering/latch across mark-equation parameter pairs; baked point marked |
| `fig_memory_ramp_fluctuation.png` | `fig_memory_ramp_fluctuation` | Memory function on current v44 default under mitogen ramp hysteresis and fluctuation low-pass buffering. |
| `fig_param_distributions.pdf` | `fig_param_distributions` | Population-of-models input draws (CyclinD1, birth-p27, EZH2) plus cycle-duration distribution |
| `fig_pool_vs_halflife.png` | `fig_pool_vs_halflife` | Shows a large pool of a fast transcript reproduces slow-reservoir coasting, supporting measured-halflife view. |
| `fig_prc2_buffering_memory.png` | `fig_prc2_buffering_memory` | PRC2 buffering-memory topologies: repressor (A) buffers dose-response only vs sustaining latch (C). |
| `fig_ramp_pop_diffgene.png` | `fig_ramp_pop_diffgene` | GNP mitogen ramp + population with differentiation tripwire; serial chain, consistent diagnostic. |
| `fig_speed_limiter.png` | `fig_speed_limiter` | Conceptual speed-limiter framework: mark=writing/dilution balance, CyclinD1-governor, bistable R-point go/no-go. |
| `fig_v44_cdmrna_ezh2_commitment.png` | `fig_v44_cdmrna_ezh2_commitment` | Supp commitment map in CyclinD1-transcript x EZH2-repression plane, current two-step-Rb model. |
| `fig_v44_cellcycle_duration.png` | `fig_v44_cellcycle_duration` | Prototype: longer cycle reduces S-phase dilution so H3K27me3 climbs and represses CyclinD1. |
| `fig_v44_differentiation_collapse.png` | `fig_v44_differentiation_collapse` | EZH2 is consequence not cause of exit; EZH2-CyclinD1 repression resists collapse |
| `fig_v44_dilution_replenishment.png` | `fig_v44_dilution_replenishment` | Conceptual: replicative dilution vs EZH2 replenishment set H3K27me3 across FAILING/SUFFICIENT/EXCESS regimes. |
| `fig_v44_ezh2_concentration.pdf` | `fig_v44_ezh2_concentration` | Contrasts dilution (current default) vs legacy concentration convention; own result current |
| `fig_v44_ezh2_knockdown.png` | `fig_v44_ezh2_knockdown` | Cell-cycle response to decreasing EZH2 expression (knockdown) vs catalytic inhibition. |
| `fig_v44_g0_bistability_confirm.png` | `fig_v44_g0_bistability_confirm` | Two-IC fixed-SHH equilibration bistability test; diagnostic consistent with monostable-G0 finding. |
| `fig_v44_h3k27_response.png` | `fig_v44_h3k27_response` | H3K27me3->CyclinD1 transfer functions, dilution shift, de-repression timing intuition plots. |
| `fig_v44_me1me2me3_chain.png` | `fig_v44_me1me2me3_chain` | Serial me0->me1->me2->me3 chain (baked default) kinetics module vs lumped mark; not wired to core. |
| `fig_v44_mitogen_paramspace.png` | `fig_v44_mitogen_paramspace` | 2D mitogen-threshold planes vs EZH2/mark feedback strength and p16 brake |
| `fig_v44_mk_accumulation_decline.pdf` | `fig_v44_mk_accumulation_decline` | Drivers of H3K27me3 accumulation vs decline over divisions (writer/eraser/turnover/dilution sweeps) |
| `fig_v44_mk_parameter_atlas.png` | `fig_v44_mk_parameter_atlas` | 28-parameter atlas of H3K27me3 plateau colored by cycle period; diagnostic sweep. |
| `fig_v44_mk_possibility_space.png` | `fig_v44_mk_possibility_space` | Maps reachable mark behaviors; latch is possibility-space only, not the calibrated model. |
| `fig_v44_mk_timing_interactions.png` | `fig_v44_mk_timing_interactions` | Does division timing gate mark plateau; buffered unless persistent-mark regime the data exclude; lever maps. |
| `fig_v44_mother_g2_ensemble.png` | `fig_v44_mother_g2_ensemble` | Negative-result diagnostic: mother-G2 p27 integrator fails to reproduce the Overton bifurcation. |
| `fig_v44_progressive_repression_timecourse.png` | `fig_v44_progressive_repression_timecourse` | Rising Hh: EZH2->H3K27me3 brake lags and progressively reins in CyclinD1 governor; consistent timecourse. |
| `fig_v44_ptch1_feedback.png` | `fig_v44_ptch1_feedback` | Ptch1 negative-feedback intact in GNP, broken in MB; Gli overshoot vs monotonic plateau |
| `fig_v44_paramsweep_entry.png` | `plot_param_exploration` | Plots cached parameter sweep: entry threshold, withdrawal, and 1-D cuts through operating point. |
| `—` | `sim_async_marker_vs_if` | Async total-protein marker (e.g. EZH2) vs measured IF distribution; abundance-noise gap tool. |
| `sim_async_pRb_vs_if.pdf` | `sim_async_pRb_vs_if` | Asynchronous baked-MB pRb distribution vs measured IF; quantifies spread from phase vs missing G0 mode. |
| `sim_cyclind1_variance_buffering.png` | `sim_cyclind1_variance_buffering` | EZH2 feedback buffering of cell-to-cell CyclinD1 variance; single-cell prediction. |
| `sim_cyclind1_variance_calibrate.png` | `sim_cyclind1_variance_calibrate` | Calibrates CyclinD1 variance-buffering to observed CV 0.70; EZH2-off widening factor. |
| `sim_ezh2_arrest_calib.png` | `sim_ezh2_arrest_calib` | Arrest-model drug calibration (HU, Palbociclib) tuned to MB55 phase proportions; model vs data. |
| `sim_ezh2_arrest_stage1.png` | `sim_ezh2_arrest_stage1` | Stage-1 arrest EdU/p27 gating fit scoring HU/Pal vs measured MB55 proportions. |
| `—` | `sim_ezh2_dilution_tension` | H3K27me3 deposition-vs-dilution tension surfaces setting CyclinD1 transcript |
| `—` | `sim_ezh2_dilution_tension_variants` | HU-axis and MB/vismo variants of EZH2 deposition-vs-dilution tension surface |
| `sim_ezh2_dmso_hu_fit.png` | `sim_ezh2_dmso_hu_fit` | EZH2 rate fit vs DMSO/HU/scRNA targets via runtime overrides; diagnostic fit. |
| `sim_ezh2_kdeez_sweep.png` | `sim_ezh2_kdeez_sweep` | Sweeps EZH2 degradation kDeEZ against withdrawal ramp and within-cycle phase-ratio honesty check. |
| `sim_ezh2_phase_fit.png` | `sim_ezh2_phase_fit` | Fit EZH2 protein-by-phase (kEZbas x kDeEZ) landscape and cycle trace vs data. |
| `fig_v44_ezh2_phaseplane.png` | `sim_ezh2_phaseplane` | Phase-plane/bifurcation of EZH2-mark/CyclinD1 feedback; monostable buffering, mitogen-gated boundary |
| `fig_v44_g0_bifurcation.png` | `sim_g0_bifurcation` | Proliferation-quiescence bifurcation from CyclinD1/p27 heterogeneity; transient-G0 fraction tuned by SHH/EZH2i. |
| `sim_gnp_developmental.png` | `sim_gnp_developmental` | GNP developmental division count with vs without mark at matched output; entry/exit timing comparison. |
| `sim_h3k27_dilution_race.png` | `sim_h3k27_dilution_race` | Phase map of H3K27me3 as EZH2-methylation vs replicative-dilution race (graded) |
| `—` | `sim_hh_division_thresholds` | Hh entry/withdrawal division thresholds, repressed vs not-repressed mark; current v44 |
| `sim_inheritance_partition.png` | `sim_inheritance_partition` | Prototype: calibrated distributions plus p27 partitioning noise flipping near-threshold fate |
| `sim_mark_dip_depth_duration.png` | `sim_mark_dip_depth_duration` | 2-D Hh dip duration x depth arrest sweep, WITH vs WITHOUT mark; mark-attributable arrest diagnostic. |
| `sim_mark_dip_duration.png` | `sim_mark_dip_duration` | Transient Hh-dip duration vs arrest/recovery, probing the mark's memory and critical dip length. |
| `sim_mark_fluctuation_response.pdf` | `sim_mark_fluctuation_response` | Frequency response of proliferation to Hh fluctuations; mark adds low-pass filtering, with vs without. |
| `sim_mark_kdeez_ramp_sweep.png` | `sim_mark_kdeez_ramp_sweep` | Sweeps EZH2 half-life; faster EZH2 weakens mark integration (entry/withdrawal/hysteresis). |
| `sim_mark_lifecycle_hysteresis.png` | `sim_mark_lifecycle_hysteresis` | Entry/steady/withdrawal hysteresis loop with vs without H3K27me3 repression; calibrated single-cell ensemble. |
| `sim_mark_memory_2d_sweep.png` | `sim_mark_memory_2d_sweep` | 2D sweep: EZH2 half-life x mark read-write re-accumulation on memory/arrest |
| `sim_mark_period_crossover.png` | `sim_mark_period_crossover` | Clean mu-only sweep of mark dilution-vs-restoration race and crossover period |
| `sim_p27_bimodality.pdf` | `sim_p27_bimodality` | Tests whether sharp R-point yields emergent bimodal G0 p27 from unimodal birth-p27 input. |
| `fig_v44_population_g0.png` | `sim_population_g0` | Population GNP vs MB transient-G0 fractions from data-grounded distribution machinery. |
| `sim_population_mitogen_ramp.png` | `sim_population_mitogen_ramp` | Heterogeneous-population mitogen ramp entry/exit thresholds with vs without mark repression |
| `sim_ramp_hysteresis_area.png` | `sim_ramp_hysteresis_area` | Mark-lag CyclinD1-vs-SHH hysteresis loop area over ramp/mark/cycle timescales; diagnostic. |
| `sim_ramp_hysteresis_phaseplot.png` | `sim_ramp_hysteresis_phaseplot` | Phase-plane orbits of CyclinD1 vs mark vs SHH, WITH vs WITHOUT mark repression (governor feedback). |
| `sim_ramp_hysteresis_trace_control.pdf` | `sim_ramp_hysteresis_trace_control` | Triangle Hh ramp with vs without mark repression; mark makes entry/exit history-dependent (hysteresis loop). |
| `sim_ramp_parameter_exploration.png` | `sim_ramp_parameter_exploration` | Entry threshold and withdrawal divisions over mark-strength x cycle-duration plane. |
| `sim_ramp_timing_exploration.png` | `sim_ramp_timing_exploration` | Fine-grid hysteresis gap vs mark turnover, EZH2 t1/2, cycle period across ramp speeds |
| `sim_ramp_traces_diagnostic.png` | `sim_ramp_traces_diagnostic` | Full mitogen-ramp traces across mark-strength x cycle-duration corners; mechanism diagnostic |
| `sim_singlecell_withdrawal_g0.pdf` | `sim_singlecell_withdrawal_g0` | Single-cell ensemble through Hh withdrawal; transient-G0 and interval lengthening with/without mark |
| `—` | `sim_speed_limiter` | Self-consistent speed-limiter toy model; mark as proliferation-coupled rate-sensor governor. |
| `fig_ws1_transient_g0.png` | `sim_ws1_transient_g0` | MB transient G0 on baked dynamic-CDKI: CdP21 buffer + high INK4 + birth-p27 dwell/re-entry. |
| `sim_ws2_mitogen_g1.png` | `sim_ws2_mitogen_g1` | Mitogen dose sets G1 duration in GNP vs MB via p27 buffering; baked dynamic-CDKI, consistent. |

---

## Compute / calibration helpers (no standalone figure) (32)
Write caches / optimizer output read by the figure scripts, or are prototypes / levers / smoke tests.

| Script | Role |
|---|---|
| `sim_buffering_parameters` | One-at-a-time buffering sensitivity sweep; writes sim_buffering_parameters_cache.pkl, no figure. |
| `sim_cdc20_lever` | Compute-only lever test: CHK1-slowed Cdc20/CycA destruction raising S-phase EZH2; prints metrics |
| `sim_cdc25a_lever` | Compute-only lever test of Cdc25A checkpoint on S-progression; no standalone figure |
| `sim_cdk46i_dose_response` | Compute feeder: CDK4/6i modeled as graded dose-response in GNP vs MB; writes JSON, no figure. |
| `sim_cdk46i_washout` | CDK4/6i washout cycling-fraction kinetics; writes JSON cache, no standalone figure |
| `sim_checkpoint_diagnostics` | v43 intra-S checkpoint diagnostics (HU effect on period/S/EZH2); diagnostics helper |
| `sim_diff_gene` | Bivalent differentiation gene latch model; writes sim_diff_gene_cache.pkl, no standalone figure. |
| `sim_ezh2_bursting` | Transcriptional-bursting EZH2 noise; writes model json cache, no standalone figure; t1/2 ~42h current. |
| `sim_ezh2_calib` | DMSO-only bursting-amplitude calibration feeder; prints one parseable CALIB line, no figure |
| `sim_ezh2_joint_fit` | Joint 4-parameter EZH2 rate optimizer fit to DMSO/HU protein+mRNA; identifiability feeder. |
| `sim_ezh2_prc2_kinetics` | Feeder: writes ezh2_prc2_kinetics_results.json consumed by fig_ezh2_prc2_kinetics |
| `sim_ezh2_repression_mitogen` | Overnight compute feeder writing _cache.pkl read by fig_ezh2_repression_mitogen; full v44 default. |
| `sim_ezh2_treatments` | MB population EZH2 under DMSO/CDK4-6i/HU; dumps JSON cache, no standalone figure. |
| `sim_ezh2_treatments_smoke` | Smoke-test grid for EZH2-by-phase across four arrest drugs vs model. |
| `sim_ezh2_variance_decomp` | Decomposes EZH2 snapshot variance into cell-to-cell vs within-cell; writes json, no figure. |
| `sim_g1_lengthening_buffer` | Reduced graded-G1-lengthening prototype; writes _cache.npz, figure done by fig_g1_lengthening_buffer. |
| `sim_hu_dose_response` | Diagnostic: sweeps HU dose x CHK1-Pb testing longer-S->more-EZH2; prints ratios only, no figure. |
| `sim_latch_consequences` | Baked-vs-latch regime experiments; caches sim_latch_consequences_cache.pkl, no standalone figure |
| `sim_mark_cellcycle_recovery` | Computes mark S-dilution recovery vs cycle length Tc; writes pickle cache. |
| `sim_mark_memory_cooptimize` | Optimizer feeder maximizing withdrawal-memory mark while pinning CyclinD1 MB/GNP to 5.07; writes JSON. |
| `sim_mark_persistence_sweep` | Overnight parameter sweep for longer mark persistence while keeping MB<GNP; no figure |
| `sim_matched_withdrawal` | Transcript-matched withdrawal population sim; caches npz read by fig_matched_withdrawal.py. |
| `sim_memory_ramp_fluctuation` | Ramp-hysteresis and fluctuation memory on current default; writes _cache.npz, figure is separate script. |
| `sim_population_mitogen` | Population mitogen dose-response feeder; writes sim_population_mitogen_cache.pkl, figure made separately. |
| `sim_prc2_buffering_memory` | Writes prc2_buffering_memory cache npz; figure produced by separate fig_ script |
| `sim_prc2_writer_dilution_states` | Toy writing-vs-dilution two-state model feeder; writes _cache.npz read by fig_prc2_writer_dilution_states. |
| `sim_punctate_bistable_proto` | Prototype exploring bistable latched mark with punctate bursts; prints stats, no figure |
| `sim_punctate_meanfield_proto` | Reduced 1-locus stochastic prototype testing punctate vs continuous transcription eviction |
| `sim_ramp_pop_diffgene` | Compute feeder writing sim_ramp_pop_diffgene_cache.pkl read by fig_ramp_pop_diffgene. |
| `sim_realistic_mitogen` | Evidence-driven Shh input + low-pass decomposition feeder; writes _cache.pkl read by fig_realistic_mitogen. |
| `sim_recalibrate_hill` | Re-calibration sweep reporting all EZH2 targets; optimizer/calibration feeder, no figure |
| `sim_stress_ezh2` | Calibrates replication-stress to EZH2-transcription coupling (v43 k_EZH2_stress sweep); no figure. |

---

## Archived — stale narrative or broken (18; kept for provenance, NOT current)

**Moved to `simulations/archive/` (2026-10-17)** — out of the active set, git-tracked (reversible). See `simulations/archive/README.md`. These describe a **superseded model state** (old EZH2 t½ ~77 h, 11–12 h CyclinD1 mRNA reservoir, EZH2-concentration convention, bistable / "lock-in" G0, old validation scores, the v42 architecture) or read data files that no longer exist. **Do not cite as current.**

| Script | Why archived |
|---|---|
| `fig_param_biplot` | _broken_ — Reads sensitivity.json/param_meta.json from a stale scratchpad path (different session UUID) that no longer exists; figure present but regen fails. |
| `fig_param_heatmap` | _broken_ — Reads a specific /private/tmp scratchpad path (allparams.json, param_cats.json) that is gone. |
| `fig_param_settings` | _broken_ — Reads a hardcoded scratchpad sensitivity.json (and param_meta.json) that does not exist; fails at import. |
| `fig_cyclind1_memory` | _stale_ — Main narrative is the superseded slow-transcript reservoir; current model = measured 3.4h mRNA, memory = large reservoir draining past threshold. |
| `fig_g1_lengthening_withdrawal` | _stale_ — Old score 28/29 (current 30/32); latent not-baked hypothesis presented with stale score. |
| `fig_prc2_writer_dilution_states` | _stale_ — Main narrative is bistable latch/hysteresis memory; current G0 is monostable/reversible, EZH2 lock-in retired. |
| `fig_realistic_mitogen` | _stale_ — Main narrative = slow CyclinD1 mRNA reservoir as the low-pass/memory carrier; superseded by measured 3.4h mRNA + large-reservoir-draining memory. |
| `fig_v42_architecture` | _stale_ — Depicts the v42 architecture; superseded by the v44 3-tier structural map; not current model structure. |
| `fig_v44_cellcycle_exit` | _stale_ — Main narrative = 'epigenetically locking the arrested state' / lock-in; current model G0 arrest is MONOSTABLE/reversible, EZH2 lock-in retired. |
| `fig_v44_g0_mark_regime` | _stale_ — Operating point/uncertainty box uses kDeEZ=1.5e-4 (~77h); current baked EZH2 t1/2 ~42h (kDeEZ=2.77e-4) |
| `fig_v44_mb_g0_birthp27` | _stale_ — Birth-p27 crutch retired; current MB transient G0 comes from CDK4/6-alone hyper-P escape in the dynamic-CDKI stack. |
| `fig_v44_mitogen_ramp_onoff` | _stale_ — docstring states old validation score '25/27 parameters' as current (now 30/32) |
| `sim_ezh2_compensation` | _stale_ — uses superseded build_model_v42_mycn with MYCN amplification; MYCN is elevated-not-amplified in current model |
| `sim_ezh2_data_constrained` | _stale_ — Uses build_model_v42_mycn (superseded by v44); presents 1.13x gap/MYCN->EZH2 as open though MYCN_expr is now baked. |
| `sim_ezh2_feedback_impact` | _stale_ — Uses superseded build_model_v42_mycn and the dead K_EZH2_repression knob; predates v44 baked model. |
| `sim_mark_battery` | _stale_ — Memory = CyclinD1 reservoir draining, not slow mark discharge; uses fast t1/2 10h / 25/27 framing. |
| `sim_mark_memory_chromatin_sweep` | _stale_ — built around superseded fast EZH2 t1/2 10h and old 25/27 score as the data-supported value |
| `sim_ramp_hysteresis_trace` | _stale_ — Memory framed as the H3K27me3 mark driving a hysteresis loop; current model: memory = CyclinD1 reservoir, dilution is a rate-sensor, G0 monostable (no latch/hysteresis). |

> **v45 figures** (the parallel stochastic-commitment successor) remain archived; **v44 is the working model**.
