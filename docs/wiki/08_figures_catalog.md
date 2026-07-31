# Figures Catalog

> The figure suite is the visible record of everything this project *tried*. It spans two eras — a handful of **manuscript panels** that ship with Chahin et al., and dozens of **exploratory figures** that interrogated the model as JP's hypotheses evolved (mRNA-reservoir memory, EZH2 governor, transient G0, two-cyclin, dynamic CDKI). This page organizes them by theme so a reader can find "the figure that shows X." For how these figures map to the model's history see [Evolution](03_evolution_timeline.md); for the mechanisms behind them see [Mechanisms Explored](04_mechanisms_explored.md) and [Things Tried and Abandoned](05_things_tried_and_abandoned.md); for the numbers they're checked against see [Calibration and Validation](06_calibration_and_validation.md).

All figure scripts live in `Ezh2_CyclinD1_Sphase/simulations/`. Each writes a `.png` (preview) and a `.pdf` (vector, for the paper). The canonical manifest is `docs/FIGURES.md`; this page is a themed cross-index of it plus the newer (post-manifest) figures. Model time is real-time minutes plotted as hours (`/60`, `UPH=62.8` in newer scripts).

---

## How the suite is regenerated: FAST vs SLOW tiers

`simulations/regenerate_figures.sh` is the canonical "re-run everything after a model change" script. It has **two tiers**:

- **FAST tier** — auto-globs *every* `simulations/fig_v44_*.py` (per-condition, deterministic single-cell or algebraic; seconds to a couple of minutes each). New `fig_v44_*` scripts are picked up automatically with no edit to the manifest.
- **SLOW tier** — ensemble / N-cell heterogeneity figures that do **not** match the `fig_v44_*` glob and take minutes to hours. They are listed **explicitly** in the manifest and run with `--fresh` so they recompute against the changed model instead of silently reusing a stale `.npz`/`.pkl` cache. The current SLOW list: `sim_g0_bifurcation` (~10 min), `sim_population_g0 --n=100` (~12 min), `sim_ezh2_phaseplane` (~35 min), `sim_ezh2_dilution_tension` (~10 min), `sim_ezh2_dilution_tension_variants` (~13 min), `sim_ezh2_mitogen_buffering` (~20 min), `sim_mitogen_withdrawal` (~2.5 h), `sim_vismo_withdrawal` (~2.5 h), `fig_v44_ezh2i_rescue_kinetics` (~15 min).

> **Caveat (roadrunner):** `reset()` resets *species* but not *parameters*, so every sweep call must re-establish a full baseline or values leak (a documented gotcha that corrupted an early tornado). Many caches (`sweep_param_exploration`, the ramp-hysteresis family) are **not** in the auto-regen and must be refreshed by hand.

Many newer figures split **compute** (a `sim_*.py` that writes a cache) from **plot** (a `fig_*.py` that reads it). Where that split exists it's noted below.

The unified harness `simulations/validate_v44.py` defines the conditions, the division counter, the p27 phase classifier, and the experimental-target checks that the figure scripts import. Convenience knobs: CDK4/6i = `kPhRbCd=0`; serum-starve = `k_Cd_translation=0`; HHi/vismo = `GDC0449=1`; EZH2i = `EZH2i=1`; MB context = `SHH 0.5, Ptch1_copy_number 0.1, MYCN_amplification 2.8`.

---

## 1. Manuscript figures (the paper panels)

These are the model-derived panels in Chahin et al. See also [Paper model figure map](08_figures_catalog.md) cross-referenced in memory `ezh2-paper-model-figure-map`.

| Panel | Script | What it shows |
|---|---|---|
| **Fig. 5A–C** | `fig_v44_fig5_rescue.py` | Ptch+/− MB / +HHi (arrest) / +HHi+Ezh2i (rescue) Cyclin B traces — the central rescue result. |
| **Fig. 5C (pop.)** | `fig_v44_fig5_population.py` | Population cycling fraction (N=120, CyclinD1/p27 heterogeneity): vismo↓, EZH2i restores, CDK4/6i no rescue. |
| **Supp. Fig. 7A** | `fig_v44_wiring.py` | Model architecture / wiring diagram (pure schematic, no sim). |
| **Supp. Fig. 7B–M** | `fig_v44_validation.py` | Validation grid: Cyclin B traces, CyclinD1/Mycn/Gli1 model-vs-data bars, EZH2, division counts, time courses. |
| **Supp. Fig. 7x** | `fig_v44_ptch1_feedback.py` | Gli→Ptch1 negative feedback and its disruption in MB (`Ptch1_copy_number` = functional fraction, 1.0 GNP / 0.1 MB). |
| **Supp. Fig. 8A–H** | `fig_v44_feedback.py` | EZH2→CyclinD1 feedback impact ±feedback (Cyclin B, CyclinD1 mRNA, EZH2 protein, phase). Feedback OFF = `f0_mk=1.0`. |
| **Supp. Fig. 8I,J** | `fig_v44_fig5_rescue.py` | MB+CDK4/6i and MB+CDK4/6i+Ezh2i (no rescue) — emitted alongside Fig. 5. |
| **Conditions → param space** | `fig_v44_condition_paramspace.py` | Input-knob matrix per condition + CyclinD1×p16 plane with cycle/arrest boundary; shows MB+HHi arresting at the same CyclinD1 a GNP cycles at (p16 raises the threshold; EZH2i crosses back). |
| **push/pull schematic** | hand-authored `fig_v44_pushpull_schematic.svg` | Balance-sheet capstone: green drives up, red pushes down; EZH2→H3K27me3⊣CyclinD1 leaky-rheostat; three emergent-behavior callouts. |

Drop-in Results/Methods/captions: `simulations/PAPER_MODEL_SECTIONS_v44.md`.

---

## 2. Validation & cell-type figures

| Script | What it shows | Supports |
|---|---|---|
| `fig_v44_validation.py` | The full model-vs-data validation grid (see manuscript table). | Supp. Fig. 7 |
| `fig_celltype_validation.py` | **Cell-type-split validation** (GNP ~17.2h / MB ~23.5h core): left = parameter differences (only biological identities differ; cycle-length split emerges from MB CDKi via `k_mu_cki`); middle/right = per-cell-type validation forest; bottom = MB↔GNP cross-ratios. Baked two-cyclin recal, 30/32. | validation-targets spec |
| `fig_v44_ckoo_feedback_strength.py` | Is the ~3× CyclinD1 supported by the cKO? Recombination correction of the bulk fold + feedback-strength sweep. | cKO scoping |

---

## 3. Mechanism / architecture schematics

| Script | What it shows |
|---|---|
| `fig_v44_wiring.py` | Full v44 wiring diagram (Supp. Fig. 7A). |
| `fig_v42_architecture.py` | The older v42 architecture schematic (MYCN–EZH2–CyclinD1–HH). |
| `fig_v44_gli_kdm6b_eraser.py` | Detail schematic + equations of the Gli1→Kdm6b/Jmjd3→H3K27me3 **eraser** arm; makes explicit that Kdm6b is quasi-static ∝ Gli1 (lumped flux `k_jmjd3_gli·Gli1·Mk`). Pure schematic. |
| `fig_axis_map.py` | Editorial schematic of the model's **~4 emergent behavioral axes** (proliferation↔quiescence, cycle speed, etc.) with each parameter placed on the axis it moves; shows CKI hitting two axes and `decouple_commit` severing the speed→threshold leak. |
| `fig_speed_limiter.py` | JP's **speed-limiter framework**: mark = writing/dilution balance → a max sustainable division rate; CyclinD1-as-target = governor; division decision bistable, mark graded; GNP/MB regime map. (Model in `sim_speed_limiter.py`.) |

---

## 4. EZH2 / H3K27me3 mark figures

The largest theme — the mark's kinetics, its repression of CyclinD1, and whether it can carry memory.

### Mark response & repression
| Script | What it shows |
|---|---|
| `fig_v44_h3k27_response.py` | **Two figures.** `_response`: static transfer functions — leaky-Hill readout `R(Mk)`, methylation dose-response `Mk_ss(EZH2)`, the throttle, and per-parameter small-multiples (`del_mk`, `K_mk`, `k_w_mk` dominate). `_timing`: full-model kinetics — de-repression t½≈15h, dilution sawtooth, mean-Mk-vs-period, entry timing on a Hh ramp. |
| `fig_v44_progressive_repression.py` | Core consequence: EZH2 is mitogen-responsive → **progressive** H3K27me3 repression. Surviving fraction R falls 80%→35% across the Hh range; EZH2 raises the proliferation threshold from ≤0.08 to SHH≈0.35. |
| `fig_v44_progressive_repression_timecourse.py` | Same, as a real-time trace under gradually rising Hh (the slow brake lags the gas). |
| `fig_v44_replicative_dilution.py` | Replicative dilution (S-phase halving `Mk=0.5·Mk`) limits EZH2 repression; dilution ON vs OFF across a mitogen sweep — the gap widens with mitogen (faster cycle → more dilution). |
| `fig_v44_cellcycle_duration.py` | Cell-cycle duration = the dilution frequency (1/period) sets the brake; longer cycle → more mark; developmental slowing → progressive repression toward exit. |
| `fig_ezh2_repression_mitogen.py` | EZH2⊣CyclinD1 ON vs OFF → the mitogen requirement (GNP & MB), ramp hysteresis, and a 2-D search for a regime with both the mitogen-shift AND a ChIP-correct mark. (memory `ezh2-repression-mitogen-requirement`) |
| `fig_governor_noise_filter.py` | **DEFINITIVE governor figure** (post-`k_Cd_tx_basal ×0.5` bake): compares ON / OFF (`f0_prc2=1`) / EZH2i. Proliferation now Shh-dependent for all (all arrest at SHH=0); EZH2 raises the Shh threshold ~15× and buffers the dose-response, matching cKO. (memory `ezh2-governor-reservoir-remeasure`) |
| `fig_buffering_parameters.py` | One-at-a-time sensitivity of the EZH2/H3K27me3 buffering of CyclinD1 + the summary equation. |

### Mark kinetics: chain, dilution, accumulation
| Script | What it shows |
|---|---|
| `fig_v44_me1me2me3_chain.py` | Serial me0→me1→me2→me3 chain (promoted to default) vs the lumped single-Mk mark: the slow me2→me3 step adds a lag + a cell-cycle-period dependence the lumped mark lacks. (memory `v44-serial-methylation-chain-promoted`) |
| `fig_ezh2_prc2_kinetics.py` | Kinetics of EZH2i (catalytic, complex stays) vs EED KO (complex removed); the mark decay = the de-repression timescale; historical→current model progression. |
| `fig_v44_mk_accumulation_decline.py` | What drives accumulation vs decline of Mk over divisions: writer/read-write↑ accumulate, eraser/turnover/dilution↑ decline; writer×eraser phase diagram. (memory `v44-mk-accumulation-decline-sweep`) |
| `fig_v44_mk_parameter_atlas.py` | 28-panel **atlas** of the Mk plateau, one per parameter, colored by cycle period (division-timing-mediated accumulation made visible). |
| `fig_v44_mk_timing_interactions.py` | Does division timing gate the plateau? Plateau is buffered against period unless the mark is persistent (strong-memory regime the data exclude) + pairwise interaction maps. |
| `fig_mark_cellcycle_recovery.py` | Cell-cycle duration sets H3K27me3 recovery from S-phase dilution (the missing buffering term). |
| `fig_prc2_writer_dilution_states.py` | Two proliferation states from EZH2-writing vs division-dilution balance; ramp hysteresis; division-coupled dilution makes memory less sticky than constant turnover. |
| `fig_ccnd1_halflife.py` | **CCND1 transcript half-life vs the reservoir assumption** (RNAdecayCafe atlas): measured CCND1 t½ ≈1.5–2.5h vs the model's baked ~11h reservoir — the reservoir is ~5–7× too slow. (memory `ccnd1-halflife-data-vs-reservoir`) |

### EZH2 protein-by-phase fits (calibration diagnostics)
| Script | What it shows |
|---|---|
| `sim_ezh2_kdeez_sweep.py` | EZH2 protein t½ (`kDeEZ`) sweep vs the Hh-withdrawal ramp and the within-cycle phase ratio — a single constitutive `kDeEZ` can't give both a modest within-cycle contrast AND a big G0 drop. |
| `sim_ezh2_phase_fit.py` | `kEZbas × kDeEZ` fit vs measured phase ratios — raising `kEZbas`→~0.004 fixes both protein (1.48×) and mRNA (~2×) ratios; decouples the ratio from `kDeEZ`. (memory `ezh2-baseline-transcription-fix`) |
| `sim_ezh2_dmso_hu_fit.py` | Finer fit vs MB55 DMSO+HU IF + scRNA-seq; documents that Palbociclib & RO3306 collapse the model to a 2N low-E2f state (a model-file limitation). |
| `sim_ezh2_joint_fit.py` | Joint 4-parameter (`kEZbas×kEZE2f×kDeEZ×kDeEZm`) optimization; finds `kDeEZm`→~0.0016 (EZH2 mRNA ~12× more stable); structural miss: data protein is buffered relative to its mRNA. |
| `sim_ezh2_arrest_calib.py` | Arrest-model calibration vs MB55 phase proportions; documents the two structural gaps (active-replication S-gate; no stable G1-arrest). |
| `sim_ezh2_phase_validation.py`, `sim_ezh2_treatments_smoke.py`, `sim_ezh2_data_constrained.py`, `sim_ezh2_compensation.py`, `sim_ezh2_prc2_kinetics.py` | Supporting EZH2 phase/treatment validation, data-constrained fits, compensation checks, and PRC2 kinetics compute. |
| `fig_v44_ezh2_concentration.py` | EZH2 concentration convention: smooth turnover-reset carryover vs the legacy amount-halving sawtooth (division-neutral). (memory `v44-ezh2-concentration-promoted`, later reversed) |
| `fig_v44_ezh2_knockdown.py` | Response to decreasing EZH2 **expression** (`kTlEZ`↓) vs catalytic inhibition — separates period / division-rate / transient-G0 senses of "more rapid." |

### The mark as memory: latch vs buffer
| Script | What it shows |
|---|---|
| `fig_latch_dichotomy.py` | Didactic monostable (v44 baked, high `a0`) vs bistable (low `a0` + cooperative `n_rw=3`) for the reduced mark equation. |
| `fig_v44_mk_possibility_space.py` | Possibility space: calibrated architecture is monostable at all params; cooperative read-write + suppressed nucleation gives a bistable latch. (memory `v44-mk-possibility-space`) |
| `fig_latch_consequences.py` | The full v44 model in the latch regime: latch reachable but arrest-gated, never spontaneous, and breaks calibration (GNP loses Hh-dependence). (memory `v44-latch-regime-consequences`) |
| `fig_memory_parameter_pairs.py` | Memory phase diagrams over pairs of mark-equation params; two faces of memory — passive persistence (τ = ln2/δ_eff) vs latch (bistability width). |
| `fig_prc2_buffering_memory.py` | Two topologies for "PRC2 on the CyclinD1 axis": A (represses CyclinD1) buffers the dose-response only; C (double-negative, sustains proliferation) is needed for withdrawal memory. (memory `prc2-buffering-memory-topologies`) |
| `fig_pool_vs_halflife.py` | Does a big CyclinD1 mRNA pool substitute for slow degradation? Half-life is pool-independent; a big pool buys time-above-threshold logarithmically; too much pool (>~6×) makes it Shh-independent. |

### G0 mark dynamics (is arrest bistable?)
| Script | What it shows |
|---|---|
| `fig_v44_g0_mark_dynamics.py` | **Consolidated**: mark plateaus in arrest (writer E2f-gated, collapses); a race between writer persistence and mark turnover; monostable/reversible from two initial conditions. |
| `fig_v44_g0_mark_regime.py` | Phase diagram (EZH2 t½ × mark turnover): RELAX (reversible) vs ACCUMULATE (ratchet); operating point in RELAX near the boundary. |
| `fig_v44_g0_bistability_confirm.py` | Definitive two-initial-condition equilibration → **monostable**; the mark does not create a bistable G0 lock. (memory `v44-g0-arrest-monostable`) |
| `fig_v44_g0_hysteresis.py` | **SUPERSEDED** — false-positived bistability (a slow-transient artifact); kept as a cautionary record. |

---

## 5. G0 / G1 / transient-quiescence figures

| Script | What it shows |
|---|---|
| `fig_v44_transient_g0_definition.py` | **Defines** transient G0 (p27-high dwell that resolves) vs immediate G1 vs sustained arrest; CyclinD1×birth-p27 Fan–Meyer phase diagram. p27 (not pRb) is the marker. (memory `v44-transient-g0-finding`) |
| `fig_v44_cycle_anatomy.py` | One-cycle anatomy: transient G0 = the p27-high window (phospho-Rb saturates early). |
| `fig_v44_g0_mechanism.py` | CyclinD1 sets transient-G0 duration (p27 marker); EZH2i collapses it; arrest below threshold. |
| `fig_v44_g0_bifurcation.py` (`sim_g0_bifurcation.py`) | Proliferation–quiescence bifurcation from CyclinD1 / birth-p27 (`P21_div`) heterogeneity (N=140). |
| `fig_v44_two_step_rpoint.py` | Two-step Rb R-point hand-off (mono-P Rbm → hyper-P pRb); restores hyper-P Rb as a faithful G0 marker. (memory `v44-transient-g0-reframe...`) |
| `fig_v44_mb_g0_birthp27.py` | MB p27-high/G0-rich via the MB-specific birth-p27 (GNP 0% vs MB ~21% G0), fold-safe. (memory `mb-g0-regression-and-p27-reconception`) |
| `fig_v44_population_g0.py` (`sim_population_g0.py`) | GNP-vs-MB population transient-G0 + quiescence-onset-vs-mitogen split + Spencer ~95–98% sister concordance. Flow-comparable count-fraction metric (GNP ~9% / MB ~18%). |
| `fig_v44_mother_g2_ensemble.py` | **Negative result**: the mother-G2 p27 integrator does NOT reproduce the Overton graded→binary bifurcation; documents the two design flaws + the fix. (memory `v44-mother-g2-negative-result`) |
| `sim_ws1_transient_g0.py` → **`fig_ws1_transient_g0`** | **WS1 (JP 2026-07-30, on the baked dynamic-CDKI model):** why MB shows transient G0 and GNP doesn't — commitment is the `CdP21` buffer; MB's high INK4 + birth-p27 → weak buffering → free p27 stays high → pre-S pRb-low dwell → re-entry. |
| `sim_ws2_mitogen_g1.py` → **`fig_ws2_mitogen_g1`** | **WS2 (JP 2026-07-30):** how mitogen (Shh) dose sets **G1 duration** in GNP vs MB — higher Cd → more p27 buffered → CDK2 fires early → short G1; as Shh falls, G1 lengthens graded then arrests; MB sits at longer G1 and arrests at a higher Shh floor. |
| `fig_g1_lengthening_buffer.py` (`sim_g1_lengthening_buffer.py`) | Graded G1/cycle lengthening (no transient G0) as the EZH2/CyclinD1 buffer; current model holds Tc constant then abruptly arrests. (memory `gnp-g1-lengthening-buffering`) |
| `fig_g1_lengthening_withdrawal.py` | Withdrawal memory via graded G1-lengthening (lever = latent `mu_eff` slowdown, `mu_min_frac`/`K_g1len`); validation-compatible; JP's EZH2-reduction driver breaks validation. (memory `withdrawal-graded-g1-lengthening`) |
| `probe_transient_g0_v2.py`, `analyze_cyclind1_g0.py`, `probe_two_step.py` | Transient-G0 / CyclinD1 / two-step-Rb probes (no figure). |

---

## 6. Mitogen dose, buffering & noise-filtering

| Script | What it shows |
|---|---|
| `fig_v44_cellcycle_entry.py` | GNP cell-cycle **entry** under rising Hh, ±EZH2: EZH2-active commits later (t≈106h, SHH≈0.59) vs EZH2i (t≈50h, SHH≈0.22) — EZH2 gates entry to a genuine Hh threshold. |
| `fig_v44_cellcycle_exit.py` | Linear Hh **withdrawal**: EZH2 exits earlier/at higher Hh (SHH≈0.54) vs EZH2i (SHH≈0.18); GNP arrest reversible (EZH2 falls with exit). |
| `fig_v44_mitogen_ramp_onoff.py` | Single-cell staircase ramp-ON→sustain→ramp-OFF ±mark: entry threshold + exit collapse. |
| `fig_v44_mitogen_ramp_onoff_population.py` | Same as an ensemble (calibrated draws) → distributions of entry/exit thresholds. |
| `fig_v44_mitogen_sensitivity.py` | Mitogen sensitivity in Gli1 units: GNP=SHH, MB=vismodegib; ±feedback flattens the dose-response and raises the threshold. |
| `fig_v44_mitogen_paramspace.py` | 2-D proliferation planes: SHH × feedback strength `f0_mk`, and SHH × p16; threshold rises with both. |
| `fig_v44_condition_paramspace.py` | (Manuscript) conditions → input params + CyclinD1×p16 arrest boundary. |
| `fig_v44_dilution_tension.pdf` (`sim_ezh2_dilution_tension.py`) | Methylation (S-phase + mitogen) vs replicative dilution → CyclinD1 transcript, over (mitogen × S-phase). |
| `fig_v44_dilution_tension_variants.pdf` (`sim_ezh2_dilution_tension_variants.py`) | HU-axis (replication stress) + MB/vismo variants of the tension surface. |
| `fig_v44_mitogen_buffering.pdf` (`sim_ezh2_mitogen_buffering.py`) | EZH2 buffers CyclinD1 vs mitogen: sensitivity/onset, **fluctuation frequency-response** (Berry–Howard low-pass), partial withdrawal, R-point commitment memory. |
| `fig_realistic_mitogen.py` (`sim_realistic_mitogen.py`) | **Biologically realistic Shh input** (tonic + slow OU drift + fast OU noise, cilium duty cycle) + Gli decomposition + demonstration that the CyclinD1 mRNA reservoir low-pass-filters it. (memory `v44-realistic-mitogen-input`) |
| `fig_cyclind1_memory.py` | **CyclinD1 = the Hh memory carrier** — the memory is the cumulative CyclinD1 **transcript** reservoir (slow ~12h), protein is the fast readout; EZH2 = incoherent feedforward that coasts ~1 cycle after Hh-off. Baked default. (memory `v44-cyclind1-hh-memory-carrier`) |
| `fig_population_mitogen.py` (`sim_population_mitogen.py`) | Many-cell population mitogen dose-response, EZH2⊣CyclinD1 ON vs OFF, GNP & MB. |
| `sim_population_mitogen_ramp.py` → `sim_population_mitogen_ramp.pdf` | Population Hh ramp up/down: the mark raises median entry (0.32→0.54) and exit (0.30→0.46) thresholds — symmetric gatekeeper. |
| `fig_memory_ramp_fluctuation.py` (`sim_memory_ramp_fluctuation.py`) | The memory function on the current default (Cdh1-EZH2 + G1-lengthening + Skp2) under ramp (hysteresis) and fluctuation (low-pass). |
| `fig_buffering_parameters.py` (`sim_buffering_parameters.py`) | Which parameters support buffering (static B=1−R + dynamic fluctuation attenuation). |

### High-resolution parameter sweeps (not in auto-regen)
| Script | What it shows |
|---|---|
| `sweep_param_exploration.py` → `plot_param_exploration.py` | **High-res** 50×50×10 sweep → `fig_v44_paramsweep_{entry,withdrawal,cumprolif,cuts}.pdf`; entry threshold, withdrawal divisions, cumulative proliferation (2^divisions), 1-D cuts. Cached (~2h). |
| `fig_v44_param_exploration.py` | Coarse predecessor (EZH2 responsiveness × inhibition depth × Hh decline rate). Superseded by the above. |
| `sim_ramp_parameter_exploration.py` | What shifts entry/exit thresholds: H3K27me3 strength × cycle duration; entry dominated by intrinsic features, exit by duration/S-phase. |
| `sim_ramp_timing_exploration.py` | Timing sweep: mark/EZH2 kinetics × ramp speed → entry/exit hysteresis; documents the noisy threshold metric. |

---

## 7. EZH2 / mark ramp-hysteresis & fluctuation-response family

A dense sub-suite probing whether the mark is a memory or just a filter (mostly a "the mark is a slow trend-tracker, not a fluctuation filter" conclusion).

| Script | What it shows |
|---|---|
| `sim_ramp_hysteresis_area.py` | Hysteresis as a smooth **area** of the CyclinD1-vs-SHH loop (supersedes the noisy threshold-gap readout); pre-equilibration scaled to the kinetic timescale. |
| `sim_ramp_hysteresis_trace.py` | Worked single-triangle-ramp trace — the transient/kinetic entry-exit hysteresis (entry 0.32 / exit 0.09). |
| `sim_ramp_hysteresis_trace_control.py` | Same, ±repression (`f0_mk` 0.233 vs 1.0); caches both runs. |
| `sim_ramp_hysteresis_phaseplot.py` | Phase-plane portraits (CyclinD1 vs mark): negative-feedback manifold WITH, slope flips positive WITHOUT. |
| `sim_ramp_hysteresis_mark_control.py` | Is the loop caused by the mark? Mark-ON loop ~0.33 flat across speeds; mark-OFF decays to 0 — >85% of the loop is H3K27me3. |
| `sim_ramp_traces_diagnostic.py` | Full trajectories for the 4 corners of (mark strength × cycle duration); confirms the "first-division-during-ramp" metric is unreliable. |
| `sim_mark_fluctuation_response.py` | Frequency response to Hh(t)=sin: the mark is **not** a fluctuation filter (scales CyclinD1, doesn't low-pass it); G0/commitment IS filtered but by the cell-cycle switch, not the mark. Also emits `sim_mark_fluctuation_transientg0.pdf`. |
| `sim_mark_dip_duration.py` / `sim_mark_dip_depth_duration.py` | Response to a transient Hh **dip** of varying duration — the mark makes cells far more vulnerable (50% arrest at a ~35h dip); an epigenetic memory of Hh loss. |
| `sim_mark_kdeez_ramp_sweep.py` | EZH2 t½ sweep: the mark's ramp effects turn on sharply above t½~25–40h — quantifies the tension (data pins t½~10h; memory story needs ≳25h). |
| `sim_mark_memory_2d_sweep.py` / `sim_mark_memory_chromatin_sweep.py` / `sim_mark_memory_cooptimize.py` | 2-D memory: EZH2 t½ × read-write re-accumulation (`k_w_mk`). Iso-memory contours run diagonally → the tension is **resolved** by a self-reinforcing read-write (fast writer + strong chromatin memory). |
| `sim_mark_period_crossover.py` | Clean μ-only sweep: the dilution-vs-restoration race + the crossover period (~16h); default period sits above it (restoration-dominated). |
| `sim_withdrawal_mark_dilution.py` | Hh withdrawal + the dilution race; transient G0 is *less* with the mark (converts to deep arrest). |
| `sim_mark_lifecycle_hysteresis.py` | Entry + steady + withdrawal ±repression with the upgraded machinery (N=2000); the mark's dominant effect is on entry; `--sdmu=0.22` companion. |
| `sim_mark_persistence_sweep.py`, `sim_mark_battery.py`, `sim_h3k27_dilution_race.py`, `sim_mark_cellcycle_recovery.py` | Supporting persistence sweeps, battery of mark tests, dilution-race compute, and recovery. |
| `sim_hh_division_thresholds.py` → `sim_hh_entry_threshold.pdf` / `sim_hh_withdrawal_threshold.pdf` | The Hh threshold for division, repressed vs not — entry delayed (0.32 vs 0.16), exit sooner (0.48 vs 0.15); mark narrows the proliferative window from both sides. |

---

## 8. Withdrawal & EZH2i-rescue kinetics

| Script | What it shows |
|---|---|
| `fig_v44_ezh2i_doseresponse.py` | Graded EZH2i dose-response of CyclinD1 + proliferation in GNP / MB / MB+HHi / MB+CDK4/6i — rescue threshold + no-rescue control. |
| `sim_ezh2i_population_dose.py` → `fig_v44_ezh2i_population_dose.pdf` | Population EZH2i dose-response: % of a heterogeneous tumor pushed back into cycle vs dose. |
| `fig_v44_ezh2i_rescue_kinetics.py` | **EZH2i fast-tracks CyclinD1 de-repression → rescues vismo-arrested MB** (with_h3k27_memory model); timed dosing @0/24/48/72h. (SLOW tier) |
| `fig_v44_feedback_role.py` | ±EZH2⊣CyclinD1 across mitogen: the feedback installs a mitogen-gated quiescence threshold. |
| `sim_mitogen_withdrawal.py` → `fig_v44_mitogen_withdrawal.pdf` | Sudden/gradual mitogen withdrawal × depth × ±feedback on a heterogeneous GNP ensemble; feedback gap depth-dependent, maximal at partial withdrawal. (SLOW, ~2.5h) |
| `sim_vismo_withdrawal.py` → `fig_v44_vismo_withdrawal.pdf` | MB+vismodegib (partial-withdrawal analog) × dose × ±feedback; bimodal EZH2i rescue at near-full vismo. (SLOW, ~2.5h) |
| `fig_matched_withdrawal.py` (`sim_matched_withdrawal.py`) | **Transcript-matched** ON/OFF/EZH2i started at the same Ccnd1 (~1.61, different Shh) then withdrawn at multiple rates — the governor's withdrawal effect is **rate-dependent** (fast=memory-dominated, slow=threshold-dominated). (memory `matched-withdrawal-rate-dependent`) |
| `sim_singlecell_withdrawal_g0.py` | Single-cell ensemble (N=1000, parallel) through withdrawal — recovers the transient G0 & ~2.5× period-lengthening cycle-averaging hides; `--shallow` companion. |
| `sim_ezh2_phaseplane.py` → `fig_v44_ezh2_phaseplane.pdf` | Ferrell-style phase-plane / bifurcation of the EZH2→Mk⊣CyclinD1 feedback. (SLOW, ~35 min) |
| `sim_steadystate_g0_mark.py` | Steady-state transient-G0 & arrest, repressed vs not (disentangles the withdrawal observation from its arrest confound). |
| `sim_inheritance_partition.py` | Prototype: data-calibrated abundance draws + real division-to-division inheritance with partitioning noise. (memory `v44-deterministic-multiperiodic-g0`) |
| `fig_v44_divisions_to_arrest.py` | MB's reduced vismodegib sensitivity: a Gli1 residual (Gli1_epi memory, 2.3% vs GNP 0.4% at 24h) + a resistant proliferating fraction. (memory `v44-gli1-autoregulation-memory`) |

---

## 9. Two-cyclin, CDKI, CDK4/6i & differentiation

| Script | What it shows |
|---|---|
| `fig_cdk46i_dose_response.py` (`sim_cdk46i_dose_response.py`) | CDK4/6i dose-response on a real concentration axis, fit to JP's two palbo points (1µM→16-18% pRb⁺, 5µM→~2%); Langmuir K≈1.5µM. |
| `fig_cdk46i_washout.py` (`sim_cdk46i_washout.py`) | CDK4/6i washout kinetics in MB: 1µM holds ~16-18% residual & rebounds fast; 5µM arrests fully, re-enters slowly; both reversible. |
| `fig_v44_cdk46i_norescue` | (via `fig_v44_fig5_rescue.py`) MB+CDK4/6i ±EZH2i — no rescue (Supp. Fig. 8I,J). |
| `fig_v44_differentiation_collapse.py` | Does the model support a collapse to differentiation, and is EZH2's fall cause or consequence of exit? (EZH2 is E2F-driven → consequence.) |
| `fig_diff_gene.py` (`sim_diff_gene.py`) | Poised bivalent differentiation gene → premature activation when division is too fast or under EZH2i; serial me-chain + high threshold. (memory `ezh2-two-arms-cko-scoping`, developmental arm) |
| `fig_ramp_pop_diffgene.py` (`sim_ramp_pop_diffgene.py`) | GNP mitogen ramp + population WITH the differentiation tripwire. |
| `sim_cdc20_lever.py`, `sim_cdc25a_lever.py`, `sim_checkpoint_diagnostics.py`, `sim_hu_dose_response.py`, `sim_stress_ezh2.py` | CDC20/CDC25A lever probes, checkpoint diagnostics, HU dose-response, EZH2-under-stress (S-phase/HU calibration support). |
| `sim_recalibrate_hill.py`, `sweep_mycn_hill.py` | Hill / MYCN recalibration search support (two-cyclin & MYCN elevated-expression era). |

> **Note on the CDKI era:** the newest baked default is the **dynamic-CDKI** model (INK4 p16/p18/p19 + CIP/KIP p21/p57 as species, CDK2 sequestration). The WS1/WS2 figures above run against it. See [Evolution](03_evolution_timeline.md).

---

## 10. Population / distribution-model figures

Virtual-cell populations where each cell draws a few decisive parameters from log-normal distributions. **Only the CV/shape is used, never the absolute means.** (memory `distribution-model-thread`)

| Script | What it shows |
|---|---|
| `sim_distribution_rescue.py` | Fractional EZH2i rescue across three measured-variation axes (CyclinD1 CV 0.70, birth-p27 CV 0.33, EZH2 CV ~0.55, MB55 IF); rescue phase diagram in (CyclinD1 × EZH2). |
| `sim_cyclind1_variance_buffering.py` | EZH2 buffers cell-to-cell CyclinD1 variance; prediction: EZH2i broadens the single-cell CyclinD1 CV, mainly the high tail. |
| `sim_cyclind1_variance_calibrate.py` | Calibrated variance-buffering prediction: σ*=0.73 reproduces the measured 0.70; without-EZH2 CV ~0.82, buffering factor ~1.18×. |
| `sim_p27_bimodality.py` | Bimodal G0 p27 emerges from the transient-G0/immediate-re-entry switch (~4.7× apart, matching the measured ~4×). |
| `sim_memory_population_fold.py` | Population-level fold from the memory mechanism. |
| `sim_gnp_developmental.py` → `sim_gnp_developmental.pdf` + `sim_gnp_developmental_ezh2.pdf` | The whole GNP developmental proliferative period ±mark at matched output (mean ~8 divisions); the mark re-sets the Hh operating point + sharpens the developmental shutdown; companion EZH2 panel. |

### Parameter-landscape figures (JP 2026-07-27)
| Script | What it shows |
|---|---|
| `fig_param_heatmap.py` | Heatmap of all **202** v44 parameters grouped by module, colored by log10(value); the 6 cell-type-varying params show both GNP and MB. |
| `fig_param_settings.py` | Per-parameter provenance (measured/inherited/soft/free), tested range + landed value, GNP\|MB, and a sensitivity strip (elasticity of each of the 29 targets). |
| `fig_param_biplot.py` | SVD biplot of the parameter→target sensitivity matrix — the model's major axes of variation; aligned vs antagonistic parameters. |
| `fig_param_distributions.py` | The three per-lineage input draws + the emergent cell-cycle-duration distribution (Nakashima 2015: mean 15.9h, CV 0.23). |
| `fig_v44_cd_ezh2_commitment.py` / `fig_v44_cdmrna_ezh2_commitment.py` | The G0/G1 commitment in the (CyclinD1 drive × EZH2 repression) plane (and its transcript version) — unifies withdrawal and EZH2i-rescue as crossings of the same boundary; CDK4/6i is off-plane (no rescue). |

---

## 11. Calibration / diagnostic scripts (no figure)

| Script | Purpose |
|---|---|
| `validate_v44.py` | Unified validation harness — all target checks + figure-critical behaviors (run to score the model). |
| `v44_calibrate_hh.py`, `v44_optimize_hh.py` | HH/MYCN→CyclinD1 between-condition ratio calibration. |
| `v44_calibrate_ezh2.py` | EZH2 within-cycle gradient + HU-boost calibration. |
| `v44_calibrate_final.py`, `v44_probe_sets.py` | Constrained whole-model parameter search (figure-safe guard). |
| `probe_transient_g0_v2.py`, `analyze_cyclind1_g0.py`, `probe_two_step.py` | Transient-G0 / CyclinD1 / two-step-Rb probes. |
| `sim_punctate_bistable_proto.py`, `sim_punctate_meanfield_proto.py` | Punctate-chromatin bistable / mean-field prototypes (mark-topology exploration). |

---

## 12. Archived: v45 stochastic-commitment suite

The v45 stochastic-commitment exploration (model, validator, and its 10-figure suite) is in **`archive/v45/`** (`archive/v45/README.md`). **v44 is the committed working model**; its figures are the ones catalogued above. (memory `v45-stochastic-commitment`, `v45-feedback-paramspace`.)

---

### Reading the catalog

- **Superseded / cautionary figures are kept, not deleted** (git-preserved): `fig_v44_g0_hysteresis` (false-positive bistability), `fig_v44_param_exploration` (coarse), `fig_v44_mitogen_rampdown` (removed duplicate). This is deliberate — the dead-ends are part of the record. See [Things Tried and Abandoned](05_things_tried_and_abandoned.md).
- Figures whose conclusions have flipped across the project (e.g. EZH2 concentration convention, the "EZH2-maintained lock-in" framing → monostable G0) carry in-file notes; trace them via [Evolution](03_evolution_timeline.md).
- For the underlying model structure that every `fig_v44_*` runs against, see [Model Architecture](02_model_architecture.md); for where the scripts live and how to run them, [Repos, Directories, Reproduction](09_repos_directories_reproduction.md).
