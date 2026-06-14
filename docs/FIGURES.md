# Model figures — which code generates each

All scripts live in `simulations/` and run against the calibrated **v44 single-step default**
(`src/build_model_v44_heldt.py`). Each writes a `.png` (preview) and a `.pdf` (vector, for the paper)
into `simulations/`. Regenerate any figure with, e.g.:

```bash
./venv/bin/python simulations/fig_v44_fig5_rescue.py
```

The unified validation harness `simulations/validate_v44.py` defines the conditions, the division
counter, the phase (p27) classifier, and the experimental-target checks that the figure scripts import.

---

## Paper figures (model)

| Manuscript panel(s) | Output PDF | Generating script | Notes |
|---|---|---|---|
| **Fig. 5A–C** — Ptch+/− MB / +HHi (arrest) / +HHi+Ezh2i (rescue), Cyclin B traces | `fig_v44_fig5_rescue.pdf` | `fig_v44_fig5_rescue.py` | |
| **Supp. Fig. 7A** — model architecture / wiring diagram | `fig_v44_wiring.pdf` | `fig_v44_wiring.py` | pure schematic (no simulation) |
| **Supp. Fig. 7B–M** — validation grid (Cyclin B traces; CyclinD1/Mycn/Gli1 model‑vs‑data; EZH2; division counts; time courses) | `fig_v44_validation.pdf` | `fig_v44_validation.py` | |
| **Supp. Fig. 8A–H** — EZH2→CyclinD1 feedback (±feedback: Cyclin B, CyclinD1 mRNA, EZH2 protein, phase relationship) | `fig_v44_feedback.pdf` | `fig_v44_feedback.py` | feedback OFF = `K_EZH2_repression = 1e6` |
| **Supp. Fig. 8I,J** — MB+CDK4/6i and MB+CDK4/6i+Ezh2i (no rescue) | `fig_v44_cdk46i_norescue.pdf` | `fig_v44_fig5_rescue.py` | emitted alongside Fig. 5 |
| **Supp. Fig. 7x** — Gli→Ptch1 negative feedback & its disruption in MB (Gli1 step overshoot GNP vs MB; Ptch1 mRNA induction; Smo/Gli broken‑loop) | `fig_v44_ptch1_feedback.pdf` | `fig_v44_ptch1_feedback.py` | Ptch1 = Gli target; `Ptch1_copy_number` = functional fraction (1.0 GNP / 0.1 MB) |
| **Fig. 5C** — population cycling fraction (vismo↓, Ezh2i restores; CDK4/6i no rescue) | `fig_v44_fig5_population.pdf` | `fig_v44_fig5_population.py` | N=120 ensemble, CyclinD1/p27 heterogeneity |
| **Conditions → input parameters + parameter space** (A: input-knob matrix per condition; B: CyclinD1×p16 plane with cycle/arrest boundary + each condition placed) | `fig_v44_condition_paramspace.pdf` | `fig_v44_condition_paramspace.py` | shows MB+HHi arresting at the *same* CyclinD1 a GNP cycles at, because p16 raises the threshold; EZH2i crosses it back |

Drop‑in replacement Results/Methods/captions for the manuscript: `simulations/PAPER_MODEL_SECTIONS_v44.md`.
Data sources used to parameterize the model: `docs/PARAMETERIZATION.md`.

---

## Exploratory analyses (developed while interrogating the model; not current manuscript figures)

| Output PDF | Generating script | What it shows |
|---|---|---|
| `fig_v44_ezh2i_doseresponse.pdf` | `fig_v44_ezh2i_doseresponse.py` | Graded EZH2i dose‑response of CyclinD1 (transcript+protein) and proliferation in GNP / MB / MB+HHi / MB+CDK4/6i — the rescue **threshold** and the no‑rescue control. |
| `fig_v44_ezh2i_population_dose.pdf` | `sim_ezh2i_population_dose.py` | Population EZH2i dose‑response: % of a heterogeneous (cyclin D1/p27) tumour pushed back into cycle vs dose — the "low‑dose rescues a fraction" curve. (~15 min) |
| `fig_v44_feedback_role.pdf` | `fig_v44_feedback_role.py` | WITH vs WITHOUT the EZH2 ⊣ CyclinD1 arm across mitogen: the feedback installs a mitogen‑gated quiescence threshold (cell‑cycle exit). |
| `fig_v44_ckoo_feedback_strength.pdf` | `fig_v44_ckoo_feedback_strength.py` | Is the ~3× CyclinD1 supported by the cKO? Recombination correction of the bulk fold + model feedback‑strength sweep. |
| `fig_v44_g0_mechanism.pdf` | `fig_v44_g0_mechanism.py` | CyclinD1 sets transient‑G0 duration (p27 marker); EZH2i collapses it; arrest below threshold. |
| `fig_v44_cycle_anatomy.pdf` | `fig_v44_cycle_anatomy.py` | One‑cycle anatomy: the transient G0 is the p27‑high window (phospho‑Rb saturates early — why we classify G0 by p27). |
| `fig_v44_g0_bifurcation.pdf` | `sim_g0_bifurcation.py` | Proliferation–quiescence bifurcation from cyclin D1/p27 heterogeneity; SHH and EZH2i tune the split. (~10 min) |
| `fig_v44_mitogen_sensitivity.pdf` | `fig_v44_mitogen_sensitivity.py` | **mitogen sensitivity** on the Heldt‑engine model: sweep SHH (Hh dose) ± the EZH2⊣CyclinD1 feedback. (A) CyclinD1 transcript+protein dose‑response — feedback represses & flattens it (lower sensitivity); (B) EZH2 tracks mitogen (Fig 4H); (C) the feedback **installs** the Hh threshold to cycle (without it, basal CyclinD1 lets GNP cycle even at SHH=0). Caches the sweep; `--fresh` to recompute. (~10 min) |
| `fig_v44_mitogen_paramspace.pdf` | `fig_v44_mitogen_paramspace.py` | **2‑D mitogen‑sensitivity planes**: proliferation (divisions/168 h) over (A) SHH × EZH2‑feedback strength `K_EZH2_repression` and (B) SHH × p16, with the cycle/arrest boundary drawn and the defaults marked. The Hh threshold to proliferate **rises** with both stronger feedback and more CKI. Cached grids; `--fresh` recomputes. (~20–30 min cold) |

---

## Calibration / diagnostic scripts (no figure)

| Script | Purpose |
|---|---|
| `validate_v44.py` | Unified validation harness: all experimental‑target checks + figure‑critical behaviors (run this to score the model). |
| `v44_calibrate_hh.py`, `v44_optimize_hh.py` | HH/MYCN→CyclinD1 between‑condition ratio calibration. |
| `v44_calibrate_ezh2.py` | EZH2 within‑cycle gradient + HU boost calibration. |
| `v44_calibrate_final.py`, `v44_probe_sets.py` | constrained whole‑model parameter search (figure‑safe guard). |
| `probe_transient_g0_v2.py`, `analyze_cyclind1_g0.py` | transient‑G0 vs CyclinD1 probes (p27 marker). |
| `probe_two_step.py` | two‑step‑Rb variant probes (`with_two_step_rb=True`). |

---

## v45 figures (stochastic-commitment model)

The **v45** parallel-successor model (`simulations/v45_stochastic_commitment.py`,
`docs/v45_MODEL_DESCRIPTION.md`) has its own figure suite. Each writes a `.png` (preview) + `.pdf`
(vector) into `simulations/`. They import the v45 module directly (ensemble + `condition_drive` +
`fixed_points` + `ezh2_stats`). Regenerate with, e.g., `./venv/bin/python simulations/fig_v45_rescue.py`.

| Output PDF | Generating script | Shows |
|---|---|---|
| `fig_v45_rescue.pdf` | `fig_v45_rescue.py` | **headline**: quiescence by condition (permanent arrest vs transient G0), the MYCN floor (MB+HHi does not permanently arrest), CyclinD1⊣EZH2 (EZH2i de-represses), and the MB→+HHi→+EZH2i rescue arc |
| `fig_v45_commitment_bifurcation.pdf` | `fig_v45_commitment_bifurcation.py` | the bistable CDK2-p27 toggle nullclines (OFF/separatrix/ON), fixed points vs Rb concentration (the G1 timer), and the commit-or-arrest threshold vs CyclinD1 drive |
| `fig_v45_single_cell.pdf` | `fig_v45_single_cell.py` | single-cell anatomy: a committed cycler (commit→S→G2→divide) vs a mass-capped arrested cell, with mass/RbC, CDK2/p27, DNA/G2, EZH2/CyclinD1 |
| `fig_v45_ezh2_period.pdf` | `fig_v45_ezh2_period.py` | EZH2 by phase (peaks S/G2), EZH2i de-repression of CyclinD1, the period mixture (fast cyclers + G0 tail), and MB count-fractions vs flow data |
| `fig_v45_mitogen_ramp.pdf` | `fig_v45_mitogen_ramp.py` | gradual mitogen increase→decrease, with vs without the EZH2 feedback: the feedback **raises the entry threshold ~3×** (the paper's "increases the mitogen threshold to enter cycle") and represses/buffers CyclinD1. The dynamic-withdrawal "extra divisions" buffer is **weak** in v45 — feedback cells arrest *sooner* (the higher threshold dominates), consistent with v44 (EZH2 feedback doesn't lengthen the period) |
| `fig_v45_cyclind1_decline.pdf` | `fig_v45_cyclind1_decline.py` | CyclinD1 **transcript & protein** under a linear mitogen decline, ± EZH2 feedback (the CyclinD1↔EZH2 module, v45 forms). The feedback **represses** CyclinD1 (lower absolute level) AND, because EZH2 now tracks the falling mitogen, **buffers** its decline (sub-linear, normalized panel); protein lags transcript. Buffering is modest (~2 h on a 30 h half-life), so it doesn't overcome the baseline repression |

Validation harness: `simulations/validate_v45.py` (20/20). **v44 remains the committed working model**
for the manuscript figures above; the v45 figures document the mechanistic successor.
| `fig_v45_feedback_paramspace.pdf` | `fig_v45_feedback_paramspace.py` | **parameter-space** analysis of the CyclinD1↔EZH2 feedback over loop STRENGTH (EZH2 gain) × DELAY (EZH2 stability): mitogen-sensitivity heatmap (homeostasis), overshoot/ringing heatmap, dose-response and regime time-courses. KEY: the measured CyclinD1 MB/GNP fold (5.07×) pins the feedback to a weak/mild regime (sensitivity ~0.81) — v45 sits there; strong-buffering (~36% of swept space) and oscillatory (~13%) regimes are **data-excluded** |
| `fig_v45_feedback_prediction.pdf` | `fig_v45_feedback_prediction.py` | **flipped analysis → falsifiable prediction**: the measured CyclinD1 (5.07×) + EZH2 (2.05×) folds pin a single weak-feedback strength (self-consistency, panel A), which **predicts the mitogen-entry threshold shift ~1.85×** (EZH2 raises the rShh-dose-to-cycle vs EZH2i — a separable experiment, panel B). Discriminator (panel C): the CyclinD1 dose-response *slope* distinguishes weak (steep, data) vs strong (flat) feedback. Structural note: EZH2 only lowers CyclinD1, so EZH2i always ≥ proliferation (the rescue) — "EZH2→extra divisions" is impossible |
| `fig_v45_cutoff_rate.pdf` | `fig_v45_cutoff_rate.py` | does the **rapidity of mitogen withdrawal** matter? A fast cut exposes the slow-EZH2 lag → a transient CyclinD1 **undershoot** (~7%, recovers in ~6 h = 1/kdEZ); a slow cut is quasi-static. Effect is modest (weak feedback); feedback CyclinD1 < no-feedback at every rate, gap largest for a fast cut |
| `fig_v45_vs_v44.pdf` | `fig_v45_vs_v44.py` | **v44 vs v45** behavioral comparison: MB phase fractions (v45 closer to data), period (v44 deterministic ~23 h vs v45 distribution/mixture), and the rescue (v44 single-cell hard-arrests at MB+HHi → 0 div; v45 graded ~55% cycling). See `docs/v45_MODEL_DESCRIPTION.md` §8 for the structural table |
