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
