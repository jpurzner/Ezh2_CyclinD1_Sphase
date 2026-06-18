# Model figures — which code generates each

All scripts live in `simulations/` and run against the calibrated **v44 single-step default**
(`src/build_model_v44_heldt.py`). Each writes a `.png` (preview) and a `.pdf` (vector, for the paper)
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
| `fig_v44_g0_bifurcation.pdf` | `sim_g0_bifurcation.py` | Proliferation–quiescence bifurcation from cyclin D1 / **birth-p27** (`P21_div`) heterogeneity (N=140). Low birth p27 → immediate re-entry, high → prolonged transient G0 (the ultrasensitive Fan–Meyer threshold; mapped by `probe_birthp27_dwell.py`); low CyclinD1 → arrest. SHH and EZH2i tune the split (immediate ↑, arrest ↓, dwell ↓). On the baked model. (~10 min) |
| `fig_v44_mitogen_withdrawal.pdf` | `sim_mitogen_withdrawal.py` | **Mitogen-withdrawal experiment on the heterogeneous GNP ensemble (N=16, birth-p27 + CyclinD1-setpoint draws).** Cells equilibrated cycling at high SHH, then mitogen withdrawn SUDDENLY or GRADUALLY to a range of final levels; counts post-withdrawal divisions ± the EZH2⊣CyclinD1 feedback (no-feedback = K_EZH2_repression→∞). (A/B) divisions vs final mitogen, sudden/gradual — the **feedback gap is depth-dependent, maximal at partial withdrawal (SHH 0.15–0.3)** and closes at full withdrawal (both arrest) and at baseline (both cycle); sudden > gradual. (C) per-cell distributions at SHH→0.3: sudden +fb 2.4 vs −fb 4.6; gradual +fb 3.3 vs −fb 4.5. (D) representative cell: no-feedback CyclinD1 stays high → coasts (5 div), feedback represses CyclinD1 → arrests (1 div). State-reuse + tol-retry; `--probe` / `--fresh`. (~2.5 h) |
| `fig_v44_ezh2i_rescue_kinetics.pdf` | `fig_v44_ezh2i_rescue_kinetics.py` | **EZH2i fast-tracks CyclinD1 de-repression → rescues vismo-arrested MB** (uses the calibrated `with_h3k27_memory` model — an explicit H3K27me3 mark at the CyclinD1 locus). Vismo lowers Gli-driven CyclinD1 but the MYCN floor remains; sustained high EZH2 holds the mark → stable EZH2-maintained arrest. EZH2i blocks deposition → mark cleared by demethylation + replication-dilution → CyclinD1 de-represses over ~24–48h → cycle re-entry. (A) CyclinD1 ± EZH2i; (B) the H3K27me3 mark (mechanism); (C) timed dosing @0/24/48/72h — de-represses whenever given (the fast-track); (D) rescue outcome (5–6 div vs vismo-only's 1) on-demand. NB the rescue is MB/vismo-specific: in mitogen-withdrawn GNPs EZH2 falls so the mark self-clears (arrest not EZH2-maintained). `--fresh`. (~15 min) |
| `fig_v44_vismo_withdrawal.pdf` | `sim_vismo_withdrawal.py` | **Vismodegib treatment of the heterogeneous MB ensemble — the MB analog of mitogen withdrawal.** Vismo (GDC0449) blocks Smo, removing Gli-driven CyclinD1 but leaving the MYCN floor → even full vismo is a *partial* withdrawal, so the EZH2 feedback is decisive (= the model's vismo-arrest / EZH2i-rescue). MB cells equilibrated cycling, then vismo applied SUDDENLY/GRADUALLY to a range of doses; counts post-treatment divisions ± feedback (no-fb = K_EZH2_repression→∞ ≡ EZH2i/tazemetostat). (A/B) divisions vs dose — feedback gap opens at high vismo. (C) at near-full vismo (GDC=0.95): +fb → uniform arrest (1 div); −fb (EZH2i) → **bimodal rescue** (high-CyclinD1-setpoint cells coast 4–5, low ones arrest = "rescues a fraction"). (D) representative cell: no-fb MYCN-floor CyclinD1 sustains cycling, fb represses → arrest. Capped at GDC=0.95 (GDC=1.0 = Smo singularity, spurious arrest; vismo never 100%-blocks Smo). State-reuse + tol-retry; `--probe`/`--fresh`. (~2.5 h) |
| `fig_v44_ezh2_phaseplane.pdf` | `sim_ezh2_phaseplane.py` | **Ferrell-style phase-plane / bifurcation analysis of the EZH2 ⊣ CyclinD1 feedback.** (A) EZH2–CyclinD1 nullcline portrait (GNP): Cd-nullcline analytic (Cd fast, EZH2 represses it ↓), EZH2-nullcline traced by breaking the loop (EZH2i=1, sweep drive) (↑); they meet at one stable node (monostable, homeostatic) — the negative feedback buffers Cd ~2.2× (no-feedback 4.1 → operating 1.85); fast/slow direction field + limit-cycle orbit overlaid. (B) GNP/MB/+EZH2i — mitogen drive & feedback removal move the operating point. (C) bifurcation vs CyclinD1 drive: feedback shifts the proliferation threshold +0.24 higher (mean-E2F order parameter). (D) bifurcation vs mitogen (Gli1): feedback raises the mitogen threshold for proliferation (gated quiescence). Calibrated by `probe_ezh2_phaseplane.py`; two-stage npz cache; `--fresh`. (~35 min) |
| `fig_v44_mitogen_sensitivity.pdf` | `fig_v44_mitogen_sensitivity.py` | **mitogen sensitivity in Hh-activity (Gli1) units**, with the correct knobs: (A) **GNP — SHH** is the mitogen (sweep up; vismo only extends the axis below the SHH=0 basal floor); (B) **MB — vismodegib** titrates the tonically-high Hh DOWN toward the MYCN floor (MB CKI tones on). Each: CyclinD1 vs Gli1 ± the EZH2⊣CyclinD1 feedback — the feedback represses/flattens CyclinD1 (lower mitogen sensitivity) and raises the cycling threshold (GNP +fb Gli1~0.008, MB +fb ~0.015; arrest region shaded). NB the basal floor still sits close to threshold (small sub-threshold range). Caches; `--fresh`. (~12 min) |
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

## v45 figures — ARCHIVED

The v45 stochastic-commitment exploration (model, validator, and its 10-figure suite) has been
moved to **`archive/v45/`** (see `archive/v45/README.md`). **v44 is the committed working model**;
its figures are the ones listed above.
