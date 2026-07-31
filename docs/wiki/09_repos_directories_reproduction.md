# Repositories, Directories & Reproduction

> The practical navigation guide to the whole project: the three local directories, the two GitHub repos, what lives in each folder, and exactly how to build, validate, and regenerate figures. If you are new (or returning after months away), start here, then jump to [Architecture](02_model_architecture.md), [Evolution](03_evolution_timeline.md), and [Calibration & Validation](06_calibration_and_validation.md).

---

## TL;DR — where do I go?

| I want to… | Go to |
|---|---|
| Run the current model | `Ezh2_CyclinD1_Sphase/` (branch `transient-g0-reframe`), `src/build_model_v44_heldt.py` |
| Validate the current model | `./venv/bin/python simulations/validate_v44.py` |
| Regenerate all figures | `bash simulations/regenerate_figures.sh` |
| Understand where a number came from | `docs/DATA_PROVENANCE.md`, `docs/PARAMETERIZATION.md` |
| See the published-baseline (paper) model | `src/build_model_v42_mycn.py` (frozen; tag `v42-published-baseline`) |
| Read the foundational BioModels experiments | `ezh2_cyclind1_sym/` (Gerard, Novak, Heldt, Tyson, Yao…) |
| Read the distilled "what we tried" | `~/.claude/.../memory/MEMORY.md` (60 files) |

---

## 1. The three local directories

All three live under `/Users/jpurzner/Dropbox/Q_research/py_projects/`.

### 1a. `Ezh2_CyclinD1_Sphase/` — MAIN ACTIVE REPO

- **Size / history:** ~798 MB, 233 commits.
- **Git remote:** `https://github.com/jpurzner/Ezh2_CyclinD1_Sphase.git` (GitHub repo #1, see §2).
- **Current checked-out branch:** `transient-g0-reframe` (this is where day-to-day v44 work happens — *not* `main`).
- **What it is:** the v42 → v43 → v44 lineage plus the CDKI / transient-G0 / two-cyclin / H3K27me3 work. This is the model you run today.
- **The current model:** `src/build_model_v44_heldt.py` — a real-time (minutes, no global clock) rebuild on the **Heldt et al. 2018 (PNAS)** explicit-replication cell-cycle core (`BIOMD0000000700`), with the EZH2 / Hedgehog / H3K27me3 / MYCN biology on top. Entry point `build_model_v44(...)` (signature begins at line 450).

> ⚠️ **Doc-vs-reality drift on the validation score.** `README.md` still says "22/27 targets" and describes the AUM dilution module as the default. The repo has moved well past that: the formal target spec is now **cell-type-separated with ~29–32 targets** ([validation spec](06_calibration_and_validation.md)), and the most recent commits report **28/29 → 30/32** as modules were baked in. Trust the git log and `simulations/validate_v44.py` over the README. Recent commits: `Bake dynamic-CDKI module as the DEFAULT (30/32)`, `Add CyclinD2 as a separate D-cyclin (two-cyclin model), 30/32`, `Bake the P7 GNP / SHH-MB cell-type split as the model default`.

### 1b. `ezh2_cyclind1_sym/` — FOUNDATIONAL / OLD REPO

- **Size:** ~864 MB (largest of the three; contains published-model SBML + PDFs).
- **Git remote:** `https://github.com/jpurzner/Ezh2_Hh_Ccnd1_model.git` (GitHub repo #2). Note the local **directory name (`sym`) differs from the GitHub name (`Ezh2_Hh_Ccnd1_model`)** — same repo.
- **What it is:** the v17 → v31 foundational work, built by importing published BioModels and experimenting with which cell-cycle engine to adopt. This is where the Gerard–Goldbeter engine (v42's core) was chosen. See [Origins & Foundations](01_origins_and_foundations.md).
- **Published models present (SBML/CPS/OMEX):** `Gerard2009.xml` / `Gerard2010.xml`, `Novak2004.xml` (+ `.cps`, `_no_mass_event.xml`), `Tyson2001.xml`, `Yao2008.xml`, `yao2008_rb_e2f.xml`, `swat2004_rb_e2f.xml`, `Abroudi2017.xml`, `Schwarz2018.xml`, `schoeberl2002_egf_mapk.xml`, `MODEL1812130001.4.omex`, plus the `CyclinD_EZH2_module.xml` prototype.
- **Design docs (the v17–v31 paper trail):** `BIOMODELS_REVIEW_REPORT.md`, `EZH2_CyclinD1_Model_Integration_Spec.md`, `MODEL_DEVELOPMENT_SUMMARY.md`, `MODEL_DOCUMENTATION_FOR_PAPER.md`, `SOURCES_AND_CITATIONS.md`, `MINIMAL_OSCILLATOR_OPTIONS.md`, `LIVE_IMAGING_INSIGHTS.md`, `IMPLEMENTATION_GUIDE.md`, the `V28`/`V29`/`V30`/`V31` status files, `v17_SUCCESS_SUMMARY.md`, GLI-module wireframes, and analysis scripts (`analyze_biomodels.py`, `deep_analysis.py`).
- **Source papers (PDFs):** Gérard–Goldbeter 2009, Abroudi 2017, `msb4100126.pdf`, `s41420-025-02605-7.pdf`, `s41420` etc.
- **Note:** the README says use system `python3` here (tellurium was fragile in venvs at the time). This repo is largely historical — do not do new modeling here.

### 1c. `ezh2_cyclind1_scheckpoint/` — SMALL PLANNING DIR

- **Size:** ~36 KB. Not a code repo — a planning/checkpoint stub.
- **Contents:** `checkpoint_model_plan.md` (the AMOS intra-S checkpoint + Ezh2–CyclinD1 feedback extension plan — the design doc that spawned v43/v44; identical to the copy in the main repo's `CHECKPOINT_PLAN.md`), plus two equation Word docs: `Ptch1_mRNA_equation.docx` and `Smo_active_equation.docx` (the Hedgehog-input equations JP drafted for the realistic-mitogen work — see [Mechanisms](04_mechanisms_explored.md)).
- **This directory is also the CWD anchor for the Claude Code project** whose memory system lives under `~/.claude/projects/-Users-jpurzner-...-scheckpoint/memory/` (§4).

---

## 2. The two GitHub repos

| GitHub repo | Local dir | Role |
|---|---|---|
| **`jpurzner/Ezh2_CyclinD1_Sphase`** | `Ezh2_CyclinD1_Sphase/` | **Active.** The v42→v44 S-phase/checkpoint/EZH2 work. Push here. |
| **`jpurzner/Ezh2_Hh_Ccnd1_model`** | `ezh2_cyclind1_sym/` | **Old foundation.** The v42 model submitted with the Chahin paper; frozen baseline. |

> **Stale-clone gotcha (from project memory):** a GitHub clone of the Sphase repo was found to lag far behind the local working copy. **Treat the local `Ezh2_CyclinD1_Sphase/` as the source of truth**, not a fresh `git clone`. See memory note `github-clone-is-stale-use-local.md`.

### Branches in the active repo

```
* transient-g0-reframe                    <- current working branch
  main
  feature/intra_s_checkpoint
  feature/v44-validation-figures-twostep
  v44-explicit-replication
  remotes/origin/{main, transient-g0-reframe, ezh2i-catalytic-fix,
                  structure-analysis, feature/intra_s_checkpoint,
                  feature/v44-validation-figures-twostep, v44-explicit-replication}
```

### Tags (recoverable milestones)

- `v42-published-baseline` — the frozen paper model.
- `v43-milestone` — intra-S checkpoint via phase-gated clock (12/12).
- `v44-aum-default` — AUM H3K27me3 dilution module made default.
- `v44.1-prc2-readwrite` — PRC2 read-write chain.

See [Evolution Timeline](03_evolution_timeline.md) for the full version story.

---

## 3. What lives where — the `Ezh2_CyclinD1_Sphase/` layout

```
Ezh2_CyclinD1_Sphase/
├── README.md              # overview (validation score STALE — see §1a note)
├── CHECKPOINT_PLAN.md     # the intra-S checkpoint + Ezh2 feedback extension plan
├── requirements.txt       # tellurium>=2.2, numpy, scipy, matplotlib, pandas, sympy, networkx
├── .gitignore             # excludes venv/, large caches
├── venv/                  # the pinned Python 3.9 environment (use venv/bin/python)
│
├── src/
│   ├── build_model_v44_heldt.py   # CURRENT model (~91 KB). build_model_v44(...) @ line 450
│   ├── build_model_v43_checkpoint.py  # intra-S checkpoint (phase-gated clock), 12/12
│   ├── build_model_v42_mycn.py    # FROZEN published baseline (do NOT modify; keeps its own GDC0449)
│   └── __init__.py
│
├── simulations/           # ~661 files: validators, figure scripts, probes, ensembles, caches
│   ├── validate_v44.py    # THE current validation suite (env-driven; see §5)
│   ├── validate_v43.py / validate_v42.py
│   ├── regenerate_figures.sh   # canonical figure-regen manifest (FAST glob + SLOW list)
│   ├── fig_v44_*.py       # ~75 figure scripts (rescue, phase-plane, withdrawal, g0, ...)
│   ├── sim_*.py / probe_*.py   # ensemble + mechanism experiments
│   ├── calibrate_h3k27_*.py    # H3K27me3 module calibration + parameter search
│   ├── *_results.json / *.npz  # cached sim outputs + figure caches
│   └── *.png / *.pdf      # generated figures (png preview + pdf for paper)
│
├── docs/                  # ~50 markdown design/state docs + the writeup PDF
│   ├── model_writeup.pdf / .md          # full walkthrough (every equation + all figures)
│   ├── v44_MODEL_DESCRIPTION.md, v44_status.md, v44_verification_findings.md
│   ├── v44_recalibration_and_search_results.md   # the parameter bake + wide search
│   ├── DATA_PROVENANCE.md                # audit trail: data -> param -> MEASURED/INHERITED/SOFT/FREE
│   ├── PARAMETERIZATION.md, PARAMETER_IDENTIFIABILITY_PLAN.md
│   ├── validation_targets.md            # the target compendium (companion to data/validation_targets.json)
│   ├── H3K27me3_CyclinD1_literature_review.md, h3k27me3_repression_model_state.md
│   ├── FIGURES.md                       # which script makes each figure
│   ├── (many dated exploration docs: cyclind1_hh_memory_, realistic_mitogen_input_,
│   │    speed_limiter_framework_, prc2_buffering_memory_, g1_lengthening_buffering_, ...)
│   └── wiki/            # 00_home.md (this wiki lives here)
│
├── data/
│   ├── validation_targets.json    # CANONICAL machine-readable targets (source of the doc tables)
│   ├── AvgKdegs_genes_v1.1.csv     # RNAdecayCafe transcript half-life atlas (CCND1 t1/2 work)
│   ├── RNAdecayCafe_database_v1.1.rds
│   └── skp2_timecourse.json
│
├── models_external/       # the published core, as imported
│   ├── heldt2018.ant              # Antimony source of the v44 engine
│   └── heldt2018_BIOMD700.xml     # BIOMD0000000700 SBML
│
├── archive/v45/           # the v45 stochastic-commitment successor (parallel, COMPLETE 20/20)
│   ├── README.md
│   └── fig_v45_*.py/.png/.pdf     # commitment bifurcation, feedback param-space, rescue, ...
│
├── deprecated/            # retired scripts (e.g. fig_prc2_transcription_paradox.py)
└── papers/                # source PDFs (ncomms6425.pdf, ...)
```

**Notes on the sub-directories that trip people up:**

- **`archive/v45/`** is not an old version — it is a *parallel successor* (bistable CDK2–p27 toggle, stochastic commitment) that reached 20/20 but was set aside; **v44 remains the working model**. See memory `v45-stochastic-commitment.md` and [Evolution](03_evolution_timeline.md).
- **`src/build_model_v42_mycn.py`** is byte-frozen and uses the *old* drug name `GDC0449`; v44 renamed the Hedgehog inhibitor input to **`HHi`**. Never modify v42.
- **`simulations/` is huge (~661 files)** because it accumulates result JSONs and `.npz` figure caches alongside the ~75 `fig_v44_*.py` scripts. Many SLOW figures are cached; regen scripts pass `--fresh` to force recompute (see §6).

---

## 4. The project memory system (`~/.claude/.../memory/`)

Path: `~/.claude/projects/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/memory/`

This is the **distilled "what we tried" over ~4 months** — 60 markdown files (59 notes + `MEMORY.md` index), persisted across Claude Code sessions. **Read `MEMORY.md` first** (it is a dense linked index), then the specific note you need. Highlights:

- **Baked defaults & their evidence:** `v44-four-mark-mechanisms.md`, `v44-cyclind1-hh-memory-carrier.md`, `ezh2-governor-reservoir-remeasure.md`, `v44-skp2-mitogen-dose.md`, `dynamic-cdki-module.md`, `mb-celltype-transient-g0-parameterization.md`, `mycn-elevated-expression-not-amplified.md`, `cyclind1-d2-two-cyclin-model.md`.
- **Dead-ends / negative results:** `recalibration-16h-cycle-attempt.md`, `v44-latch-regime-consequences.md`, `v44-mother-g2-negative-result.md`, `withdrawal-graded-g1-lengthening.md`, `v44-chromatin-memory-validation-bound.md`, `withdrawal-loop-overdamped-governor.md`.
- **Structural maps & gotchas:** `v44-structural-map-3tier.md`, `roadrunner-reset-param-gotcha.md`, `github-clone-is-stale-use-local.md`, `validation-targets-formal-spec.md`.

This memory is the human-readable backbone of [Mechanisms Explored](04_mechanisms_explored.md) and [Things Tried & Abandoned](05_things_tried_and_abandoned.md).

---

## 5. Reproduction — build & validate

The environment is a **pinned Python 3.9 venv** already present in the repo. Always call it explicitly (`./venv/bin/python`) rather than a system Python.

```bash
cd /Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase

# (first-time setup only; the venv is usually already present)
python3.9 -m venv venv && ./venv/bin/pip install -r requirements.txt

# validate the current model (writes simulations/validation_v44_results.json)
./venv/bin/python simulations/validate_v44.py

# inspect / print the model
./venv/bin/python src/build_model_v44_heldt.py
```

**Dependencies** (`requirements.txt`): `tellurium>=2.2` (→ Antimony + libRoadRunner, the ODE engine), `numpy>=1.24`, `scipy>=1.10` (`find_peaks` for period detection), `matplotlib>=3.7`, `pandas>=2.0`, `sympy>=1.14`, `networkx>=3.2`. Integrator is CVODE (stiff) via roadrunner — deterministic, so results reproduce numerically modulo roadrunner version.

### Environment-variable flags on `validate_v44.py`

The validator (and the builder) are driven by env vars so you can toggle modules without editing code. The current **baked defaults** are shown in parentheses (from `simulations/validate_v44.py` lines 59–71):

| Env var | Default | Effect |
|---|---|---|
| `H3K27_DILUTION` | `1` | mean-field AUM H3K27me3 dilution repression module (`0` = legacy direct repression) |
| `H3K27_CHAIN` | `1` | serial me1/me2/me3 methylation chain (`with_h3k27_chain`) |
| `TWO_STEP_RB` | `1` | two-step Rb hyperphosphorylation (G0 marker) |
| `DECOUPLE_COMMIT` | `1` | **BAKED 2026-07-27** — G0/commitment gated on CyclinD1/CDKi, not cell size (`0` = legacy size-gated) |
| `CDKI_SPECIES` | `1` | **BAKED 2026-07-30** — individual dynamic CDKI species (p16/p18/p19/p21/p57) (`0` = legacy lumped CDKI) |
| `EZH2_CONC` | `0` | legacy EZH2 concentration convention (default is the dilution/Cdh1-boost convention) |
| `MITOGEN_TRACKER` | `0` | low-pass mitogen tracker (off) |
| `MYCN_AUTOREG` | `0` | MYCN from Gli1 + bistable self-activation; deprecated/inert (SHH-MB is not MYCN-amplified) |
| `VALIDATE_PARAMS` | `{}` | JSON param overrides — the driver used by the wide parameter search |
| `PRB_G0_THR`, `P18_MB`, `P19_MB`, `P21_DIV_MB` | see file | tunable G0 / CDKI target thresholds |
| `VALIDATE_RESULTS` | (auto) | output path for the results JSON |

Example — run with the legacy lumped-CDKI model and a param override:

```bash
CDKI_SPECIES=0 VALIDATE_PARAMS='{"k_Cd_tx_basal":0.5}' \
  ./venv/bin/python simulations/validate_v44.py
```

The builder call inside the validator is:
```python
build_model_v44(with_ezh2=True, with_hh=True,
                with_h3k27_dilution=..., with_h3k27_chain=..., with_two_step_rb=...,
                decouple_commit=..., with_cdki_species=..., ...)
```

---

## 6. Reproduction — figures

```bash
bash simulations/regenerate_figures.sh    # run from repo root; calls ./venv/bin/python internally
```

`regenerate_figures.sh` is the **canonical figure manifest**. It runs in two tiers:

- **FAST** — auto-globs every `simulations/fig_v44_*.py` (per-condition figures). New `fig_v44_*` scripts are picked up automatically, no edit to the manifest needed.
- **SLOW** — an explicit list of ensemble/analysis figures that take minutes-to-hours and are **cached** (`.npz`). The manifest passes `--fresh` to force recompute against a changed model (otherwise they silently reuse stale caches). Examples and their runtimes are documented inline: `sim_g0_bifurcation.py` (~10 min, N=140), `sim_population_g0.py --n=100` (~12 min), `sim_ezh2_phaseplane.py --fresh` (~35 min), `sim_mitogen_withdrawal.py --fresh` (~2.5 h), `sim_vismo_withdrawal.py --fresh` (~2.5 h).

Each script writes `.png` (preview) + `.pdf` (paper) into `simulations/`. `docs/FIGURES.md` maps each figure to its script; the figure catalog is on [Figures Catalog](08_figures_catalog.md).

**Run figure scripts individually** the same way:
```bash
./venv/bin/python simulations/fig_v44_fig5_rescue.py
./venv/bin/python simulations/sim_g0_bifurcation.py --fresh
```

> **Cache/gotcha reminders:** (1) always run from the repo root — scripts resolve `src/` and write outputs via paths relative to root; (2) `roadrunner reset()` resets *species, not parameters* — sweeps must re-apply a full baseline each call or values leak (memory `roadrunner-reset-param-gotcha.md`); (3) cached ensemble figures need `--fresh` after any model change.

---

## 7. Reproducing the old models

**v42 published baseline** (`Ezh2_CyclinD1_Sphase/`, tag `v42-published-baseline`):
```bash
./venv/bin/python simulations/validate_v42.py     # the 10-condition drug/rescue panel
./venv/bin/python src/build_model_v42_mycn.py      # print the model
```

**Foundational BioModels work** (`ezh2_cyclind1_sym/` = GitHub `Ezh2_Hh_Ccnd1_model`): use system `python3` (README caveat: tellurium was fragile in venvs there). The published SBML models sit at the repo root (`Gerard2009.xml`, `Novak2004.xml`, `Yao2008.xml`, etc.); walkthrough notebooks (`Gerard2009_only_walkthrough.ipynb`, `Novak_Ezh2.ipynb`) and analysis scripts (`analyze_biomodels.py`, `deep_analysis.py`) reproduce the model-selection experiments. See [Origins & Foundations](01_origins_and_foundations.md) and [Data & Evidence](07_data_and_evidence.md).

---

## 8. Cross-references

- Version story and why v44 replaced v42's engine: [Evolution Timeline](03_evolution_timeline.md)
- What each module does and how they wire together: [Model Architecture](02_model_architecture.md)
- The full validation-target spec and scores: [Calibration & Validation](06_calibration_and_validation.md)
- Every figure and its generating script: [Figures Catalog](08_figures_catalog.md)
- Dead-ends catalogued from the memory notes: [Things Tried & Abandoned](05_things_tried_and_abandoned.md)
- Open items and reviewer prep: [Reviewer & Open Questions](10_reviewer_and_open_questions.md)
