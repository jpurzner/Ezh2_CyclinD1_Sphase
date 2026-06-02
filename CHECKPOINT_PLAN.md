# AMOS Cell Cycle Model Extension: Intra-S Checkpoint and Ezh2-Cyclin D1 Feedback

## Scientific question

The Chahin et al. model (v42) integrates Hh signaling, MYCN, Ezh2, and the mammalian cell cycle, and demonstrates that EZH2 inhibition can rescue HHi-treated MB cells while CDK4/6 inhibition cannot. The model in its current form treats S-phase duration as fixed and does not include an intra-S checkpoint.

Our recent experimental data show that Ezh2 protein accumulates progressively through the cell cycle, peaks in G2, and is boosted by 1.31-fold in S-phase cells under Hydroxyurea treatment (10 µM, p = 0.005). This implies that the *duration* of S-phase is a determinant of total Ezh2 protein produced per cycle, which in turn should suppress Cyclin D1 more strongly and lengthen the subsequent G0/G1.

This creates a hypothesized self-reinforcing dynamic:

> Longer S-phase → more Ezh2 accumulated → stronger Cyclin D1 repression → less CDK4/6 activity → longer G0/G1

The goal of this work is to extend the v42 model with an explicit intra-S checkpoint, parameterize it against the Hydroxyurea data, and explore the parameter space to determine whether this dynamic produces (a) cell cycle period changes proportional to S-phase extension, (b) bistability between a fast-cycling, low-Ezh2 state and a slow-cycling, high-Ezh2 state, or (c) something else entirely. This bears directly on how GNP differentiation timing is controlled, since differentiation requires a sustained drop in Ezh2.

## Repository strategy

The current `Ezh2_Hh_Ccnd1_model` repo is being submitted with the Chahin paper and should not be modified during this work. The plan:

1. Fork `jpurzner/Ezh2_Hh_Ccnd1_model` into a new repo. Suggested name: `jpurzner/cellcycle_checkpoint_ezh2` (final naming up to JP).
2. Create a branch from the current v42 state called `feature/intra_s_checkpoint`. Tag the parent commit clearly so the published model state is recoverable.
3. All new work happens in the fork. Periodically pull from the original if hotfixes to v42 occur before submission.

## Current model state to inherit

The v42 model lives at `src/build_model_v42_mycn.py` and contains 38 cell cycle species adapted from Gerard and Goldbeter 2009, plus added modules for Hh signaling, MYCN, and Ezh2. Key species relevant to the checkpoint extension:

- CycE/CDK2 (active and inactive forms)
- Pe (Cdc25E phosphatase, the functional equivalent of CDC25A acting on Cyclin E/CDK2 in this model)
- pRB, pRBp, pRBpp phosphorylation states
- Free E2F
- p27/Kip1
- Cyclin D1 mRNA and protein (central integration node)
- EZH2 mRNA and protein

The new checkpoint module will target Pe activity and CycE/CDK2 dynamics directly, since these are the canonical points where the intra-S checkpoint suppresses progression in published models.

## Implementation plan

### Phase 1: Repo setup and module identification

Fork the existing repo. Create the new branch. Confirm v42 reproduces its 16/18 validation tests in the new repo before any modification (run `simulations/validate_v42.py` to confirm).

Search the BioModels Database (ebi.ac.uk/biomodels) for cell cycle models that include intra-S checkpoint modules. Candidates to investigate, in rough order of expected relevance:

- Iwamoto, Tashima, Kimura, Okamoto 2008-2011 family (Biosystems) — Tyson-Novak backbone with ATM/ATR-CHK1/2-CDC25-CDK2 axis. **JP should verify the exact BioModels accession before download; I believe it is in the BIOMD0000000380 range but this is not confirmed.**
- Aguda and Tang 1999 — older, conceptually clean, basis for several later models.
- Branzei-Foiani family — fork stalling kinetics under HU and aphidicolin specifically.

For the first pass, extract only the CHK1-CDC25A-CDK2 axis. The ATM/CHK2 arm can be added later if double-strand-break agents (etoposide, IR) become relevant.

### Phase 2: Minimal checkpoint module in Antimony

Add the following species to v42:

- `D`: DNA damage / replication stress signal. Input function parameterized by HU dose.
- `ATR_active`: activated by D, deactivated first-order.
- `CHK1_active`: phosphorylated by ATR_active, deactivated first-order.
- Modification to existing `Pe` dynamics: active CHK1 either accelerates Pe degradation or inhibits Pe activation. Choose based on which Iwamoto formulation is cleaner.
- Origin firing rate modifier: a multiplicative term on the rate of CycE/CDK2-dependent origin firing, such that active CHK1 suppresses *new* origin firing without halting elongation of already-fired origins.

This is approximately 4 new ODEs and 8-10 new parameters. Use the Iwamoto parameter set as initial values where available; document any deviations from published values in the model construction script.

The integration point with v42 is the Pe → CycE/CDK2 axis. Do not modify the Gerard-Goldbeter cell cycle core itself; only inject the checkpoint as an external regulator of Pe.

### Phase 3: Calibration to Hydroxyurea data

Experimental targets to match:

- 10 µM HU produces a 1.31-fold increase in Ezh2 protein in S-phase cells (p = 0.005, n = 12)
- DMSO Ezh2 protein ratios relative to G0: G1 = 1.16x, S = 1.32x, G2 = 1.48x
- HU Ezh2 protein ratios relative to G0: G1 = 1.42x, S = 1.62x, G2 = 1.47x
- HU should extend S-phase duration without abolishing it (cells still complete S, just slowly)
- Cell cycle period under HU should be longer than under DMSO

Parameter sweep over the checkpoint strength parameter (origin firing suppression coefficient) and the HU input strength to find the regime that matches all targets simultaneously.

Create `simulations/validate_checkpoint.py` following the structure of `validate_v42.py`.

### Phase 4: Parameter space exploration

The core biological question: does the Ezh2-Cyclin D1 feedback combined with checkpoint-mediated S-phase extension produce qualitatively distinct cell cycle dynamics? Specific sweeps:

1. Vary checkpoint strength (α, the CHK1 inhibition coefficient on origin firing) from 0 (no checkpoint) to high (full block). Record: cell cycle period, S-phase duration, mean Ezh2 protein, mean Cyclin D1 protein, mean fraction of cells in G0.
2. At each checkpoint strength, vary the Ezh2-Cyclin D1 feedback strength (K_EZH2 in v42 nomenclature) over a range that includes the published value. Look for bistability or fold bifurcations.
3. Add a constant low-level damage signal `D_baseline` to simulate chronic replication stress (analogous to what may occur in proliferating GNPs approaching differentiation). Ask whether this produces a stable shift into a low-Cyclin D1, high-Ezh2, long-G0 regime.

Output: a heatmap or phase diagram showing cell cycle period and steady-state Ezh2 as functions of checkpoint strength and feedback strength. Identify any separatrices or fold bifurcations. Save figures to `figures/checkpoint_phase_diagram/`.

### Phase 5: Predictions for new experiments

The extended model should generate at least two falsifiable predictions:

1. Aphidicolin (replication stress without nucleotide depletion) should produce a comparable Ezh2 boost to HU per unit of S-phase extension. A different boost magnitude would imply non-checkpoint mechanisms.
2. Sub-lethal etoposide (DSB induction, activates ATM/ATR) should produce a larger Ezh2 boost than HU at equivalent S-phase extension, because the ATM arm provides additional CHK1 activation. This requires Phase 6 (ATM arm addition) before it can be tested.

## Technical environment

- Python 3, Tellurium for ODE simulation, Antimony for model definition.
- Frontenac HPC (Compute Canada, Rocky Linux) available for parameter sweeps. Claude Code v2.0.42 is pinned on Frontenac due to a regression in newer versions; use this version.
- scipy.signal.find_peaks for cell cycle period measurement (already in use in v42).
- matplotlib for visualization, pandas for parameter sweep output, seaborn for heatmaps.
- All simulation scripts should follow v42 conventions: explicit simulation conditions table at the top, validation against numerical targets, output figures saved to `figures/`.

## Starting prompt for Claude Code

Suggested first prompt to open the Claude Code session:

> "Fork the `jpurzner/Ezh2_Hh_Ccnd1_model` repo to `jpurzner/cellcycle_checkpoint_ezh2`. Verify that v42 reproduces its 16/18 validation tests in the new repo by running `simulations/validate_v42.py`. Then search BioModels for an intra-S checkpoint model from the Iwamoto / Tyson-Novak family, download the SBML, and convert the CHK1-CDC25A-CDK2 axis to Antimony format. Show me the converted module before integrating it into v42."

After that initial step, subsequent prompts can address Phase 2 (integration), Phase 3 (calibration), and Phase 4 (parameter sweeps), each as a separate session.

## Open questions to resolve as work proceeds

- Should the checkpoint module use ATR alone (HU-relevant) or both ATM and ATR from the start? Recommendation: ATR-only first pass, add ATM later only if needed for etoposide or IR predictions.
- The v42 model uses Pe for CDC25 function on Cyclin E/CDK2. Confirm that Iwamoto's CDC25A maps cleanly onto Pe, or whether a separate species is needed for clarity.
- EZH2 itself has been reported to localize to stalled forks and undergo replication-stress-dependent post-translational modification (Campbell et al. and others, EMBO J / Nat Comm). The current model treats Ezh2 as a transcriptionally regulated species only. Decide whether to add a stress-stabilization term to Ezh2 protein dynamics or document this as a known limitation.
- The 1.31-fold HU boost is modest. Check whether the model predicts this magnitude with the Iwamoto parameter set out of the box or whether tuning is required. Large deviations from published parameter sets are informative and should be flagged for biological interpretation rather than blindly corrected.
- Eventually: how do these dynamics connect to the GNP differentiation transition? Does the model predict that GNPs near differentiation are operating closer to a fold bifurcation that tips them into stable G0 / Ezh2-low state?

## File structure to create in new repo

```
cellcycle_checkpoint_ezh2/
├── README.md                                    # describe extension vs v42
├── CHECKPOINT_PLAN.md                           # this document
├── src/
│   ├── build_model_v42_mycn.py                 # inherited, do not modify
│   └── build_model_v43_checkpoint.py           # new
├── simulations/
│   ├── validate_v42.py                         # inherited
│   ├── validate_checkpoint.py                  # new
│   ├── sim_hu_dose_response.py                 # new
│   └── sim_parameter_sweep.py                  # new
├── figures/
│   ├── v42_validation/                         # inherited
│   └── checkpoint_phase_diagram/               # new
└── data/
    └── ezh2_protein_phase_treatment.csv        # experimental targets
```
