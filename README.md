# Ezh2_CyclinD1_Sphase

Modeling the **intra-S-phase checkpoint** extension of the published MYCN–EZH2–Cyclin D1–Hedgehog
cell-cycle model, applied to cerebellar granule-neuron progenitors (GNPs) and SHH-subgroup
medulloblastoma (SHH-MB).

This repository extends the published **v42** model (Chahin et al., on EZH2/Cyclin D1 negative
feedback in GNPs and medulloblastoma) to test a specific hypothesis:

> **Longer S-phase → more EZH2 → stronger Cyclin D1 repression → longer subsequent G0/G1.**

Experimental anchors driving the work: 10 µM hydroxyurea (HU) → **1.31× EZH2 in S-phase**
(p = 0.005); within-cycle DMSO EZH2 ratios vs G0 of **G1 = 1.16×, S = 1.32×, G2 = 1.48×**; and a
small HU-induced increase in the G0 fraction (delineated experimentally by Rb phosphorylation).

> ⚠️ The published model `src/build_model_v42_mycn.py` is **frozen** (byte-identical to the
> `v42-published-baseline` tag) and is never modified. All extensions are built on top of it.

---

## The three model versions

| Version | File | What it is | Status |
|---------|------|-----------|--------|
| **v42** | `src/build_model_v42_mycn.py` | Published MYCN–EZH2–CyclinD1–HH model on the Gérard–Goldbeter (2009) cell-cycle engine (ε = 150 global clock). | Frozen baseline |
| **v43** | `src/build_model_v43_checkpoint.py` | Intra-S checkpoint via a **phase-gated clock**: φ(HU) slows ε *specifically during S* (gated on Cyclin A). EZH2 transcription rewired to a Cyclin E + Cyclin A (S-window) gate. 12/12 validation; emergent 1.31× S-phase EZH2 at HU = 1.0. | Complete (tag `v43-milestone`) |
| **v44** | `src/build_model_v44_heldt.py` | **Ground-up rebuild** on the Heldt et al. 2018 (PNAS) real-time core with **explicit DNA replication**, so S-phase duration is genuinely **concentration-dependent** (set by replication fork speed). | Structure complete; calibration pending |

### Why v44 exists (the key scientific finding)

The Gérard–Goldbeter engine is a **relaxation oscillator**, and its S-phase duration turns out to
be a **structural invariant** — exhaustive testing (`simulations/diag_v44_*.py`) showed that *no*
molecular concentration (CHK1→Cdc25, Cdk1 throttling, replication-licensing gates, …) can lengthen
S; only the global time-scale ε can. The v43 "phase-gated ε" checkpoint works and matches the data,
but slowing a clock is phenomenological. To make S-phase length emerge from real replication
kinetics, v44 rebuilds the engine on **Heldt 2018** (BioModels `BIOMD0000000700`): a real-time
(minutes, **no ε**) model in which DNA is synthesized by active replication forks. Full rationale and
the literature survey are in `docs/`.

---

## v44 architecture

Built by `build_model_v44(hu=, with_ezh2=, with_hh=)` on top of the Heldt core:

1. **Heldt G1/S core** — Rb–E2F restriction point, Cyclin E/A–CDK2, p21/p53, PCNA, and **explicit
   DNA replication** (`Dna` synthesized by active replication complexes `aRc` at fork speed `kSyDna`).
2. **HU → fork speed** — hydroxyurea/dNTP depletion slows forks: `vfork(HU)` scales `kSyDna`, so S
   lengthens concentration-dependently (G1 unchanged, mitosis waits for `Dna = 1`, no arrest).
3. **CyclinB/CDK1 (MPF) mitotic switch** — Cdc25/Wee1 hysteresis; ignites once the checkpoint clears.
4. **CHK1 intra-S checkpoint** — gated on active forks (`aRc`); holds mitosis until replication completes.
5. **Mitotic APC/Cdc20 + division-reset event** — gives sustained cycling (Heldt alone is one-shot).
6. **EZH2 layer** — transcription gated on `E2f × (CycE + CycA)`; EZH2 is a stable protein diluted at
   division, so it **integrates** S-phase duration. EZH2 represses Cyclin D1 (`EZH2i` toggles the feedback).
7. **Hedgehog + MYCN → Cyclin D1** — ported from v42; SHH → Ptch1 → Smo → Gli → CyclinD1, MYCN drives
   CyclinD1 GDC-independently. Inputs: `SHH`, `GDC0449`, `Ptch1_copy_number`, `MYCN_amplification`.

### What v44 reproduces so far (qualitative; magnitudes pending calibration)

- **Concentration-dependent S-phase**: fork speed ↓ → S 5.5 → 16.5 h, G1 flat, mitosis always at
  `Dna ≈ 1.0`, no arrest (`simulations/v44_fork_doseresponse.py`).
- **EZH2 integrates S-duration**: EZH2-in-S boost 1.18× at HU = 1.0, 1.59× at HU = 2.0
  (`simulations/v44_ezh2_hypothesis.py`).
- **Weak S→next-G1 coupling** (EZH2 brake ≈ cancelled by residual CyclinA carried into the daughter)
  — consistent with the small HU-induced G0 increase and with Chao et al. 2018 (largely uncoupled phases).
- **HH/MYCN/GDC biology**: GDC0449 ablates CyclinD1 in WT (~84%) but MYCN-amplified MB retains much
  more (~55%) — the MB GDC-resistance.

---

## Repository layout

```
src/
  build_model_v42_mycn.py        # FROZEN published v42 model (do not modify)
  build_model_v43_checkpoint.py  # v43: phase-gated-ε intra-S checkpoint + EZH2 rewire
  build_model_v44_heldt.py       # v44: Heldt-based explicit-replication rebuild (current)
models_external/
  heldt2018_BIOMD700.xml         # Heldt 2018 source SBML (BIOMD0000000700)
  heldt2018.ant                  # Antimony conversion (what the v44 builder reads)
docs/
  model_dynamics.md                          # how the model sets phase durations; cycle walkthrough
  research_concentration_dependent_checkpoint.md  # literature survey -> the v44 decision
  v44_verification_findings.md               # proof GG cannot give concentration-dependent S
  v44_status.md                              # current v44 state + open items
simulations/
  validate_v42.py / validate_v43.py          # validation suites (v43: 12/12)
  calibrate_ezh2_rewire.py                   # v43 EZH2-gate calibration + phase helpers
  v44_build_test.py                          # v44 mitotic switch + cycling
  v44_fork_doseresponse.py                   # v44 fork speed -> S-phase duration
  v44_ezh2_hypothesis.py                     # v44 HU -> EZH2 -> CyclinD1 -> G1 test
  diag_v44_*.py                              # diagnostics: why GG fails; coupling-point search
  diag_phi_*.py                              # v43 phase-gated-clock calibration diagnostics
CHECKPOINT_PLAN.md                           # original extension plan
```

(The v42 paper figures and `V42_MODEL_*.md` summaries from the published release also live under
`simulations/`.)

---

## Branches & tags

- `v44-explicit-replication` — current work (the v44 rebuild). **Active branch.**
- `feature/intra_s_checkpoint` — the v43 milestone (`v43-milestone` tag) — the retreat point.
- `main` — published v42 release.
- Tags: `v42-published-baseline`, `v43-milestone`.

---

## Getting started

```bash
python -m venv venv && ./venv/bin/pip install -r requirements.txt
```

Built on Antimony/Tellurium (libRoadRunner, CVODE integrator). Inspect any model with, e.g.:

```bash
./venv/bin/python -m src.build_model_v44_heldt     # builds + prints a cycling sanity check
```

Run a v44 experiment (from the repo root):

```bash
./venv/bin/python simulations/v44_fork_doseresponse.py   # fork speed (HU) -> S-phase duration
./venv/bin/python simulations/v44_ezh2_hypothesis.py     # HU -> EZH2 -> CyclinD1 -> G1
```

---

## Status / next step

All v44 structural modules are present and the model cycles. The next phase is a **single
whole-model parameter-calibration pass** against the full experimental dataset (phase-duration
absolutes ~20 h, the HU 1.31× EZH2 anchor + G0 shift, the DMSO within-cycle EZH2 gradient,
GDC0449/MYCN CyclinD1 reductions, EZH2 mean/G0 ratios), producing a `validate_v44.py` suite.
Parameter *magnitudes* in v44 are currently placeholders pending that fit.

## Citation

Cite Gérard & Goldbeter (2009) *PNAS* for the original cell-cycle engine, Heldt, Barr, Cooper,
Bakal & Novák (2018) *PNAS* (`BIOMD0000000700`) for the v44 real-time core, and Chahin et al.
(in preparation) for the EZH2/CyclinD1/MYCN/Hedgehog integration and the intra-S-checkpoint extension.
```
