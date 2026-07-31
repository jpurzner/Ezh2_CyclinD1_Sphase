# EZH2 – CyclinD1 – Hedgehog cell-cycle model — Project Wiki

A durable map of **what we've built, why, how it evolved, and everything we've tried** (including the dead-ends).
Compiled 2026-07-31 from the two code repos, ~50 design docs, 60 curated project-memory notes, and the full
233-commit history.

## What this project is

A tellurium/roadrunner **ODE model of the EZH2 – CyclinD1 – Hedgehog network controlling the cell cycle** in
**P7 cerebellar granule-neuron precursors (GNP)** and **SHH-medulloblastoma (MB)**, built to support the
**Chahin et al.** paper ("GNP cell-cycle exit and differentiation regulated through negative feedback of Ezh2 and
Cyclin D1"). The central hypothesis the model formalizes: **EZH2 / H3K27me3 acts as a graded *governor* on CyclinD1**,
setting the mitogen (Shh) threshold for proliferation rather than an on/off switch — one cell-cycle engine, two cell
types (GNP vs MB) differing only in their biological identities (Hedgehog tone, INK4/CDKi, MYCN, CyclinD2).

## The big picture (three tiers)

1. **Mitogen input** — Shh → Ptch1/Smo → Gli → CyclinD1 (fast, slaved; optional realistic OU noise + cilium + Ptch1 adaptation).
2. **Cell-cycle oscillator core** — Heldt-2018 engine + CycB/CDK1 mitotic switch + a **two-step Rb** R-point + Skp2–p27 commitment + a **dynamic CDKI module** (INK4 p18/p19, CIP/KIP p27/p21).
3. **Slow epigenetic / growth memory** — the **CyclinD1 transcript reservoir** (the Hh memory), the serial **H3K27me3 me0→me3 chain** on Ccnd1, EZH2 as a proliferation-coupled writer, and growth-timed cycle length.

## Current state (2026-07-30)

- **Default model:** v44 + two-step Rb + serial H3K27me3 chain + two-cyclin (D1+D2) + **dynamic-CDKI module (baked)**.
- **Validation:** **30/32** targets (cell-type-separated GNP/MB/ratios); the two failures are the structural HU-S/HU-G2 pair.
- **Repo/branch:** `Ezh2_CyclinD1_Sphase`, branch `transient-g0-reframe` (GitHub `jpurzner/Ezh2_CyclinD1_Sphase`).

## How to read this wiki

New here? Read **[Origins & Foundations](01_origins_and_foundations.md)** → **[Evolution Timeline](03_evolution_timeline.md)**
→ **[Model Architecture](02_model_architecture.md)**. Trying to remember whether we already tried something? Go straight
to **[Things Tried & Abandoned](05_things_tried_and_abandoned.md)** and **[Mechanisms Explored](04_mechanisms_explored.md)**.

| Page | What's in it |
|---|---|
| [01 — Origins & Foundations](01_origins_and_foundations.md) | How the model was first assembled from published BioModels (Gerard, Novak, Abroudi, Heldt…); the v17→v31 lineage in the old `sym` repo. |
| [02 — Model Architecture](02_model_architecture.md) | The current default model: the three tiers, species, modules, key equations, load-bearing dials, build flags. |
| [03 — Evolution Timeline](03_evolution_timeline.md) | The chronological story v17→v42→v43→v44→two-step-Rb→me-chain→two-cyclin→CDKI, dated, with the validation state at each step. |
| [04 — Mechanisms Explored](04_mechanisms_explored.md) | Each biological mechanism and how it's modeled: EZH2/H3K27me3 governor, CyclinD1 reservoir, two-cyclin, MYCN, CDKI, mitogen noise, commitment/R-point, transient G0, G1-lengthening, withdrawal. |
| [05 — Things Tried & Abandoned](05_things_tried_and_abandoned.md) | The dead-ends and negative results, with **why** each was dropped (MYCN latch, D2 equal-marking, 16h recal, oscillation, discrete-G0, strong H3K27me3 memory…). |
| [06 — Calibration & Validation](06_calibration_and_validation.md) | The validation-target set, the calibration methods used over time, the current best + what fails, the structural limits. |
| [07 — Data & Evidence](07_data_and_evidence.md) | The experimental data (RNA-seq, flow, ChIP, cycle length, half-life) and literature grounding the model; measured vs assumed. |
| [08 — Figures Catalog](08_figures_catalog.md) | Every figure: script, what it shows, which concept/paper panel it supports. |
| [09 — Repos, Directories & Reproduction](09_repos_directories_reproduction.md) | Map of the three directories + two GitHub repos + how to run everything. |
| [10 — Reviewer Engagements & Open Questions](10_reviewer_and_open_questions.md) | The expert reviews that shaped the model, JP's key conceptual corrections, and the live open-questions agenda. |

*Provenance: this wiki summarizes the state as of the dates cited within each page; where a claim is about code behavior or a
parameter value, it reflects what was true when written — verify against current source before relying on a specific number.*
