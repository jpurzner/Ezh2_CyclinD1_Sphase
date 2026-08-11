# Two-compartment PRC2 / H3K27me3 module — reservoir-sustained transient-G0 amplifier (2026-08-09)

Implements JP's molecular mechanism (2026-08-09) for how heavily-H3K27me3-marked CDK6 & CyclinD1 are also
highly expressed, and how the mark can nonetheless deepen a transient G0. Flag `with_proximal_distal`
(default OFF); env `PROXIMAL_DISTAL=1`. Requires `with_cdk6_gli` + the H3K27 chain (`Mk`).

## The mechanism (JP)

- **Repression is the PRC2 COMPLEX, not H3K27me3.** The mark recruits/holds the complex (EED read-write);
  the complex represses.
- **Two compartments at the CDK6/CyclinD1 loci:**
  - **DISTAL RESERVOIR** = broad H3K27me3 (`Mk`, the existing me-chain: read-write substrate, replicative
    dilution, slow turnover) = the memory.
  - **PROXIMAL OCCUPANCY** (`P_prox`, NEW [0,1] species) = PRC2 complex at the TSS = the actual repressor.
    Loaded by the reservoir (`a_P·Mk^n_anch`) **only when Pol2 is NOT elongating** (`1 − elong_gate`), and
    evicted by active elongation (`elong_gate`) + a small basal. `elong_gate = E2f^n_el/(K_el^n_el+E2f^n_el)`
    (the proliferation/elongation state).
- **Marked-but-expressed:** in cycling cells Pol2 elongation keeps the proximal promoter cleared
  (`P_prox ≈ 0`) even though the distal reservoir is high → the gene is ON despite heavy mark.
- **Reservoir-sustained transient G0:** when a cell stops cycling (E2f collapses → Pol2 stops), the
  persistent reservoir re-loads `P_prox` → represses CDK6/CyclinD1 → the arrest deepens & persists. This does
  **not** need EZH2 (fixes the prior `PRC2_rep = EZH2×mark` collapse in arrest). It RELEASES when
  mitogen/cycling returns (elongation evicts `P_prox`) → **reversible** = transient.

## Why baseline is preserved by construction

`P_prox ≈ 0` whenever the cell is cycling (elongation gate ≈ 1 + slow P kinetics ride over transient G1
E2f-dips), so `prox_amp ≈ 1` and CDK6/CyclinD1 are untouched. Confirmed: base cdk6 5.08, Cd 7.8, 30 div/750h
across all amplifier strengths. The **default model (flag OFF) is unchanged, still 31/32.**

## Parameterization — the scenario map (`amplifier_hysteresis.py` + `amp_scenario_wf`)

The amplifier strength (`a_P` load, `f_P` = 1−max-extra-repression) tunes the arrest behavior. All settings
keep baseline; validation and transient dwell vs strength:

| scenario | a_P | f_P | verdict | dwell (h) | validation | new fails |
|---|---|---|---|---|---|---|
| B-leaky | 0.05 | 0.90 | transient | 89 | **30/32** | GNP+HHi CyclinD1 |
| C-lite | 0.08 | 0.85 | transient | 105 | **30/32** | GNP+HHi CyclinD1 |
| **C (default)** | **0.10** | **0.80** | **transient** | **148** | **30/32** | GNP+HHi CyclinD1 |
| C-mid | 0.12 | 0.75 | permanent-lock | ∞ | 30/32 | GNP+HHi CyclinD1 |
| C-deep | 0.15 | 0.65 | permanent-lock | ∞ | 30/32 | GNP+HHi CyclinD1 |
| A-mid | 0.30 | 0.45 | permanent-lock | ∞ | 29/32 | + Palbo EZH2 |
| A-strong | 0.50 | 0.30 | permanent-lock | ∞ | 28/32 | + Palbo EZH2, MB+HHi CyclinD1 |

- **Reversible TRANSIENT G0** (the biologically-sensible regime, flag default): `a_P ≤ ~0.11` → a tunable
  ~90–150 h dwell (≈4–6 missed cycles = JP's "over a couple of divisions"), then re-entry. **30/32.**
- **PERMANENT-LOCK cliff** at `a_P ≈ 0.11` (sharp): the reservoir self-sustains (P_prox high → drive low →
  E2f low → elongation off → P_prox stays high) = a bistable latch. At high strength it additionally breaks
  **Palbo-EZH2** and **MB+HHi CyclinD1** → ≤29/32 = the known chromatin-latch exclusion, reached from a new
  direction. So the deep permanent latch is data-excluded; the **transient** regime is not.

## The one new cost (a prediction, not just a failure)

Every amplifier setting drops **CyclinD1 GNP+HHi/GNP** from 0.157 (bulk target) to ~0.08, because Hedgehog
withdrawal (GNP+HHi) IS a sustained arrest → the reservoir re-loads the proximal promoter → CyclinD1 is
further repressed. This is mechanistically consistent (a Hh-withdrawn GNP exiting cycle should have very low
CyclinD1) and reads as a **testable prediction**: CyclinD1 in Hedgehog-withdrawn GNPs is lower than the 0.157
bulk estimate (single-cell / sorted measurement would settle it). The other miss (MB HU G2) is the
pre-existing baseline structural near-miss, unrelated.

## Open scenario (JP: G0 behavior is unknown; the experiment isn't done)

The map above assumes the reservoir re-loads the proximal promoter in arrest. The competing G0 unknown JP
named — EZH2 declines to a floor vs the mark has more time to accumulate — is captured by `a_P`/`n_anch`
(how strongly the reservoir sustains `P_prox` independent of EZH2). The **discriminator experiment**: sort
Atoh1+ G0 vs cycling MB cells and CUT&RUN/ChIP H3K27me3 + PRC2 at the CDK6/CyclinD1 **proximal promoter** —
does proximal PRC2 occupancy rise in G0 (reservoir-sustained, scenario A/C) or fall (EZH2-gated, scenario
B)? Temporarily arresting GNPs isn't feasible; the MB arrest experiment is the one to do.

## Files
- `src/build_model_v44_heldt.py` — `with_proximal_distal` block (P_prox species, elong_gate, prox_amp).
- `simulations/amplifier_hysteresis.py` — transient-vs-permanent dwell probe.
- `simulations/amp_scenario_wf.js` — scenario sweep (validation + dwell).
- `simulations/fig_v44_proximal_distal.py` — the review figure (4 panels).
- Nothing baked into the default; `with_proximal_distal` is a flag, default OFF, and JP decides the scenario.
