# v44 status — explicit-replication rebuild on the Heldt 2018 core

Branch `v44-explicit-replication`. Retreat point: tag `v43-milestone`.

## What v44 is

A real-time (minutes, **no `eps`**) cell-cycle model built on Heldt, Barr, Cooper, Bakal &
Novák 2018 (PNAS; BIOMD0000000700), with **explicit DNA replication** (`Dna` synthesized by
active forks `aRc` at fork speed `kSyDna`). On top of Heldt's one-shot G1/S core we added:

- **CycB/CDK1 (MPF) mitotic switch** — Cdc25/Wee1 hysteresis (MPF activates Cdc25 & inhibits
  Wee1; basal `a25/aWee` ignite the switch once the checkpoint clears).
- **CHK1 intra-S checkpoint** gated on active forks (`aRc`) — holds mitosis until replication
  completes (`Dna→1`, forks disassemble).
- **Mitotic APC/Cdc20** + **division-reset event** → sustained cycling.
- **HU → fork speed** (`vfork(HU)` scales `kSyDna`) — hydroxyurea/dNTP depletion slows forks.
- **EZH2 layer** — `Cd` (CyclinD) made dynamic and EZH2-repressible; EZH2 transcription gated
  on `E2f × (CycE+CycA)` (the Heldt analog of v43's `E2F × (Me+Ma)`); EZH2 is a stable protein
  reset by **dilution at division**, so it *integrates* S-phase duration.

Builder: `src/build_model_v44_heldt.py` (`build_model_v44(hu=, with_ezh2=)`).

## What works (the rebuild's primary goal — achieved)

**Concentration-dependent S-phase**, which the Gérard–Goldbeter oscillator structurally could
not produce (`simulations/v44_fork_doseresponse.py`, `v44_ezh2_hypothesis.py`):

| HU | S-phase | G1 (origin-firing marker) | mitoses | Dna@mitosis | EZH2-in-S boost |
|----|---------|---------------------------|---------|-------------|-----------------|
| 0  | 5.5 h   | 2.4 h | 15 | 1.00 | 1.00 |
| 1.0| 8.2 h   | 2.2 h | 11 | 1.00 | **1.18** |
| 2.0| 16.4 h  | 1.9 h | 6  | 1.00 | 1.59 |

- S-phase lengthens mechanistically with HU (slower forks), **G1 stays flat**, mitosis always
  waits for `Dna = 1.0` (perfect checkpoint), cells keep dividing (no arrest).
- **EZH2 integrates S-duration**: EZH2-in-S boost 1.18× at HU=1.0, 1.59× at HU=2.0 (vs the
  experimental 1.31× at 10 µM — close; tunable via EZH2 turnover/gate).
- EZH2→CyclinD feedback is a real tonic brake (G1 with feedback ON 2.38 h vs OFF 1.99 h).

## HH/MYCN module (added) — all structural modules now present

`Cd` (CyclinD) is now driven by the full Hedgehog/MYCN module ported from v42
(`with_hh=True`): SHH→Ptch1→Smo→Gli→CyclinD1, MYCN→CyclinD1 (GDC-independent), EZH2 represses
CyclinD1 transcription. HH species are initialized near the proliferating steady state so the
cell does not false-quiesce during the startup transient. Inputs: `SHH`, `GDC0449`,
`Ptch1_copy_number`, `MYCN_amplification`. Condition responses (qualitatively correct,
magnitudes pending calibration):

| condition | Cd | note |
|---|---|---|
| baseline (SHH=0.5) | 0.70 | — |
| SHH=0 | 0.38 | HH off (MYCN floor) |
| GDC0449=1 (WT) | 0.11 | ~84% CyclinD1 reduction |
| MYCN_amp=2.8 (MB) | 1.20 | MYCN drives CyclinD1 |
| GDC + MYCN_amp (MB+GDC) | 0.66 | only ~45% reduction → **MYCN buffers GDC (MB GDC-resistance)** |

Builder: `build_model_v44(hu=, with_ezh2=, with_hh=)`. All three configs compile and cycle;
the HU→S dose-response is unchanged by adding HH (S 5.6→16.5 h, EZH2-in-S boost 1.18× at HU=1.0).

**STRUCTURE COMPLETE.** Next phase per plan: compile all experimental data and do a single
whole-model parameter-calibration pass (not piecemeal).

## Calibration progress (vs docs/Ezh2_CcnD1_model_targets.md)

**LOCKED — EZH2 / HU core targets** (`simulations/v44_calibrate_ezh2.py`, MB context):
- HU 10µM → EZH2(protein)-in-S / DMSO = **1.29×** (target 1.31×, alt-text 1.22×) ✓
- EZH2 **transcript** within-cycle gradient: S/G0 = 2.0, G2/G0 = 2.05 (Section D MB range 1.8–2.5) ✓
- Monotone gradient G0 < G1 < S < G2 ✓
- Calibrated `kDeEZ = 0.0003` (stable EZH2 integrates S-duration); boost is decoupled from the
  within-cycle gradient (tunable independently).
- Division reset starts the daughter in a true G1 (low CyclinA/E) — biologically correct, and
  lengthened G1 (preS 21% → 33%).

**STRUCTURAL REDESIGN #1 — real G2 phase (DONE).** The mitotic switch was reworked so CyclinB-CDK1
is synthesized (inactive, Tyr15-P) only as replication completes (`g2gate` on `Dna`), building
during G2 instead of being pre-stocked during S and flipping instantly. Result:
- **G2 phase: 2% → ~20%** (a genuine multi-hour G2); period 9 → ~11–13 h (toward 22 h).
- HU redistribution now strong & correct: S 50→58→74%, **G2 20→14→4%** across HU 0→2 (Section F:
  S up, G2 down), divisions 23→18→11 (no arrest).
- Absolute durations now realistic: S ≈ 6 h, G2 ≈ 2.5 h (only G1 remains short).
- Protein EZH2 gradient improved (G2/G0 1.99 → 1.89); transcript gradient still 2.0; HU boost 1.23.
- Cycles robustly in BOTH GNP (`MYCN_amp=1`) and MB (`MYCN_amp=2.8`); division reset moderated
  (`Ca_div=0.30`) so GNP (lower drive) does not false-quiesce. HH/GDC biology intact (GNP+GDC −87%,
  MB+GDC −49%).

**STILL OPEN — G1 too short.** G1 remains ~8% (target cycling ~58%): the Heldt bistable restriction
point caps G1 at ~30–35% of the cycle in a robustly-cycling regime — below `kPhRbCd ≈ 0.15` the cell
either quiesces or the HU response goes pathological. A longer G1 needs a different G1 model (e.g. a
cell-growth/size-gated restriction point or an explicit G1 timer), not parameter tuning. This is the
binding constraint on the absolute proportions (S/G1 ratio inverted vs data).

**DEFERRED — G0 as a separate quiescent pool** (structural redesign #2) + absolute period (~22h).

**Original deferred note — absolute MB phase proportions** (Section F: G1 58 / S 21 / G2 21 among cycling).
The current model is S-dominated (S ≈ 63%, G1 ≈ 11%, G2 ≈ 2%) and cannot reach these by parameter
tuning — three structural walls (mapped in `v44_calibrate_phases.py`):
- **G1** capped ~33–40%: below `kPhRbCd ≈ 0.14` the bistable restriction point either quiesces or
  makes the HU response pathological (S decreases under HU).
- **G2** intrinsically ~2%: mitosis fires immediately after replication (the checkpoint holds a
  preMPF reservoir that flips instantly). Needs a redesigned mitotic switch (CyclinB-CDK1 building
  *after* S) for a genuine multi-hour G2.
- **G0**: the within-cycle low-pRb window is a proxy; the experimental ~25% G0 is a separate
  quiescent/exited pool. This also explains the steep **protein** EZH2 gradient (model G2/G0 ≈ 1.9
  vs Section E 1.48×): the model's brief post-mitotic G0 reference is artificially low.

Next structural task: (a) mitotic-switch redesign for a real G2, (b) model G0 as a separate
quiescent pool. Then absolute proportions + the protein gradient + absolute period (~22 h) become
fittable. The HU/EZH2 mechanism (the project's core claim) is already captured.

## What is weak / open (honest assessment)

1. **The "longer S → longer *next* G1" loop is weak.** G1 does not lengthen with HU (the
   apparent increase was a measurement artifact of a fork-dependent S-onset threshold; with a
   fork-independent marker, true G1 is flat/slightly shorter). The EZH2 brake on G1 (~0.4 h) is
   roughly **cancelled by the residual CycA carried into the daughter**, which shortens G1.
   This matches Chao et al. 2018 (human cell-cycle phases are largely uncoupled/memoryless) —
   so it may be a *correct* prediction, not just a calibration gap. Making the loop dominant
   would require giving the EZH2–CyclinD–G1 axis more leverage (e.g. CyclinD as the dominant Rb
   kinase, or EZH2 acting on more than `Cd`).
2. **EZH2-in-S boost is 1.18× at HU=1.0 vs 1.31× target** — a calibration gap (slow EZH2 turnover
   / stronger gate would close it).
3. **HH/MYCN module not yet ported** — `Cd` is driven by a placeholder `mitogen = 1.0`. The full
   SHH/Gli/MYCN→CyclinD input (and the GDC0449 between-condition invariants) remain to attach.
4. Phase-duration absolutes are Heldt-default (cycle ~8.5 h at HU=0); not yet calibrated to the
   GNP/MB system's ~20 h.

## Files
- `src/build_model_v44_heldt.py` — builder
- `models_external/heldt2018{.ant,_BIOMD700.xml}` — Heldt core
- `simulations/v44_build_test.py` — mitotic switch + cycling
- `simulations/v44_fork_doseresponse.py` — fork-speed → S-phase
- `simulations/v44_ezh2_hypothesis.py` — HU → EZH2 → CyclinD → G1 test
- `simulations/diag_v44_*.py` — the diagnostics proving GG could not do this (motivation)
