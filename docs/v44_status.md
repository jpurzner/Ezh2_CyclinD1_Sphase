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
