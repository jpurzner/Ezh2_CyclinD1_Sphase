# v44 structure analysis

Executing the dependency-structure exploration plan — reading v44 as a *machine* instead of 180 dials.
Run any script with `./venv/bin/python simulations/struct/<name>.py` (each writes a `.json` + `.png`).

## The answer: v44 is a three-tier machine

| tier | timescale | contents | role |
|---|---|---|---|
| **FAST · input** | ~1 min | whole Hedgehog cascade + CyclinD1 | *slaved* → collapses to `Cd* = F(SHH, Ptch1, EZH2, mark)`; feed-forward periphery |
| **MEDIUM · engine** | min–1 h | one 35-species feedback core (Rb-E2f R-point, Skp2-p27 toggle, MPF-Cdc20 oscillator, replication) | runs the cycle |
| **SLOW · memory** | hours | mass (33 h), Gli1_epi (21 h), EZH2 (14 h), me3 mark (7 h) | *modulates* the oscillator — **the fold/G0 story lives here** |

## The five analyses

- **`E_conservation.py`** — 3 conserved pools (total Rb=5, Cdh1=1, Gli=0.6); true dim **44**. Rb conserved even across division ⇒ mass-gate ≡ Rb-dilution is an exact identity. (Corrects the plan's guess of 5 pools.)
- **`A_loops.py`** — signed influence graph: 203 edges, 0 sign-flips. Feedback **core = one 35-species SCC** (10 modules); Hedgehog + mass are feed-forward. Loops L1/L2/N1/N2 verified with signs. mass-homeostat (N5) closes through the *division event*, not continuous dynamics.
- **`B_timescales.py`** — τ = 1/|Jᵢᵢ|: ~5 orders of magnitude. Slow memory drivers (mass, EZH2, mark) vs fast slaved Hedgehog+CyclinD1 (~1 min).
- **`C_dynamics.py`** — phase portraits: MPF-Cdc20 relaxation oscillator, Rb-E2f R-point jump, Skp2-p27 toggle, EZH2 as a ~4%-amplitude slow buffer; the commitment threshold is a **sharp switch** (~166 %/unit birth-p27).
- **`F_knockouts.py`** — load-bearing ranking: 6 knobs break the cycle (E2f, CyclinE, p27-clearance, CyclinD1, CDK4/6, growth); EZH2 feedback breaks only 7 (separable modulator); **Skp2 & Cdc20 arms redundant** (break 0); **`K_EZH2_repression` dynamically DEAD** in the chain default (break 0).

Full write-up + figures: the "v44 Structure — Findings" artifact (companion to the exploration plan).
