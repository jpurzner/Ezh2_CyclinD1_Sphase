# archive/v45 — the stochastic-commitment prototype (set aside)

**v45 was an exploration; v44 is the committed working model** (`src/build_model_v44_heldt.py`,
`docs/v44_MODEL_DESCRIPTION.md`). This folder preserves the v45 work so it isn't lost, but it is **not
maintained** and is **not part of the active model/figure suite**.

## What v45 was

A parallel-successor prototype that replaced v44's cell-cycle engine with a reduced kernel: a **bistable
CDK2–p27 commitment toggle** + **Rb-dilution G1 timer** + **mitogen-gated CyclinD bootstrap** +
**mass-capped quiescence** + a **dynamic, mitogen-tracking EZH2** species, run as a **stochastic ensemble**
of single cells. It validated 20/20 against the same targets (`validate_v45.py`) and produced the figure
suite here (rescue, commitment bifurcation, single-cell anatomy, EZH2/period, mitogen ramp, CyclinD1
decline, feedback parameter-space, the flipped falsifiable prediction, cutoff-rate, and a v44-vs-v45
comparison).

## Why it was set aside (the lesson)

v45 demonstrated useful ideas — commitment as a bifurcation, intrinsic noise giving the fractional
arrest/rescue — but its reduced kernel **sacrificed Heldt's explicit DNA replication / variable S-phase**,
which is the substrate for the paper's central mechanism ("S-duration integrates EZH2", Fig 4G). And on
inspection, **v44 already contains the good parts**: its Skp2–p27 feedforward *is* the bistable CDK2–p27
commitment, and its growth/size gate *is* algebraically the Zatulovskiy Rb-dilution timer (on Heldt's
conserved Rb, `[Rb]=tRb/mass ∝ 1/mass`). So porting v45's mechanism into v44 would have been cosmetic for
real recalibration risk. v44 keeps the mechanistic S-phase **and** the bistable, Rb-dilution-grounded
commitment. See `docs/v44_MODEL_DESCRIPTION.md` §E and the git history around this decision.

## Running these scripts

The path constants assume the original `simulations/` location (they reach `src/` and import
`simulations.validate_v44`). To run from `archive/v45/`, the `sys.path` depth would need adjusting. Treat
this as a frozen snapshot rather than a live target.
