# EZH2 ⊣ CyclinD1 repression and the mitogen requirement — large GNP/MB study (2026-07-20)

*JP: large simulation of GNP + MB over multiple divisions, GNP mitogen ramp up/down, varying parameters; emphasize
the EZH2→CyclinD1 repression ON vs OFF, remembering that with the repression ON you need MORE mitogen than with it
OFF. "I feel like we are almost there." Files: `sim_ezh2_repression_mitogen.py`, `fig_ezh2_repression_mitogen.py`.
Repression ON = f0_prc2 baked; OFF = f0_prc2=1.0 (zero CyclinD1 repression).*

## Your emphasis, confirmed
- **GNP: repression ON needs 2.6× more mitogen than OFF** (proliferation threshold Shh 0.104 vs 0.040). Exactly the
  behaviour you described — the mark raises the commitment threshold (Fig A).
- **MB does NOT shift** (1.0×): its constitutive Gli keeps CyclinD1 above threshold at every tested Shh, so it
  overpowers the repression (Fig B). MB proliferates regardless; the repression only gates the mitogen-limited GNP.
- Multi-division trajectories (Fig D): ~13 GNP divisions over 300 h, mark stays high through the cell-cycle dilution
  dips. GNP ramp up/down shows the proliferation hysteresis (Fig C).

## The "almost there" tension — and its resolution
A 2D parameter search exposed a clean **anti-correlation** with the two obvious knobs (f0_prc2 = repression floor,
mark de-saturation): the *saturated* mark gives the big mitogen shift but goes the **wrong** way with mitogen
(slope < 0), while *de-saturating* the mark fixes the direction (mark peaks at proliferation, ChIP-seq-correct) but
**kills the shift** (→1.0×). You can't get both from those two knobs.

**The missing knob is `K_prc2`** (the PRC2-Hill threshold — how low a mark still represses hard). Adding it to the
sweep **resolves the tension**: with **lower K_prc2 (~0.002–0.004) + a de-saturated mark**, the model gets **BOTH**
(Fig E, ★ cells):
- a real mitogen-requirement shift (repression ON needs **~3–10× more mitogen**), AND
- a mark that **tracks proliferation UP** (slope +0.02 to +0.04 — maximal at max division, matching your ChIP-seq;
  Fig F).

The mechanism: de-saturating makes the mark *sub-saturating* so it tracks EZH2/proliferation (the right direction);
lowering K_prc2 makes the repression strong even at that lower mark level, so the threshold shift survives. The
sweet spot is roughly **K_prc2 ≈ 0.003–0.004, a_rw_prc2 ≈ 0.0005, del_mk ≈ 0.0024** (shift ~3×, slope +0.02) — a
biologically moderate requirement-shift with the correct mark direction. Much lower K_prc2 (≤0.001) over-shoots the
shift to ~10× (the cell needs almost saturating mitogen).

## What this means / next
- The model **can** reproduce your two co-requirements (repression raises the mitogen requirement AND H3K27me3 peaks
  at max division) — they are not incompatible; they just need `K_prc2` moved, not only the mark strength. That is
  very likely the "almost there."
- **Caveat (not yet checked):** lowering K_prc2 + de-saturating the mark also shifts the EZH2i de-repression, mark
  persistence, and MB/GNP fold. Whether a "both" point is compatible with the full 28/29 validation needs a **joint
  re-optimization** with K_prc2 + the mark knobs in the search and the mark-slope added as a target — that is the
  concrete next step to fold this into the working model (I did NOT bake it; the default is unchanged).
- Retired `fig_prc2_transcription_paradox` (its global-up/local-down separation contradicted the ChIP-seq).
