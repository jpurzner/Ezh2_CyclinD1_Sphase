# Buffering-memory phase diagrams over parameter pairs (2026-07-21)

*JP: "generate multiple parameter-pair 2D plots to demonstrate changes to memory based on these parameters within
the equation." Self-contained REDUCED equation (no full model needed). `simulations/fig_memory_parameter_pairs.py`
→ `fig_memory_parameter_pairs.png/.pdf`.*

## The equation
```
    dM/dt = λ·W0·(a0_prc2 + a_rw_prc2·M^n_rw)·(1 − M)  −  δ_eff·M
    δ_eff = del_mk + k_jmjd3_gli·Gli1 + ln2/Tc + k_h33·tx        (the three-term effective decay)
```
W0 (writing scale ∝ EZH2) is calibrated so the baked config sits at the saturated mark M*≈0.96 (⇒ W0=18.7,
baseline δ_eff=0.00386, **baseline τ = 2.86 h** — matches the full model's ~2.8 h mark half-life). ★ = baked v44.

## Two FACES of memory (the two rows)

**Row 1 — LATCH memory (bistability = the mark self-sustains ⇒ remembers ON/OFF through mitogen change).**
Metric = hysteresis width in mitogen λ (0 = monostable = no latch). Governed by the WRITE shape (a0, a_rw, n_rw):
- **The latch needs suppressed nucleation `a0` + cooperative read-write `n_rw≥2`.** Bistability appears only for
  **a0 ≲ 1×10⁻⁵** (~300× BELOW the baked a0=2.9×10⁻³) — panels A, B, D. **v44 is deep in the monostable zone**
  (★ far above the boundary): the baked mark has *dose-response* buffering only, NO latch (consistent with
  `v44-g0-arrest-monostable`, `prc2-buffering-memory-topologies`).
- **Cooperativity is the enabler** (panel B): the bistable region GROWS with n_rw — at n_rw=1 it barely exists,
  at n_rw≥3 it extends up to a0~10⁻⁴ with wide hysteresis. Panel C (a_rw × n_rw at a0=10⁻⁵) shows the latch needs
  BOTH read-write strength and cooperativity.
- **Fast cycling fights the latch** (panel D): short Tc (heavy dilution) raises δ_eff, shrinking the bistable region.

**Row 2 — PASSIVE memory τ = ln2/δ_eff (how long a transient/withdrawal is remembered).**
Set purely by the DECAY params. Ranges ~1–11 h across the panels:
- **Eraser (panels E, F, at Gli1=0.5 = high-mitogen/MB): memory COLLAPSES as k_jmjd3_gli rises** — the mitogen-gated
  loss. In high-Gli MB the mark is short-lived (~1–3 h); in low-Gli differentiating/withdrawn cells (Gli→0) the
  eraser term vanishes and τ extends (this is the mitogen-gating that gives ChIP MB<GNP + the withdrawal story).
- **Cell-cycle dilution (panels F, G): slow cycling ⇒ longer memory.** τ rises with Tc (less ln2/Tc dilution) — the
  cell-cycle-coupled memory. Fast-cycling GNP/MB dilute the mark; a lengthened/withdrawn cell retains it.
- **H3.3 turnover (panels G, H): raising k_h33 shortens τ**, and (because it is tx-gated) does so specifically where
  Ccnd1 is transcribed — a second, replication-independent memory-erasing channel.

## One-line summary
v44's memory today is **passive + monostable**: a low-pass buffer with half-life τ = ln2/δ_eff (~2.9 h at the
operating point), tuned by the eraser, dilution (Tc), and H3.3 — NOT a bistable latch. A true self-sustaining
memory latch is reachable in the equation but only by **suppressing nucleation a0 ~300× AND keeping n_rw≥2** — a
regime the saturated-mark calibration and the 5.07× fold wall (`v44-chromatin-memory-validation-bound`) currently
exclude. The plots show exactly how far the operating point sits from that boundary.
