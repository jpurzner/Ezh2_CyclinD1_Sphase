# Overnight work — 2026-07-24 (mRNA-reservoir memory + optimized noisy Hh input)

*JP: "the cumulative reservoir of cyclin D1 transcript is probably the memory mechanism... the protein is degraded
for the majority of cell cycle but when not degraded it forms more rapidly when there is lots of Cyclin D1 RNA, Gli1
increases the pool, Ezh2 decreases the pool... this was the memory we were trying to explore before" then "once the
model is rebuilt can you test it and then try an optimize your noisy Hh input, take your time."*

## 1. CyclinD1 memory re-baked as the mRNA RESERVOIR (not the protein) — validation 28/29
The Hh memory now lives in the **CyclinD1 transcript reservoir** (slow mRNA), with the **protein fast** (a labile
readout). Replaced the earlier protein-slowing bake.
- **In `_ts_bake` (`src/build_model_v44_heldt.py`):** `k_Cd_mRNA_deg` 0.8→0.001 (mRNA τ ≈ 16 h reservoir); the three
  transcription rates `k_Cd_tx_{basal,Gli_max,MYCN}` scaled ×1/800 to **hold both mRNA and protein levels** (so fold
  ratios and the g_prc2/H3.3 terms reading Cd_mRNA are unchanged). Protein stays fast (`k_Cd_deg=1.0` default).
- **EZH2 = incoherent feedforward (IFFL):** Gli raises the transcript AND raises EZH2 (which represses it); when
  mitogen falls, EZH2 (slow, τ≈14 h) de-represses with a lag → transcript "stays up longer."
- **Validation 28/29** (HU-G2 the lone structural miss, unchanged); levels held (Cd≈1.4–1.5, Cd_mRNA≈1.8, ~1.06
  ÷/day); after Hh-off the reservoir coasts ~1 more cycle (mRNA t½ ≈ 13 h ≈ the paper's ~11–12 h).
- Figure `fig_cyclind1_memory.py` rewritten (6 panels): A two-tier (slow mRNA → fast protein); B previous-cycle
  restriction point (coasting vs mRNA-fast/slaved); **C EZH2 IFFL buffers the dose-response** (rep factor 0.26→0.44 as
  mitogen falls → transcript held ~1.7× higher than a frozen-repression control); D reservoir "stays up" after the
  input drops (Gli resets <1 h, reservoir t½ ≈ 10 h, protein follows the reservoir); E cilium duty cycle; F robustness.
- Memory: `[[v44-cyclind1-hh-memory-carrier]]`; doc `docs/cyclind1_hh_memory_2026-07-21.md`.

## 2. Optimized the noisy Hh input — with a key correction to the low-pass story
Made the input distribution/temporal properties evidence-driven (`docs/hh-signalling-timing.md`) and characterized the
pathway's actual filter.
- **Input:** `L(t) = tonic(0.45) + slow OU drift (τ≈40 h, passband) + fast OU noise (τ≈0.5 h, stopband)`, re-centred so
  mean = tonic; CV ≈ 0.22. No pulse train. Cilium duty cycle (Cb-gated, ON-fraction ≈ 0.86); weak ~20 h adaptation
  (Ptch1_slow, kon=0.008/w=1.5 → Gli −35%, strength a scanned free parameter). All OFF by default → validation 28/29.
- **KEY CORRECTION:** the transfer function (sinusoidal SHH + lock-in) shows **v44 Gli is QUASI-STATIC** — ligand→Gli
  gain is **flat (~0.15) at every period** → Gli *tracks* the ligand and *carries* the fast noise; it does NOT
  low-pass. The earlier "Gli low-passes the noise" claim was wrong.
- **The low-pass filter is the CyclinD1 mRNA RESERVOIR** (τ = 15.9 h; the baked memory) + EZH2 IFFL. Measured
  ligand→Cd_mRNA transfer = a clean **first-order low-pass** (matches |H| = 1/√(1+(2πτ/period)²), τ≈16 h): fast noise
  (period < 3 h) attenuated **~12×** (HF amplitude) / up to **~36×** (Bode gain). This is the review's "temporal
  averaging by slow mRNA turnover" (§7) — the noisy-input story now ties directly to the memory mechanism.
- Figure `fig_realistic_mitogen.py` rewritten (6 panels): A input; B Gli decomposition (Gli tracks the noise);
  **C reservoir removes the fast noise (Gli 4.7% → Cd_mRNA 0.4% HF)**; **D Bode transfer function** (Gli flat vs
  Cd_mRNA first-order low-pass, measured + theory, cell-cycle contamination flagged); E cilium duty cycle;
  F heterogeneity (per-cell Gli noisy vs Cd_mRNA reservoir-smoothed).
- Memory: `[[v44-realistic-mitogen-input]]`; doc `docs/realistic_mitogen_input_2026-07-21.md`.

## The unified picture
Noisy paracrine Shh → quasi-static Gli that carries the noise → the slow **CyclinD1 mRNA reservoir** (the Hh memory)
low-pass-filters it, buffering the proliferation decision; EZH2's incoherent feedforward adds further buffering. Same
slow-transcript mechanism explains both the **previous-cycle restriction point** (Ho, Tsai & Stearns 2020) and the
**robustness to ligand noise** — one memory, two consequences.

## Open / deferred (scope-limited)
- Adaptation strength = the central GNP unknown (GNPs proliferate for days) → free parameter.
- Withdrawal/low-mitogen → differentiation trigger is MISSING (adaptation makes the cell quiescent, not
  commit-then-differentiate). Deferred per scope-limiting.
