# Biologically realistic Shh mitogen input + Gli-signal decomposition (2026-07-21, updated 07-24)

*JP: "finish a biologically realistic mitogen input where the distribution and temporal properties are somewhat
evidence-driven, then show what it looks like with a decomposition of what is making the Gli signal look that way";
then "optimize your noisy Hh input." Grounded in `docs/hh-signalling-timing.md`. `simulations/sim_realistic_mitogen.py`
+ `fig_realistic_mitogen.py`. All new model features are OFF by default → validation unchanged (28/29).*

## The input model  `Shh_effective(t) = L(t) · G_cilium(cell-cycle) · A(t)`
Three evidence-driven components (from the review):
1. **Ligand L(t)** — paracrine, Purkinje-derived, **tonic / slowly-varying**, with two OU components: a **slow drift**
   (τ≈40 h, the "real" slow variation, PASSBAND) + a **fast noise** (τ≈0.5 h, sub-hour, STOPBAND — the thing the
   pathway should filter). NOT a fast pulse train (that is ERK, not Hh). The components are **re-centred so the mean =
   tonic (0.45)** — a controlled proliferating setpoint — giving CV ≈ 0.22 (slow 0.16 + fast 0.14). The GNP ligand
   fluctuation timescale/amplitude is **unmeasured** → amplitudes are a stated assumption; the point is the *timescale
   separation* (fast in the stopband, slow in the passband).
2. **Cilium duty cycle G_cilium** — ~94% of EGL GNP progenitors ciliated but cilia **resorb pre-mitotically**, so Shh
   reception is phase-gated ON in G1/S, OFF through G2/M, restored after division. Model term
   `cilium_gate := 1 − w_cilium·Cb^n/(K_cil^n+Cb^n)` on SHH→Ptch1 reception (gated by CyclinB = G2/M). ON-fraction ≈
   0.86 (K_cil=0.70). `w_cilium=0` = off (default).
3. **Adaptation A(t)** — the ~20 h temporal adaptation (Ptch1 feedback + Gli2-down + GliA-destabilization). Model term
   `Ptch1_slow` (slow, Gli²-driven, adds to Smo inhibition). Tuned **WEAK** here (`k_ptch1_slow_on=0.008,
   w_ptch1_slow=1.5` → Gli drops ~35% over ~20 h, sub-lethal) because **GNPs proliferate for days under sustained Shh**
   → strength is the central GNP unknown, a **scanned free parameter**, exactly as the review recommends.

## KEY correction (2026-07-24): Gli is quasi-static — the LOW-PASS filter is the CyclinD1 mRNA RESERVOIR
Measuring the pathway's transfer function (sinusoidal SHH, lock-in at the drive frequency) shows **v44 Gli is
quasi-static**: the ligand→Gli gain is **FLAT (~0.15) at every period** from 0.5 h to 160 h. So Gli does NOT low-pass —
it **tracks the ligand and carries the fast noise** (consistent with the review: Gli protein τ ≈ 40 min–2 h « the
~22 h cycle; no intrinsic Gli oscillation). The earlier figure claim that "Gli low-passes the noise" was WRONG.

The low-pass filter that actually buffers the **proliferation decision** is **downstream**: the slow **CyclinD1 mRNA
RESERVOIR** (τ = 15.9 h from `k_Cd_mRNA_deg=0.001`; the baked Hh memory) plus the EZH2 incoherent feedforward. The
measured transfer function ligand→Cd_mRNA is a clean **first-order low-pass** (matches |H| = 1/√(1+(2πτ/period)²),
τ≈16 h): fast fluctuations (period < 3 h) attenuated **~12×** in high-frequency amplitude and up to **~36×** in Bode
gain (fast vs slow limit), while the slow drift passes. This IS the review's "temporal averaging by slow mRNA
turnover" (§7) and it ties the noisy-input story directly to `[[v44-cyclind1-hh-memory-carrier]]`.

## What the figure shows (6 panels)
- **A · Input:** paracrine Shh = tonic + slow drift (τ≈40 h) + fast noise (τ≈0.5 h); CV ≈ 0.22.
- **B · Gli decomposed** (add one component at a time): (1) tonic → **flat** Gli; (2) + ligand noise → Gli **TRACKS**
  it (quasi-static, noisy); (3) + cilium duty cycle → **G2/M dips**; (4) + weak ~20 h adaptation → slow **setpoint
  drop** (~35%). Gli carries the noise — it does not filter.
- **C · The reservoir removes the fast ligand noise:** high-frequency (period < 3 h) component of Gli (amp ~4.7%) vs
  CyclinD1 mRNA (amp ~0.4%) → **~12× attenuation** (isolated from the slow cell-cycle oscillation, which dominates
  Cd_mRNA's *total* variance and is NOT noise).
- **D · Transfer function (Bode):** ligand→Gli is flat (quasi-static); ligand→CyclinD1 mRNA is a first-order low-pass
  (τ≈16 h theory curve overlaid on the measured points). The cell-cycle period (~22.8 h) is shaded — the lock-in is
  contaminated there (the 24 h outlier).
- **E · Cilium duty cycle:** `cilium_gate` drops through each G2/M (tracking CyclinB), restored at division;
  ON-fraction ≈ 0.86.
- **F · Cell-to-cell heterogeneity:** per-cell Gli (noisy) + per-cell CyclinD1 mRNA (reservoir-smoothed) — different
  ligand realizations + per-cell CyclinD1/EZH2 draws; precision emerges downstream (the review's point).

## Scope note
This finishes the **mitogen input** and its characterization. The two new terms (cilium gate, slow adaptation) are OFF
by default so the calibrated model is untouched; they are turned on only for this figure. Adaptation strength is
deliberately a free parameter per the review's central caveat. The unresolved withdrawal→differentiation trigger
(adaptation makes the cell quiescent, not commit-then-differentiate) is deferred (scope-limited).
