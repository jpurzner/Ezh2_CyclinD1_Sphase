# Withdrawal memory via graded G1-lengthening — overnight exploration (2026-07-24)

*JP (going to sleep): "the withdrawal memory would require the reduction of ezh2 so the mitogen would need to extend
g0, the pool has to keep them dividing long enough to lengthen g1 but not to lead to exit — can you explore parameter
combinations that would allow this." Grounded in: CyclinD1 is 11th-percentile expressed, CyclinD2 top 0.3% (bulk
RNA-seq) → the large-pool premise is real. `simulations/sim/fig_g1_lengthening_withdrawal.py`. NOT baked.*

## The lever was already in the model (latent)
G1 is **growth-limited** in v44 (committed cells wait for mass = M_size). The growth rate already has a
CyclinD1-dependence:
```
mu_eff := mu*(mu_min_frac + (1 - mu_min_frac)*Cd^n_g1len/(K_g1len^n_g1len + Cd^n_g1len))
Growth: => mass; Cell*mu_eff*mass
```
so growth slows as CyclinD1 falls → G1 lengthens. It's baked with `K_g1len = 0.6` (below the operating CyclinD1
~1.5), so it only engages near arrest → the current **abrupt** behavior. NB kPhRbCd/K_CdRb (CyclinD1→Rb) do NOT set
Tc — growth does — so those are the wrong knobs.

## Result — graded G1-lengthening (validation-compatible)
Tuning **`mu_min_frac` 0.4→0.1, `K_g1len` 0.6→0.8** (n_g1len=4, kPhRbCd 0.35 unchanged): the steady period lengthens
**smoothly 24 h (Shh 0.5) → 84 h (Shh 0.08), then arrests at Shh 0** — Shh-dependent. In a gradual withdrawal the
division intervals **ramp up** (25→32 h in the transient; the steady dose-response is larger) vs the current model's
even-then-stop. **Validation-neutral: 28/29** (GNP period 24 h passes, CyclinD1 MB/GNP 4.76, folds intact). This is
JP's "lengthen G1 but not exit," realizing the goal of `gnp-g1-lengthening-buffering`.

Caveats: rate-dependent (a faster withdrawal outpaces the lengthening cycle → shorter ramp, as in
`matched-withdrawal-rate-dependent`); pushing kPhRbCd→0.55 gives bigger ramps but tips into Shh-INDEPENDENCE.

## The EZH2-reduction driver — NOT supported as-is
JP's proposed driver is that **EZH2 reduces** during withdrawal → de-represses the CyclinD1 pool → sustains division.
The model does **not** do this: EZH2 stays **flat** during withdrawal (panel C) because it turns over slowly (τ~14 h)
and tracks the **cycle**, not the mitogen — it "banks." Making it reduce requires faster EZH2 turnover, which
**breaks validation (23/29)** (CyclinD1 folds, period 30 h, Palbo) because `kDeEZ` is pinned by the HU-S in-S boost +
G0-withdrawal data. (An earlier apparent "more divisions" from fast turnover was an artifact of an uncompensated
*lower* EZH2 level, not the reduction dynamics.) So the graded phenotype comes from **growth-slowdown**, not EZH2
reduction; JP's specific EZH2 driver is an unresolved tension (would need a joint re-optimization, and may fight the
HU-S accumulation data).

## Recommendation for JP
- The **growth-slowdown route** (mu_min_frac 0.4→0.1, K_g1len 0.6→0.8) cleanly gives graded G1-lengthening and is
  validation-compatible — ready to bake if you want the withdrawal memory as graded lengthening.
- The **EZH2-reduction route** (your proposed mechanism) needs a joint re-opt of kDeEZ/kTlEZ against the CyclinD1
  folds + period + the HU-S data before it's viable. Decision pending.
