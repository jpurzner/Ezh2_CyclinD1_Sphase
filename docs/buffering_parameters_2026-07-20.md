# Parameters that support buffering + the summary equation (2026-07-20)

*JP: determine the parameters that support buffering, generate parameter plots + a summary equation. Files:
`sim_buffering_parameters.py`, `fig_buffering_parameters.py`. GNP, current v44 default. Buffering = the EZH2/H3K27me3
suppression of CyclinD1. Metrics: static B = 1 − R with R = Cd_ON/Cd_OFF (repression depth); dynamic = mark low-pass
under an Shh oscillation.*

## Which parameters support buffering
- **`f0_prc2` (repression floor) — the clean, monotonic knob.** R tracks it almost exactly (R ≈ f0 when the mark is
  engaged): buffering fell 0.82 → 0.20 as f0 went 0.05 → 0.80. It sets the **ceiling**: max buffering = 1 − f0.
- **`K_prc2` (repression sensitivity) — the engagement knob.** Lower K → the repression engages at a lower mark
  occupancy → more buffering (buffering 0.83 at K=0.0005, falling as K rises). This was the knob that resolved the
  "both" trade-off last night.
- **`n_prc2`** sharpens the engagement.
- **Mark-LEVEL knobs (`a_rw_prc2`, `del_mk`, `kEZbas_Cd`, `kEZE2f`, `kDeEZ`)** push the PRC2 occupancy P up (more
  EZH2/mark = more repression), but the P ↔ CyclinD1 loop is **BISTABLE**, so these don't move buffering smoothly —
  they **SWITCH** it between a low-repression branch and a latched high-repression branch (the non-monotonic panels).
- **`del_mk` also sets the DYNAMIC buffering** (the mark's low-pass corner frequency): slower mark → deeper filtering
  of a given fluctuation frequency (fluc-atten rises with `kDeEZ` and falls with slower `del_mk`).

## Summary equation (corrected 2026-07-20 to include the CELL-CYCLE / S-phase-dilution term)

The original version below missed the **cell-cycle duration Tc** (JP). H3K27me3 is HALVED at S-phase (replication)
and re-methylated over the rest of the cycle, so the mark's decay has TWO parts — intrinsic turnover `del_mk` plus a
**replication-dilution rate `ln2/Tc`** — and the steady mark is set by the writing/decay balance:
```
    P*  ≈  w / (del_mk + k_jmjd3_gli·Gli1 + ln2/Tc)     w = writing rate (∝ EZH2·PRC2, proliferation-coupled)
    δ_eff = del_mk + k_jmjd3_gli·Gli1 + ln2/Tc          effective mark decay (THREE terms)
```
(2026-07-21: the middle term — the mitogen-gated Gli→JMJD3/KDM6B **eraser** — was re-enabled nerfed and belongs
in the denominator; it is the mitogen-dependent loss channel that gives ChIP MB<GNP. See `v44-four-mark-mechanisms`.)
Validated: Mk ≈ w/(del_mk+ln2/Tc) fits with corr 0.95 (`sim_mark_cellcycle_recovery.py`, `fig_mark_cellcycle_recovery.py`).
The two decay terms are **comparable** (del_mk ≈ 0.028/h vs ln2/Tc ≈ 0.032/h at Tc=22 h); **dilution dominates at fast
cycling** (short Tc), turnover at slow cycling, crossover **Tc ≈ ln2/del_mk**. Since **Tc is mitogen-dependent** (the
graded G1-lengthening: subthreshold mitogen → longer Tc), the mark level P* — and thus the buffering — is coupled to
the cell-cycle-length response. **At max proliferation Tc is short (heavy dilution), so the mark peaks at max division
only if the writing w outpaces ln2/Tc** (the writer–dilution balance behind the ChIP-seq).

**Static buffering** (dose-response flattening; mitogen-requirement shift):
```
    B = 1 − R(P*),   R(P*) = f0_prc2 + (1 − f0_prc2) / (1 + (P*/K_prc2)^n_prc2),   P* = w/(del_mk + ln2/Tc)
```
- Ceiling: **max B = 1 − f0_prc2** (when the mark saturates the repression, P* ≫ K).
- Engages once **P* ≳ K_prc2**; n_prc2 sets the sharpness. Mitogen-requirement shift = **1/R** (validated: GNP shift ≈ 1/R).

**Dynamic buffering** (fluctuation low-pass; the mark is a first-order filter, corner **ω_c = δ_eff = del_mk + ln2/Tc**):
```
    atten(ω) = ω_c / √(ω_c² + ω²),   ω_c = del_mk + ln2/Tc
```
- The mark filters on the **CELL-CYCLE timescale** (not just del_mk): fluctuations slower than ~Tc are buffered; faster pass through to CyclinD1.

**Combined:**
```
    B(ω) ≈ [1 − R(P*)] · ω_c / √(ω_c² + ω²),   P* = w/(del_mk + ln2/Tc),   ω_c = del_mk + ln2/Tc,   Tc = Tc(mitogen)
```
- Slow changes (ω→0, withdrawal / dose-response): full static buffer 1 − R(P*).
- Fast fluctuations: absorbed by the slow mark, attenuated ∝ ω_c/ω.

## Caveat
The bistability of the P↔CyclinD1 loop means the mark-level knobs give switch-like (not graded) control, and R can
jump between ~0.1 (latched) and ~1.0 (off) with small changes near the bifurcation. The clean, dial-able buffering
controls are **f0_prc2 (ceiling) and K_prc2 (engagement)**; the mark-level knobs set which side of the switch you are on.
