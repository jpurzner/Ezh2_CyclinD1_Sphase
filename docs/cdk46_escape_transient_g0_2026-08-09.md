# Reinstating MB transient G0 — the CDK4/6-alone escape (2026-08-09)

**JP:** "MB cells DO have transient G0 arrest; there is huge CyclinD1/CDK6 which should pull the cells out
of G0." This doc lays out the options and the recommended one.

## Why deep MB arrests are currently PERMANENT (the problem)

The permanent lock is the **two-step Rb R-point trap** (see `permanent_lock_decomp.py`): once E2f collapses,
it is trapped on mono-phosphorylated Rb; CyclinD–CDK4/6 can only *mono*-phosphorylate Rb, and the
hyper-phosphorylation that releases E2f needs CyclinE/A–CDK2, which needs E2f — a deadlock. So in the model
the **huge CyclinD1/CDK6 CANNOT pull cells out** (restoring CDK6 fully does not release the lock), and deep
MB arrests go permanent — which is wrong for the transient Atoh1+ G0.

## The options

1. **CDK4/6-alone Rb-hyper-P escape (Yang 2020 eLife 44571 / Chung 2019)** — *recommended.* Let
   CyclinD–CDK4/6 *weakly* hyper-phosphorylate Rb (`with_cd_hyper_escape`, `w_cd_hyper`). The huge MB
   CyclinD1/CDK6 then pulls the arrested cell back out once it rebuilds → **reversible transient G0**. Directly
   implements JP's insight; biologically grounded; **31/32**; tunable + self-bounded.
2. **Two-compartment amplifier** (`with_proximal_distal`) — *ruled out for this.* It *represses* CyclinD1/CDK6
   in arrest (the opposite direction) → inert (safe regime) or a permanent lock (excluded). Does not give
   transient G0.
3. **Fluctuating / periodic CDKI (p27)** — makes the *arrest itself* transient (p27 pulses); re-entry still
   needs option 1 (or the drive can't pull out). A frequency knob on top of option 1.
4. **Birth-p27 / CDKI-variance tail** — sets *which* MB cells arrest (the rare high-CDKI tail,
   `mb-atoh1-g0-rare-cdki-tail`); combined with option 1 for reversible re-entry → the population fraction.

The escape (1) is the **re-entry** mechanism that makes any arrest transient via the high drive; a brake
(3/4) sets the arrest frequency. Option 1 is foundational and is the direct answer to JP's question.

## Result — the escape reinstates MB transient G0 at 31/32

Fixed `with_cd_hyper_escape` for the integrated model (`_dact` now matches the option-B drive + cdk6-sinked
INK4 brake; added `CD_HYPER_ESCAPE` env flag). At fixed MB CDKI, divisions in 60000 min (period h):

| w_cd_hyper | validation | 1× | 4× | 6× | 8× |
|---|---|---|---|---|---|
| OFF | 31/32 | cycle 24h | slow 52h | **permanent** | permanent |
| 0.05 | **31/32** | 24h | 52h | **transient 70h** | permanent |
| 0.10 | **31/32** | 24h | 52h | transient 70h | permanent |
| 0.20 | **31/32** | 24h | 52h | 71h | **transient 89h** |
| 0.40 | 30/32 | 24h | 52h | 71h | 89h *(breaks GNP-SHH arrest)* |

- The escape converts the **permanent** deep-CDKI arrest into **multi-periodic transient-G0 cycling**:
  CyclinD1 accumulates in the G0 dwell → pulls the cell out → it divides → re-arrests → re-enters. (fig B)
- **Drive-graded dwell**: deeper arrest → longer G0 dwell (52 → 70 → 89 h). Higher CyclinD1/CDK6 (MB) pulls
  cells out faster — exactly JP's picture.
- **Validation-neutral (31/32)** up to `w_cd_hyper ≈ 0.2` (only new-vs-baseline miss is the pre-existing MB
  HU G2). **Better than the amplifier (30/32).** Self-bounded: `> 0.2` lets the escape fire without mitogen
  → breaks the GNP-SHH=0 arrest, so the data caps the strength.
- Baseline cycling (24h) + graded G1 (52h) untouched — the escape only rescues the deep arrests.

## Recommendation

**Option 1 (CDK4/6-alone escape) is the mechanism for MB transient G0**, and it is the opposite of the
amplifier: the amplifier *represses* the drive (→ permanent lock, excluded), while the escape lets the drive
*pull cells out* (→ transient G0, 31/32). Nothing baked — flag default OFF; JP to decide whether to bake at
`w_cd_hyper ≈ 0.1` (rescues 6× MB arrests, comfortably within the GNP-SHH bound). Files:
`escape_probe.py`, `escape_eval_wf.js`, `fig_v44_cdk46_escape.py`, `permanent_lock_decomp.py`.
