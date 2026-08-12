# Permanent cell-cycle exit vs transient G0: H3K27me3 dilution × CyclinD1/CDK6 drive (2026-08-09)

Overnight exploration (JP): "how H3K27me3 accumulation vs loss through replicative dilution, and the
CyclinD1/CDK6 drive difference between GNPs and MBs, determine permanent cell-cycle exit vs transient G0."
On the baked default (integrated governor + CDK4/6-alone escape, 31/32).

## Headline

Three factors act, on **different things**:

1. **CyclinD1/CDK6 DRIVE sets the FATE.** MB's high drive creates a wide **transient-G0 escape regime**;
   GNP's low drive has **none** — it goes straight from cycling to **permanent exit**.
2. **The H3K27me3 mark's static repression TUNES the fate boundary** (~2× CDKI): the mark lowers the
   effective drive, pushing cells toward permanent exit ~2× earlier.
3. **Replicative dilution of the mark is a CYCLING-RATE sensor, NOT a fate switch.** It sets the mark (and
   cycle period) in *cycling* cells, but does **not** move the permanent/transient boundary — because that
   fate is decided in the *arrested* state, where there is **no replication to dilute**.

## The data

### Drive sets the fate (recon_permanent_transient.py) — CDKI drive-reduction axis

| CDKI | 1× | 2× | 3× | 4× | 5× | 6× | 7× | 8× |
|---|---|---|---|---|---|---|---|---|
| **GNP** | cyc | cyc | cyc | cyc | **PERM** | perm | perm | perm |
| **MB** | cyc | cyc | **TRANS** | trans | trans | trans | trans | **PERM** |

GNP (low CyclinD1 ~3 / CDK6 ~1) has **no transient regime** — cycling → permanent at 5×. MB (CyclinD1 ~7.8 /
CDK6 ~5) has a **wide transient-G0 regime** (3–7×, the CDK4/6-alone escape pulls it back out; period
43→70 h) before permanent at 8×. The **Hedgehog axis** is even starker: GNP permanently exits at SHH=0
(drive collapses), while **MB keeps cycling even at SHH=0** (MYCN floor + CDK6 reservoir = vismo-resistance).

### The mark tunes the boundary (mechanism_mark_dilution.py)

Turning the H3K27me3/PRC2 repression OFF (`f0_prc2=f0_cdk6=1`) shifts the permanent-exit onset **~2× CDKI
later**: GNP 5→7×, MB 8→10×. So the mark's *static* repression lowers effective drive and brings cells ~2×
closer to permanent exit.

### Dilution is a rate-sensor, not a fate switch (mechanism_mark_dilution.py)

Varying `k_dil_S` (dilution ×2 / ×0.5 / OFF) leaves the permanent/transient boundary **completely
unchanged** (GNP 5×, MB 8× in every case). The mark is near-saturated in cycling cells (Mk ≈ 0.94–0.96), so
dilution only nudges it; and the *fate* is decided in the arrested state, where replication (hence dilution)
has stopped. Dilution's real role — reducing the mark as the cycle speeds up — governs the **cycling
proliferation rate** (the graded slowdown / rate-sensor), a distinct phenomenon from the fate.

### Why the fate is drive-decided — the arrested state

In a fully-arrested cell the mark **drops** (Mk 0.95→0.85–0.88) because the EZH2 writer is E2f-gated and
collapses when the cell stops. What remains is a two-step-Rb R-point trap; the *only* way out is the
CDK4/6-alone escape, which needs enough CyclinD1/CDK6. MB retains it (CDK6 ~5–9 in the transient state →
escapes); GNP does not (CDK6 ~1 → permanent). This is why the same perturbation gives GNP permanent exit and
MB transient G0.

## Interpretation (answering JP)

- The **permanent-vs-transient dichotomy is a DRIVE phenomenon**, not a mark phenomenon: MB transiently
  arrests because its huge CyclinD1/CDK6 escapes the R-point trap; GNP permanently exits because its low
  drive cannot. The mark **modulates** this (~2× CDKI) by adding repression, and would push a borderline
  cell over into permanent exit.
- **Replicative dilution of H3K27me3 does not set the fate** in this model — it sets the *cycling rate*
  (the mark-as-speed-limiter). The two things JP named therefore operate at different stages: dilution in
  the cycling/slowing phase (rate control), drive in the arrested phase (fate).
- **Tested the unsaturated regime (writer ×0.15):** the cycling mark becomes strongly dilution-dependent
  (Mk 0.38 at dil ×0.5 vs 0.15 at dil ×2) — yet the permanent/transient boundary is **still identical**
  (MB permanent at 10× for both). So dilution does **not** feed the fate *even when the mark is unsaturated
  and highly dilution-sensitive in cycling*. The reason is structural: the fate is decided in the arrested
  state, where replication (hence dilution) has stopped and the mark re-equilibrates to its
  writer-vs-eraser value. **Replicative dilution is decoupled from the fate regardless of mark saturation.**
  (Reducing the writer did shift the boundary 8×→10× — the mark's *static* repression tuning again, since
  less writer = less repression = higher effective drive.)

## Discriminator experiment

The model predicts: (1) GNP has no transient-G0 escape regime (deep arrests are permanent = differentiation);
MB does (deep arrests re-enter). (2) Blocking H3K27me3 (EZH2i) shifts the permanent-exit threshold ~2× —
i.e. EZH2i lets cells tolerate more CDKI/less mitogen before permanent exit. (3) Changing the proliferation
rate (hence dilution) changes the mark level but not whether a given arrest is permanent vs transient.

## Files
`recon_permanent_transient.py`, `mechanism_mark_dilution.py`, `fig_v44_permanent_transient.py`.
