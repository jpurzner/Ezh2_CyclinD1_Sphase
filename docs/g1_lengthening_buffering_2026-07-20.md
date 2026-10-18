# Graded G1 lengthening (no transient G0) as the EZH2/CyclinD1 buffering mechanism (2026-07-20)

*JP: EZH2 has a mitogen-independent BASELINE plus an E2F/S-phase component, so it responds to mitogen SLOWLY — it
falls only when the cycle spends less time in S, i.e. when the cycle LENGTHENS. GNPs have NO transient G0, so the
slowdown must be graded G1 lengthening, not a discrete arrest. JP directive: "GNPs likely have the capacity to
lengthen the cell cycle without using a transient G0 — introduce this in the model. Not required to be frequent
(unperturbed WT shows little lengthening)." Grounded in JP's developmental-kinetics research review (Contestabile
2009 GNP: Tc 16.25 h, G1 7.63 h, S 7.11 h, G2 1.18 h, M 0.33 h; S-phase index ~40%).*

Files: `simulations/sim_g1_lengthening_buffer.py`, `fig_g1_lengthening_buffer.py`. EXPLORATION — default untouched.

## The buffering logic
- **CyclinD1 tracks mitogen FAST** (withdrawal: half-falls in ~1.4 h in the prototype).
- **EZH2 = baseline + S-phase(E2F) component**, so it tracks mitogen only through the CELL CYCLE. It does NOT drop
  when CyclinD1 drops; it drops only when the cell does fewer S-phases, i.e. when the cycle lengthens (prototype:
  EZH2 half-falls in ~14.8 h). The current bake already has this structure (kEZbas baseline + kEZE2f·E2f·(CycE/CycA);
  measured EZH2 baseline ≈ 1.57, rising to ≈ 2.78 with mitogen — never zero, matching JP).
- So EZH2 (and the H3K27me3 it writes) **low-pass-filter mitogen**, and the filter time-constant IS the
  cell-cycle-length response. A fast dip that doesn't lengthen the cycle passes EZH2 by; a sustained subthreshold
  drop lengthens the cycle and pulls EZH2 down slowly.

## The gap in the current v44 model (measured this session)
Under graded SHH the current default holds **Tc ≈ constant (~20–22 h)** across the whole proliferative range, then
**abruptly ARRESTS** (0 divisions) below SHH ≈ 0.08. That abrupt arrest is the transient-G0 behaviour JP says does
not occur in GNPs. **The model has no graded G1-lengthening capacity** — so it cannot produce the slow EZH2 response
by lengthening; it can only keep EZH2 high (cycling) or drop it fast (arrest). This is the thing to change.

## Prototype (GNP-grounded)
`Tc(CyclinD1)` lengthens G1 smoothly as CyclinD1 falls (TG1 from 7.65 h up, arrest only as D→0), so Tc runs **17 h
(saturating) → 68 h (deep subthreshold)** continuously — no discrete G0 (Fig A, vs the current model's flat-then-cliff).
Withdrawal then shows the **timescale separation** that is the buffering: CyclinD1 ½-fall 1.4 h vs EZH2 ½-fall 14.8 h
(Fig B). Through a transient 60 h subthreshold dip, graded lengthening retains EZH2 (dips only to 2.32, near
baseline) **better than an abrupt G0** (dips to 1.62) — the graded slowdown is a gentler, more buffered response than
a discrete arrest (Fig C).

## The proposed v44 change (for review)
Replace the sharp commitment-threshold arrest with a **mitogen/CyclinD1-graded G1 duration**: below the current
commitment point, cells should CYCLE SLOWLY (long G1), not arrest — arrest only at very low mitogen. Concretely, make
the G1→S progression rate (Rb-phosphorylation / commitment / origin-licensing) a graded, saturating function of
CyclinD1 so Tc rises smoothly as mitogen falls, with the GNP Tc_min ≈ 16.25 h at saturation. This softens the
bistable commitment into a graded slowdown and gives EZH2 its slow (cycle-length-gated) response — the enabling
mechanism for the buffering — WITHOUT a transient-G0 pool. Latent under normal WT (little lengthening at high
mitogen), engaged under subthreshold mitogen. Needs re-validation against the phase-fraction targets (S% 43.7%,
G2+M short) using the GNP numbers.

## Complementary piece (mark direction / ChIP-seq)
With a stable mark + substantial EZH2 baseline + division-coupled dilution, the *steady-state* mark comes out
slightly HIGHER in slow cells (dilution falls faster than EZH2 when cycling slows). So matching the ChIP-seq
(H3K27me3 maximal at max division, declining with differentiation) needs the **active differentiation erasure** from
[[prc2-buffering-memory-topologies]] — G1 lengthening supplies the *buffering* (slow EZH2), the erasure supplies the
*ChIP-seq direction*. They are complementary, not alternatives.

## Discriminating experiments (from JP's research review, highest value)
1. **Extend cumulative EdU past Tc (≥36–48 h at P2)** to MEASURE the growth fraction rather than assume it — the
   reported G1-lengthening effect is exactly what a growing non-cycling fraction would counterfeit (all published GNP
   G1 numbers derive G1 by subtraction under GF=1). If the plateau < 1 and lower when perturbed, graded G1 lengthening
   is real; if not, the buffering should instead be built on a variable-duration CDK2low pre-G1 interval.
2. **Phospho-Rb (S807/811) IF on P7 EGL** for the committed fraction in fixed tissue (Moser 2018).
3. **DHB CDK2 sensor in GNPs under graded Shh** — the CDK2low fraction + residence-time distribution is the exact
   input this module needs.

Caveat (JP's review): the GF-unbiased measurement (Legué clone S-phase index) RISES toward differentiation (46.9% vs
34.6%), arguing that GNPs speed up in final divisions rather than progressively lengthening G1 — so the lengthening
capacity should be a mitogen-driven option, not the developmental default.
