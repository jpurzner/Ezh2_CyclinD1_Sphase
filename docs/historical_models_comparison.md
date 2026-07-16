# How the EZH2 → CyclinD1 repression model evolved (historical vs current)

*2026-07-16. Companion to `simulations/sim_ezh2_prc2_kinetics.py` / `fig_ezh2_prc2_kinetics.png` (panel D). All
CyclinD1 folds below are the model output on EZH2i in an MB cell (CyclinD1 relative to pre-EZH2i baseline, at 72h).
The kinetics of EZH2/PRC2 inhibition on CyclinD1 are the sharpest place the modeling generations diverge, so this
is the cleanest lens on "how historical models differ from the current one."*

## The foundation (unchanged across generations)
v44 = **Heldt et al. PNAS 2018** cell-cycle core (Cdk1/APC/Cdc20 + Rb/E2F restriction switch, conserved Rb/Cdh1/Gli
pools) + v44 additions: an explicit **mitotic switch** (Cdc25/Wee1 hysteresis), **HU→fork-speed** S-phase coupling,
a **Hedgehog/MYCN** module driving CyclinD1 transcription, cell-growth-gated commitment, and the **EZH2/H3K27me3**
epigenetic layer. What changed generation-to-generation is *how EZH2/PRC2/H3K27me3 represses CyclinD1* — everything
else (the engine) is shared, so behavioural differences isolate to this layer.

## The four generations of CyclinD1 repression

| Gen | flag | repressor | EZH2i CyclinD1 fold | de-repression kinetics |
|----|------|-----------|--------------------|------------------------|
| **1. Instant (direct-EZH2)** | `with_h3k27_dilution=False` | EZH2 level, algebraic `K/(K+EZH2·(1−EZH2i))` | **6.2× (step)** | INSTANT — EZH2i collapses repression the same timestep |
| **2. Explicit-mark memory** | `with_h3k27_memory=True` | the mark `H3K27_Cd` directly | **5.6× (slow)** | slow (~72h), set by mark demethylation `k_demeth` |
| **3. PRC2-occupancy (lumped)** | `with_h3k27_chain=False` | **PRC2 occupancy** `f0+(1−f0)/(1+(PRC2/K)^n)` | **1.25×** | fast + small (complex stays) |
| **4. Serial me-chain (current)** | *default* | PRC2 occupancy; mark = me3 via me0→me1→me2→me3 | **1.7×** | gradual (~12h) + me2→me3 lag |

## The three conceptual shifts (why each move was made)

1. **Instant → a real timescale (Gen 1 → 2).** The original model repressed CyclinD1 by the EZH2 *level*, so an
   EZH2 inhibitor de-repressed instantly and hugely. This is the "old wrong instant mechanism" — and it is why the
   original EZH2i-fold validation target (2.2×) was mis-set: it was fit to an instantaneous full collapse. Making
   the **H3K27me3 mark** (not EZH2 directly) the repressor gives de-repression a **real, mark-decay-limited
   timescale** (EZH2 stops writing → the mark decays over hours → CyclinD1 rises), matching the experimental
   observation that EZH2i de-repression is gradual, not a step.

2. **Mark-as-repressor → PRC2-occupancy-as-repressor (Gen 2 → 3), and the catalytic-inhibitor correction.** JP's
   EED-KO-vs-catalytic-dead-KO data settled the mechanism: it is **PRC2 complex OCCUPANCY** that represses
   transcription, not the H3K27me3 mark itself. So (i) a SAM-competitive **EZH2 inhibitor is CATALYTIC** — it blocks
   the SET-domain (writing) but the PRC2 complex stays bound, so repression falls only as the mark decays and is
   **partial** (an occupancy floor `a0` remains); (ii) an **EED KO removes the complex** → fuller de-repression
   (panel A: EZH2i 1.7× vs EED KO 2.25×). The mark became a **read-write amplifier** (EED reads me3 → recruits more
   PRC2) rather than the repressor. This is also what keeps CyclinD1 responsive to EZH2 over-expression/inhibition
   (all repression is EZH2/PRC2-mediated) while letting the Gli→Jmjd3/Kdm6b eraser keep the *mark* low in MB (ChIP
   MB<GNP) even as occupancy stays high via high EZH2.

3. **Lumped mark → serial me1/me2/me3 chain (Gen 3 → 4).** Resolving me0→me1→me2→me3 with a **slow, rate-limiting
   me2→me3** step gives me3 an accumulation **lag**, making it a division-timing readout (fast cycling dilutes me3
   out before the slow step finishes). Cosmetic for settled-state validation but adds the mechanistic origin of
   "persistence" and a de-novo lag.

## Behavioural upshot (panel D)
- **Instant (Gen 1)** overshoots wildly and has no kinetics — the model artefact that motivated the whole rebuild.
- **Memory (Gen 2)** has kinetics but an over-large fold (the mark is too strong a lone repressor).
- **PRC2/chain (Gen 3–4, current)** give a **gradual, PARTIAL** de-repression (~1.25–1.7×) because the complex
  stays bound under catalytic EZH2i — the data-consistent behaviour (EZH2i < EED KO; de-repression is not a step).
- **The mark AMOUNT matters** (panel C): more H3K27me3 read-write = more baseline repression = a larger, still-
  gradual de-repression when lifted. In MB the Gli→Jmjd3 eraser keeps the mark (and thus the fold) modest.

## Other (non-repression) model evolution, for completeness
- **Single-step → two-step Rb** (mono-P by CyclinD-CDK4/6 → hyper-P by CyclinE/A-CDK2): required for the
  dose-response and pRb-marker work (`with_two_step_rb`, default).
- **EZH2 division-halving → concentration convention** (`with_ezh2_conc`): soluble EZH2 = concentration, preserved
  at division.
- **Transient-G0 reframe + CDK4/6i dose-response** (2026-07): CDK4/6-activity-gated p27 clearance, pRb-Ser807/811
  G0 marker, and CDK4/6i modelled as a graded dose-response (not a resistant clone) — see `docs/transient_g0_synthesis.md`.
