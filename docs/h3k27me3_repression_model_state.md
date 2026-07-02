# H3K27me3 repression of CyclinD1 — model state & outstanding items
*Checkpoint after the replicative-dilution work. Captures the model variants, the calibrated leaky
dilution module, validation status, key findings, and what remains — so nothing is dropped.*

---

## 1. Three repression mechanisms exist; the AUM module is now the DEFAULT

| Build flag | Repression of CyclinD1 | Role | Validation |
|---|---|---|---|
| **`with_h3k27_dilution=True` (DEFAULT)** | **mean-field AUM module** (Purzner spec + lit review) — see §2 | the mechanistic, biology-faithful model, now the default | **22/27** |
| `with_h3k27_dilution=False` | `K_rep/(K_rep + EZH2·(1−EZH2i))` — EZH2 level represses directly, instantaneous | legacy direct-repression model | **22/27** |
| `with_h3k27_memory=True` | mark `H3K27_Cd` tracks EZH2, `K_rep/(K_rep+H3K27_Cd)`; EZH2i → mark decays over ~24–48 h | de-repression KINETICS variant; used by `fig_v44_ezh2i_rescue_kinetics` (takes priority over dilution) | 22/27 |

Promoted June 2026: `build_model_v44()` now defaults `with_h3k27_dilution=True`. `validate_v44.py` tests it
by default (legacy via `H3K27_DILUTION=0`). All registry figures were regenerated on the AUM model; the
~12 feedback figures were remapped from the legacy `K_EZH2_repression` knob to `f0_mk` ("feedback OFF" =
`f0_mk=1.0`), and `sim_ezh2_phaseplane` was reworked into the (H3K27me3-mark Mk, CyclinD1) plane.

---

## 2. The mean-field AUM module (the current best mechanistic model)

State `Mk ∈ [0,1]` = H3K27me3 occupancy over the ~7 kb / ~35-nucleosome bivalent *Ccnd1* domain. This
is the **mean-field reduction of the Berry–Dean–Howard A/U/M per-nucleosome model**: one state (the
U↔M arm with read-write); the "active/A" (H3K27ac) state's antagonism is supplied not by a separate
acetyl species but by the **transcription→PRC2 reciprocal arm** below. Well-mixed mean-field is valid
at N≈35 because PRC2 read-write recruitment is non-local (Dodd 2007).

- **Continuous methylation** (one PRC2 catalytic-rate term, all factors multiply):
  `Mk_methylation = EZH2·(1−EZH2i)·(k_w·Mk + k0)·(1−Mk)·(1 − g·Cd_mRNA^p/(K_tx^p+Cd_mRNA^p))`
  - autocatalytic read-write `k_w·Mk` (EED reads me3) + de-novo nucleation floor `k0` (reseeds a halved
    locus); both on substrate `(1−Mk)`.
  - **EZH2-scaled** (the whole term ∝ EZH2): all PRC2 catalysis is EZH2-dependent, so when a cell exits
    the cycle and EZH2 falls (E2F-coupled), methylation weakens and the arrested mark does NOT run away —
    this is what keeps vismo-arrested CyclinD1 in band (was over-silenced when only read-write was scaled).
  - **blocked by EZH2i** (`1−EZH2i`): all EZH2 repression flows through the mark.
  - **transcription→PRC2 reciprocal arm** `(1 − g·tx-Hill)`: ongoing *Ccnd1* transcription (proxied by
    `Cd_mRNA` = nascent output) evicts PRC2 / inhibits its HMTase (lit-review §3.3). This is a
    double-negative ⇒ **positive feedback**: the bivalent, sub-saturating M_ss now EMERGES from the
    antagonism balance (it is not imposed), the switch is sharpened, and vismo arrest is self-reinforcing.
    `g` is the deliberately WEAK eviction strength — the measured 5.07× fold pins it to the weak regime.
    `g=0` recovers the prior read-write-only module exactly.
- **Turnover** `δ·Mk` (demethylase + histone exchange), **lit-review-paced**: `δ=0.0007` gives a
  de-repression t½ ≈ 16 h (was 9.9 h), inside the lit-review τ_restore band (10–18 h) and far less
  turnover-dominated, so de-repression is now genuinely dilution-paced.
- **Discrete**: replicative **halving at early S** (`Dna>0.05` event, once/cycle) → the emergent cycle
  period sets the dilution frequency (faster cycle ⇒ more dilution).
- **Readout (LEAKY)**: `R(Mk) = f0 + (1−f0)/(1+(Mk/K)^n)` multiplying the full Gli+MYCN drive — throttles
  initiation/burst frequency with a residual floor **f0** (loaded Ser5P Pol II → leaky firing); it
  DAMPENS, does not lock out. Mark inherited through division = the memory.

**Calibrated params:** `k_w=0.00160, k0=0.000258, δ=0.00070, K=0.305, n=4.15, f0=0.233,
g=0.30, K_tx=5.0, p=2.0`. (Calibration: `simulations/calibrate_h3k27_aum.py` — slowest turnover that
keeps all four CyclinD1 targets in band, at HHi=1.0 with validate's stiff-retry.)

**Behavior:** GNP & MB cycle, GNP+HHi & MB+HHi arrest (residual transcript; vismo removes the Gli DRIVE,
not via lockout), **MB+HHi+EZH2i RESCUES** (6 div), GNP+EZH2i de-represses. MB/GNP fold ~3.9, EZH2i fold
~1.7, both in band. Three emergent behaviors verified:
- **EZH2i de-repression**: CyclinD1 ↑ to ~1.8× by 24 h, plateaus ~48 h; mark t½ ≈ 18 h (matches the
  lit-review ~24 h and JP's "clear by 24 h, max by 48 h").
- **Self-reinforcing arrest** (the reciprocal arm's signature): after vismo, Mk *rises* 0.36→0.59 as
  transcription collapses and eviction relaxes → durable silencing.
- **Replicative-dilution phenotype**: the S-phase halving holds Mk 17 % (GNP) / 25 % (MB) below its
  no-dilution ceiling — read-write racing dilution.

---

## 3. Key findings

- **The leaky biology fixed the fold.** A sharp lockout Hill capped MB/GNP at ~4; the graded leaky
  repressor reaches target. The biologically-correct model is also the better-fitting one.
- **The reciprocal (nascent-RNA eviction) arm makes the bivalent state emerge.** Without it, M_ss is just
  set by read-write/turnover; with it, MB transcribes hard → keeps PRC2 partly evicted → poised/bivalent
  despite high EZH2, and vismo (transcription off) flips it to consolidated silencing. It also made the
  5.07× fold REACHABLE (previously capped ~4×).
- **EZH2-scaling all PRC2 catalysis** (not just read-write) was required: slow turnover otherwise lets the
  mark accumulate in non-dividing arrested cells and over-silences them. Tying methylation to the
  (E2F-coupled, cycle-falling) EZH2 level reconciles slow turnover with the in-band vismo residual.
- **Dilution limits EZH2's repressive reach** — each S-phase halving raises R; a faster cycle dilutes
  more (the endogenous-period loop, basis of the T_cc phenotype).
- **The EZH2i rescue is MB/vismo-SPECIFIC**: MB sustains high EZH2 → EZH2-maintained arrest → EZH2i
  rescues; mitogen-withdrawn GNP loses EZH2. See [[v44-h3k27-memory-derepression]].
- **HHi=1.0 needs no fix** under the leaky model — the residual floor keeps CyclinD1 off the stiff zero.

---

## 4. Outstanding (nothing here is a dilution-module regression)

**Pre-existing validation fails (the 5 keeping BOTH default and dilution at 22/27):**
1. **EZH2 transcript S/G0** (default 5.7, dilution 6.4; target 1.8–2.5) — within-cycle EZH2 gradient
   too steep. The longest-standing open item.
2. **MB HU S-fold (4.95) & G2-fold (0.84)** — HU fork-coupling needs re-tuning for the short S-phase.
3. **MYCN GNP+HHi (0.91 vs 0.78)** — cascade target.
4. **MB G2+M duration (3.6 vs 2.5 h)** — the accepted G2 sacrifice from the bake.
- The **expanded wide-search** (add S/G0, HU folds, MYCN-GNP+HHi as scored targets) would make a high
  score a strict 27-target improvement.

**Decisions / next steps:**
- **DONE — promoted**: `with_h3k27_dilution` is now the default repression; all registry figures were
  regenerated on it (feedback figures remapped to `f0_mk`; phase-plane reworked to the (Mk, CyclinD1) plane).
- **Reconcile the two H3K27 models** (memory vs dilution): the dilution model subsumes the memory model's
  de-repression role, but `fig_v44_ezh2i_rescue_kinetics` still uses the memory variant (via
  `with_h3k27_memory=True`, which now takes priority over the default dilution). Re-pointing it to the AUM
  model is still open.
- **Spec §8 payoff (PARTLY done):** the dilution-on/off causal test confirms the halving lowers Mk by
  17–25 % (above). Still to run as a figure: the full **T_cc sweep** (a CLEAN period knob that does not
  also move EZH2 — mitogen titration muddles it because EZH2 rises with dose) and **EZH2i × CDK4/6i
  synergy** (CDK4/6i lengthens T_cc → fewer dilutions → mark restored → differentiation).

**Figures:** the current suite regenerates on the DEFAULT model now via `simulations/regenerate_figures.sh`
(unaffected by the gated dilution work).

---

## 5. Files
- Module: `src/build_model_v44_heldt.py` (`with_h3k27_dilution` = the AUM module, `with_h3k27_memory` flags).
- Evidence base: `docs/H3K27me3_CyclinD1_literature_review.md` (the lit review that set the AUM design +
  kinetic anchors).
- Calibration: `simulations/calibrate_h3k27_aum.py` (the AUM δ/k_w sweep + t½, current),
  `calibrate_h3k27_dilution.py` (older candidate eval), `calibrate_h3k27_dilution_search.py` (6-param).
- Probes: `probe_h3k27_dilution.py` (sawtooth), `probe_derepression_kinetics.py` (memory kinetics).
- Figure: `fig_v44_ezh2i_rescue_kinetics.py` (uses the memory model).
- Validation of AUM model: `H3K27_DILUTION=1 ./venv/bin/python simulations/validate_v44.py`.
