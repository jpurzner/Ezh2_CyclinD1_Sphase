# H3K27me3 repression of CyclinD1 — model state & outstanding items
*Checkpoint after the replicative-dilution work. Captures the model variants, the calibrated leaky
dilution module, validation status, key findings, and what remains — so nothing is dropped.*

---

## 1. Three repression mechanisms now exist (all gated; default is unchanged)

| Build flag | Repression of CyclinD1 | Role | Validation |
|---|---|---|---|
| **default** (`build_model_v44()`) | `K_rep/(K_rep + EZH2·(1−EZH2i))` — EZH2 level represses directly, instantaneous | the production model | **22/27** |
| `with_h3k27_memory=True` | mark `H3K27_Cd` tracks EZH2, `K_rep/(K_rep+H3K27_Cd)`; EZH2i → mark decays (k_demeth+dilution) over ~24–48 h | de-repression KINETICS only (mark tracks EZH2 at steady state → validation-preserving); used by `fig_v44_ezh2i_rescue_kinetics` | 22/27 (steady states identical) |
| `with_h3k27_dilution=True` | **explicit replicative-dilution module** (Purzner spec) — see §2 | the mechanistic, biology-faithful model | **22/27** (parity; same 5 fails as default) |

Default flag is OFF on all → the committed/validated production model is untouched. `validate_v44.py`
runs the dilution model via `H3K27_DILUTION=1` env var.

---

## 2. The replicative-dilution module (the current best mechanistic model)

State `Mk ∈ [0,1]` = H3K27me3 occupancy over the ~7 kb *Ccnd1* domain.
- **Continuous**: autocatalytic read-write `k_w·EZH2·Mk` (EED reads me3) + de-novo nucleation floor
  `k0` (lets a halved locus reseed), both on `(1−Mk)` and **both blocked by EZH2i** (all EZH2
  repression flows through the mark); minus demethylation/turnover `δ·Mk`.
- **Discrete**: replicative **halving at early S** (`Dna>0.05` event, once/cycle) → the emergent cycle
  period sets the dilution frequency (the loop the spec closes; faster cycle ⇒ more dilution).
- **Readout (LEAKY, per JP)**: `R(Mk) = f0 + (1−f0)/(1+(Mk/K)^n)` multiplying the full Gli+MYCN drive.
  R=1 at Mk=0; saturates at residual **f0** at high mark — H3K27me3 DAMPENS transcription (impedes
  elongation, Pol II stays), it does NOT lock out. Mark inherited through division = the memory.

**Calibrated params (search, 6/6 guards + 6/6 soft):** `k_w=0.00261, k0=0.000258, δ=0.00117,
K=0.305, n=4.15, f0=0.233`.

**Behavior:** GNP & MB cycle (Cd ~1.9–2.1 / 8.3, periods ~22 h), GNP+HHi & MB+HHi arrest (with
RESIDUAL transcript — MB+HHi Cd ~1.1, arrests because vismo removed the Gli DRIVE, not via lockout),
**MB+HHi+EZH2i RESCUES** (4 div), GNP+EZH2i de-represses. MB/GNP fold ~4.0–4.4 (within validate's 35%
tol of 5.07; search hit 5.06), EZH2i fold 2.17. Operating repression is partial everywhere (GNP ~46%
transcription, MB ~34%, MB+HHi ~60%) = "lots of transcript even when present".

---

## 3. Key findings

- **The leaky biology fixed the fold.** A sharp lockout Hill capped MB/GNP at ~4 and hit only 5/6
  soft; the graded leaky repressor reaches the target (6/6 soft). The biologically-correct model is
  also the better-fitting one.
- **Dilution limits EZH2's repressive reach** — each S-phase halving raises R; a faster cycle dilutes
  more. This is the spec's endogenous-period loop and is the mechanistic basis of the T_cc phenotype.
- **The EZH2i rescue is MB/vismo-SPECIFIC** (from the kinetics work): MB sustains high EZH2 → stable
  EZH2-maintained arrest → EZH2i rescues; mitogen-withdrawn GNP loses EZH2 so its arrest is not
  EZH2-maintained. See [[v44-h3k27-memory-derepression]].
- **GDC=1.0 needs no fix** under the leaky model — the residual floor keeps CyclinD1 off the stiff zero.

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
- **Canonical vs research module**: is `with_h3k27_dilution` promoted to the default repression (then
  regenerate all ~24 figures on it), or does it stay an additional model? It is at validation parity,
  so promotion is viable. NOT yet decided.
- **Reconcile the two H3K27 models** (memory vs dilution) if dilution becomes canonical — the dilution
  model subsumes the memory model's de-repression role; re-point `fig_v44_ezh2i_rescue_kinetics`.
- **Spec §8 payoff (NOT yet run):** the **T_cc sweep** (critical cycle length where restoration can't
  keep pace and Mk ratchets below K — *the* replicative-dilution phenotype curve) and **EZH2i ×
  CDK4/6i synergy** (CDK4/6i lengthens T_cc → fewer dilutions → mark restored → differentiation).

**Figures:** the current suite regenerates on the DEFAULT model now via `simulations/regenerate_figures.sh`
(unaffected by the gated dilution work).

---

## 5. Files
- Module: `src/build_model_v44_heldt.py` (`with_h3k27_dilution`, `with_h3k27_memory` flags).
- Calibration: `simulations/calibrate_h3k27_dilution.py` (candidate eval),
  `calibrate_h3k27_dilution_search.py` (6-param search), `calibrate_h3k27_kinetics.py` (memory model).
- Probes: `probe_h3k27_dilution.py` (sawtooth), `probe_derepression_kinetics.py` (memory kinetics).
- Figure: `fig_v44_ezh2i_rescue_kinetics.py` (uses the memory model).
- Validation of dilution model: `H3K27_DILUTION=1 ./venv/bin/python simulations/validate_v44.py`.
