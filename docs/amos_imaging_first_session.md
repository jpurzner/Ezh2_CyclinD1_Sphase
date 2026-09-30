# AMOS imaging: first session

Read `amos/dsa/MODEL_AGENT_BRIEF.md`, then ran `amos.imaging_model_targets()` and followed up.
2026-09-30.

**Headline: the data blocker is gone.** The HCI/IF data I had listed as blocking is available
programmatically, with per-cell access, replicate structure, gate definitions and explicit trust
flags. That is better than the flat export I asked for.

---

## 1. `imaging_model_targets()` — as advertised

12 targets: **2 reproduced, 5 close, 3 mismatch, 2 not addressable** — exactly the brief's stated
counts.

**Reproduced to 0.0%:** `ezh2cko_g0` (46.5) and `ezh2cko_g0_wt` (53.3), AAI, n ≈ 20k each.

**Gate verified sound.** The brief's free check: pRb median in `G0` = 317.9 vs `diploid` = 315.4,
**0.79% apart** — inside the ~1% criterion. G1 = 1094 and G2 = 1439 sit 3.4–4.5× above that
background, so pRb genuinely separates the states it gates. Values match the brief's 318/315/1094/1439
exactly.

## 2. The G0 discrepancy is a definition difference, not a measurement disagreement

The note on `mb_g0_frac` suspected the manuscript folds the `aneuploid` class into G0. **It does.**

| quantity | substrate | target | rel diff |
|---|---|---|---|
| **G1 alone** | 43.48 | 43.4 | **+0.2%** |
| G0 alone | 20.29 | 24.8 | −18.2% |
| **G0 + aneuploid** | **25.75** | 24.8 | **+3.8%** |
| 2N = G0+G1 | 63.77 | 68.2 | −6.5% |
| **2N + aneuploid** | **69.23** | 68.2 | **+1.5%** |

`aneuploid` = 5.46%. G1 matching to +0.2% shows the gate is placed correctly; the entire gap is the
one class. **Two of the five "close" targets become effectively reproduced once the convention is
fixed.**

### But this creates a real problem for the ODE model

`aneuploid` is DNA-high / pRb-negative — doublets or genuine aneuploidy. **v44 has no compartment
for it.** So comparing the model's 2N fraction against a target that silently includes 5.5% of
cells the model cannot represent is comparing different quantities.

v44 currently produces 68.7% for `MB 2N (G0+G1) count%` against a target of 68.2. If the target
includes aneuploid (69.2), the model is fine. If it should exclude them (63.8), the model is 7.7%
high. **Which convention the target uses needs deciding, and it changes the score.**

## 3. Two of v44's 28 validation targets have no traceable source

This is the finding with the most consequence, and it revises how the "28/28" headline should be
read.

### `EZH2 protein G2/G0 = 1.48` — unsourced, and only passing on a wide tolerance

AMOS: *"no extracted experiment measures Ezh2 intensity under a gate that resolves G2. Ezh2 is
measured in `gnp_ezh2_shh_dose` (a CONFLICTED experiment) and `ezh2cko_gnp_h3k27me3`. PI must name
the canonical experiment."* The brief is explicit: **do not fit to it from imaging.**

Separately, in `validate_v44.py:325` this target carries a tolerance of **0.45**:

```python
check("EZH2 protein G2/G0 (1.48)", ez0['G2']/ez0['G0'], 1.48, 0.45)
```

The model produces **1.009** — a **32% miss** — and is scored as a pass only because the tolerance
is 45% wide. So this target is both **untraceable** and **not actually being met**.

### `MB HU S fold = 1.36` — no HU arm exists in the imaging data

AMOS: *"no hydroxyurea arm exists in any of the 31 recipes. The HU experiment is not in
`imaging_tables` — PI to locate it."* The model passes at 1.365 (tolerance 0.40), but against a
number with no extracted source.

Note the brief says this one *does* now have a source but is deliberately unscored, because the
substrate will not claim an S-phase fraction from PCNA. Either way it is not currently traceable.

### Tolerances worth knowing

Default is 0.20, but a cluster of targets run much wider:

| target | tolerance |
|---|---|
| EZH2 protein G2/G0 | **0.45** |
| EZH2 transcript S/G0 | 0.45 |
| EZH2 Palbo mRNA drop | 0.42 |
| MB S count% | 0.42 |
| EZH2 G0/cycling | 0.40 |
| MB HU S / G2 fold | 0.40 |
| HU EZH2-in-S boost | 0.35 |
| MB 2N count% | 0.30 |
| MB G2+M duration | 0.55 |

Wide tolerances are not automatically wrong for a model calibrated to ratios rather than absolute
concentrations. But **"28/28 passed" should be stated with the tolerances and the two unsourced
targets attached**, not as a bare score. I have been repeating the bare number and should stop.

## 4. The three "mismatch" targets are a dose question, not a model failure

| target | manuscript | substrate | |
|---|---|---|---|
| `mb_cdk46i_arrest` | 17.0 | **13.5** @ 5 µM | |
| `mb_cdk46i_arrest_1uM` | 17.0 | **39.2** @ 1 µM | |
| `mb_vismo_arrest` | 25.0 | 34.6 | |

The target of 17.0 sits **between** the two palbociclib arms. Neither dose reproduces it, and the
substrate reports both rather than picking one. Same pattern on the rescue: UNC1999 at 1 µM gives
69.9% (target 75, close) while **5 µM gives no rescue at all**.

So these are not disagreements about measurement — they are the manuscript target not specifying a
dose, with the two available doses bracketing it. Both AMOS notes say the same thing: *PI to
confirm which dose is canonical.*

## 5. What this changes for the noise work

**Now unblocked:** per-cell distributions, which is what every predicted CV needed.
`imaging_distribution` and `imaging_fit_distribution` give the spread *and* fit its shape by MLE
with a KS statistic, plus Gaussian mixtures ranked by BIC. The brief notes pRb is gamma with
CV 0.35 in G1 but skew-normal with CV 0.13 in G0 — negative states sit at background, positive ones
are biology. That distinction matters for comparing predicted CVs to measured ones: **a CV computed
across a state that is really two populations is meaningless**, and the BIC mixture test says which.

**Still needed, and unchanged:** an absolute scale (Ω) and measured T_cc. Intensities are per-plate
arbitrary units by design, so imaging cannot supply Ω.

**One caution carried forward:** 10 experiments have unresolved channel conflicts, including a
suspected pRb↔NeuN swap in the MB-washout family. ACD — the workhorse here — is *not* among them.

---

## Pending decisions for the PI

Not urgent, but they gate specific numbers. Recorded so they are not rediscovered:

1. **Palbociclib canonical dose** — 1 µM or 5 µM? The target 17.0 lies between them.
2. **UNC1999 canonical dose** — 1 µM rescues, 5 µM does not.
3. **`EZH2 G2/G0 = 1.48`** — name the source experiment, or retire the target. It is currently
   unsourced *and* the model misses it by 32%.
4. **`MB HU S fold = 1.36`** — locate the HU imaging experiment, or score it from the dose-response
   rather than the single fold.
5. **Does the manuscript's 2N / G0 include the aneuploid class?** It appears to. If so the model
   target should say so, and the model has no compartment for those cells.
