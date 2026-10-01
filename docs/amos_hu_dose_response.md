# The HU dose series vs v44: the target is met, the mechanism is not

Follow-up after the imaging agent located the HU data. I pulled the numbers directly from AMOS
rather than working from the relayed summary. `AC_MB55_HU_doseEsc_Feb2026`, MB55, 14 doses over
5 log units, 251,277 cells, 4 replicates per dose.

**Two findings. One is a measurement artefact in the data; the other is a structural miss in the
model, and it is the more important of the two.**

---

## 1. The imaging agent's mechanism is correct, and the data confirms it in absolute units

Their account — HU stalls forks at the start of S, so cells are biologically in S but their DNA
content has barely advanced past 2N, and a DNA-content gate scores them 2N — is right, and this
gate makes it quantitative: its 2N/S boundary is **"convention 2N × 1.25"**, so **a cell with less
than 25% of its genome replicated is invisible as S.**

Worth adding: in *relative* terms the 2N row looks flat (97–105% of control), which is how I first
misread it. In **absolute** percentage points it is not flat at all, and conservation balances:

| dose | 2N | S | G2 | aneu | ΔS+ΔG2 | Δ2N |
|---|---|---|---|---|---|---|
| 0 | 85.95 | 8.34 | 4.02 | 1.70 | — | — |
| 100 µM | 88.97 | 6.36 | 2.79 | 1.89 | −3.21 | **+3.03** |
| 10 mM | 90.36 | 4.57 | 3.13 | 1.94 | −4.65 | **+4.41** |

2N accumulates monotonically and absorbs the S and G2 losses to within 0.2 pp. **v44's S-entry
block is qualitatively validated** — HU does hold cells at 2N, exactly as `fire_gate_HU` intends.

## 2. The 10 µM point is confounded by object count, not DNA content

The 10 µM condition breaks the monotone trend in *every* state at once — 2N down 2.26 pp while S,
G2 **and** aneuploid all rise. I first hypothesised a shifted DNA-content calibration on that plate.
**The data refutes that:** the S-state DNA distribution is essentially identical across conditions
(mean 3.65–3.70 × 10⁷, p25/p75 unchanged). Not a staining or exposure shift.

The actual confound is **object count**. Totals vary 0.78×–1.69× across the series, and the classes
do not scale together — at 10 µM: 2N ×1.64, S ×1.93, G2 ×1.76, aneuploid ×2.53.

In the low-dose plateau (≤10 µM), where HU should do little, the phase proportions are predicted by
object count rather than dose:

| | Pearson r | p |
|---|---|---|
| **S% vs total objects** | **+0.823** | 0.0019 |
| aneuploid% vs total objects | +0.839 | 0.0013 |
| 2N% vs total objects | −0.805 | 0.0028 |

The three highest-count conditions (0.1, 0.3, 10 µM) average S = 9.00% against 7.89% for the rest —
a 14% inflation from counting alone, which is about the entire size of the "S rises at 10 µM"
signal. HU does not increase cell number at 24 h, so a 1.69× object count points at debris or
segmentation, consistent with the brief's own caveat that these plates are **under-filtered for
nuclear size**.

There is also a genuine dose trend (S% vs log dose, r = −0.814, p = 0.0004). Both effects are
present; the trend is real, the point-to-point scatter largely is not.

## 3. The model reproduces the target — at a dose where the data does not support it

With the validation's exact MB preset, `hu = dose_µM/10`:

| dose | model S fold | measured | |
|---|---|---|---|
| 0 | 1.00 | 1.00 | |
| 1 µM | 1.05 | 0.91 | |
| 5 µM | 1.44 | 0.85 | |
| **10 µM** | **1.37** | **1.15** | ← the validation target is 1.36 |
| 100 µM | **0.00** | 0.76 | |
| 1 mM | **0.00** | 0.77 | |
| 10 mM | **0.00** | 0.55 | |

The model hits 1.37 at 10 µM, matching `validate_v44`'s reported 1.365 — so the target is
self-consistent with the model's own HU scale. But **no dose in the data gives 1.36.** The closest
is 1.15, at the single most count-confounded condition; its neighbours give 0.85 and 0.76.

So `mb_hu_s_fold = 1.36` is not supported by this readout at any dose. **AMOS's
`not_addressable` call is correct**, now for two independent reasons: the gate cannot see
early-stalled S cells (their point), and the one dose that looks supportive is a counting artefact
(this analysis).

## 4. The structural miss: v44 over-arrests by two orders of magnitude

This is the finding that matters for the model.

**Above ~50 µM the model has no S-phase cells at all** — 100% G1, S = 0.00. The data has S at
**76% of control at 100 µM and still 55% at 10 mM**, a thousand-fold higher dose.

The cause is explicit in the builder:

```
KmHU_fire = 0.25;   hHU_fire = 2;
fire_gate_HU := KmHU_fire^hHU_fire/(KmHU_fire^hHU_fire + HU^hHU_fire);
```

| dose | HU (model) | `fire_gate_HU` | `vfork` |
|---|---|---|---|
| 10 µM | 1.0 | 0.059 | 0.107 |
| 100 µM | 10 | **0.0006** | 0.100 |
| 10 mM | 1000 | **~0** | 0.100 |

Origin firing is suppressed to 0.06% by 100 µM, so no cell can enter S and the population drains
into G1. `vmin_fork = 0.1` was intended to keep S finite — but it cannot, because the *entry* gate
has already closed.

**This is not explainable by the gate limitation.** Cells scored S have DNA content between 1.25×
and ~1.7× 2N, i.e. they have measurably replicated 25–70% of their genome. 4.57% of cells are doing
that at 10 mM. The model says zero cells are replicating at all.

### The fix, and the data to do it with

`KmHU_fire` and `hHU_fire` were "calibrated to MB55 HU phase folds" — i.e. to the single fold that
§3 shows is unsupported. A 14-point, 5-log dose series now exists. Replacing one fold with a
dose-response fit would constrain both the IC50 and the Hill coefficient, and would have to be done
against **absolute percentage points** with the object-count confound handled, not against
normalised folds.

Expect the fit to want a **higher `KmHU_fire`** and a **shallower `hHU_fire`**, plus a floor on
origin firing so S never goes identically to zero — the biological reality being that even saturating
HU leaves some replication running.

---

## Recommendation

1. **Retire `mb_hu_s_fold = 1.36` as a scored target.** It is met by the model and unsupported by
   the data. Keeping it rewards the wrong behaviour.
2. **Replace it with a dose-response target** over the 5-log series, fitted on absolute percentage
   points, with the three high-count conditions down-weighted or the counts regressed out.
3. **Recalibrate `KmHU_fire` / `hHU_fire`** against that series. The current values make the model
   arrest completely at doses where the data still has a quarter to three-quarters of its S-phase
   population.
4. **Flag the object-count confound to the imaging side.** Totals ranging 0.78×–1.69× with
   class-dependent scaling will distort any proportion computed from these plates, not just this
   target. The brief's nuclear-size under-filtering caveat is the likely cause.

## Caveat on comparability

v44's `classify()` defines S as a **BrdU/EdU analogue** (`vfork*aRc > 0.02`), deliberately excluding
HU-stalled cells. The AMOS gate is DNA-content based. These are different observables and should not
be compared naively at intermediate doses. The §4 finding survives that objection — zero is zero on
any readout — but the §3 dose-by-dose comparison should be read with it in mind. If a BrdU or EdU
arm exists anywhere in the imaging set, it would be the correct comparator and is worth looking for.
