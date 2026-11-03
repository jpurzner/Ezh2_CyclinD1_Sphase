# Stochastic chromatin layer on the two-compartment amplifier (2026-10-22)

Tests whether the all-or-none behaviour reported in
`two_compartment_prc2_2026-08-09.md` survives making the P_prox transition
stochastic. Script: `simulations/stochastic_amplifier_population.py`.
Figures: `fig_v44_stochastic_amplifier`, `fig_v44_stoch_fine_N12`,
`fig_v44_stoch_fine_N4` (+ .json).

---

## 1. BUG: `population_amplifier.py`'s "amplifier OFF" arm never turns it off

```python
CONDS = [('amplifier OFF', None, ...), ('transient (C)', dict(a_P=0.10, f_P=0.80), ...), ...]
def make(amp):
    P = dict(BP)
    if amp: P.update(amp)          # amp=None -> NO override
    rr = te.loada(build_model_v44(..., with_proximal_distal=True, params=P))
```

`amp=None` leaves the Antimony defaults `a_P=0.10, f_P=0.80` in place, and the
flag is hard-coded True. **"amplifier OFF" is bit-identical to "transient (C)".**
That is why the script reported them as giving identical dwells, and it is the
sole basis for this conclusion in the design doc:

> "Amplifier OFF and transient (C) give IDENTICAL dwells at every excursion —
> the validation-safe amplifier does **not** prolong the transient-G0 dwell."

That conclusion does not hold. With the amplifier **genuinely** off (`a_P=0`),
same protocol, 5x CDKI excursion:

| condition | dwell |
|---|---|
| genuinely OFF (a_P=0) | **59.8 h** |
| transient C (a_P=0.10, f_P=0.80) | **85 h** |

So the validation-safe amplifier prolongs the transient-G0 dwell by **~42%**,
not 0%. **Fix:** `('amplifier OFF', dict(a_P=0.0), ...)`, or build that arm with
`with_proximal_distal=False`.

(Instrumentation check: with noise disabled, the build used here reproduces the
existing model to 3 decimals — transient C 84/85/85 h at 4x/5x/6x, permanent A
LOCK, P_prox identical. The discrepancy is the CONDS entry, not the model.)

## 2. f_P, not a_P, is the controlling axis

Deterministic lock map, 5x excursion:

```
         a_P=0.10   0.15    0.20    0.30    0.50
f_P=0.80     85h     95h    104h    132h    LOCK
f_P=0.70    LOCK    LOCK    LOCK    LOCK    LOCK
f_P<=0.60   LOCK    LOCK    LOCK    LOCK    LOCK
```

Sweeping a_P alone at f_P=0.80 never reaches the cliff (0% locked out to
a_P=0.14). The documented "sharp boundary at a_P~0.11" only holds along the
scenario path where f_P co-varies. The real boundary at a_P=0.10 is a **step
between f_P 0.760 and 0.765** — a 0.7% parameter change flips permanent lock to
a finite 148 h dwell. The BAKED default (f_P=0.80) sits ~5% from permanent lock.

Note the deterministic dwell at f_P=0.765 is **148 h** — the "~148 h" figure in
the scenario table sits exactly on the knife-edge, one grid step from lock.

## 3. The stochastic layer

CLE on P_prox for a promoter with `Nsites` proximal nucleosome sites:

```
dP_prox/dt = lam_on - lam_off + sqrt((lam_on + lam_off)/Nsites) * eta(t)
```

Nsites -> inf recovers the current model exactly. Amplitude **calibrated, not
tuned**: scaled to the exact binomial stationary variance P(1-P)/Nsites at one
operating point, then verified at a second independent one (different rate,
different N, different P*) to 3.6%.

`eta(t)` is a band-limited Gaussian process (random-phase sum of 48 cosines,
periods 40-1500 min), not white noise and not a jump process. **This is forced,
not preferred:** v44 cannot be restarted during deep arrest — chunked
integration dies at t~6250 min with CV_CONV_FAILURE, which reproduces with
`with_proximal_distal=False` and is unfixable by any tolerance or integrator
setting. That rules out chunk+SSA and timed SBML events, both of which force a
CVODE restart. `population_amplifier.py` survives only because it restarts once,
at t=22000, in a safe state.

## 4. Result: the locked fraction IS graded — but by ESCAPE, not commitment

% of cells permanently locked, a_P=0.10, 5x excursion, 24 cells/point:

| f_P | deterministic | Nsites=12 | Nsites=4 |
|---|---|---|---|
| 0.745 | LOCK | 100% | **20.8%** |
| 0.750 | LOCK | 95.8% | 4.2% |
| 0.755 | LOCK | 75.0% | 0% |
| 0.760 | LOCK | 4.2% | 0% |
| 0.765 | 148 h | 0% | 0% |
| 0.770 | 109 h | 0% | 0% |

- A **tunable intermediate locked fraction exists** (96/75/4%) — the Bintu
  signature: cells are individually all-or-none, the parameter sets the fraction.
- But it lies **entirely inside the deterministic lock region**. It is not a new
  middle regime on the safe side; it is the lock boundary smeared into the lock.
- **More noise does not widen the window — it shifts it.** N=4 (3x the variance
  of N=12) moves the whole curve to lower f_P and *reduces* locking everywhere.
  So the mechanism is **noise-induced escape from a metastable lock** (Kramers),
  not stochastic commitment into it. The locked state needs P_prox held high
  continuously; any downward excursion that lets Pol II re-engage is
  self-amplifying, so the escape route is much wider than the entry route.
- On the safe side (f_P >= 0.765) stochastic and deterministic are
  indistinguishable in locked fraction (both 0%); noise shortens dwells modestly
  (N=4: 90-110 h vs deterministic 109-148 h).

## 5. What this does and does not change

**Survives:** "no robust moderately-prolonged transient dwell in between" on the
validation-safe side. Noise does not manufacture a tunable dwell at f_P >= 0.765.

**Changes:** (a) the safe amplifier is *not* inert — it prolongs the dwell ~42%
(bug above); (b) a tunable locked *fraction* does exist, in the region previously
written off as permanent-lock. A parameter set that locks 75% of cells is a
different validation proposition from one that locks 100%, so the blanket
exclusion of the deep latch deserves re-checking against the stochastic fraction
rather than the deterministic verdict.

**Direction of the correction:** stronger chromatin noise makes locking *harder*,
so a real promoter with few proximal sites would reach permanent lock less
readily than the deterministic model says. That strengthens the exclusion of the
deep latch rather than weakening it.

## 6. Caveats

- Band-limited noise: the window's absolute width depends on the correlation
  time chosen (40-1500 min). The amplitude is calibrated; the bandwidth is a
  modelling choice.
- 24 cells/point => roughly +/-10% on a 75% fraction. Sigmoid shape is solid,
  midpoint good to maybe +/-0.003 in f_P.
- One excursion magnitude (5x) and one a_P (0.10).
- `Nsites` is unmeasured. Since position and width both depend on it, treat the
  existence of the graded window as the result and the numbers as illustrative.
- The discriminator experiment from the original doc is unchanged and would also
  constrain Nsites: sorted G0 vs cycling MB cells, CUT&RUN H3K27me3 + PRC2 at the
  CDK6/CyclinD1 proximal promoter. Single-cell rather than bulk would give the
  locked *fraction* directly.
