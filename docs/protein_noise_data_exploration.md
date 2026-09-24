# What the downloaded data actually say, and what is still missing

Third document in the set. `protein_noise_plan.md` is the reconnaissance, `protein_noise_fixes.md`
is the tested change list, this one is the exploration of the datasets themselves.

Everything below is reproducible from `data/priors/raw/` (provenance in `data/priors/SOURCES.md`).

> **Two corrections to `protein_noise_plan.md` are folded in below and flagged where they occur.**
> The headline one: **EZH2 is slow but it is *not* quiet.** §1.3 of that document said CV ≈ 0.26;
> the correct figure is **≈ 0.89**, and the error was mine — see §3.

---

## 1. Schwanhäusser's half-lives ARE growth-corrected. T4's open question is settled.

The dilution audit could not classify Schwanhäusser's values without reading the paper's Methods,
so the CSV recorded `growth_corrected = unknown`. It can be settled from the data instead.

**The test.** Schwanhäusser measured cycling NIH3T3. Mathieson measured primary cells that are
post-mitotic, so their half-lives contain no dilution term. If Schwanhäusser's numbers were *not*
growth-corrected they would carry an extra removal rate µ = ln2/T_cc ≈ 0.029–0.035 /h relative to
Mathieson. Dilution is only detectable where degradation is slow, so stratify by stability.

Mouse-vs-mouse (NIH3T3 vs mouse neurons), 3613 paired genes:

| neuron t½ bin | n | median t½ Schwanh. | median t½ neuron | median (1/t_S − 1/t_M) | as a fraction of µ |
|---|---|---|---|---|---|
| 10–25 h | 318 | 15.9 | 20.7 | +0.0116 /h | 0.33× |
| 25–50 h | 1213 | 33.1 | 38.6 | +0.0028 /h | 0.08× |
| 50–100 h | 1667 | 70.8 | 68.0 | −0.0003 /h | −0.01× |
| 100–300 h | 400 | 122.0 | 116.2 | +0.0002 /h | 0.01× |

The gap **shrinks to zero exactly where dilution would be largest** — the opposite of what a
baked-in µ would do. Sharpest form, on the 413 most stable proteins (neuron t½ > 100 h):

| hypothesis | predicted median Schwanhäusser t½ |
|---|---|
| Schwanhäusser = neuron degradation **+ dilution** | 23.1 h |
| Schwanhäusser = **pure degradation** | 116.9 h |
| **observed** | **121.3 h** |

**Conclusion: Schwanhäusser protein half-lives are growth-corrected — class (b).** So `-µ·x` is a
genuinely *missing* term for any parameter derived from them, not a double-count. The residual
+0.012 /h in the 10–25 h bin is the wrong shape for dilution and is more likely real biology
(fast-turnover proteins differ between dividing fibroblasts and post-mitotic neurons).

Caveat worth keeping: Mathieson's own cell types span 4× in median half-life (NK cells 170 h vs
neurons 47 h), so cross-cell-type biology is a large confounder. The stratified result is robust
to it; the unstratified one is not.

---

## 2. Our species are drawn from the unstable tail of the proteome

Schwanhäusser genome-wide (mouse NIH3T3), and where v44's species sit:

| | p5 | p25 | **p50** | p75 | p95 |
|---|---|---|---|---|---|
| protein copies/cell | 1,156 | 8,141 | **42,775** | 240,097 | 2,700,630 |
| mRNA copies/cell | 3.9 | 9.6 | **17.7** | 36.0 | 129.7 |
| protein t½ (h) | 7.6 | 23.6 | **47.3** | 88.0 | 215.3 |
| mRNA t½ (h) | 4.5 | 7.3 | **10.0** | 13.3 | 20.2 |

| gene | protein t½ | **%ile** | mRNA t½ | %ile | copies | %ile |
|---|---|---|---|---|---|---|
| **Ccnd1** | **1.3 h** | **0** | 14.0 h | 78 | 89,631 | 61 |
| Ccna2 | 4.3 h | 1 | 8.5 h | 36 | 21,381 | 38 |
| Skp2 | 3.4 h | 1 | **2.8 h** | **0** | 3,402 | 14 |
| Cdk4 | 7.0 h | 4 | 10.5 h | 54 | 298,984 | 78 |
| Chek1 | 7.7 h | 5 | 5.9 h | 14 | 9,371 | 27 |
| Ezh2 | 22.0 h | 23 | 9.8 h | 49 | 2,156 | 10 |
| Cdkn2c (p18) | **189.5 h** | **94** | 4.8 h | 7 | 4,356 | 16 |

**`Ccnd1` protein is at the 0th percentile — effectively the least stable protein in the dataset.**
Every cell-cycle regulator here sits at the 0th–5th percentile for protein stability. The brief's
§2.1 "typical mammalian values" table (median t½ 46 h) is therefore the **wrong reference class
for almost every node in this model** and should not be used to set defaults.

The two exceptions run the other way: `Cdkn2c` (p18) at 189.5 h is removal-by-dilution-only, which
is consistent with v44 treating it as a constant parameter; and `Ezh2` at 22 h is the most stable
*dynamic* species in the model.

---

## 3. Correction: EZH2 is slow, but it is **not** quiet

### The error

`analyze_noise_budget.py` originally took `⟨m⟩` from Schwanhäusser when available and `B` from
Larsson. **Those are not the same unit system**, and the mix does not cancel.

Larsson is allele-resolved scRNA-seq: its counts are *detected* molecules. Calibrating against
Schwanhäusser over 3080 shared genes gives a median ratio of **4.8×**, i.e. a UMI capture
efficiency **p ≈ 0.21**. Binomial thinning of a Poisson-Beta leaves `kon` and `koff` untouched and
scales `ksyn` by p, so

```
B_obs = p·B_true        ⟨m⟩_obs = p·⟨m⟩_true        burst FREQUENCY is unaffected
⇒  (1 + B_true)/⟨m⟩_true  ==  (p + B_obs)/⟨m⟩_obs
```

The ratio is **capture-invariant — but only within one dataset.** Pairing an absolute `⟨m⟩` with a
deflated `B` understates the mRNA noise term by ≈1/p. Measured on our gene set: **4.4× to 15.6×**.

Verified internally consistent: Larsson's own `f·B` reproduces its reported mean expression
(median ratio 1.05), so `f`, `B` and `⟨m⟩` are one system.

`analyze_noise_budget.py` is fixed — `⟨m⟩` now always comes from Larsson alongside `B`, with the
`(p + B)/⟨m⟩` form; Schwanhäusser's absolute copy number is kept only for the `1/⟨p⟩` Poisson term,
which needs absolute counts and is negligible anyway above ~10³ copies.

### What changes

| gene | CV as published in the plan doc | **corrected CV** | τ (unchanged) |
|---|---|---|---|
| **Ezh2** | 0.26 | **0.89** | 14.2 h |
| Ccna2 | 0.52 | **1.87** | 12.3 h |
| Skp2 | 0.39 | **1.00** | 4.1 h |
| Suz12 | 0.29 | **1.00** | 13.2 h |
| Ccnd1 | 0.91 | **0.87** | 20.2 h |

Only `Ccnd1` is essentially unchanged — because its `⟨m⟩` already came from Larsson, so it was
self-consistent by accident. Correlation times are unaffected: they depend only on half-lives.

### Does this overturn the recommendation?

**No, but it changes the reason.** The plan said "EZH2 doesn't need bursts because it is slow and
quiet". Drop the second half. EZH2 is a **large-amplitude, slow** fluctuator.

The picture that emerges is coherent and rather striking. With `⟨m⟩_true ≈ 30`, `B_true ≈ 43` and
`f ≈ 0.055/h` per cell, **EZH2 receives roughly one large transcriptional burst per cell cycle**,
and because both its mRNA (9.8 h) and its protein (22 h) are long-lived, that burst structure is
preserved as large, slowly-varying cell-to-cell heterogeneity.

What decides drift-vs-bursts is the **timescale**, not the amplitude. At τ ≈ 14 h ≈ one generation
a cell keeps its EZH2 rank for the whole cycle, which is exactly frozen heterogeneity. So the
recommendation stands: **option 2 (per-cell sampled parameters) or a slow OU for EZH2, not a PDMP.**
But it should be sampled with **CV ≈ 0.4–0.9, not 0.26** — a substantially wider spread of the
epigenetic brake across the population than the plan doc implied.

---

## 4. The biggest finding: almost nothing fluctuates fast

Computing the full budget genome-wide over the 3298 genes with complete data (T_cc = 16 h):

| | p5 | p25 | **p50** | p75 | p95 |
|---|---|---|---|---|---|
| predicted CV_p | 0.39 | 0.52 | **0.66** | 0.87 | 1.61 |
| correlation time τ (h) | 11.1 | 15.9 | **18.6** | 21.1 | 29.0 |

```
fraction with tau >  8 h  (half a cycle) : 99%
fraction with tau > 16 h  (a full cycle) : 74%
fraction with tau <  2 h  (averages out) :  0%
```

**The brief's central dilemma is largely answered by the priors.** §1 of the hand-off asks whether
a snapshot reflects frozen heterogeneity or fast fluctuation, and treats the two as equally
plausible. For essentially the whole measured proteome the answer is **frozen heterogeneity**:
protein lifetimes are simply long enough that nothing fluctuates fast. This is a strong prior in
favour of the brief's **option 2**, and it means option 3 (drift) should be parameterised with
long correlation times by default rather than short ones.

It also reframes the noise work's premise. If nothing moves on a sub-cycle timescale, spontaneous
basin-hopping driven by *fast* fluctuation is rare, and cell-to-cell outcome differences are set
mostly at birth. That is a testable, and quite different, model of G0 exit.

### The one exception in our set — and the brief missed it

| gene | CV | %ile | τ (h) | %ile | reading |
|---|---|---|---|---|---|
| Ccnd1 | 0.87 | 74 | 20.2 | 67 | slow + loud → frozen heterogeneity |
| Ezh2 | 0.89 | 77 | 14.2 | 16 | slow + loud → frozen heterogeneity |
| Ccna2 | 1.87 | 97 | 12.3 | 8 | slow + loud |
| **Skp2** | **1.00** | **83** | **4.1** | **0** | **fast + loud → genuine fluctuation** |
| Chek1 | 1.62 | 95 | 8.5 | 1 | slow + loud |
| Cdk4 | 0.48 | 17 | 15.1 | 21 | slow + quiet |

**`Skp2` is the only species in the model whose noise is genuinely fast** — τ = 4.1 h, the 0th
percentile genome-wide — because it is the rare gene with *both* a short mRNA (2.8 h, 0th
percentile) and a short protein (3.4 h, 1st percentile) half-life. And it is loud (CV ≈ 1.0).

That matters because **Skp2 sits directly on the Skp2–p27 bistable commitment switch** that
decides G0 exit in the Heldt core. The brief nominates `Ccnd1`, `Cdkn1a` and `Mycn` for explicit
bursts on the grounds of low copy number and short protein half-life; on the priors, **`Skp2` is
the strongest candidate for fast-fluctuation-driven basin escape and is not on the list.**

Recommend adding Skp2 to the T6 candidate set — arguably ahead of `Ccnd1`, whose τ ≈ 20 h makes it
a slow variable whatever its amplitude.

---

## 5. Burst kinetics: `B` transfers between cell types, `f` does not

Larsson gives mouse primary fibroblasts (7870 genes) and mouse ES cells (7309), 5970 shared.

| | median B | median f [per δ] | median ⟨m⟩ [per allele] |
|---|---|---|---|
| fibroblast | 3.05 | 0.413 | 1.05 |
| ESC | 3.50 | 0.370 | 1.02 |

Across the 5970 shared genes:

| | median ESC/fibroblast ratio | log–log correlation |
|---|---|---|
| burst **size** B | 1.11 | **r = +0.50** |
| burst **frequency** f | 1.01 | **r = +0.13** |

No systematic shift either way, but **per-gene burst frequency is essentially uncorrelated between
cell types** while burst size partially transfers. Practical consequence: `B` borrowed from
fibroblasts is defensible for GNPs; **`f` is not**, and `f` is what sets burst timing and therefore
the tails that drive basin escape. Since GNPs are proliferating progenitors, ESC is arguably the
closer reference than fibroblast — but with r = +0.13 the choice barely constrains `f` anyway.

Burst size percentiles for our genes (fibroblast):

| gene | B | %ile | note |
|---|---|---|---|
| **Cdkn1a** | **19.04** | **96** | confirms it as the top burst candidate |
| Ccna2 | 18.66 | 96 | |
| Cdk1 | 19.85 | 96 | |
| Ccne1 | 13.57 | 93 | |
| **Ccnd1** | **6.76** | **80** | bursty, not extreme |
| **Ezh2** | **2.77** | **45** | completely typical — another reason not to burst it |
| Cdkn1b | 0.99 | 1 | essentially non-bursty |
| Mycn | 18.50 (ESC only) | 94 | no fibroblast value |

---

## 5b. Sigal 2006 confirms the brief's extrinsic-noise guess — with one subtlety

Fetched after the first pass (`data/priors/raw/sigal2006_SuppNotes.pdf`, parsed to
`data/priors/sigal2006_mixing_times.csv`). Suppl. Tables 1 and 2 give per-protein **mixing time**
τm and CV for 20 CD-tagged H1299 clones, by two estimators.

| estimator | n | τm median [range] (generations) | CV median [range] |
|---|---|---|---|
| protein ranks (Table 1) | 20 | 1.60 [0.8 – 2.6] | 0.20 [0.12 – 0.28] |
| protein levels (Table 2) | 20 | 1.85 [0.9 – 4.3] | 0.20 [0.12 – 0.28] |
| **pooled** | **40** | **1.75 [0.8 – 4.3]** | **0.20 [0.12 – 0.28]** |

The brief's §3.1 specifies a shared global OU with "correlation time ≈ 1–2 × T_cc; CV ≈ 0.2".
**Both numbers are confirmed by measurement**, which is a genuinely reassuring result — that part
of the plan can be taken as parameterised rather than assumed.

### The subtlety: do not set τ_ext = 1.75 generations directly

The supplement derives τm = ln2/α, where **α is the combined degradation-plus-dilution rate**, and
states explicitly that "for proteins not appreciably degraded, α reflects only dilution by cell
growth, and the memory time is τm = 1 cell generation". Memory longer than one generation requires
noise transmitted from upstream, or feedback — "for stable proteins τm ~ 2".

So the measured 1.75 generations is the **total observed memory, which already contains the
dilution-driven part that v44 models explicitly via µ.** Setting the extrinsic OU's correlation
time to 1.75 generations *on top of* the model's own dilution double-counts the same physics and
would overstate the memory.

The correct calibration is to choose τ_ext so that the *combined* simulated autocorrelation decays
to half in ≈1.75 generations, not to set τ_ext to that value. The brief's phrasing invites the
error. Note also the supplement's remark that **τm does not depend on the noise amplitude** — so
amplitude (CV ≈ 0.2) and timescale can be calibrated independently.

---

## 6. What is missing

### Blocking — the work cannot be validated without these

1. **The HCI protein distributions.** Everything above is *prediction*. The brief's T8 leaves the
   loader as a placeholder ("data format to be supplied"). Until these arrive there is nothing to
   test the predicted CVs against, and no way to separate biological spread from measurement
   noise. This is the single biggest gap.
2. **An absolute scale (Ω).** `docs/DATA_PROVENANCE.md` states the model's species are in
   **arbitrary units**, calibrated to ratios and never to absolute concentrations. T2 (Gillespie)
   and T6 (discrete bursts) both need molecule *counts*. Minimum needed: GNP cell/nuclear volume,
   or an absolute copy number for at least one species, to anchor Ω. Without it the Poisson term
   `1/⟨p⟩` and the burst discretisation have no defined magnitude.
3. **Measured T_cc for GNP and MB.** µ = ln2/T_cc scales the entire dilution analysis and the
   extrinsic correlation time. The model produces ~22–23 h; whether that matches the real system
   should be an input, not an output.

### High value — would resolve a specific open question

4. ***Ccnd1* mRNA half-life in GNP or a neural progenitor.** Still the 12× lever on CyclinD1's
   correlation time (F2). Schwanhäusser says 14.0 h; the ActD literature says <1 h; after fix F1
   the model implies 52 min.
5. ***Ezh2* mRNA half-life in GNP.** Would resolve the F3 conflict — the model needs ~1.2 h to
   reproduce the phase gradients, Schwanhäusser measures 9.85 h in NIH3T3. Also moves EZH2's CV
   between 0.42 and 0.89.
6. **Burst kinetics in a neural progenitor.** Because `f` does not transfer between cell types
   (§5), the burst *frequency* for every gene is currently unconstrained for GNPs.
7. **Turnover for `Gli1`, `Gli3`, `Ptch1`, `Smo`, `Ccnb1`, `Fbxo5`, `Cdkn2d`** — absent from all
   three datasets. The Hedgehog arm has no measured turnover at all.

### From the brief's own source list, not yet fetched

8. ~~**Sigal 2006**~~ — **obtained**, see §5b. τm = 1.75 generations, CV = 0.20, both confirming
   the brief's assumed values.
9. **Eden 2011** (bleach-chase) — **blocked, needs institutional access.** science.org returns 403
   on the article, the `doi/suppl/` listing and every `downloadSupplement` pattern; not open
   access, no PMC record, and the 2012 *Nat Protoc* follow-up is closed too. Target: Supporting
   Online Material, Table S1 (~100 proteins, human H1299, bleach-chase). Its value is as an
   **arbiter** of the "Schwanhäusser runs long" question — an independent method in living cells —
   not as a per-gene resource; it is too small and too human to fill the Hedgehog-arm gaps.
10. **McShane 2016** — non-exponential, age-dependent degradation, and free vs complexed pools.
    Directly relevant: v44 carries many complexes (`CeP21`, `CaP21`, `RbE2f`, `E1C1`, `iPcna`,
    `iRc`) and currently degrades bound and free pools at the same rate.
11. **Zatulovskiy 2020** (Rb dilution), **Padovan-Merhar 2015** (burst frequency vs cell size,
    replication dosage compensation), **Lanz 2022** (proteome scaling), **Beentjes 2020** (exact
    solutions — the analytic test cases for T6's acceptance criteria).

### Conceptual gaps no download fixes

12. **Extrinsic noise for this system.** Magnitude and correlation time are taken from human cell
    lines. No mouse-progenitor equivalent exists that I found.
13. **Protein–protein correlation structure.** Gaudet 2012's point is that correlations between
    sampled proteins matter. Nothing measures them here; the single shared global OU is the only
    handle, and it imposes a rank-1 correlation structure by construction.
14. **Measurement-noise floor for the HCI.** Needed to subtract before comparing CVs. Usually
    estimated from a stained-but-uniform control or from technical replicates.

---

## Reproducing

```bash
python3 simulations/build_turnover_priors.py     # regenerate the CSV from data/priors/raw/
python3 simulations/analyze_noise_budget.py      # corrected, self-consistent CV / tau tables
```
