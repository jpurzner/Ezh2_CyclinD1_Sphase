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

## 7. Eden 2011 and McShane 2016 (supplied 2026-09-24): the "Schwanhäusser runs long" question is settled

Both were blocked to automated download and were supplied from institutional access. Eden's
supplementary tables could not be exported, but **Table S4 is machine-readable inside the SOM PDF**
and was parsed from it — 100 proteins, matching the paper's stated count
(`data/priors/eden2011_bleachchase_halflives.csv`).

### 7.1 Schwanhäusser reads ~1.74× long, and we now know why

McShane measured pulse-chase turnover **in mouse NIH3T3 — the same line as Schwanhäusser** — so
this is a clean method-vs-method comparison with cell type held fixed. On 2549 shared genes:

| | median t½ |
|---|---|
| McShane (pulse-chase, 1-state fit) | 22.9 h |
| Schwanhäusser (SILAC) | 47.0 h |
| **median ratio S/M** | **1.74×** |

Schwanhäusser is longer for **76%** of genes (log–log r = +0.54). Decisive, though, is the
breakdown by McShane's degradation class:

| class | n | median McShane | median Schwanhäusser | ratio |
|---|---|---|---|---|
| **ED** (exponential) | 1138 | 32.9 h | 46.6 h | **1.24×** |
| **NED** (non-exponential) | 319 | 15.2 h | 56.0 h | **3.41×** |
| UN (undetermined) | 1092 | 20.6 h | 45.8 h | 2.07× |

**The bias is concentrated almost entirely in the non-exponential proteins.** That is exactly what
you expect from fitting one exponential to two-pool kinetics: the fast pool is gone before the
first useful timepoint, so the fit tracks the stable pool and reports a long half-life. So
"Schwanhäusser runs long" is confirmed *and* explained — it is not a global calibration offset but
a model-misspecification artefact affecting ~10–40% of proteins.

**Practical rule:** for an ED protein, Schwanhäusser is fine (1.24×). For an NED protein it should
not be used at all.

Eden provides an independent third reading, in living human H1299: median half-life **8.2 h**
(range 0.8–22.5). But note the convention difference in §7.3 before comparing it to anything.

**Eden covers none of our 33 genes** — as anticipated when it was requested, it is a global
arbiter, not a per-gene resource. That is the answer it was wanted for, and it delivered it.

### 7.2 Non-exponential degradation is real for this model, and EZH2 resolves cleanly

McShane classifies **11% of NIH3T3 and 10% of RPE proteins as NED** — confirming the brief's
"~10%". For NED proteins the two-state split is enormous: median `kA` = 0.419/h (t½ 99 min) decaying
to `kB` = 0.0038/h (t½ 181 h), a median **138×** ratio.

Among v44's species, **`Cdk4`, `Pcna` and `Ccnb1` are classified NED** (their fitted `kB` is at or
near zero, i.e. the assembled pool is effectively stable on the experiment's timescale). v44
currently degrades free and complexed pools at the **same** rate, so this is a structural
simplification, not a parameter choice — which is precisely why McShane was worth getting.

**EZH2 resolves the earlier ambiguity.** In NIH3T3 it is classified ED but its two-state fit still
splits sharply: `kA` = 1.465/h (t½ 28 min, unassembled) versus `kB` = 0.0569/h (t½ 12.2 h,
assembled). Which value should v44's single `EZH2` species match?

| candidate | t½ | τ with dilution (T_cc=16 h) | v44 kDeEZ / bio |
|---|---|---|---|
| Schwanhäusser 1-exp | 22.0 h | 13.4 h | 1.1× |
| **McShane assembled pool (`kB`)** | **12.2 h** | **10.0 h** | **1.4×** |
| McShane 1-state blend | 5.0 h | 5.5 h | 2.6× |
| McShane fast pool (`kA`) | 0.5 h | 0.7 h | 21.5× |

EZH2 in the model is the **catalytic subunit of assembled PRC2** — it writes H3K27me3 — so the
assembled-pool rate `kB` is the correct comparator, not the 1-exponential blend that averages in
unassembled protein destined for rapid destruction. v44's `kDeEZ` (τ = 14.3 h) sits within 1.1–1.4×
of both defensible values. **The constant is well-placed**; the earlier claim that it matched
Schwanhäusser to 0.9× was right in magnitude but for a partly coincidental reason (that comparison
was mRNA-limited, not protein-limited).

The same pattern holds for the other PRC2 subunits — `Suz12` 21.2 h (Schwanhäusser) vs 4.9 h
(McShane 1-state) vs 6.9 h (`kB`); `Eed` 37.1 / 9.0 / 26.6 h.

### 7.3 Four datasets, three different dilution conventions — do not pool them

This is now the main provenance hazard in the priors table, and each row carries its convention:

| source | measures | class | `growth_corrected` |
|---|---|---|---|
| **Mathieson 2018** | degradation only (post-mitotic cells, no dilution exists) | (b) | yes |
| **Schwanhäusser 2011** | degradation only (established empirically in §1) | (b) | yes |
| **Eden 2011** | degradation **+ dilution** (SOM Box S1 defines α as the sum) | **(a)** | **no** |
| **McShane 2016** | pulse-chase in cycling cells; convention not yet verified | (c) | unknown |

Subtracting µ from Eden before comparing is mandatory. At H1299's ~21 h doubling, Eden's median
8.2 h corresponds to a degradation-only half-life of ~12.5 h.

### 7.4 Correction: the extrinsic noise component is real and dominant

§5b cautioned that Sigal's τm = ln2/(degradation + dilution) already contains the dilution v44
models explicitly, so setting τ_ext = 1.75 generations would double-count. **Eden and Sigal used
the same lab's same CD-tagged H1299 clones**, so that can now be tested rather than assumed.

Sigal's theory predicts τm (in hours) = ln2/α = the bleach-chase half-life. On the 6 proteins
present in both:

| protein | Eden t½ (h) | Sigal τm (gen) | τm in h @ T_cc = 21 h | ratio |
|---|---|---|---|---|
| CIRBP | 6.9 | 0.90 | 18.9 | 2.74 |
| DDX46 | 7.9 | 1.55 | 32.6 | 4.12 |
| DDX5 | 8.7 | 0.85 | 17.9 | 2.05 |
| H2AFV | 13.9 | 1.80 | 37.8 | 2.72 |
| HMGA1 | 11.1 | 1.95 | 41.0 | 3.69 |
| LMNA | 12.8 | 1.60 | 33.6 | 2.62 |
| | | | **median** | **2.73** |

The observed memory is **~2.7× longer than the removal rate alone predicts** (all six between 2.0
and 4.1; none near 1.0). Reading it the other way, the theory would only hold if T_cc were 7.7 h,
and H1299 doubles in ~20–24 h.

**So the extrinsic component is not a bookkeeping artefact of dilution — it is the dominant
contributor to protein memory**, roughly tripling the correlation time beyond what degradation and
dilution give. My §5b caution was too strong: an extrinsic OU with τ ≈ 1.5–2 generations is
justified by measurement, and the brief's §3.1 specification was right.

The calibration target is now concrete: **with the model's own removal rates in place, the added
extrinsic drift should extend each species' autocorrelation half-time by ≈2.7×.** n = 6, so treat
the factor as order-2–3 rather than precisely 2.73.

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
9. ~~**Eden 2011**~~ and ~~**McShane 2016**~~ — **both obtained** (supplied 2026-09-24), see §7.
   Together they settle the "Schwanhäusser runs long" question: 1.74× in the same cell line, with
   the bias concentrated in non-exponential proteins.
10. **Zatulovskiy 2020** (Rb dilution), **Padovan-Merhar 2015** (burst frequency vs cell size,
    replication dosage compensation), **Lanz 2022** (proteome scaling), **Li, Bickel & Biggin 2014**
    (the Schwanhäusser reanalysis), **Beentjes 2020** (exact solutions — the analytic test cases for
    T6's acceptance criteria). **All reachable without help** — the first four via PMC author
    manuscripts, Beentjes via arXiv:1912.04337. Mine to fetch, see
    `protein_noise_data_requests.md` §C.

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
