# Bursting × mRNA stability × chromatin: what exists, and what it says about PRC2

Answer to: *is there a dataset linking bursting to mRNA stability, and ideally to a ChIP-seq
signature?* Yes to both, and the chromatin one is directly about PRC2.

---

## 1. The short answer

| need | dataset | cell type | status |
|---|---|---|---|
| burst kinetics ← **ChIP-seq occupancy** | **Ochiai 2020** *Sci Adv* 6:eaaz6699 | mouse ESC | obtained (PMC7299619, OA) |
| **mRNA stability** genome-wide | **Herzog 2017** *Nat Methods* 14:1198 (SLAM-seq) | mouse ESC | obtained (PMC5712218, OA) |
| burst kinetics ← promoter **sequence** | Larsson 2019 Suppl. Table 2 | mouse fibroblast | obtained |
| burst kinetics | Larsson 2019 Suppl. Table 3 | mouse ESC | already held |
| synthesis **and** degradation, cell-cycle resolved | Battich 2020 *Science* 367:1151 (scEU-seq) | RPE1 | **blocked** — no PMC, 403 |

**Ochiai, Herzog and Larsson Table S3 are all mouse ES cells**, so they join gene-by-gene. That
matters more than it sounds — see §4.

---

## 2. Ochiai 2020 is the paper you are describing, and its headline result is about PRC2

They inferred burst kinetics genome-wide from allele-resolved scRNA-seq, then regressed them
against **ChIP-seq occupancy computed separately at the promoter and the gene body** for a large
panel of DNA-binding proteins. The central finding, verbatim:

> *"Promoter localization of PRC2 subunits (EZH2, SUZ12, and JARID2) correlated **inversely with
> burst frequency**, while they correlated **positively with burst size** and normalized intrinsic
> noise."*

and the paper's own summary:

> *"transcriptional bursting kinetics was regulated by a combination of promoter- and gene
> body–binding proteins, including the polycomb repressive complex 2 and transcription elongation
> factors."*

### Why this matters for v44 specifically

The hand-off brief's §3.2 proposes that "H3K27me3 enters as a burst-frequency (k_on) modulator".
**Ochiai measured that, genome-wide, and the direction is confirmed — but the specification is
incomplete in a way that matters.** PRC2 does two things at once:

| | direction | effect on mean | effect on noise |
|---|---|---|---|
| burst frequency `k_on` | **down** | ↓ | ↑ |
| burst size `B` | **up** | ↑ | ↑ |

The two partly **cancel for the mean** but **add for the noise**. Modelling only the frequency arm
(as the brief specifies) would therefore get the repression of mean CyclinD1 roughly right while
systematically **understating** the noise that PRC2 introduces — which is precisely the quantity
the whole exercise exists to estimate. Recommend the H3K27me3 coupling act on `k_on` **and** `B`,
in opposite directions.

There is a useful mechanistic reading: PRC2 makes a gene fire *less often but bigger* — a rarer,
lumpier transcriptional output at a similar mean. For a slowly-varying downstream protein that is
a substantial increase in cell-to-cell spread at little cost in mean expression.

### Two caveats, both important

**It is context-dependent, not a universal law.** Their SUZ12 knockouts went both ways: intrinsic
noise and burst size of *Dnmt3l* and *Dnmt3b* fell, but *Peg3*'s rose; burst frequency of *Dnmt3l*
rose while *Peg3*'s fell. Their own conclusion is that "promoter localization of PRC2-related
factors influences bursting properties in a **gene-specific** manner". So this supports a
*mechanism and a prior direction*, not a transferable coefficient.

**They perturbed the complex, not the catalysis.** SUZ12 knockout removes PRC2 itself. The
correlations are with PRC2 *promoter occupancy*, not with H3K27me3 level per se. That maps onto a
distinction this repo already makes — the recent "EZH2i mechanism fix: catalytic inhibitor blocks
writing only, not the PRC2 complex". Under tazemetostat PRC2 stays bound, so Ochiai's
occupancy-driven effect on bursting would be **largely preserved**, while the mark decays. If that
is right, EZH2i should de-repress the *mean* without removing the PRC2-imposed burstiness — a
sharp, testable prediction that falls straight out of the existing model structure.

### What is deposited, and what is not

The supplementary tables are the **allele-resolved count matrices** (Tables S1/S2, ~119 MB each
uncompressed), not the fitted per-gene kinetics or the ChIP regression coefficients. To get
per-gene `k_on`/`B` with feature associations you would either re-infer from S1/S2 (the
`txburst` route, same as Larsson) or pull the values off the figures. Raw data: **GEO GSE132593**.
Their ChIP-seq came from public mESC compendia, so the same pipeline could be pointed at GNP/MB
ChIP-seq — which is the interesting version for you.

---

## 3. Sequence explains burst size a little and burst frequency not at all

Larsson's own Table S2 (three OLS regressions on core promoter elements) quantifies this:

| response | R² | strongest term |
|---|---|---|
| **burst size** | **0.047** | **TATA +0.602** (t = 4.62, p < 0.001) |
| **burst frequency** | **0.008** | nothing significant (TATA p = 0.21) |
| mean expression | 0.014 | gene length −0.107 |

So a TATA box makes bursts bigger, and **promoter sequence says essentially nothing about how
often a gene fires**. Larsson attribute frequency control to enhancers instead.

This closes a loop with the earlier finding in `protein_noise_data_exploration.md` §5, that burst
*size* partially transfers between cell types (ESC/fibroblast log–log r = +0.50) while burst
*frequency* does not (r = +0.13). Three independent lines agree:

- size is **promoter-intrinsic** — partly sequence-encoded, partly chromatin;
- frequency is **cell-state-determined** — enhancer- and chromatin-driven, hence not transferable.

**Practical consequence: ChIP-seq is the only route to burst frequency.** Sequence cannot give it
and another cell type cannot lend it. That is the strongest argument for using your own ChIP-seq.

---

## 4. The mESC join puts burst frequency into real units for the first time

Larsson reports burst frequency per **mRNA-degradation-rate** (`per δ`), not per hour — so it has
been unusable as a rate until now. Herzog supplies the missing denominator **in the same cell
type**, so `f[1/h] = f[per δ] · ln2 / t½(mRNA)`. 3927 genes have both.

| gene | t½ mRNA (h) | f (per δ) | **f (1/h)** | **h per burst** | B | ⟨m⟩/allele |
|---|---|---|---|---|---|---|
| **Cdkn1a** | 1.41 | 2.387 | 1.177 | **0.8** | 9.83 | 23.00 |
| Cdc25a | 1.40 | 1.142 | 0.565 | 1.8 | 6.35 | 4.64 |
| Cdt1 | 2.67 | 2.319 | 0.603 | 1.7 | 8.22 | 10.48 |
| **Ezh2** | 3.71 | 1.766 | 0.330 | **3.0** | 3.12 | 5.19 |
| Trp53 | 5.99 | 3.134 | 0.363 | 2.8 | 3.41 | 8.64 |
| E2f1 | 3.38 | 0.973 | 0.199 | 5.0 | 2.99 | 2.50 |
| Ccne1 | 4.01 | 0.969 | 0.168 | 6.0 | 14.88 | 10.65 |
| **Ccnd1** | 3.37 | 0.677 | 0.139 | **7.2** | 5.68 | 3.54 |
| Cdk4 | 8.04 | 1.633 | 0.141 | 7.1 | 43.22 | 25.56 |
| Gli2 | 1.90 | 0.365 | 0.133 | 7.5 | 7.22 | 1.71 |
| Mycn | 2.00 | 0.271 | 0.094 | 10.6 | 18.50 | 3.95 |
| Ccna2 | 4.82 | 0.461 | 0.066 | 15.1 | 30.30 | 6.53 |
| **Skp2** | 2.23 | 0.147 | 0.046 | **21.9** | 12.76 | 1.79 |

Per allele; a diploid cell fires twice as often. A GNP cycle is ~22 h.

**`Skp2` fires roughly once per cell cycle per allele, in bursts of ~13 transcripts.** Combined
with its short mRNA (2.2 h) and short protein (3.4 h) half-lives — the reason it was flagged in
§4 of the exploration doc as the only genuinely fast-fluctuating species in the model — this is a
node that spends most of the cycle decaying from a rare, large input. That is the textbook regime
for burst-driven basin escape, and it sits on the Skp2–p27 commitment switch. It strengthens the
case for putting Skp2 on the T6 list.

---

## 5. Herzog also collapses the biggest open parameter

`Ccnd1` mRNA half-life has been the project's largest single uncertainty — a 12× lever on
CyclinD1's noise correlation time (`protein_noise_fixes.md` §F2):

| source | *Ccnd1* mRNA t½ | method |
|---|---|---|
| Schwanhäusser 2011 | 14.0 h | SILAC + transcription-rate inference, NIH3T3 |
| **Herzog 2017** | **3.37 h** | **SLAM-seq metabolic labelling, mESC** |
| v44 after the `/60` fix (F1) | 0.87 h | model constant |

SLAM-seq is the better measurement — metabolic labelling, no transcription block, no
actinomycin-D artefact. And Schwanhäusser runs long here too, by the same sort of margin as for
protein: **median ratio 2.31× on 3071 paired genes** (Herzog median 4.2 h vs Schwanhäusser 10.1 h).

So the 12× spread narrows to roughly 4× (0.87 h in the model vs 3.37 h measured), and the right
default is **~3.4 h, not 14 h and not 52 min**. At t½ = 3.4 h, CyclinD1's noise correlation time is
≈ 5 h — sub-cycle but not fast, i.e. intermediate between the two extremes we have been hedging
between.

`Ezh2` moves the same way: **3.71 h** (Herzog) vs 9.85 h (Schwanhäusser) vs the ~1.2 h v44 needs to
reproduce its phase gradients. The F3 conflict narrows from 8.4× to 3.1× but does not disappear.

**And it fills most of the Hedgehog gap.** Genes with no mRNA half-life in any previously held
dataset now have one: `Gli1` 3.05 h, `Gli2` 1.90 h, `Mycn` 2.00 h, `Cdkn1a` 1.41 h, `Ccne1` 4.01 h,
`Ccnb1` 3.09 h, `E2f1` 3.38 h, `Trp53` 5.99 h, `Cdc20` 4.40 h, `Cdt1`, `Dtl`, `Wee1`, `Cdc25a`,
`Fbxo5`. Still missing: `Gli3`, `Ptch1`, `Smo`, `Rb1`, `Cdkn1b`.

*Caveat: mESC is not a GNP. These are progenitor cells rather than fibroblasts, which is closer,
but the numbers are still borrowed. mRNA half-lives transfer better between cell types than burst
frequency does (§3) — but not perfectly.*

---

## 6. What to do with your own ChIP-seq

If you have GNP/MB ChIP-seq, the Ochiai design is directly reusable and is the one analysis that
could fill the burst-frequency gap properly:

1. Compute ChIP RPM at **promoter** and **gene body separately** — Ochiai's key methodological
   point is that the two carry different information (PRC2 at the promoter, elongation factors in
   the body).
2. For the marks you already have — H3K27me3 above all, ideally with EZH2/SUZ12 occupancy — regress
   against burst kinetics. Absent GNP burst data, use the Ochiai/Larsson mESC coefficients as the
   prior and your ChIP as the per-gene predictor.
3. That yields a **per-gene, chromatin-conditioned `k_on` and `B` for the Hedgehog arm**, where no
   external dataset has anything — and it parameterises the H3K27me3 → bursting coupling from your
   own data rather than by assumption.

The strongest version, if the ChIP-seq spans conditions: **±vismodegib and ±EZH2i**. PRC2 occupancy
changes between those conditions are exactly the modulation the model posits, and having them
measured would turn the H3K27me3 → burst coupling from a hypothesis into a fitted relationship.

---

## 7. Still blocked

**Battich et al. 2020** *Science* 367:1151, scEU-seq. No PMC record; science.org 403s. It is the
only dataset giving **synthesis and degradation rates simultaneously, resolved across the cell
cycle** — and it was done in **RPE1**, the line Heldt 2018 was fitted in. That cell-cycle
resolution is unique and directly relevant to a model whose whole point is cycle-dependent
regulation. Added to `docs/protein_noise_data_requests.md`.
