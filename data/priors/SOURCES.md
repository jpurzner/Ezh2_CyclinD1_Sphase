# Prior-data provenance (protein/mRNA turnover and transcriptional bursting)

Downloaded 2026-09-24. Raw files are in `raw/`; `turnover_burst_priors.csv` is generated from
them by `simulations/build_turnover_priors.py` (re-run it to regenerate; it never imputes).

| file | source | sha256 (first 16) | bytes | contents |
|---|---|---|---|---|
| `schwanhausser2011_TableS3.xls` | Schwanhausser et al. 2011 Nature 473:337, Suppl. Table 3 | `5343f73b4014a21c…` | 3,394,048 | mouse NIH3T3; protein & mRNA copy number, protein & mRNA half-life, transcription rate (vsr), translation rate (ksp); 5029 rows |
| `larsson2019_TableS1_fibroblast_burst.xlsx` | Larsson et al. 2019 Nature 565:251, Suppl. Table 1 | `157eb2815d64cd9f…` | 1,295,813 | mouse primary fibroblast, allele-resolved (C57/CAST); burst frequency, burst size, mean expression; 6962 genes |
| `larsson2019_TableS3_ESC_burst.xlsx` | Larsson et al. 2019 Nature 565:251, Suppl. Table 3 | `e7b790a78c6f8540…` | 1,041,762 | mouse embryonic stem cells, allele-resolved; same columns; 6179 genes |
| `mathieson2018_SupplData1.xlsx` | Mathieson et al. 2018 Nat Commun 9:689, Suppl. Data 1 | `76c5951eeaeb4c19…` | 12,327 | time-point / experimental design legend only (no per-gene values) |
| `mathieson2018_SupplData2.xlsx` | Mathieson et al. 2018 Nat Commun 9:689, Suppl. Data 2 | `e0503f49230680cb…` | 3,338,308 | protein half-lives in B cells, NK cells, monocytes, hepatocytes and MOUSE NEURONS; 8805 high-quality rows |
| `mathieson2018_SupplData3.xlsx` | Mathieson et al. 2018 Nat Commun 9:689, Suppl. Data 3 | `e1e619c98b966c35…` | 20,706 | proteasome subunit half-lives + exchangeability classification |
| `sigal2006_SuppNotes.pdf` | Sigal et al. 2006 Nature 444:643, Suppl. Information | `9ecd5e3cdd3117ec…` | 212,999 | Suppl. Tables 1-2: per-protein mixing time tau_m (generations) and CV, 20 CD-tagged H1299 clones; the extrinsic-noise prior |

## URLs

* `schwanhausser2011_TableS3.xls`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature10098/MediaObjects/41586_2011_BFnature10098_MOESM304_ESM.xls>
* `larsson2019_TableS1_fibroblast_burst.xlsx`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0836-1/MediaObjects/41586_2018_836_MOESM3_ESM.xlsx>
* `larsson2019_TableS3_ESC_burst.xlsx`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0836-1/MediaObjects/41586_2018_836_MOESM5_ESM.xlsx>
* `mathieson2018_SupplData1.xlsx`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM4_ESM.xlsx>
* `mathieson2018_SupplData2.xlsx`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM5_ESM.xlsx>
* `mathieson2018_SupplData3.xlsx`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM6_ESM.xlsx>
* `sigal2006_SuppNotes.pdf`
  <https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature05316/MediaObjects/41586_2006_BFnature05316_MOESM1_ESM.pdf>
  Parsed to `data/priors/sigal2006_mixing_times.csv` (40 rows = Suppl. Tables 1 and 2).

## Caveats that the CSV carries per row, restated

### Schwanhausser: this is the corrigendum version, and the half-lives run long

The hosted `MOESM304` file has an internal last-saved date of **2012-11-07**, i.e. it is the
file swapped in for the **2013 corrigendum** (Nature 495:126, doi:10.1038/nature11848), not the
2011 original. Good — but note that Li, Bickel & Biggin 2014 (PeerJ 2:e270) re-analysed the same
data and argued protein copy numbers and translation rates are several-fold higher.

More importantly for this project, the mRNA half-lives are systematically long. `Ccnd1` is
reported at **14.0 h** (replicates 10.5 and 18.7), where actinomycin-D chase experiments
usually put *Ccnd1* mRNA under an hour. The CSV records 14.0 h with its source; it does not
endorse it. `simulations/analyze_noise_budget.py` sweeps this value instead of choosing one,
because it is the single number that decides whether CyclinD1 noise is fast or frozen.

**`growth_corrected` is recorded as `unknown` for every Schwanhausser row.** The Methods are
said to describe a growth correction; that has *not* been re-read and verified here. Until it
is, treat these as class (c) for the T4 dilution audit — do not add `-mu*x` on top of them.

### Larsson: burst frequency is not in per-hour units, and is per allele

The telegraph-model parameters (`kon`, `koff`, `ksyn`) are inferred in units of the **mRNA
degradation rate**, so the reported burst frequency is dimensionless (`per delta`). Convert with

```
f[1/h] = f[/delta] * ln2 / t_half_mRNA[h]
```

which means **the burst frequency inherits all the uncertainty in the mRNA half-life above.**
Burst size `B` (= ksyn/koff) is dimensionless and is *not* affected — which is why `B` is the
more trustworthy half of this dataset.

Values are **per allele**. A diploid cell has ~2x the burst frequency and ~2x `<m>`; the CSV
stores the per-allele number and the analysis script doubles it.

### Mathieson: post-mitotic, therefore dilution-free

The `Mouse Neurons` columns are primary post-mitotic neurons. They do not divide, so there is
no dilution term and the fitted half-life **is** the degradation half-life. These are the only
unambiguous **class (b)** values available, and they are marked `growth_corrected = yes`.

This matters most for `Rb1` (**18.9 h** in neurons). At T_cc = 16 h, dilution (mu = 0.043/h)
and degradation (0.037/h) are the same size for Rb — which is exactly the Zatulovskiy 2020
size-dilution premise, now with a number attached.

### Eden 2011 — NOT obtained

Eden et al. 2011 *Science* 331:764, "Proteome half-life dynamics in living human cells"
(doi:10.1126/science.1199784) is the fourth requested turnover dataset and **could not be
downloaded**. science.org returns HTTP 403 (Cloudflare) on the article, the `doi/suppl/` listing
and every `downloadSupplement` pattern tried; the paper is not open access, has no PMC record, and
the follow-up protocol (Eden et al. 2012 *Nat Protoc*, PMID 22461069) is closed as well.

It needs institutional access. What it would add: **bleach-chase protein half-lives for ~100
proteins in living human H1299 cells** — a fourth, methodologically independent turnover estimate.
Its value here is as an *arbiter*, not as a per-gene resource: it is small, human, and a cancer
line, so it will not cover our missing Hedgehog-arm genes. But it measures half-lives in *living*
cells by a completely different route from SILAC pulse-chase, so it is the cleanest available
check on whether Schwanhäusser's half-lives run systematically long (see the Ccnd1 14.0 h vs <1 h
conflict). Target: Supporting Online Material, Table S1.
