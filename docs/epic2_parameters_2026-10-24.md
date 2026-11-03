# Why epic2 keeps needing to be re-tuned

2026-10-24. Read of epic2 0.0.54 source (`biocore-ntnu/epic2`), with its
background model reimplemented locally and evaluated on our actual libraries.

## How epic2 decides what a peak is

Everything flows from one number:

```
lambda = total_chip_reads * bin_size / effective_genome_length
```

From `lambda` epic2 derives, in `src/statistics.pyx`:

- **`island_enriched_threshold` k** — the smallest integer k with
  Poisson(lambda) survival <= 0.20. A bin needs >= k reads to be *eligible*.
- **window score** = `-log(Poisson.pmf(count, lambda))`, and **0 if count < lambda**.
- **`single_gap_factor`** = P(X < k), the chance a bin is ineligible.
- **`gap_contribution`** = `1 + sum_{i=1..gaps} single_gap_factor^i`.
- **`score_threshold`** — set so the expected number of background islands is
  `E_VALUE`, hard-coded to 1000.

Islands then grow across bins while the distance to the next *eligible* bin is
`<= gaps_allowed * bin_size + 2`, accumulating score only from eligible bins,
and are kept if the total score clears `score_threshold`. Only afterwards are
ChIP counts compared to input and BH-FDR applied. So there are two independent
filters, and `--false-discovery-rate-cutoff` touches only the second.

## Finding 1: mm39 is not in epic2, so our calls silently used hg19

epic2 ships effective-genome-fraction tables for 90 genomes. **mm39 is not one
of them** — only mm9 and mm10. So `--genome mm39` is impossible, which is
presumably why our scripts pass `--chromsizes` instead.

But in `genome_info.pyx::egl_and_chromsizes`, `--chromsizes` *without* `-egf`
takes this branch:

```python
elif have_chromsizes:
    chromsizes = create_genome_size_dict_custom_genome(args["chromsizes"])
    egf = get_effective_genome_fraction(args["genome"], read_length)  # genome defaults to "hg19"
    egl = egf * sum(chromsizes.values())
```

`--genome` defaults to **hg19**, so every one of our mm39 calls used **hg19's**
effective-genome fraction against the mouse chromosome total:

| | fraction at readlength 100 |
|---|---|
| hg19 (what was used) | 0.9020 |
| mm10 (what we want) | 0.8881 |

A 1.6% error — small, and not the cause of anything dramatic, but it is
silently wrong in every call we have made and it is free to fix:
**always pass `-egf 0.8881`**.

## Finding 2: the eligibility threshold is an integer step function of depth

With mm39 (2.7282 Gb total) and mouse egf, effective genome length is 2.416 Gb:

| library | reads | bin | lambda | min tags/bin | P(X>=k) |
|---|---|---|---|---|---|
| Scott rep1 | 96,428,049 | 100 | 3.992 | **7** | 0.110 |
| Scott rep1 | 96,428,049 | 200 | 7.984 | 11 | 0.182 |
| Scott rep2 | 8,808,239 | 100 | 0.365 | **2** | 0.052 |
| EL GNP1 | 4,221,855 | 100 | 0.175 | **1** | 0.160 |

At identical settings, a Scott bin needs 7 reads to count and an EL bin needs
1. k is an integer, so it jumps:

| depth | lambda | k |
|---|---|---|
| 2M | 0.083 | 1 |
| 8M | 0.331 | 2 |
| 32M | 1.325 | 3 |
| 64M | 2.649 | 5 |
| 96M | 3.974 | 7 |
| 128M | 5.299 | 8 |

**This is the main reason tuning never transfers between samples.** Two
libraries differing 2x in depth land on different integers and behave
differently at the same parameters. The quantisation also bites hardest at low
depth: at lambda = 0.175 the achievable survival values are 0.160 (k=1) or
0.013 (k=2), so the nominal 0.20 window p-value cannot be met and a shallow
library is implicitly held to a much stricter threshold than a deep one.

## Finding 3: the gap correction saturates, the merging does not

`gap_contribution` is a geometric series in `single_gap_factor` = P(X < k),
so it converges to `1 / (1 - P(X<k))`. For Scott rep1 at bin 100, P(X<k) =
0.8902 and the ceiling is 9.105:

| gaps | merge distance | gap_contribution | % of ceiling |
|---|---|---|---|
| 1 | 102 bp | 1.890 | 20.8% |
| 3 | 302 bp | 3.388 | 37.2% |
| 6 | 602 bp | 5.072 | 55.7% |
| 10 | 1,002 bp | 6.573 | 72.2% |
| 20 | 2,002 bp | 8.314 | **91.3%** |
| 50 | 5,002 bp | 9.081 | **99.7%** |
| 100 | 10,002 bp | 9.105 | 100.0% |

The *statistical* correction for permitting gaps is essentially spent by
gaps ~20. The *physical* merge distance keeps growing linearly for ever. So
past gap 20 epic2 joins islands further and further apart without raising the
score bar to compensate.

That is exactly what our sweep measured — coverage climbing with no plateau:
3.36 / 5.62 / 9.00 / 12.95 / 18.65% of the genome at gaps 1 / 3 / 6 / 10 / 20.
**There is no plateau to tune against, and looking for one is why parameter
choice feels arbitrary.** The gap has to be set against an expected coverage
instead.

## What to do

1. **Always pass `-egf 0.8881`** for mm39. Without it we are using hg19's.
2. **Equalise depth before calling** if calls are to be compared across
   samples. Downsampling to a common depth fixes lambda, which fixes k, which
   is the only way identical parameters mean the same thing in two libraries.
   This matters more than any gap choice.
3. **Keep gaps <= 20**, and prefer 3-6. Beyond 20 the model has stopped
   compensating and the result is merging, not detection.
4. **Pick the gap against an expected genome coverage**, not by looking for
   stability. For GNP H3K27me3 against the project's ~5% figure that is
   gap 3 (5.62%, mean domain 2.65 kb).
5. Remember `-fdr` only filters the ChIP-vs-input stage. The island-formation
   stage is governed by the hard-coded E_VALUE of 1000 and the 0.20 window
   p-value, neither of which is exposed on the command line.

## Not explained by any of this

Re-running the stored Scott command verbatim gives 39,256 peaks / 9.00% where
the stored file has 1,195 / 0.11%. No parameter in this document accounts for
a 78x difference from identical inputs, which supports the conclusion that the
stored file came from a failed run rather than a bad choice of settings.
