# Validation targets — species & provenance annotation (pre-re-optimization audit, 2026-08-01)

**Why:** before re-optimizing we must be sure each data-driven target compares the **same molecular species**
the experiment measured (mRNA vs protein vs activity vs count) — fitting a model *protein* readout to an *mRNA*
measurement (or conflating GNP↔MB, published↔unpublished) corrupts the fit. This table pairs each target's
**measured species/assay** (from `data/validation_targets.json` `quantity`+`source`) against the **model species
the harness actually reads** (`simulations/validate_v44.py`), and gives a match verdict.

**Model species inventory** (which genes have BOTH transcript+protein vs a SINGLE lumped species):
`CyclinD1 = Cd_mRNA + Cd` · `EZH2 = EZH2m + EZH2` · `Gli1 = Gli1_mRNA + Gli1` — a choice exists;
`CyclinD2 = Cd2` · `MYCN = MYCN` · `Skp2 = Skp2` — single species (no transcript variant).

Verdict key: **✓ match** · **✗ mismatch (fixable — model has the right species, harness reads the other)** ·
**≈ single-species proxy (model has no transcript/protein split; mRNA fold assumed ≈ species fold)** ·
**⚠ provenance/cell-type/value issue (species aside)**.

## Level & ratio targets (species matters most here)

| Target | Cell type | Measured species / assay / source | Model reads | Verdict |
|---|---|---|---|---|
| CyclinD1 GNP+HHi/GNP (0.157) | P7 GNP | Ccnd1 **mRNA** — bulk RNA-seq §A | `Cd_mRNA` | ✓ |
| CyclinD1 MB+HHi/MB (0.144) | MB | Ccnd1 **mRNA** — RNA-seq §A | `Cd_mRNA` | ✓ but ⚠ **value conflict 0.144 vs 0.40** (compendium §A) |
| CyclinD1 MB/GNP (5.07) | MB↔GNP | Ccnd1 **mRNA** — Fig 4I RNA-seq | `Cd_mRNA` | ✓ (⚠ 5.07 Fig4I vs 7.58 raw — normalization fork wider than its 0.20 band) |
| EZH2i CyclinD1 fold GNP (2.2) | P7 GNP | Ccnd1 **mRNA** — Fig 3C qPCR (Tazemetostat) | `Cd_mRNA` | ✓ (genetic cKO 2.46× is a *separate* Fig 3A RNA-seq point) |
| **EZH2 MB/GNP (2.05)** | MB↔GNP | EZH2 **mRNA** — Fig 4J RNA-seq | **`EZH2` (protein)** | **✗ FIXABLE → should read `EZH2m`.** JSON already flags this. **Re-opt-critical.** |
| **Gli1 MB/GNP (6.9)** | MB↔GNP | Gli1 **mRNA** — raw RNA-seq (no fig pinned) | **`Gli1` (protein)** | **✗ FIXABLE → should read `Gli1_mRNA`** |
| **Gli1 GNP+HHi reduction (0.99)** | P7 GNP | Gli1 **mRNA** — RNA-seq §A | **`Gli1` (protein)** | **✗ FIXABLE → `Gli1_mRNA`** (low impact: >99% either way) |
| MYCN GNP+HHi/GNP (0.78) | P7 GNP | Mycn **mRNA** — RNA-seq §A | `MYCN` (single) | ≈ proxy |
| MYCN MB+HHi/MB (0.86) | MB | Mycn **mRNA** — RNA-seq §A | `MYCN` | ≈ proxy |
| MYCN MB/GNP (2.80) | MB↔GNP | Mycn **mRNA** — RNA-seq §A | `MYCN` | ≈ proxy |
| Skp2 MB/GNP (2.3) | MB↔GNP | Skp2 **mRNA** — scRNA + MB table (**unpublished**) | `Skp2` (single) | ≈ proxy; ⚠ unpublished, no fig |
| CyclinD2 MB/GNP (2.39) | MB↔GNP | Ccnd2 **mRNA** — Chahin RNA-seq | `Cd2` (single) | ≈ proxy |
| CyclinD2 GNP+HHi/GNP (0.66) | P7 GNP | Ccnd2 **mRNA** — RNA-seq | `Cd2` | ≈ proxy |
| CyclinD2 MB+HHi/MB (0.59) | MB | Ccnd2 **mRNA** — RNA-seq (trend, padj 0.076) | `Cd2` | ≈ proxy |
| EZH2 G0/cycling (0.6) | P7 GNP | EZH2 — **§D** (source says **transcript** ~0.5, pooled GNP+MB) | `EZH2` (protein) | **⚠ species+value: harness 0.6 protein vs source 0.5 transcript** |
| EZH2 transcript S/G0 (2.0) | MB | EZH2 **mRNA** — scRNA Fig 4A/4B | `EZH2m` | ✓ (⚠ pools GNP+MB, G1/S-vs-G0 bin) |
| EZH2 protein G2/G0 (1.48) | MB | EZH2 **protein** — Supp 6D IF | `EZH2` | ✓ |
| EZH2 Palbo mRNA drop (0.44) | MB | EZH2 **mRNA** — Fig 4F qPCR | `EZH2m`, **GNP** cond | ✓ species; **⚠ CELL-TYPE: measured MB, harness scores GNP+CDK4/6i** |
| HU EZH2-in-S boost (1.31) | MB | EZH2 **protein** — Fig 4G IF | `EZH2` (S vs S) | ✓ (⚠ alt-text conflict 1.31 vs 1.22) |

## Behavioral / timing / phase targets (no molecular-species issue; provenance notes)

| Target | Cell type | Measure / source | Model reads | Note |
|---|---|---|---|---|
| GNP+SHH cycles; GNP−SHH / +HHi / serum-starve arrest | P7 GNP | divisions — Wechsler-Reya '99 / inferred / definitional | MPF-peak division count | qualitative ✓ |
| MB cycles; MB+CDK4/6i arrest; +EZH2i no-rescue | MB | divisions — Chahin / **JP culture** / **model prediction** | division count | ⚠ CDK4/6i "arrest" is really 16–18% reversible pRb⁺ (JP culture); no-rescue is a **prediction**, fitting it is circular |
| Period GNP (22 h) | P7 GNP | cell-cycle length — Nakashima 15.9h + Contestabile 16.25h | MPF period | **⚠ 22 h is a MISCITATION; real ~16 h** |
| MB 2N (G0+G1) 68.2% ; MB S 15.7% | MB | count-fraction — **cell-culture microscopy** (mislabeled "flow") §F | age-weighted count-frac | ⚠ microscopy not flow; G0/G1 unresolvable (2N pool only) |
| MB G2+M 2.5 h | MB | duration — Fujita (direct) | duration×Tc | ✓ duration↔duration |
| MB HU S fold 1.36 ; HU G2 fold 0.23 | MB | count-fraction fold — microscopy §F2 | classify fold | ⚠ the structural HU pair; microscopy |

## Conflation risks, ranked by impact on the re-optimization

1. **EZH2 MB/GNP scored on protein, measured as mRNA (Fig 4J RNA-seq).** *Directly bears on the un-map re-opt* — this was one of the targets I tuned `Kez_cd` against. If it should be `EZH2m`, the whole EZH2-fold behaviour (and the Kez_cd value) must be re-derived on the transcript. **Fix before re-opt.**
2. **Gli1 folds scored on protein, measured as mRNA.** Model has `Gli1_mRNA`; switch the readout. Affects Gli1 MB/GNP (6.9), which the CyclinD1-drive re-opt also moves.
3. **Single-species proxies (MYCN, Skp2, CyclinD2) fit to mRNA.** For steady-state folds this is usually acceptable *if* per-transcript translation is constant GNP↔MB — but the p21/miR-17-92 case shows that assumption can fail in exactly this system. Decide per gene: accept as transcript-proxy, or add a transcript species.
4. **EZH2 Palbo drop: measured MB, scored GNP.** Cell-type mismatch (species OK).
5. **EZH2 G0/cycling: transcript 0.5 (pooled) vs protein 0.6 (deep G0).** Pick the intended species+value.
6. **Value/provenance:** CyclinD1 MB+HHi 0.144 vs 0.40; HU-EZH2-S 1.31 vs 1.22; Period 22 vs 16; MB phase = microscopy; Gli1/Skp2 figures unpinned; Skp2 unpublished.

## Resolutions (JP 2026-08-01) — applied to the harness

- **EZH2 MB/GNP** → RNA-seq transcript → harness now scores **`EZH2m`** (was `EZH2` protein). ✅ fixed.
- **Gli1 folds** → RNA-seq transcript → harness now scores **`Gli1_mRNA`** (was `Gli1` protein). ✅ fixed.
- **MYCN, Skp2, CyclinD2** → bulk RNA-seq transcript; **accepted as single-species transcript-proxies** (no
  transcript species added). Skp2 = JP unpublished; all are JP's high-quality bulk RNA-seq. ✅ recorded.
- **EZH2 Palbo drop** → measured in **MB** → harness now scores `MB+CDK4/6i` / `MB` (was GNP). ✅ fixed.
- **CyclinD1 MB+HHi = 0.144** confirmed from raw counts (3971/27590); the 0.40 was wrong. ✅ (already 0.144).
- **HU EZH2-in-S = 1.31** (use the higher of 1.31/1.22). ✅ (already 1.31).
- **Period GNP → 16 h** adopted (Nakashima 15.9 + Contestabile 16.25; 22h was a miscitation) → target now 16h;
  **requires a mu re-tune** (model runs ~22h). ✅ target changed; recal pending.
- **MB phase fractions = cell-culture microscopy** (not flow) → labels fixed. ✅
- **EZH2 dual modality:** EZH2 has BOTH scRNA (transcript) AND IF-inferred protein data. Mapping now explicit
  (see below); wiki Data & Evidence updated. The G0/cycling (0.6) target is the **IF-protein / deep-G0** measure
  (kept on `EZH2`), distinct from the scRNA within-cycle transcript (S/G0, on `EZH2m`).

### ⚠ OPEN — the one JP call still needed: CyclinD1 MB/GNP = 5.07 (Fig 4I) vs **7.58 (raw DESeq2 counts)**
JP's raw table gives Ccnd1 MB/GNP = 27590/3639 = **7.58**; every other fold from that table matches the existing
targets exactly, so the raw table is the source of truth for the rest — but MB/GNP is the one place it diverges
from the Fig 4I-normalized 5.07 the harness targets. This is *the* load-bearing, tightest-band target and the
re-opt's sticking point, so which value we fit to matters a lot. **Held at 5.07 pending JP's call (5.07 vs 7.58).**

### EZH2 target → model-species → assay map (the dual-modality clarification)
| EZH2 target | model species | assay |
|---|---|---|
| EZH2 MB/GNP (2.05) | `EZH2m` transcript | bulk RNA-seq (Fig 4J) |
| EZH2 transcript S/G0 (2.0) | `EZH2m` transcript | scRNA-seq (Fig 4A/4B) |
| EZH2 Palbo drop (0.44, MB) | `EZH2m` transcript | qPCR (Fig 4F) |
| EZH2 protein G2/G0 (1.48) | `EZH2` protein | IF (Supp 6D) |
| HU EZH2-in-S (1.31) | `EZH2` protein | IF (Fig 4G) |
| EZH2 G0/cycling (0.6) | `EZH2` protein | IF-inferred, deep-G0 serum-starve |
