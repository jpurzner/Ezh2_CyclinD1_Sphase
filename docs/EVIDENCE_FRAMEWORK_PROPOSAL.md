# Proposal: an evidence framework for target provenance

**Status:** proposal, nothing implemented. Written 2026-07-25.

The problem: a target value like "CyclinD1 MB+HHi/MB = 0.144" currently exists as four
unlinked copies in four files, and the thing that made it 0.144 — a specific RNA-seq
comparison, in a specific cell type, with a specific normalization — exists only as an
English sentence that gets re-typed each time it is quoted. When the sentence is re-typed
imperfectly, nothing notices. That is the whole failure mode, and everything below is
aimed at it.

---

## 1. Where information is actually leaking

Six distinct leaks, roughly in order of how much damage each is doing.

### 1.1 The harness is the source of truth, and it is a pile of literals

`simulations/validate_v44.py` scores the model against numbers typed directly into
`check()` calls (lines 263–334). `data/validation_targets.json` was written to be the
"single source of truth" — but **no code reads it**:

```
$ grep -rn "validation_targets" --include=*.py .   # → no matches
```

So the JSON is a *description* of the harness written alongside it, and descriptions drift
from what they describe. The drift is already visible and already documented in the JSON's
own `OPEN_ISSUES`: the period target is 22.0 in the harness and 16.25 in the spec;
`ezh2_palbo_drop` is measured in MB and scored on GNP; CyclinD1 MB+HHi is 0.144 in two
places and 0.40 in a third. Those were caught by a manual three-source audit. The next ones
will not be, because manual audits do not run on every commit.

### 1.2 Evidence has no identity

Sources are free-text strings:

> `"source": "Chahin Fig 4F (qPCR, MB: EZH2 0.444x, -55.6%, E2f-gated). CELL-TYPE MISMATCH: ..."`

This is genuinely informative prose, and it is unqueryable. There is no way to ask *"what
else does Fig 4F support?"* or *"if Fig 4F gets re-quantified during revision, which targets
and which parameters move?"* Fig 4F is in fact cited by at least three separate targets
under three differently-worded strings. A re-quantification would have to be caught by
memory.

### 1.3 Primary lab measurements enter the repo through conversation

This is the most serious one. Several load-bearing numbers are tagged like this:

- `mb_cdk46i_arrest` → *"JP experimental (culture, provided in-session)"*
- `mb_vismo_ezh2i_rescue` → *"JP experimental (culture): pRb+ falls to ~1/4 ... EZH2i rescues to ~3/4"*
- `cd_g01_cv` → *"JP flow"*
- `skp2_mb_gnp` → *"Purzner lab, unpublished; NOT in compendium/provenance audit"*

These are real primary observations with no primary record anywhere in the repo. There is
no cell line, no N, no replicate structure, no date, no gating definition, no raw file
pointer. The 16–18% pRb+ figure and the 100/25/75 rescue are structurally load-bearing —
the rescue contrast is arguably the central claim the model exists to explain — and they
survive only as sentences in derived documents and in chat transcripts. If a reviewer asks
"what was N?", the answer is currently not recoverable from the repo.

### 1.4 External literature has no bibliography

~40 distinct author-year strings appear across the docs (`Heldt 2018` ×11, `Luo 2008` ×8,
`Moser 2018` ×7, `Nakashima 2015` ×5 …). There is no `.bib`, no DOI or PMID anywhere, and
the same work appears under variant strings (`Heldt 2018` and `Heldt et al. 2018`;
`Contestabile 2009` and `Contestabile 2013` — which may be two papers or one mis-typed).
`*.pdf` is gitignored (correct, for copyright), so `papers/` holds exactly one file. The
literature the model inherits from is, as far as the repo is concerned, hearsay.

### 1.5 `docs/` is a flat pile with two incompatible kinds of file in it

45 markdown files, undifferentiated. Some are durable references maintained over time
(`DATA_PROVENANCE.md`, `PARAMETERIZATION.md`, `v44_MODEL_DESCRIPTION.md`). Most are dated
exploration notes that were true on the day they were written
(`buffering_parameters_2026-07-20.md`, `overnight_redo_2026-07-21.md`). Nothing marks which
is which, and there is no index, so a reader cannot tell whether a statement is current
doctrine or a thought someone had in July. Findings in the dated notes never get promoted
into the durable set — they just age.

### 1.6 Results are not tied to the model that produced them

`simulations/` holds 636 files: source, 19 result JSONs, `.log` files, and hundreds of
figures interleaved. A results JSON does not record which model version, which git SHA, or
which parameter set produced it, so a figure in the manuscript cannot be traced back to a
reproducible run. (Also: there is a 360 KB stray binary file `s` at the repo root — an
accidental shell redirect. Unrelated, but delete it.)

---

## 2. The proposal

One idea, three layers, one hard rule.

**The idea:** separate *what was observed* from *what the model is scored against*. Right
now these are fused. `0.157`, `"86% reduction"`, and `"compendium rounds to 0.14"` are three
representations of one observation living inside one string in one target. Split them: an
**evidence record** holds the observation exactly as reported, and a **target** holds a
pointer to that evidence plus the explicit arithmetic that turns it into a scoreable number.

**The hard rule:** *evidence is cited by ID, never by prose, and the harness reads the
target file rather than restating it.* Everything else is bookkeeping around that rule.

### Layer 1 — Evidence (`evidence/`)

One file per observation, Markdown with YAML frontmatter. Structured fields for the things
you want to query; free prose body for the caveats, which is where most of the value is and
which JSON handles badly.

Two namespaces:

- `evidence/lab/EXP-####-slug.md` — Purzner-lab measurements, published or not
- `evidence/lit/LIT-author-year.md` — external published work

Worked example, the messiest current target:

```markdown
---
id: EXP-0012
title: pRb+ cycling fraction in MB under palbociclib, and after washout
cell_type: MB
modality: immunofluorescence
assay: pRb (Ser807/811, D20B12) positive fraction
cell_line: TBD                      # ← unfilled fields are the point; they show the gap
date_measured: TBD
n_replicates: TBD
reported_by: JP, in-session 2026-07-25
publication_status: unpublished
raw_data: TBD                       # path or "not in repo"
observations:
  - quantity: pRb+ fraction, palbociclib vs DMSO
    value: [16, 18]
    unit: "% of DMSO"
  - quantity: pRb+ fraction, 72h palbociclib then washout
    value: 50
    unit: "% of DMSO"
corroborated_by: [LIT-cooksangar-2017]
---

Palbociclib does **not** produce a clean arrest in MB culture. A residual 16–18% of cells
remain pRb+, and the arrest is reversible: after 72 h of treatment plus washout the pRb+
fraction rebounds to ~50%.

Dose is the relevant axis, not a resistant subclone — see
[cdk46i_resistance_research](cdk46i_resistance_research.md).

**Open:** cell line, N, and the gate threshold for "pRb+" are not recorded. Needed before
this can go in a manuscript.
```

Note what this buys immediately. The `TBD` fields are not a defect of the format — they are
the format doing its job, making a four-item gap visible in a number that is currently
load-bearing and looks fully specified when written as a prose string.

For `LIT-*` records the frontmatter carries `doi`, `pmid`, `bibtex_key`, and `local_pdf`
(a path under the gitignored `papers/`), so the repo stays self-sufficient without
redistributing copyrighted PDFs.

### Layer 2 — Targets (`data/validation_targets.json`, upgraded)

Keep the file and its structure — it is good, and the cell-type separation in it is exactly
right. Change three things.

**Replace prose `source` with evidence IDs plus an explicit derivation.** Before:

```json
"source": "Chahin §F1 microscopy (raw 13.7% renormalized 15.7% ÷0.874). Reliable S-duration
           anchor is cumulative-BrdU (Ts~3h), NOT the count gate."
```

After:

```json
"evidence": ["EXP-0004"],
"derivation": "raw S 13.7% ÷ 0.874 (unclassified-fraction renormalization)",
"notes": "Reliable S-duration anchor is cumulative BrdU (Ts~3h), not the count gate."
```

The observation (13.7%) now lives in EXP-0004 and appears exactly once. The target
(15.7) shows its arithmetic. If the renormalization convention changes, one edit.

**Promote `OPEN_ISSUES` out of `_meta` into a real register** (`data/conflicts.json`, or a
top-level key). The nine entries there are among the most valuable content in the repo — a
list of every place the evidence genuinely disagrees with itself — and burying them in a
metadata block of a file nothing reads is how they will rot. Each conflict gets an ID, the
two-or-more evidence IDs in tension, the currently adopted resolution, and what
observation would settle it. The period conflict (Nakashima ~22–23 h vs Contestabile
16.25 h) is a live example: it is not a mistake to be fixed, it is a state to be tracked.

**Give every target a `status` that the harness respects,** which it can once it reads the
file: `active` (scored), `harness_mismatch` (recorded but the check is a known-coarse
proxy), `documented_nonformal` (real data, deliberately unscored), `proposed`.

### Layer 3 — The harness reads Layer 2

Delete the literals from `validate_v44.py`. It loads `validation_targets.json`, iterates
`status == "active"`, and scores. Writing the value in a second place stops being possible,
so §1.1's entire class of bug disappears rather than being audited for.

This also makes the harness output self-describing: a failure can print the target ID, the
evidence IDs behind it, and any open conflict touching it, instead of just a name and two
numbers.

### The linter (`tools/check_evidence.py`)

Small script, run in CI or as a pre-commit hook. Enforces the invariants the manual
three-source audit was doing by hand:

1. Every `evidence:` ID referenced by a target exists.
2. Every evidence record is referenced by at least one target, conflict, or doc — or is
   explicitly marked `unused: true` with a reason.
3. Every `LIT-*` has a DOI or PMID, and a matching key in `refs/refs.bib`.
4. Every target has `cell_type`, non-empty `evidence`, and — if `status: active` — a
   `tolerance`.
5. **Cell-type coherence:** a target scored on condition `GNP+*` may not cite evidence whose
   `cell_type` is `MB` without an explicit `cross_type_justification`. This single rule
   catches the `ezh2_palbo_drop` bug (OPEN_ISSUES #5) automatically and forever.
6. Warn on author-year strings in `docs/` that have no corresponding `LIT-*` record.

Rule 5 is the one worth building first. The GNP/MB conflation is the error this project is
most prone to, because the two cell types share most of their biology and every quantity has
a version in each.

### Layer 4 — Runs (lighter touch)

Lower priority, but: each simulation result JSON should carry a `_provenance` block with
git SHA, model builder version, parameter set, and the SHA-256 of `validation_targets.json`
at run time. Twelve lines in a shared helper. It makes "which run produced Fig 5C" answerable,
and it means a stale figure can be *detected* rather than assumed current.

### Docs reorganization

Split the flat pile by lifetime, which is the property that actually distinguishes them:

- `docs/reference/` — durable, maintained, expected to be correct now
  (`DATA_PROVENANCE`, `PARAMETERIZATION`, `v44_MODEL_DESCRIPTION`, `validation_targets`)
- `docs/notes/YYYY-MM-DD-slug.md` — dated exploration, append-only, **never edited after
  the fact**, carrying a header saying so

A note that produces a durable conclusion doesn't get edited; its conclusion gets promoted
into `reference/` or minted as an evidence record, and the note stays as the record of how
you got there. `docs/README.md` indexes both. This is the cheapest of all the changes and
removes the "is this still true?" question that currently attaches to all 45 files.

---

## 3. Suggested order

Sequenced by risk reduction per hour, not by architectural tidiness.

| # | Step | Why first | Rough effort |
|---|---|---|---|
| 1 | Mint `EXP-*` records for the unpublished JP measurements | These exist **only in chat transcripts**. Everything else in this document is recoverable from files; this is not. | 1–2 h |
| 2 | Point `validate_v44.py` at `validation_targets.json` | Kills the drift class outright. | 2–3 h |
| 3 | Convert target `source` strings → `evidence` IDs | Makes "what depends on Fig 4F" answerable. | 3–4 h |
| 4 | `refs/refs.bib` + `LIT-*` for the ~15 load-bearing citations | Needed for the manuscript regardless. | 2–3 h |
| 5 | `tools/check_evidence.py`, starting with the cell-type rule | Automates the audit that just caught six real bugs. | 2 h |
| 6 | Promote `OPEN_ISSUES` → conflicts register | They are already written; this is a move, not new work. | 1 h |
| 7 | `docs/reference` vs `docs/notes` split + index | Cheap, mostly mechanical. | 1 h |
| 8 | Run `_provenance` blocks | Nice-to-have until manuscript figures are final. | 1 h |

Steps 1 and 2 are independent and capture most of the benefit. If only one thing happens,
it should be step 1 — a repo reorganization can be redone later, but a measurement that was
only ever spoken aloud cannot be recovered once the transcript is gone.

## 4. What this deliberately does not do

No database, no ontology, no external tool (Zotero/DVC/electronic lab notebook). Everything
stays as plain files in git, because git already provides the two things that matter most
here and that a database would take away: `git blame` on a target value, and a commit
message explaining why it changed. The 2026-07-25 audit corrections are a good example —
they are currently legible only as a commit message, and that is genuinely the right place
for them.

The format is also deliberately tolerant of incompleteness. `TBD` in an evidence record is
more useful than no record, because it converts an invisible gap into a visible one.
