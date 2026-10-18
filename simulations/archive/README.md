# Archived figures (2026-10-17)

These figure/sim scripts were moved out of the active `simulations/` set during the figure-registry
reconciliation. They are **NOT current** — do not cite them as the model's behaviour. Kept for provenance
(reversible: they are git-tracked, just relocated).

They fall into two buckets:

**Stale narrative** — the figure's main claim presents a superseded model state as current:
- old EZH2 protein t½ ~77 h (`kDeEZ`=1.5e-4); the baked value is now ~42 h
- the 11–12 h CyclinD1 mRNA reservoir / "memory = slow CyclinD1 transcript"; memory is now the large
  CyclinD1 **reservoir draining** at the measured 3.4 h half-life
- the EZH2-concentration (division-halving) convention; the default is dilution + APC/Cdh1
- bistable / "lock-in" / latched G0; the baked G0 is monostable / reversible
- old validation scores (25/27, 28/29, 31/32) stated as current; now 30/32
- the v42 architecture (`build_model_v42_mycn`); superseded by v44

**Broken** — read cache/data files (sweep JSONs) that no longer exist and cannot be regenerated here.

Per-script reasons are in `docs/FIGURES.md` → "Archived" section. To restore one, `git mv` it back to
`simulations/` and update its narrative to the current baked model (`docs/FIGURES.md` → Current status).
