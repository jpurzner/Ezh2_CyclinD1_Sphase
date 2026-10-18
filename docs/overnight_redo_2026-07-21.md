# Overnight redo — GNP mitogen ramp + distribution-population (2026-07-21)

*JP (before sleep): "redo the GNP mitogen ramp up and down experiments as well as the cell simulations with drawing
from the distributions." Both re-run on the CURRENT default model (the four baked mark mechanisms: nerfed eraser,
cooperative read-write, continuous S-dilution, H3.3 turnover), plus a bonus that folds in this session's
speed-limiter governor + poised differentiation gene.*

## 1. Mitogen ramp up/down + fluctuation — `fig_memory_ramp_fluctuation.png`
Re-ran `sim_memory_ramp_fluctuation.py` on the current model.
- **CyclinD1 / proliferation hysteresis** under a ramp: the down-branch differs from the up-branch — division
  persists to lower Shh on the way down (memory). The loop needs the mark (mark-off control has none).
- **H3K27me3 is the deepest-buffered slow variable** (fluctuation CV ≈ 0.02–0.05 vs CyclinD1 ≈ 0.69–0.75); EZH2 is
  intermediate (CV ≈ 0.20). The mark low-pass-filters mitogen fluctuations; CyclinD1 tracks.
- (The mark is still saturated ~0.98, so its hysteresis/CV is small — the known de-saturation limitation.)

## 2. Distribution-population dose-response — `fig_population_mitogen.png`
Re-ran `sim_population_mitogen.py` (120 cells, calibrated per-cell distributions: CyclinD1 CV 0.77, EZH2 CV 0.53,
birth-p27 median 0.46) on the current model.
- The single-cell sharp threshold becomes a **graded population dose-response**.
- **EZH2⊣CyclinD1 repression shifts the whole GNP population threshold 3.37×** (50%-proliferating Shh: ON 0.101 vs
  OFF 0.030); **MB overpowers it (1.0×)** — constitutive Gli keeps MB above threshold at all Shh.
- Higher CyclinD1-setpoint draw → lower per-cell mitogen threshold (the cd_scale distribution generates the threshold
  distribution); EZH2 rises with mitogen (wide IQR), mark flat-saturated.

## 3. BONUS — ramp + population through the speed-limiter, with the differentiation tripwire — `fig_ramp_pop_diffgene.png`
New (`sim_ramp_pop_diffgene.py`), `with_diff_gene=True`, tying in the session's framework.
- **A · Physiological SHH ramp (0.2→1.0→0.2):** CyclinD1 governed by its mark; the poised differentiation gene stays
  **OFF throughout** (max Diff = 0.00). Robust — steady-state at every Shh (0.5–1.2) keeps Tc≈22 h, me3@diff≈0.85,
  gene poised. (Caveat: a ramp starting from *deep arrest* Shh=0.05 tripped it via a re-entry-from-quiescence burst —
  an over-dilution transient; excluded by starting above commitment.)
- **B · Overdrive ramp (growth-rate μ up/down):** pushing division fast collapses me3 at the diff locus → the gene
  **prematurely activates** → latches. The speed limit, in a ramp.
- **C · Population:** premature-activation fraction rises with growth drive (**5% at normal μ → 60% at 3× μ**); a tail
  trips even at normal drive (the weak-writer / high-setpoint cells).
- **D · Per-cell:** the cells whose me3@diff falls below the activation threshold fire; strong-writer cells keep me3 up
  and stay poised even while dividing fast.

## Bottom line
On the current model the two requested experiments reproduce the earlier findings (proliferation hysteresis; 3.4×
population threshold-shift; MB overpowers). The bonus shows the whole picture in one place: **CyclinD1 is governed by
its mark (sustainable proliferation through a mitogen ramp), the differentiation gene is a robustly-poised tripwire,
and premature activation is a drive-dependent, per-cell event** — exactly the speed-limiter framework.
Nothing baked beyond the already-approved defaults; the diff gene is off by default (validation 28/29).
