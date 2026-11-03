"""Regenerate the two IF-derived inputs that fig_ezh2_bursting.py needs, from the AMOS IF arm
(JP 2026-11-02 — closes the data-provenance gap flagged in docs/FIGURES.md: these two files previously had
no in-repo generator, so the figure could not be rebuilt from a clean checkout).

Source: AMOS imaging experiment AC_MB55_CC_Ezh2_Apr2026 (MB55 EZH2-GFP/Ki67/PCNA/DNA, 24 h treatments).
Uses the IDENTICAL phase gate as extract_if_ezh2_targets.py (S by DNA-intermediate cell_state, G2 by 4N,
2N split G0/G1 by a DMSO-fit Ki67 GMM boundary; aneuploid dropped) so the two IF exports stay consistent.

Writes (both COMMITTED, so the figure regenerates anywhere):
  data/ezh2_hist.json        {'data': {condition: {phase: [log10 EZH2 per cell]}}}  -> per-phase medians
  data/if_EZH2_target.json   {'raw_log10_sample': [...]}  pooled DMSO sample (bimodal, incl. the dim lump)

Run with the AMOS python (not the project venv):
  /Users/jpurzner/miniconda3/bin/python3 simulations/extract_if_ezh2_bursting_inputs.py
"""
import sys, os, json, numpy as np
sys.path.insert(0, '/Users/jpurzner/Dropbox/Q_research/py_projects/amos/toolkit')
import amos
from sklearn.mixture import GaussianMixture

EXP = 'AC_MB55_CC_Ezh2_Apr2026'
CONDS = ['DMSO_0_24h', 'Pal_1_24h', 'HU_10_24h', 'HU_100_24h', 'RO3306_5_24h']
PHASES = ['G0', 'G1', 'S', 'G2']
CAP_PHASE = 3000      # per condition x phase (medians are robust; keeps the committed file small)
N_POOL = 8000         # pooled DMSO sample for the panel-B histogram
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(0)


def gmm_boundary(x):
    """2-comp GMM on 1-D x; return the posterior-0.5 crossing between the two means."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    g = GaussianMixture(2, random_state=0, n_init=3).fit(x.reshape(-1, 1))
    mu = g.means_.ravel(); lo, hi = sorted(mu)
    grid = np.linspace(lo, hi, 2000).reshape(-1, 1)
    post = g.predict_proba(grid)[:, np.argmax(mu)]
    return float(grid[np.argmin(np.abs(post - 0.5)), 0])


cells = {c: amos.imaging_cells(EXP, condition=c, max_cells=50000) for c in CONDS}
d0 = cells['DMSO_0_24h']
thr_ki = gmm_boundary(np.log10(d0['Ki67'].clip(1)))
thr_ez = gmm_boundary(np.log10(d0['Ezh2'].clip(1)))
print(f"DMSO-fit thresholds: log10Ki67 {thr_ki:.3f} | log10Ezh2 dim/bright {thr_ez:.3f}")


def phase_of(df):
    cs = df['cell_state'].values
    ki_hi = np.log10(df['Ki67'].clip(1)).values > thr_ki
    ph = np.where(cs == 'G2', 'G2', np.where(cs == 'S', 'S', np.where(ki_hi, 'G1', 'G0')))
    ph[cs == 'aneuploid'] = 'aneuploid'
    return ph


def sub(x, n):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return x if len(x) <= n else rng.choice(x, n, replace=False)


hist = {'experiment': EXP, 'source': 'AMOS imaging_cells; gate matches extract_if_ezh2_targets.py',
        'thresholds': {'log10Ki67': thr_ki, 'log10Ezh2_dimbright': thr_ez}, 'data': {}}
for c in CONDS:
    df = cells[c]; ph = phase_of(df)
    lez = np.log10(df['Ezh2'].clip(1)).values
    hist['data'][c] = {}
    for p in PHASES:
        x = sub(lez[ph == p], CAP_PHASE)
        hist['data'][c][p] = np.round(x, 4).tolist()
    ns = {p: len(hist['data'][c][p]) for p in PHASES}
    meds = {p: (round(float(np.median(hist['data'][c][p])), 3) if ns[p] else None) for p in PHASES}
    print(f"  {c:14s} n/phase {ns}  median log10 {meds}")

# pooled DMSO sample (all phases, INCLUDING the dim lump — panel B shows the full bimodal distribution)
lez0 = np.log10(d0['Ezh2'].clip(1)).values
pool = sub(lez0, N_POOL)
g1 = np.array(hist['data']['DMSO_0_24h']['G1'])
tgt = {'experiment': EXP, 'raw_log10_sample': np.round(pool, 4).tolist(),
       'n_pooled_available': int(np.isfinite(lez0).sum()),
       'dim_fraction': float(np.mean(lez0 <= thr_ez)),
       'g1_median_log10': float(np.median(g1)),
       'g1_log10_sd': float(g1.std(ddof=1)),
       'g1_cv_linear': float(np.std(10 ** g1, ddof=1) / np.mean(10 ** g1))}
print(f"  pooled DMSO sample n={len(pool)} (of {tgt['n_pooled_available']}), dim frac {tgt['dim_fraction']:.2f}, "
      f"G1 linear CV {tgt['g1_cv_linear']:.2f}")

for name, obj in [('ezh2_hist.json', hist), ('if_EZH2_target.json', tgt)]:
    p = os.path.join(ROOT, 'data', name)
    json.dump(obj, open(p, 'w'))
    print(f"wrote {p}  ({os.path.getsize(p)/1024:.0f} KB)")
