"""Persist the EZH2 per-cell-cycle-phase DISTRIBUTION targets (JP 2026-10-17) for the distribution test.
PROTEIN: IF (AMOS AC_MB55_CC_Ezh2_Apr2026, MB55 EZH2-GFP). Proper 4-phase gate from per-cell features
(G0/G1 by Ki67 since pRb absent; S by PCNA replication foci; G2 by 4N DNA; aneuploid dropped). Thresholds
fit by 2-comp GMM on the DMSO vehicle control and applied to every condition (fixed boundaries). The
EZH2-low 'dim lump' (contaminant lineage, JP) is split off by a 2-comp GMM on log10(Ezh2); targets use the
BRIGHT population. Per condition x phase we store the log10 histogram (SHAPE) + CV/log10_sd/median/skew/n
(MAGNITUDE). TRANSCRIPT: scRNA per-phase FOLDS (Gershon MB & GNP) - all scRNA supports.
Writes data/if_ezh2_phase_targets.json (committed). Run with the amos python:
  /Users/jpurzner/miniconda3/bin/python3 simulations/extract_if_ezh2_targets.py
"""
import sys, os, json, numpy as np
sys.path.insert(0, '/Users/jpurzner/Dropbox/Q_research/py_projects/amos/toolkit')
import amos
from sklearn.mixture import GaussianMixture

EXP = 'AC_MB55_CC_Ezh2_Apr2026'
CONDS = ['DMSO_0_24h', 'Pal_1_24h', 'HU_10_24h', 'HU_100_24h', 'RO3306_5_24h']
LABEL = {'DMSO_0_24h': 'DMSO', 'Pal_1_24h': 'CDK4/6i', 'HU_10_24h': 'HU_10',
         'HU_100_24h': 'HU_100', 'RO3306_5_24h': 'RO3306_G2block'}
BINS = np.linspace(2.3, 5.0, 46)


def gmm_boundary(x):
    """2-comp GMM on 1-D x; return the crossing (posterior-0.5 boundary between the two means)."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    g = GaussianMixture(2, random_state=0, n_init=3).fit(x.reshape(-1, 1))
    mu = g.means_.ravel(); lo, hi = sorted(mu)
    grid = np.linspace(lo, hi, 2000).reshape(-1, 1)
    post = g.predict_proba(grid)[:, np.argmax(mu)]            # P(high component)
    return float(grid[np.argmin(np.abs(post - 0.5)), 0])


cells = {c: amos.imaging_cells(EXP, condition=c, max_cells=50000) for c in CONDS}
d0 = cells['DMSO_0_24h']
thr_ki = gmm_boundary(np.log10(d0['Ki67'].clip(1)))           # G0 (low) vs G1+ (high) Ki67
thr_foci = gmm_boundary(d0['PCNA_foci_fraction'])             # S = high replication-foci fraction
thr_ez = gmm_boundary(np.log10(d0['Ezh2'].clip(1)))           # dim lump vs bright
print(f"thresholds (DMSO-fit): log10Ki67 {thr_ki:.3f} | PCNA_foci_frac {thr_foci:.3f} | log10Ezh2(dim/bright) {thr_ez:.3f}")


def phase_of(df):
    # S by DNA-intermediate (cell_state), NOT PCNA_foci_fraction (over-calls: foci punctate in all cycling
    # cells). G2 by 4N DNA. G0/G1 split the 2N by Ki67 (pRb absent in this experiment). aneuploid dropped.
    cs = df['cell_state'].values
    g2 = cs == 'G2'; s = cs == 'S'
    ki_hi = np.log10(df['Ki67'].clip(1)).values > thr_ki
    ph = np.where(g2, 'G2', np.where(s, 'S', np.where(ki_hi, 'G1', 'G0')))
    ph[cs == 'aneuploid'] = 'aneuploid'
    return ph


out = {'experiment': EXP, 'thresholds': {'log10Ki67': thr_ki, 'pcna_foci_frac': thr_foci, 'log10Ezh2_dimbright': thr_ez},
       'bins_log10': BINS.round(4).tolist(), 'protein': {}}
for c in CONDS:
    df = cells[c]; ph = phase_of(df)
    lez = np.log10(df['Ezh2'].clip(1)).values
    # FULL population per phase (the S-phase boost is G0-low -> cycling-high; the EZH2-low cells are the
    # REAL G0-low fraction, NOT excluded). dim_fraction kept only to characterize how G0-enriched the low mode is.
    cyc = ph != 'aneuploid'
    rec = {'n_total': int(len(df)), 'dim_fraction': float(np.mean((lez <= thr_ez) & cyc)), 'phases': {}}
    for p in ['G0', 'G1', 'S', 'G2']:
        m = (ph == p)
        x = lez[m]; x = x[np.isfinite(x)]
        if len(x) < 20:
            rec['phases'][p] = {'n': int(len(x))}; continue
        h, _ = np.histogram(x, bins=BINS, density=True)
        rec['phases'][p] = {'n': int(len(x)), 'median_log10': float(np.median(x)), 'log10_sd': float(x.std(ddof=1)),
                            'cv': float((10 ** x).std(ddof=1) / (10 ** x).mean()), 'skew': float(((x - x.mean()) ** 3).mean() / x.std() ** 3),
                            'dim_frac_in_phase': float(np.mean(x <= thr_ez)), 'hist_density': h.round(5).tolist()}
    out['protein'][LABEL[c]] = rec
    fr = {p: round(float(np.mean(ph == p) / max(cyc.mean(), 1e-9)), 3) for p in ['G0', 'G1', 'S', 'G2']}
    print(f"{LABEL[c]:16s} n={len(df):6d} dim={rec['dim_fraction']:.2f} phase-frac {fr}")

# transcript per-phase folds (scRNA)
out['transcript_fold'] = {}
for ds in ['Gershon_MB', 'Gershon_GNP']:
    r = amos.sc_gene_by_phase_fold('Ezh2', ds)
    out['transcript_fold'][ds] = {'S_vs_G0': r['fold_S_vs_G0'], 'G2M_vs_G0': r['fold_G2M_vs_G0']}
print("transcript folds:", out['transcript_fold'])

# per-phase PROTEIN medians -> fold-over-cycle (bright, DMSO) for the trajectory target
dm = out['protein']['DMSO']['phases']
if all(k in dm and 'median_log10' in dm[k] for k in ['G0', 'S', 'G2']):
    g0 = dm['G0']['median_log10']
    out['protein_fold_DMSO'] = {'S_vs_G0': float(10 ** (dm['S']['median_log10'] - g0)),
                                'G2_vs_G0': float(10 ** (dm['G2']['median_log10'] - g0)),
                                'G1_vs_G0': float(10 ** (dm['G1']['median_log10'] - g0)) if 'median_log10' in dm.get('G1', {}) else None}
    print("protein fold-over-cycle (DMSO bright):", out['protein_fold_DMSO'])

OUT = '/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/data/if_ezh2_phase_targets.json'
os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(out, open(OUT, 'w'), indent=1)
print('saved', OUT)
