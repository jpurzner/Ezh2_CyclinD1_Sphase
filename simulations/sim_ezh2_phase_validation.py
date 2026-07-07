"""Model vs data: EZH2 protein by cell-cycle phase, across drug treatments (MB55 IF, Apr 2026).

The model's EZH2 is an E2f-gated, stable protein that INTEGRATES S-phase and is HALVED at division.
The single-cell IF data test each of those. Scale-free: only shapes / CVs / phase-ratios are used
(absolute fluorescence a.u. are arbitrary) -- EZH2 by phase is normalized to G1.

We sample an asynchronous cycling MB population over time, bin each sample by phase, and compare the
model's EZH2-by-phase to the data, for four conditions:
  DMSO        -- baseline: expect EZH2 rise G0->G1->S, peaking S (data key finding)
  HU (10 uM)  -- fork slowdown lengthens S: expect S-phase EZH2 BOOST (~1.31x; EZH2 integrates S)
  Palbociclib -- CDK4/6i: expect lower cycling EZH2 (E2f-gated synthesis off)
  RO3306      -- CDK1i, blocks mitosis: no division -> no EZH2 halving -> G2 EZH2 stays high (data:
                 removes the ~19% low-EZH2 post-mitotic G2 subpopulation)

Phase calls (model analogs of the pRb/DAPI/PCNA gates): S = active forks (aRc>0.05); G2 = 4N (Dna>0.9)
& post-S; G0/G1 = 2N & not-S, split by E2f (commitment ~ the data's pRb 'dividing vs not').

Run:  ./venv/bin/python simulations/sim_ezh2_phase_validation.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/sim_ezh2_phase_validation_cache.npz'
N = 8
KTL0 = 0.801
CD_SDLOG = 0.633
T_END, SETTLE = 16000, 6000
DRUGS = {'DMSO': {}, 'HU': {'HU': 1.0}, 'Palbociclib': {'kPhRbCd': 0.0}, 'RO3306': {'k25': 0.0}}

# --- DATA (scale-free): DMSO EZH2 by phase (÷G1) and CV; drug S-boost etc. ---
DATA_RATIO_G1 = {'G0': 21.705 / 26.251, 'G1': 1.0, 'S': 29.545 / 26.251, 'G2': 30.031 / 26.251}
DATA_CV = {'G0': 0.37, 'G1': 0.61, 'S': 0.62, 'G2': 0.46}
DATA_S_BOOST_HU = 38.784 / 29.545        # 1.31x
DATA_DRUG_G2 = {'DMSO': 1.0, 'HU': 31.49 / 30.03, 'Palbociclib': 17.84 / 30.03, 'RO3306': 36.09 / 30.03}

M = build_model_v44(with_ezh2=True, with_hh=True)
MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 400000)
except Exception: pass
PHASES = ['G0', 'G1', 'S', 'G2']


def trace(cd_scale, drug):
    _R.reset()
    for k, v in MB.items(): _R[k] = v
    _R['k_Cd_translation'] = float(KTL0 * cd_scale)
    for k, v in drug.items(): _R[k] = v
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            d = _R.simulate(0, T_END, int(T_END / 2), selections=['time', 'EZH2', 'Dna', 'aRc', 'E2f', 'MPF'])
            m = d['time'] >= SETTLE
            return d['EZH2'][m], d['Dna'][m], d['aRc'][m], d['E2f'][m]
        except Exception:
            continue
    return None


def classify(dna, arc, e2f, e2f_thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05
    inG2 = (~inS) & (dna > 0.9)
    in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'
    ph[in2N & (e2f >= e2f_thr)] = 'G1'; ph[in2N & (e2f < e2f_thr)] = 'G0'
    return ph


if FRESH or not os.path.exists(CACHE):
    rng = np.random.default_rng(7)
    cd_scale = np.clip(np.exp(rng.normal(0.0, CD_SDLOG, N)), 0.3, 4.0)
    pooled = {}     # drug -> dict(EZH2, Dna, aRc, E2f)
    for drug, dp in DRUGS.items():
        EZ, DN, AR, E2 = [], [], [], []
        for i in range(N):
            r = trace(cd_scale[i], dp)
            if r is None:
                continue
            ez, dn, ar, e2 = r
            EZ.append(ez); DN.append(dn); AR.append(ar); E2.append(e2)
        pooled[drug] = dict(EZH2=np.concatenate(EZ), Dna=np.concatenate(DN), aRc=np.concatenate(AR), E2f=np.concatenate(E2))
        print(f"  {drug:12s} samples {len(pooled[drug]['EZH2'])}", flush=True)
    # E2f threshold from DMSO 2N cells (median)
    dm = pooled['DMSO']; twoN = (dm['aRc'] <= 0.05) & (dm['Dna'] <= 0.9)
    e2f_thr = float(np.median(dm['E2f'][twoN]))
    np.savez(CACHE, e2f_thr=e2f_thr,
             **{f'{d}_{k}': pooled[d][k] for d in DRUGS for k in ['EZH2', 'Dna', 'aRc', 'E2f']})
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
e2f_thr = float(z['e2f_thr'])
pooled = {d: {k: z[f'{d}_{k}'] for k in ['EZH2', 'Dna', 'aRc', 'E2f']} for d in DRUGS}
cv = lambda x: np.std(x) / np.mean(x) if len(x) else np.nan

# per-phase EZH2 (median + CV) per drug
byphase = {}
for d in DRUGS:
    ph = classify(pooled[d]['Dna'], pooled[d]['aRc'], pooled[d]['E2f'], e2f_thr)
    byphase[d] = {p: pooled[d]['EZH2'][ph == p] for p in PHASES}

dm = byphase['DMSO']
g1med = np.median(dm['G1']) if len(dm['G1']) else np.nan
model_ratio = {p: (np.median(dm[p]) / g1med if len(dm[p]) else np.nan) for p in PHASES}
model_cv = {p: cv(dm[p]) for p in PHASES}
# HU S-boost (model)
hu_S = np.median(byphase['HU']['S']) if len(byphase['HU']['S']) else np.nan
dm_S = np.median(dm['S'])
model_S_boost_HU = hu_S / dm_S
# drug G2 EZH2 (model, ÷DMSO G2)
dmg2 = np.median(dm['G2'])
model_drug_g2 = {d: (np.median(byphase[d]['G2']) / dmg2 if len(byphase[d]['G2']) else np.nan) for d in DRUGS}
print(f"\nmodel EZH2 ÷G1: {[f'{p} {model_ratio[p]:.2f}' for p in PHASES]}")
print(f"model CV: {[f'{p} {model_cv[p]:.2f}' for p in PHASES]}")
print(f"HU S-boost: model {model_S_boost_HU:.2f}  data {DATA_S_BOOST_HU:.2f}")
print(f"drug G2 ÷DMSO: model {[(d, round(model_drug_g2[d],2)) for d in DRUGS]}")

fig, ax = plt.subplots(1, 4, figsize=(19, 4.8))
x = np.arange(len(PHASES))

# (A) EZH2 by phase, model vs data (÷G1)
ax[0].plot(x, [DATA_RATIO_G1[p] for p in PHASES], '-o', color='#c0392b', lw=2.4, ms=8, label='data (MB55 IF)')
ax[0].plot(x, [model_ratio[p] for p in PHASES], '-s', color='#2c7fb8', lw=2.4, ms=8, label='model')
ax[0].set_xticks(x); ax[0].set_xticklabels(PHASES); ax[0].axhline(1, color='k', ls=':', lw=0.8, alpha=0.5)
ax[0].set_ylabel('EZH2  (÷ G1)'); ax[0].set_title('(A) EZH2 by phase — peaks in S\nmodel vs data (scale-free)', fontweight='bold'); ax[0].legend(fontsize=9); ax[0].grid(alpha=0.15)

# (B) EZH2 CV by phase
w = 0.38
ax[1].bar(x - w/2, [DATA_CV[p] for p in PHASES], w, color='#c0392b', alpha=0.8, label='data')
ax[1].bar(x + w/2, [model_cv[p] for p in PHASES], w, color='#2c7fb8', alpha=0.8, label='model')
ax[1].set_xticks(x); ax[1].set_xticklabels(PHASES); ax[1].set_ylabel('EZH2 CV (SD/mean)')
ax[1].set_title('(B) EZH2 variability by phase\nmodel vs data', fontweight='bold'); ax[1].legend(fontsize=9); ax[1].grid(alpha=0.15, axis='y')

# (C) HU S-phase boost (EZH2 integrates S)
ax[2].bar(['data', 'model'], [DATA_S_BOOST_HU, model_S_boost_HU], color=['#c0392b', '#2c7fb8'], alpha=0.85, edgecolor='k')
ax[2].axhline(1, color='k', ls=':', lw=0.8)
for i, v in enumerate([DATA_S_BOOST_HU, model_S_boost_HU]):
    ax[2].text(i, v + 0.02, f'{v:.2f}×', ha='center', fontsize=11, fontweight='bold')
ax[2].set_ylabel('S-phase EZH2  (HU ÷ DMSO)'); ax[2].set_title('(C) HU boosts S-phase EZH2\n(EZH2 integrates S-phase duration)', fontweight='bold'); ax[2].grid(alpha=0.15, axis='y')

# (D) drug effect on G2 EZH2 (÷DMSO): HU↑, Palbo↓, RO3306↑ (no mitotic halving)
dd = list(DRUGS)
ax[3].plot(range(len(dd)), [DATA_DRUG_G2[d] for d in dd], '-o', color='#c0392b', lw=2.2, ms=8, label='data')
ax[3].plot(range(len(dd)), [model_drug_g2[d] for d in dd], '-s', color='#2c7fb8', lw=2.2, ms=8, label='model')
ax[3].axhline(1, color='k', ls=':', lw=0.8, alpha=0.5)
ax[3].set_xticks(range(len(dd))); ax[3].set_xticklabels(dd, rotation=20, ha='right', fontsize=9)
ax[3].set_ylabel('G2 EZH2  (÷ DMSO)'); ax[3].set_title('(D) Drug effects on G2 EZH2\nRO3306 (block mitosis) → no halving → ↑', fontweight='bold'); ax[3].legend(fontsize=9); ax[3].grid(alpha=0.15)

fig.suptitle('EZH2 protein by cell-cycle phase & drug — model vs MB55 single-cell data (scale-free): S-phase peak, HU S-boost, RO3306 removes mitotic halving',
             fontsize=12.5, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('simulations/sim_ezh2_phase_validation.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ezh2_phase_validation.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ezh2_phase_validation.png')
