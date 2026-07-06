"""Distribution model (focused population-of-models) — fractional EZH2i rescue.

Each virtual cell draws a small set of DECISIVE parameters from distributions and is run through the
v44 engine. Goal: reproduce that EZH2i rescues a SUBSET of vismo-arrested MB cells, and place the
outcome on a phase diagram whose axes map onto the single-cell data (CyclinD1, p27, EZH2).

DRAWS (focused; widths data-grounded where available):
  cd_scale  ~ LogNormal(0, 0.633)   -> k_Cd_translation = 0.801 * cd_scale
              *** DATA-GROUNDED: measured CyclinD1 (GFP) G0/1 sdlog=0.633, CV=0.70 (TMP=0, MB neurons) ***
  p21_div   ~ LogNormal(ln 0.85, 0.30) -> P21_div  (birth p27; placeholder CV~0.3, lit range; awaiting IF)
  kez_scale ~ LogNormal(0, 0.30)    -> kEZE2f = 0.022 * kez_scale
              (EZH2 responsiveness; placeholder -- to be replaced by the mCherry-Ezh2 / TMP series)

Per cell, 3 short runs in the MB context (MYCN-amp, Ptch1-LOH, MB CKI tone):
  cyc   : baseline, no vismo (HHi=0)      -> is it cycling? + its CyclinD1 setpoint (for the dist. check)
  arr   : + vismodegib (HHi=0.9)          -> arrested?
  resc  : + vismodegib + EZH2i            -> re-enters? (rescued = arr arrested AND resc cycling)

OUTPUTS:
  (A) model CyclinD1 setpoint distribution vs the measured log-normal (CV check).
  (B) rescue phase diagram in (CyclinD1 setpoint x birth p27), coloured by outcome + the commitment line.
  (C) rescued-fraction summary.

Run:  ./venv/bin/python simulations/sim_distribution_rescue.py [--pilot] [--fresh]
Cache: sim_distribution_rescue_cache.npz  (pilot N=24 ~ a few min; full N=200 later).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

PILOT = '--pilot' in sys.argv or True   # default to pilot for now
FRESH = '--fresh' in sys.argv
N = 24 if PILOT else 200
CACHE = 'simulations/sim_distribution_rescue_cache.npz'
KTL0, KEZ0 = 0.801, 0.022
CD_SDLOG_DATA = 0.633                    # measured CyclinD1 G0/1 sdlog (CV 0.70)
T_END, SETTLE = 12000, 6000

M = build_model_v44(with_ezh2=True, with_hh=True)
MB = dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 300000)
except Exception: pass


def cell(ktl, p21d, kez, hhi, ezh2i):
    _R.reset()
    for k, v in MB.items(): _R[k] = v
    _R['HHi'] = hhi; _R['EZH2i'] = ezh2i
    _R['k_Cd_translation'] = float(ktl); _R['P21_div'] = float(p21d); _R['kEZE2f'] = float(kez)
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            d = _R.simulate(0, T_END, int(T_END / 0.5), selections=['time', 'MPF', 'Cd', 'EZH2'])
            m = d['time'] >= SETTLE
            pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
            hrs = (d['time'][m][-1] - d['time'][m][0]) / 60.0
            return dict(cyc=len(pk) >= 2, rate=len(pk) / hrs * 168.0,
                        cd=float(np.nanmean(d['Cd'][m])), ez=float(np.nanmean(d['EZH2'][m])))
        except Exception:
            continue
    return dict(cyc=False, rate=np.nan, cd=np.nan, ez=np.nan)


if FRESH or not os.path.exists(CACHE):
    rng = np.random.default_rng(7)
    cd_scale = np.clip(np.exp(rng.normal(0.0, CD_SDLOG_DATA, N)), 0.2, 6.0)     # DATA
    p21_div = np.clip(np.exp(rng.normal(np.log(0.85), 0.30, N)), 0.4, 2.5)       # placeholder
    kez_scale = np.clip(np.exp(rng.normal(0.0, 0.30, N)), 0.3, 3.0)              # placeholder
    rows = []
    for i in range(N):
        ktl, kez = KTL0 * cd_scale[i], KEZ0 * kez_scale[i]
        cyc = cell(ktl, p21_div[i], kez, hhi=0.0, ezh2i=0)     # baseline, no vismo
        arr = cell(ktl, p21_div[i], kez, hhi=0.9, ezh2i=0)     # + vismo
        res = cell(ktl, p21_div[i], kez, hhi=0.9, ezh2i=1)     # + vismo + EZH2i
        rescued = (not arr['cyc']) and res['cyc']
        rows.append((cd_scale[i], p21_div[i], kez_scale[i], cyc['cyc'], cyc['cd'], cyc['ez'],
                     arr['cyc'], res['cyc'], rescued))
        print(f"  cell {i+1:2d}/{N}  cd_scale {cd_scale[i]:.2f} p27 {p21_div[i]:.2f}  base-cyc {cyc['cyc']}  vismo-arr {not arr['cyc']}  +EZH2i-cyc {res['cyc']}  {'RESCUED' if rescued else ''}", flush=True)
    A = np.array(rows, float)
    np.savez(CACHE, A=A, cd_sdlog=CD_SDLOG_DATA, N=N)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
A = z['A']
cd_scale, p21_div, kez_scale = A[:, 0], A[:, 1], A[:, 2]
base_cyc, base_cd, base_ez = A[:, 3].astype(bool), A[:, 4], A[:, 5]
arr_cyc, res_cyc, rescued = A[:, 6].astype(bool), A[:, 7].astype(bool), A[:, 8].astype(bool)
arrested = ~arr_cyc
n_arr = arrested.sum(); n_res = rescued.sum()
print(f"\nN={len(A)}  baseline cycling {base_cyc.sum()}  vismo-arrested {n_arr}  rescued by EZH2i {n_res}"
      f"  → rescued fraction of arrested = {100*n_res/max(n_arr,1):.0f}%")

fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))

# (A) CyclinD1 setpoint distribution: model vs measured log-normal (shape/CV check)
cd = base_cd[np.isfinite(base_cd) & base_cyc]
cv_model = np.std(cd) / np.mean(cd) if len(cd) else np.nan
ax[0].hist(cd / np.median(cd), bins=12, density=True, color='#b0c4de', edgecolor='k', alpha=0.8, label=f'model (CV {cv_model:.2f})')
xx = np.linspace(0.15, 3.0, 200)
sig = 0.633; pdf = 1 / (xx * sig * np.sqrt(2 * np.pi)) * np.exp(-(np.log(xx) - (-sig**2/2))**2 / (2 * sig**2))  # median-normalized LN
ax[0].plot(xx, pdf, color='#c0392b', lw=2.3, label='measured G0/1\nLN σ=0.633 (CV 0.70)')
ax[0].set_xlabel('CyclinD1 setpoint (÷ median)'); ax[0].set_ylabel('density')
ax[0].set_title('(A) CyclinD1 distribution:\nmodel vs measured (CV check)', fontweight='bold'); ax[0].legend(fontsize=8)

# (B) rescue phase diagram in (CyclinD1 setpoint × birth p27) — mutually-exclusive vismo outcomes
for mask, c, lab, mk in [(arr_cyc, '#2ecc71', 'cycles despite vismo', 'o'),
                         (arrested & res_cyc, '#e67e22', 'rescued by EZH2i', 's'),
                         (arrested & ~res_cyc, '#b03a2e', 'stays arrested', 'X')]:
    if mask.any():
        ax[1].scatter(cd_scale[mask], p21_div[mask], c=c, label=lab, s=60, marker=mk, edgecolor='k', linewidth=0.6)
ax[1].set_xlabel('CyclinD1 setpoint  (cd_scale; DATA CV 0.70)'); ax[1].set_ylabel('birth p27  (P21_div)')
ax[1].set_title('(B) Who gets rescued?\nCyclinD1 × p27 plane', fontweight='bold'); ax[1].legend(fontsize=8, loc='upper right'); ax[1].grid(alpha=0.15)

# (C) rescued-fraction summary
ax[2].bar(['cycles\n(vismo-resistant)', 'rescued\nby EZH2i', 'stays\narrested'],
          [base_cyc.sum(), n_res, (arrested & ~rescued).sum()],
          color=['#2ecc71', '#e67e22', '#b03a2e'], edgecolor='k')
ax[2].set_ylabel('# cells'); ax[2].set_title(f'(C) Outcomes (N={len(A)})\nrescued = {100*n_res/max(n_arr,1):.0f}% of arrested', fontweight='bold')
for i, v in enumerate([base_cyc.sum(), n_res, (arrested & ~rescued).sum()]):
    ax[2].text(i, v + 0.2, str(int(v)), ha='center', fontsize=10)

fig.suptitle(f'Distribution model (PILOT, N={len(A)}) — fractional EZH2i rescue of vismo-arrested MB; CyclinD1 axis data-grounded (σ=0.633)',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_distribution_rescue.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_distribution_rescue.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_distribution_rescue.png')
