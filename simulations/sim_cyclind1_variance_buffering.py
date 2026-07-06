"""Does EZH2-mediated repression BUFFER cell-to-cell CyclinD1 variance?

Scale-free test (only the CV — variance relative to mean — is meaningful; absolute fluorescence units
are arbitrary). A population of MB cells is given a heterogeneous CyclinD1 SETPOINT drawn log-normal
(the measured shape). Each cell's steady CyclinD1 is measured WITH EZH2->H3K27me3 repression (feedback
ON, default) and WITHOUT it (feedback OFF = EZH2i / no H3K27me3 repression). We compare the resulting
CyclinD1 distributions and their CVs.

Hypothesis: a high-CyclinD1 cell makes more EZH2 -> more H3K27me3 -> more repression, clawing back its
excess. So the feedback should COMPRESS the population CV (buffer the variance), and removing it (EZH2i)
should WIDEN it back toward the input setpoint spread -- a prediction testable in single-cell CyclinD1
imaging +/- EZH2 inhibition.

Only the shape (log-normal) and the CV are used; distributions are shown normalized to their mean.
Run:  ./venv/bin/python simulations/sim_cyclind1_variance_buffering.py [--pilot] [--fresh]
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

PILOT = not ('--full' in sys.argv)
FRESH = '--fresh' in sys.argv
N = 40 if PILOT else 200
CACHE = 'simulations/sim_cyclind1_variance_buffering_cache.npz'
KTL0 = 0.801
SIG_IN = 0.633            # input CyclinD1-setpoint spread (log-normal SHAPE; only the CV ratio is used)
T_END, SETTLE = 12000, 6000

M = build_model_v44(with_ezh2=True, with_hh=True)
MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 300000)
except Exception: pass


def cell(cd_scale, ezh2i):
    _R.reset()
    for k, v in MB.items(): _R[k] = v
    _R['EZH2i'] = ezh2i; _R['k_Cd_translation'] = float(KTL0 * cd_scale)
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            d = _R.simulate(0, T_END, int(T_END / 0.5), selections=['time', 'MPF', 'Cd'])
            m = d['time'] >= SETTLE
            pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
            return len(pk) >= 2, float(np.nanmean(d['Cd'][m]))
        except Exception:
            continue
    return False, np.nan


if FRESH or not os.path.exists(CACHE):
    rng = np.random.default_rng(7)
    cd_scale = np.clip(np.exp(rng.normal(0.0, SIG_IN, N)), 0.15, 8.0)
    on = np.full(N, np.nan); off = np.full(N, np.nan); cyc_on = np.zeros(N, bool); cyc_off = np.zeros(N, bool)
    for i in range(N):
        cyc_on[i], on[i] = cell(cd_scale[i], 0)      # feedback ON (EZH2->H3K27me3 represses)
        cyc_off[i], off[i] = cell(cd_scale[i], 1)    # feedback OFF (EZH2i: no repression)
        print(f"  cell {i+1:2d}/{N}  cd_scale {cd_scale[i]:.2f}  Cd_on {on[i]:.2f} (cyc {cyc_on[i]})  Cd_off {off[i]:.2f} (cyc {cyc_off[i]})", flush=True)
    np.savez(CACHE, cd_scale=cd_scale, on=on, off=off, cyc_on=cyc_on, cyc_off=cyc_off, sig_in=SIG_IN, N=N)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
cd_scale, on, off, cyc_on, cyc_off = z['cd_scale'], z['on'], z['off'], z['cyc_on'], z['cyc_off']
ok = cyc_on & cyc_off & np.isfinite(on) & np.isfinite(off)      # cells cycling in BOTH conditions
on_o, off_o = on[ok], off[ok]
cv = lambda x: np.std(x) / np.mean(x)
cv_on, cv_off = cv(on_o), cv(off_o)
buffer = cv_off / cv_on
print(f"\nN={ok.sum()} cycling in both.  CyclinD1 CV: feedback ON {cv_on:.3f}  OFF (EZH2i) {cv_off:.3f}"
      f"   → buffering factor CV_off/CV_on = {buffer:.2f}×")

C_ON, C_OFF = '#7a5aa6', '#c0392b'
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

# (A) CyclinD1 distributions, normalized to mean (scale-free)
b = np.linspace(0, 3.0, 16)
ax[0].hist(off_o / np.mean(off_o), bins=b, density=True, color=C_OFF, alpha=0.5, edgecolor='k', label=f'no repression (EZH2i)\nCV {cv_off:.2f}')
ax[0].hist(on_o / np.mean(on_o), bins=b, density=True, color=C_ON, alpha=0.55, edgecolor='k', label=f'EZH2 repression ON\nCV {cv_on:.2f}')
ax[0].axvline(1.0, color='k', ls=':', lw=1)
ax[0].set_xlabel('CyclinD1  (÷ mean)'); ax[0].set_ylabel('density')
ax[0].set_title('(A) CyclinD1 distribution\ncompresses under EZH2 repression', fontweight='bold'); ax[0].legend(fontsize=8.5)

# (B) per-cell: feedback pulls high cells down more (the buffering mechanism)
mx = max(off_o.max(), on_o.max()) * 1.05
ax[1].plot([0, mx], [0, mx], 'k:', lw=1, label='y = x (no buffering)')
ax[1].scatter(off_o, on_o, s=45, c=cd_scale[ok], cmap='viridis', edgecolor='k', linewidth=0.4)
ax[1].set_xlabel('CyclinD1 without repression (EZH2i)'); ax[1].set_ylabel('CyclinD1 with EZH2 repression')
ax[1].set_title('(B) Per cell: high-setpoint cells are\npulled down more → variance shrinks', fontweight='bold')
ax[1].legend(fontsize=8.5, loc='upper left'); ax[1].grid(alpha=0.15)
cb = fig.colorbar(ax[1].collections[0], ax=ax[1], label='CyclinD1 setpoint (cd_scale)', fraction=0.046, pad=0.03)

# (C) CV comparison
ax[2].bar(['EZH2 repression\nON', 'no repression\n(EZH2i)'], [cv_on, cv_off], color=[C_ON, C_OFF], edgecolor='k')
ax[2].set_ylabel('CyclinD1 coefficient of variation (SD/mean)')
ax[2].set_title(f'(C) EZH2 buffers CyclinD1 variance\nCV {cv_on:.2f} → {cv_off:.2f}  ({buffer:.2f}× wider without)', fontweight='bold')
for i, v in enumerate([cv_on, cv_off]):
    ax[2].text(i, v + 0.01, f'{v:.2f}', ha='center', fontsize=11, fontweight='bold')
ax[2].grid(alpha=0.15, axis='y')

fig.suptitle(f'EZH2-mediated repression buffers cell-to-cell CyclinD1 variance (PILOT, N={ok.sum()}; scale-free, log-normal setpoint)',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_cyclind1_variance_buffering.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_cyclind1_variance_buffering.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_cyclind1_variance_buffering.png')
