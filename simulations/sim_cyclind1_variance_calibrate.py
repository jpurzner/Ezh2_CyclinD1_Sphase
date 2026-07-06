"""Calibrate the CyclinD1 variance-buffering prediction to the measured CV.

Only the CV (variance/mean) is meaningful (fluorescence a.u. are arbitrary). We want the CALIBRATED
statement: given the OBSERVED single-cell CyclinD1 CV (0.70, with EZH2 present), how much wider does the
model predict the distribution becomes WITHOUT EZH2 repression?

Method (efficient): the map from CyclinD1 SETPOINT (cd_scale) to steady CyclinD1 is a deterministic,
nonlinear transfer function, measured ONCE on a grid, with repression ON (default) and OFF (EZH2i). For
any input log-normal spread sigma we then push samples cd_scale=exp(sigma*z) through the interpolated
transfer functions and read off CV_on(sigma) and CV_off(sigma) -- no re-simulation. We solve for sigma*
where CV_on = 0.70 (matches the data), and report CV_off and the buffering factor there.

Run:  ./venv/bin/python simulations/sim_cyclind1_variance_calibrate.py [--fresh]
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

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/sim_cyclind1_variance_calibrate_cache.npz'
KTL0 = 0.801
CV_DATA = 0.70                     # measured single-cell CyclinD1 CV (with EZH2 present), G0/1
T_END, SETTLE = 12000, 6000
CDSCALE = np.geomspace(0.12, 12.0, 44)    # transfer-function grid

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
            d = _R.simulate(0, T_END, int(T_END / 0.5), selections=['time', 'Cd'])
            return float(np.nanmean(d['Cd'][d['time'] >= SETTLE]))
        except Exception:
            continue
    return np.nan


if FRESH or not os.path.exists(CACHE):
    cd_on = np.full(len(CDSCALE), np.nan); cd_off = np.full(len(CDSCALE), np.nan)
    for i, s in enumerate(CDSCALE):
        cd_on[i] = cell(s, 0); cd_off[i] = cell(s, 1)
        print(f"  {i+1:2d}/{len(CDSCALE)}  cd_scale {s:6.2f}  Cd_on {cd_on[i]:7.2f}  Cd_off {cd_off[i]:7.2f}", flush=True)
    np.savez(CACHE, CDSCALE=CDSCALE, cd_on=cd_on, cd_off=cd_off)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
CDSCALE, cd_on, cd_off = z['CDSCALE'], z['cd_on'], z['cd_off']
good = np.isfinite(cd_on) & np.isfinite(cd_off)
CDSCALE, cd_on, cd_off = CDSCALE[good], cd_on[good], cd_off[good]

rng = np.random.default_rng(7)
zz = rng.normal(0, 1, 40000)
SIG = np.linspace(0.30, 1.40, 45)
cv = lambda x: np.std(x) / np.mean(x)
CV_on, CV_off = [], []
for s in SIG:
    cs = np.clip(np.exp(s * zz), CDSCALE.min(), CDSCALE.max())
    CV_on.append(cv(np.interp(cs, CDSCALE, cd_on)))
    CV_off.append(cv(np.interp(cs, CDSCALE, cd_off)))
CV_on, CV_off = np.array(CV_on), np.array(CV_off)

# sigma* where CV_on == CV_DATA
if CV_on.max() >= CV_DATA >= CV_on.min():
    sig_star = float(np.interp(CV_DATA, CV_on, SIG))
else:
    sig_star = SIG[np.argmin(np.abs(CV_on - CV_DATA))]
cs = np.clip(np.exp(sig_star * zz), CDSCALE.min(), CDSCALE.max())
cvon_star = cv(np.interp(cs, CDSCALE, cd_on)); cvoff_star = cv(np.interp(cs, CDSCALE, cd_off))
ratio_star = cvoff_star / cvon_star
print(f"\nsigma* (model CyclinD1 CV_on = {CV_DATA}) = {sig_star:.3f}")
print(f"  at sigma*:  CV_on {cvon_star:.3f}   CV_off (EZH2i) {cvoff_star:.3f}   buffering {ratio_star:.2f}x")

C_ON, C_OFF = '#7a5aa6', '#c0392b'
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

# (A) transfer function
ax[0].plot(CDSCALE, cd_off, '-o', color=C_OFF, ms=4, lw=2, label='no repression (EZH2i)')
ax[0].plot(CDSCALE, cd_on, '-o', color=C_ON, ms=4, lw=2, label='EZH2 repression ON')
ax[0].plot(CDSCALE, KTL0 * CDSCALE * (cd_off[0] / (KTL0 * CDSCALE[0])), 'k:', lw=1, alpha=0.5, label='linear ref')
ax[0].set_xscale('log'); ax[0].set_yscale('log')
ax[0].set_xlabel('CyclinD1 setpoint (cd_scale)'); ax[0].set_ylabel('steady CyclinD1')
ax[0].set_title('(A) Transfer function\n(EZH2 caps the high end)', fontweight='bold'); ax[0].legend(fontsize=8.5); ax[0].grid(alpha=0.15, which='both')

# (B) CV vs input sigma
ax[1].plot(SIG, CV_off, color=C_OFF, lw=2.4, label='without EZH2 (EZH2i)')
ax[1].plot(SIG, CV_on, color=C_ON, lw=2.4, label='with EZH2 repression')
ax[1].axhline(CV_DATA, color='k', ls='--', lw=1.3); ax[1].text(0.32, CV_DATA + 0.015, f'measured CV {CV_DATA}', fontsize=8.5)
ax[1].axvline(sig_star, color='#555', ls=':', lw=1.2)
ax[1].plot(sig_star, cvon_star, 'o', color=C_ON, ms=9, mec='k'); ax[1].plot(sig_star, cvoff_star, 'o', color=C_OFF, ms=9, mec='k')
ax[1].annotate(f'σ* = {sig_star:.2f}\nCV_off = {cvoff_star:.2f}', xy=(sig_star, cvoff_star), xytext=(sig_star + 0.06, cvoff_star + 0.08),
               fontsize=9, arrowprops=dict(arrowstyle='->'))
ax[1].set_xlabel('input CyclinD1-setpoint spread  σ (log-normal)'); ax[1].set_ylabel('CyclinD1 CV (SD/mean)')
ax[1].set_title('(B) Calibrate: σ* matches the measured\nwith-EZH2 CV → predicts the without-EZH2 CV', fontweight='bold')
ax[1].legend(fontsize=8.5, loc='upper left'); ax[1].grid(alpha=0.15)

# (C) buffering ratio vs sigma
ax[2].plot(SIG, CV_off / CV_on, color='#2c7fb8', lw=2.4)
ax[2].axvline(sig_star, color='#555', ls=':', lw=1.2)
ax[2].plot(sig_star, ratio_star, 'o', color='#2c7fb8', ms=10, mec='k')
ax[2].annotate(f'at σ*: {ratio_star:.2f}×\n(EZH2i widens CyclinD1\nCV {cvon_star:.2f} → {cvoff_star:.2f})',
               xy=(sig_star, ratio_star), xytext=(sig_star - 0.55, ratio_star + 0.02), fontsize=9,
               arrowprops=dict(arrowstyle='->'))
ax[2].set_xlabel('input setpoint spread  σ'); ax[2].set_ylabel('buffering factor  CV_off / CV_on')
ax[2].set_title('(C) Predicted buffering vs setpoint spread\n(modest, ~1.1–1.24×; weakly spread-dependent)', fontweight='bold'); ax[2].grid(alpha=0.15)

fig.suptitle(f'Calibrated CyclinD1 variance-buffering: at σ* where model reproduces the measured CV {CV_DATA}, EZH2i is predicted to widen CyclinD1 CV to {cvoff_star:.2f} ({ratio_star:.2f}×)',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_cyclind1_variance_calibrate.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_cyclind1_variance_calibrate.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_cyclind1_variance_calibrate.png')
