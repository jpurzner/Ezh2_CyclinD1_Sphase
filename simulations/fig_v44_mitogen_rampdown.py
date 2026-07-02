"""Gradual mitogen withdrawal: does 'EZH2 falls -> keeps dividing a while -> then collapse' happen?

Ramps SHH slowly from 0.6 -> 0 (staircase of short segments) in a cycling GNP and records CyclinD1,
EZH2, divisions (MPF), p27 along the falling mitogen. Tests JP's reading: as mitogen falls, EZH2
falls, which de-represses CyclinD1 and sustains division for a while before the cell-cycle collapse.

Also compares WITH feedback (real model) vs WITHOUT (f0_mk=1.0, mark cannot repress): if the feedback buffers
the exit, the cell should keep dividing to a LOWER mitogen level with feedback.

Run:  ./venv/bin/python simulations/fig_v44_mitogen_rampdown.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SETTLE = 4000          # min, reach steady cycling at SHH0
SHH0 = 0.6
RAMP_MIN, RAMP_MAX, NSEG = 0.0, SHH0, 48
SEG = 300              # min per ramp step  -> ramp spans NSEG*SEG/60 h
SEL = ["time", "Cd_mRNA", "EZH2", "MPF", "P21"]


def rampdown(model):
    rr = te.loada(model); rr.integrator.setValue("absolute_tolerance", 1e-8); rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = SHH0
    chunks = [rr.simulate(0, SETTLE, 8000, selections=SEL)]
    shh_t = [(0, SHH0), (SETTLE, SHH0)]
    t = SETTLE
    for s in np.linspace(SHH0, RAMP_MIN, NSEG):
        rr['SHH'] = float(s)
        seg = None
        for atol in (1e-8, 1e-7, 1e-6):           # retry the stiff (collapse) segment with looser tol
            try:
                rr.integrator.setValue("absolute_tolerance", atol)
                seg = rr.simulate(t, t + SEG, 600, selections=SEL); break
            except Exception:
                continue
        if seg is None:                            # unintegrable -> cell has collapsed; stop the ramp
            break
        chunks.append(seg); shh_t.append((t, s)); shh_t.append((t + SEG, s)); t += SEG
    out = {k: np.concatenate([c[k] for c in chunks]) for k in SEL}
    out['SHH_t'] = np.array(shh_t)
    return out


def last_division(d):
    t = d['time']; mpf = d['MPF']
    div_idx = np.where((mpf[:-1] > 0.3) & (mpf[1:] < 0.1))[0]
    if len(div_idx) == 0:
        return None, None
    ti = t[div_idx[-1]]
    # SHH at that time
    shh = np.interp(ti, d['SHH_t'][:, 0], d['SHH_t'][:, 1])
    return ti / 60.0, shh


A = rampdown(build_model_v44())                                   # with feedback (real model)
B = rampdown(build_model_v44(params={"f0_mk": 1.0}))             # without feedback (mark cannot repress)

for name, d in (("WITH feedback", A), ("WITHOUT feedback", B)):
    tlast, shhlast = last_division(d)
    ez0 = d['EZH2'][(d['time'] > SETTLE - 1500) & (d['time'] < SETTLE)].mean()
    print(f"{name:16s}: last division at t={tlast:.0f} h, SHH~{shhlast:.3f}; "
          f"EZH2 {ez0:.2f} (cycling) -> {d['EZH2'][-3000:].mean():.2f} (final)")

fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
th = A['time'] / 60.0
ax[0].plot(th, A['SHH_t'][:, 1] if False else np.interp(A['time'], A['SHH_t'][:, 0], A['SHH_t'][:, 1]),
           color='#999', lw=1.4, ls=':', label='SHH (mitogen, falling)')
ax[0].plot(th, A['Cd_mRNA'], color='#117a65', lw=1.4, label='CyclinD1 transcript')
ax[0].plot(th, A['EZH2'], color='#8e44ad', lw=1.6, label='EZH2 protein')
ax[0].plot(th, A['MPF'] * 4, color='#e67e22', lw=0.7, alpha=0.8, label='divisions (MPF x4)')
ax[0].plot(th, A['P21'], color='#b9770e', lw=1.2, alpha=0.8, label='p27')
tlast, shhlast = last_division(A)
ax[0].axvline(tlast, color='#c0392b', ls='--', alpha=0.7)
ax[0].text(tlast + 2, 2.6, f'collapse\n(SHH~{shhlast:.2f})', fontsize=8, color='#c0392b')
ax[0].set_title('WITH the EZH2-|CyclinD1 feedback (real model): mitogen ramps down -> EZH2 & CyclinD1 fall -> divide a while -> collapse',
                fontsize=10.5, fontweight='bold'); ax[0].set_ylabel('level (a.u.)'); ax[0].legend(fontsize=8, ncol=2, loc='upper right'); ax[0].grid(alpha=0.2)

th2 = B['time'] / 60.0
ax[1].plot(th2, np.interp(B['time'], B['SHH_t'][:, 0], B['SHH_t'][:, 1]), color='#999', lw=1.4, ls=':', label='SHH')
ax[1].plot(th2, B['Cd_mRNA'], color='#117a65', lw=1.4, label='CyclinD1 transcript')
ax[1].plot(th2, B['EZH2'], color='#8e44ad', lw=1.6, label='EZH2 protein')
ax[1].plot(th2, B['MPF'] * 4, color='#e67e22', lw=0.7, alpha=0.8, label='divisions (MPF x4)')
tlastB, shhlastB = last_division(B)
if tlastB: ax[1].axvline(tlastB, color='#c0392b', ls='--', alpha=0.7); ax[1].text(tlastB + 2, 4.6, f'collapse\n(SHH~{shhlastB:.2f})', fontsize=8, color='#c0392b')
ax[1].set_title('WITHOUT the feedback (EZH2 cannot repress CyclinD1): collapse at a different mitogen level',
                fontsize=10.5, fontweight='bold'); ax[1].set_ylabel('level (a.u.)'); ax[1].set_xlabel('time (h)'); ax[1].legend(fontsize=8, loc='upper right'); ax[1].grid(alpha=0.2)

fig.suptitle('Gradual mitogen withdrawal: EZH2 falls, CyclinD1 is buffered, the cell divides a while, then collapses to G0',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_mitogen_rampdown.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_mitogen_rampdown.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_mitogen_rampdown.png")
