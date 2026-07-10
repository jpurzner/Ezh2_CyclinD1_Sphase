"""The core consequence of the model: EZH2 is mitogen-responsive, so added Hedgehog raises both the
CyclinD1 DRIVE (via Gli1) and the H3K27me3 BRAKE (via EZH2) -> the CyclinD1 transcript is progressively
repressed, and a threshold amount of Hh is required to initiate proliferation. "Push the gas and engage
the brake."

For each Hh (SHH) level, steady state in the GNP context: Gli1 (Hh activity), EZH2, H3K27me3 mark Mk,
CyclinD1 transcript (Cd_mRNA), and the division rate WITH vs WITHOUT EZH2 repression (f0_prc2=1). The
un-repressed drive is recovered self-consistently as Cd_mRNA / R(Mk), R = f0 + (1-f0)/(1+(Mk/K)^n).

Panels:
  A  The brake tracks the gas: Gli1 (Hh activity) and EZH2 + H3K27me3 both rise with Hh.
  B  Progressive repression: CyclinD1 raw drive (steep) vs actual transcript (flattened); shaded gap
     = amount repressed by H3K27me3, growing with Hh.
  C  Surviving fraction R = transcript/drive falls with Hh — the progressive repression, quantified.
  D  Hh needed to INITIATE proliferation: division rate vs Hh, WITH vs WITHOUT EZH2 repression;
     EZH2 raises the Hh threshold (a mitogen gatekeeper).

Run:  ./venv/bin/python simulations/fig_v44_progressive_repression.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

M = build_model_v44(with_ezh2=True, with_hh=True)
_rr = te.loada(M)
F0, KMK, NMK = _rr['f0_prc2'], _rr['K_prc2'], _rr['n_prc2']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SEL = ['time', 'Cd_mRNA', 'Mk', 'EZH2', 'Gli_act', 'Gli1', 'MPF']


def run(shh, ezh2=True, settle=6000):
    rr = te.loada(M); rr.integrator.setValue('relative_tolerance', 1e-6); rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    rr['SHH'] = shh
    if not ezh2:
        rr['f0_prc2'] = 1.0
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, 11000, 44000, selections=SEL); m = r['time'] >= settle; t = r['time'][m]
            pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / (t[1] - t[0])))
            d = {k: float(r[k][m].mean()) for k in ['Cd_mRNA', 'Mk', 'EZH2', 'Gli_act', 'Gli1']}
            d['div'] = len(pk) / ((t[-1] - t[0]) / 60) * 168.0
            return d
        except Exception:
            continue
    return {k: np.nan for k in ['Cd_mRNA', 'Mk', 'EZH2', 'Gli_act', 'Gli1', 'div']}


SHH = np.round(np.concatenate([np.linspace(0.08, 0.5, 12), np.linspace(0.6, 2.5, 10)]), 3)
on = [run(s, True) for s in SHH]
off = [run(s, False) for s in SHH]
cdm = np.array([d['Cd_mRNA'] for d in on]); mk = np.array([d['Mk'] for d in on])
ez = np.array([d['EZH2'] for d in on]); gli1 = np.array([d['Gli_act'] + d['Gli1'] for d in on])
div_on = np.array([d['div'] for d in on]); div_off = np.array([d['div'] for d in off])
R = F0 + (1 - F0) / (1 + (mk / KMK) ** NMK)
drive = cdm / R
cyc = SHH >= 0.34   # cycling range, for the decomposition panels B/C
for s, r_, d1, d0 in zip(SHH, R, div_on, div_off):
    print(f"SHH={s:.2f}  R={r_:.2f}  div(EZH2)={d1:.1f}  div(noEZH2)={d0:.1f}")


def onset(shh, div):
    i = np.where(np.nan_to_num(div) > 0.5)[0]
    return float(shh[i[0]]) if len(i) else np.nan


thr_on, thr_off = onset(SHH, div_on), onset(SHH, div_off)
C_GLI, C_EZ, C_CD, C_DRIVE = '#e67e22', '#8e44ad', '#117a65', '#c0392b'
fig, ax = plt.subplots(2, 2, figsize=(13.5, 10))

# ---- A ----
a = ax[0, 0]; a2 = a.twinx()
a.plot(SHH, gli1, '-o', color=C_GLI, lw=2.4, ms=3.5, label='Gli1 (Hh activity — the gas)')
a2.plot(SHH, ez, '-^', color=C_EZ, lw=2.2, ms=3.5, label='EZH2 (the brake)')
a2.plot(SHH, mk, '-s', color='#5b2c6f', lw=1.8, ms=3, alpha=0.8, label='H3K27me3 mark')
a.set_xlabel('Hedgehog (SHH dose)'); a.set_ylabel('Gli1 (Hh activity)', color=C_GLI); a2.set_ylabel('EZH2 / H3K27me3', color=C_EZ)
a.set_title('(A) The brake is mitogen-responsive\nadded Hh raises Gli1 AND EZH2→H3K27me3', fontweight='bold', fontsize=11)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left'); a.grid(alpha=0.15)

# ---- B ----
b = ax[0, 1]
b.plot(SHH[cyc], drive[cyc], '--', color=C_DRIVE, lw=2.2, label='CyclinD1 drive (Gli+MYCN, un-repressed)')
b.plot(SHH[cyc], cdm[cyc], '-o', color=C_CD, lw=2.6, ms=3.5, label='CyclinD1 transcript (actual)')
b.fill_between(SHH[cyc], cdm[cyc], drive[cyc], color=C_EZ, alpha=0.18, label='repressed by H3K27me3')
b.set_xlabel('Hedgehog (SHH dose)'); b.set_ylabel('CyclinD1 transcript (a.u.)')
b.set_title('(B) Progressive repression of CyclinD1\nmore Hh → more drive, but a growing H3K27me3 brake', fontweight='bold', fontsize=11)
b.legend(fontsize=8, loc='upper left'); b.grid(alpha=0.15)

# ---- C ----
c = ax[1, 0]
c.plot(SHH[cyc], R[cyc], '-o', color=C_EZ, lw=2.6, ms=3.5)
c.fill_between(SHH[cyc], R[cyc], 1.0, color=C_EZ, alpha=0.12)
c.axhline(1.0, color='k', ls=':', alpha=0.4); c.text(SHH[-1], 1.0, ' no repression', fontsize=7.5, va='bottom', ha='right', color='#555')
c.set_ylim(0, 1.05)
c.set_xlabel('Hedgehog (SHH dose)'); c.set_ylabel('surviving fraction  R = transcript / drive')
c.set_title('(C) The repression is PROGRESSIVE\nfraction of CyclinD1 drive escaping H3K27me3 falls with Hh', fontweight='bold', fontsize=11)
c.annotate(f'{R[cyc][0]*100:.0f}% → {R[-1]*100:.0f}% survive', xy=(SHH[cyc][2], R[cyc][2]),
           xytext=(0.45, 0.55), textcoords='axes fraction', fontsize=9.5, color=C_EZ, fontweight='bold',
           arrowprops=dict(arrowstyle='->', color=C_EZ))
c.grid(alpha=0.15)

# ---- D ----
d = ax[1, 1]
d.plot(SHH, div_off, '--s', color=C_DRIVE, lw=2.0, ms=3.5, label='without EZH2 repression (f0=1)')
d.plot(SHH, div_on, '-o', color=C_CD, lw=2.4, ms=3.5, label='with EZH2 repression (default)')
for thr, col, lab in [(thr_off, C_DRIVE, f'{thr_off:.2f}'), (thr_on, C_CD, f'{thr_on:.2f}')]:
    if thr == thr:
        d.axvline(thr, color=col, ls=':', lw=1.6, alpha=0.8)
d.annotate('', xy=(thr_on, 0.6), xytext=(thr_off, 0.6), arrowprops=dict(arrowstyle='<->', color='k', lw=1.6))
d.text(0.62, 2.6, f'EZH2 raises the Hh threshold\nto proliferate:  ≤{thr_off:.2f} → {thr_on:.2f}', ha='left', fontsize=9, fontweight='bold')
d.set_xlabel('Hedgehog (SHH dose)'); d.set_ylabel('proliferation (divisions / 168 h)')
d.set_title('(D) Hh needed to INITIATE proliferation\nEZH2 repression is a mitogen gatekeeper (higher Hh threshold)', fontweight='bold', fontsize=11)
d.legend(fontsize=8, loc='lower right'); d.grid(alpha=0.15)

fig.suptitle('EZH2 is mitogen-responsive → progressive H3K27me3 repression of CyclinD1: adding Hh drives Gli1↑ and CyclinD1↑, but the brake grows with it and sets the Hh threshold to proliferate',
             fontsize=11.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_progressive_repression.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_progressive_repression.pdf', bbox_inches='tight')
plt.close()
print(f"\nR: {R[cyc][0]:.2f} → {R[-1]:.2f};  drive {drive[cyc][-1]/drive[cyc][0]:.1f}x vs transcript {cdm[cyc][-1]/cdm[cyc][0]:.1f}x;  Hh threshold {thr_off:.2f} (no EZH2) → {thr_on:.2f} (EZH2)")
print("Saved: fig_v44_progressive_repression.png / .pdf")
