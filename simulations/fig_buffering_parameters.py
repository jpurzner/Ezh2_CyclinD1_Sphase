"""Figure: which parameters support BUFFERING? One-at-a-time sensitivity of the EZH2/H3K27me3 buffering of
CyclinD1 (GNP), + the summary equation. Reads sim_buffering_parameters_cache.pkl. (JP 2026-07-20)
Buffering (static) = 1 - R, R = Cd_ON/Cd_OFF (repression depth). fluc_atten = mark low-pass (dynamic).
Run: ./venv/bin/python simulations/fig_buffering_parameters.py"""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_buffering_parameters_cache.pkl'), 'rb'))
SHH = R['SHH']; cd_off = R['cd_off']; si = int(np.argmin(np.abs(SHH - 0.5)))
BASE = R['BASE']
c_buf, c_flu, c_slope = '#0e6251', '#c0392b', '#6c3483'

PANELS = [('f0_prc2', 'repression FLOOR — the clean knob', True),
          ('K_prc2', 'repression SENSITIVITY (mark→K)', True),
          ('n_prc2', 'Hill steepness', False),
          ('a_rw_prc2', 'mark read-write / level', False),
          ('del_mk', 'mark turnover (timescale + level)', False),
          ('kEZbas_Cd', 'mitogen-dose EZH2', False),
          ('kEZE2f', 'cycle EZH2', False),
          ('kDeEZ', 'EZH2 turnover', False)]

fig, ax = plt.subplots(3, 3, figsize=(16.5, 12.4))
fig.subplots_adjust(left=0.05, right=0.965, top=0.92, bottom=0.06, hspace=0.42, wspace=0.42)

for idx, (pname, subt, clean) in enumerate(PANELS):
    a = ax[idx // 3, idx % 3]
    vals = sorted(R[pname]); x = np.array([v for v, m in vals])
    buf = np.array([1 - m['cd'][si] / (cd_off[si] + 1e-9) for v, m in vals])
    flu = np.array([m['fluc'] for v, m in vals])
    a.plot(x, buf, '-o', color=c_buf, lw=2.2, ms=6, label='buffering (1−R)')
    a.axvline(BASE.get(pname, np.nan), color='#888', ls=':', lw=1.3)
    a.set_xlabel(pname); a.set_ylabel('buffering  (1 − R)', color=c_buf); a.tick_params(axis='y', labelcolor=c_buf)
    a.set_ylim(-0.05, 1.0)
    if pname in ('K_prc2', 'a_rw_prc2', 'del_mk', 'kEZbas_Cd', 'kEZE2f', 'kDeEZ'):
        a.set_xscale('log')
    a2 = a.twinx()
    a2.plot(x, flu, '-s', color=c_flu, lw=1.6, ms=4, alpha=0.8, label='fluc atten')
    a2.set_ylabel('fluc atten', color=c_flu, fontsize=8); a2.tick_params(axis='y', labelcolor=c_flu, labelsize=7); a2.set_ylim(0, 0.13)
    tag = '  (monotonic ✓)' if clean else '  (BISTABLE — mark latches)'
    a.set_title(f'{chr(65+idx)}  {pname}\n{subt}{tag}', fontsize=9.5, fontweight='bold')
    a.grid(alpha=0.2)

# ---- summary panel: the equation (no overlapping inset) ----
s = ax[2, 2]; s.axis('off')
s.text(0.5, 0.99, 'SUMMARY EQUATION', ha='center', va='top', fontsize=12, fontweight='bold', transform=s.transAxes)
s.text(0.0, 0.93,
       "MARK level (cell-cycle-coupled):\n"
       "   P* = w / (del_mk + ln2/Tc)\n"
       "   w = writing (∝ EZH2·PRC2);  Tc = cycle len\n"
       " • S-phase HALVES the mark → dilution ln2/Tc\n"
       " • del_mk≈ln2/Tc at Tc≈22h; dilution wins fast Tc\n"
       " • Tc is mitogen-dependent (G1-lengthening)\n\n"
       "STATIC buffering (dose-resp; mitogen shift):\n"
       "   B = 1 − R(P*),  R = f0+(1−f0)/(1+(P*/K)ⁿ)\n"
       " • ceiling max B = 1−f0_prc2;  engages P*≳K_prc2\n"
       " • mitogen shift = 1/R (ON needs 1/R× more)\n\n"
       "DYNAMIC (fluctuation low-pass):\n"
       "   atten(ω)=ωc/√(ωc²+ω²), ωc = del_mk+ln2/Tc\n"
       "COMBINED: B(ω) ≈ [1−R(P*)]·ωc/√(ωc²+ω²)",
       ha='left', va='top', fontsize=8.4, transform=s.transAxes, family='monospace')
s.text(0.0, 0.045,
       "Buffering knobs: f0_prc2↓ & K_prc2↓ (clean); mark-level\n"
       "(a_rw↑, del_mk↓, kEZbas_Cd↑, kEZE2f↑, kDeEZ↓, Tc↑) raise\n"
       "P* but SWITCH (bistable). Tc couples it to the cell cycle.",
       ha='left', va='top', fontsize=7.6, transform=s.transAxes, style='italic', color='#555')

fig.suptitle('Parameters that support the EZH2/H3K27me3 buffering of CyclinD1 — one-at-a-time sensitivity + summary equation',
             fontsize=13, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_buffering_parameters.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
