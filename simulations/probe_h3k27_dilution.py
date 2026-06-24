"""Sawtooth sanity check (design spec sec 8.1) for the replicative-dilution H3K27me3 module.
Does Mk halve at early-S and recover (one-division margin)? What M_ss / CyclinD1 / divisions emerge
for GNP vs MB? This tells us the operating regime before any fold-recalibration.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_dilution=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
print("built OK with_h3k27_dilution=True")

CONDS = {'GNP': dict(SHH=0.5, MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002),
         'MB': dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)}

fig, ax = plt.subplots(2, 2, figsize=(15, 9))
for col, (name, c) in enumerate(CONDS.items()):
    rr.reset()
    for k, v in c.items():
        rr[k] = v
    rr['EZH2i'] = 0; rr['HHi'] = 0; rr['Mk'] = 0.2     # start derepressed (cycling)
    try:
        r = rr.simulate(0, 16000, 32000, selections=['time', 'Mk', 'Cd', 'Dna', 'MPF', 'EZH2'])
    except Exception as e:
        print(f'{name}: integration FAILED: {e}'); continue
    t = r['time'] / 60.0; m = r['time'] >= 9000
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=400)
    mk = r['Mk'][m]
    print(f"{name}: Mk mean {mk.mean():.3f} range [{mk.min():.3f},{mk.max():.3f}] | "
          f"Cd mean {r['Cd'][m].mean():.2f} | EZH2 mean {r['EZH2'][m].mean():.2f} | divisions {len(pk)} "
          f"(period {np.mean(np.diff(t[m][pk]))*1 if len(pk)>1 else float('nan'):.1f}h)")
    # plot last ~5 cycles
    w = t >= t[-1] - 120
    ax[0, col].plot(t[w], r['Mk'][w], '-', color='#7d3c98', lw=1.5, label='Mk (H3K27me3)')
    ax[0, col].axhline(rr['K_mk'], color='r', ls='--', alpha=0.6, label=f"K={rr['K_mk']:.2f}")
    ax[0, col].plot(t[w], r['Dna'][w], '-', color='#aaaaaa', lw=0.8, alpha=0.7, label='Dna (S-phase)')
    ax[0, col].set_title(f'{name}: H3K27me3 sawtooth', fontweight='bold'); ax[0, col].legend(fontsize=8); ax[0, col].set_ylabel('Mk / Dna')
    ax[1, col].plot(t[w], r['Cd'][w], '-', color='#27ae60', lw=1.5, label='CyclinD1')
    ax[1, col].plot(t[w], r['MPF'][w], '-', color='#c0392b', lw=0.8, alpha=0.7, label='MPF (division)')
    ax[1, col].set_title(f'{name}: CyclinD1 & cycle', fontweight='bold'); ax[1, col].legend(fontsize=8)
    ax[1, col].set_xlabel('time (h)'); ax[1, col].set_ylabel('CyclinD1 / MPF')

plt.tight_layout(); plt.savefig('simulations/probe_h3k27_dilution.png', dpi=150, bbox_inches='tight'); plt.close()
print('Saved: probe_h3k27_dilution.png')
