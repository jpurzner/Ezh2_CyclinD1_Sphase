"""v44 -- mitogen sensitivity expressed in Hh activity (Gli1), probing BELOW the division threshold.

Mitogen = Hh pathway activity, read out as Gli1 transcript (the natural readout). We titrate Hh
ACTIVITY from the floor up by sweeping the Smo inhibitor GDC0449 (vismodegib) from 1 (Smo fully
blocked -> Gli_act=0 -> Hh-independent CyclinD1 floor) down to 0 (full Hh), at saturating SHH=1.5, in a
GNP context. Sweeping GDC (not SHH) is what removes the basal Smo leak and takes CyclinD1 to the true
floor, so BOTH the EZH2-feedback and no-feedback cases have a sub-threshold (arrest) regime and we can
read where each crosses the division threshold.

  (A) CyclinD1 transcript & protein vs Gli1, +/- feedback -- the dose-response from the floor up.
  (B) proliferation (divisions/168h) vs Gli1, +/- feedback -- the sub-threshold regime and the
      feedback's threshold SHIFT. GNP (physiological, GDC=0) and MB (Gli1~0.19, off-scale) are marked.

Note GNP Hh activity is intrinsically low (Gli1 ~0-0.06) and MB is ~3-7x higher -- absolute Gli1 values
are model units; the point is the threshold position relative to the floor.

Run:  ./venv/bin/python simulations/fig_v44_mitogen_sensitivity.py   (~8 min cold; caches)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from simulations.validate_v44 import _new_rr, count_divisions

XSEL = ["time", "MPF", "Cd", "Cd_mRNA", "Gli1_mRNA", "Gli1", "Gli_act", "EZH2"]
T_END, N_PTS = 9000, 12000
WIN_H = (T_END - 4000) / 60.0
GDC = np.array([1.0, 0.95, 0.9, 0.84, 0.77, 0.68, 0.57, 0.44, 0.3, 0.15, 0.0])   # Hh suppression: 1=floor, 0=full Hh
_rr = _new_rr()


def run(gdc, ezh2i, shh=1.5, ptch=1.0, mycn=1.0):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        _rr.reset()
        try: _rr.integrator.setValue("maximum_num_steps", 300000)
        except Exception: pass
        _rr['SHH'] = shh; _rr['Ptch1_copy_number'] = ptch; _rr['GDC0449'] = gdc
        _rr['EZH2i'] = ezh2i; _rr['MYCN_amplification'] = mycn
        _rr.integrator.setValue("absolute_tolerance", atol)
        try:
            r = _rr.simulate(0, T_END, N_PTS, selections=XSEL)
            m = r['time'] >= 4000
            nd, _, _ = count_divisions(r)
            return {s: float(np.mean(r[s][m])) for s in XSEL[1:]} | {'div': nd / WIN_H * 168.0}
        except Exception:
            continue
    return {s: np.nan for s in XSEL[1:]} | {'div': np.nan}


def sweep(ezh2i):
    rows = [run(g, ezh2i) for g in GDC]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mitogen_sensitivity_cache.npz')
if os.path.exists(_CACHE) and '--fresh' not in sys.argv:
    z = np.load(_CACHE); FB = {k[3:]: z[k] for k in z.files if k.startswith('fb_')}
    NF = {k[3:]: z[k] for k in z.files if k.startswith('nf_')}; MB = {k[3:]: float(z[k]) for k in z.files if k.startswith('mb_')}
    print('loaded cache; --fresh to recompute')
else:
    print('sweeping GDC WITH feedback ...'); FB = sweep(0.0)
    print('sweeping GDC WITHOUT feedback ...'); NF = sweep(1.0)
    print('MB reference ...'); _mb = run(0.0, 0.0, shh=0.5, ptch=0.1, mycn=2.8); MB = _mb
    np.savez(_CACHE, **{f'fb_{k}': v for k, v in FB.items()}, **{f'nf_{k}': v for k, v in NF.items()},
             **{f'mb_{k}': v for k, v in MB.items()})

GFB, GNF = FB['Gli1_mRNA'], NF['Gli1_mRNA']
CFB, CNF = '#762A83', '#1B7837'
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

# (A) CyclinD1 transcript & protein vs Gli1
axA = axes[0]
axA.plot(GNF, NF['Cd_mRNA'], '--', color=CNF, lw=1.8, label='transcript, no feedback')
axA.plot(GFB, FB['Cd_mRNA'], '-', color=CFB, lw=2.2, label='transcript, EZH2 feedback')
axA.plot(GNF, NF['Cd'], '--o', color=CNF, lw=1.2, ms=3, alpha=0.7, label='protein, no feedback')
axA.plot(GFB, FB['Cd'], '-o', color=CFB, lw=1.2, ms=3, alpha=0.7, label='protein, EZH2 feedback')
axA.axvspan(0, GFB[GFB > 0].min() if (GFB > 0).any() else 0, color='#EEE', alpha=0.5)
axA.set_xlabel('Hh activity = Gli1 transcript (a.u.)'); axA.set_ylabel('CyclinD1 (transcript / protein, a.u.)')
axA.set_title('A  CyclinD1 vs Hh activity (Gli1)', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=7.5, loc='upper left')
axA.text(0.002, 0.2, 'floor\n(GDC=1,\nGli_act=0)', fontsize=6.8, color='#666')

# (B) proliferation vs Gli1 -> sub-threshold regime + threshold shift
axB = axes[1]
axB.plot(GNF, NF['div'], '--s', color=CNF, lw=2, ms=5, label='no feedback (EZH2i)')
axB.plot(GFB, FB['div'], '-o', color=CFB, lw=2.2, ms=5, label='EZH2 feedback')
def first_cross(G, d, level=1.0):
    a = np.where(np.asarray(d) > level)[0]
    if not len(a):
        return np.nan
    i = a[0]
    return G[0] if i == 0 else 0.5 * (G[i - 1] + G[i])
tFB, tNF = first_cross(GFB, FB['div']), first_cross(GNF, NF['div'])
for t, c in [(tNF, CNF), (tFB, CFB)]:
    if not np.isnan(t):
        axB.axvline(t, color=c, ls=':', lw=1.3, alpha=0.8)
axB.axvspan(0, np.nanmax([tFB, tNF]), color='#F2EEF6', alpha=0.5)
axB.set_xlabel('Hh activity = Gli1 transcript (a.u.)'); axB.set_ylabel('proliferation (divisions / 168 h)')
axB.set_title('B  Threshold sits AT the floor (too little room)', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=8, loc='center right')
axB.annotate(f'both thresholds (Gli1 {tNF:.3f} / {tFB:.3f})\nsit right against the floor (Gli1=0):\nno usable sub-threshold range to\ncompare feedback vs no-feedback\n-> lower basal CyclinD1 to open it up',
             xy=(max(tFB, tNF), 0.5 * np.nanmax(FB['div'])), xytext=(0.012, 0.30 * np.nanmax(FB['div'])),
             fontsize=6.8, color='#B2182B', arrowprops=dict(arrowstyle='->', color='#B2182B', lw=1))
# mark physiological GNP (full Hh) and MB
axB.annotate('GNP\n(full Hh)', xy=(GFB[-1], FB['div'][-1]), xytext=(GFB[-1] - 0.012, np.nanmax(FB['div']) * 0.8),
             fontsize=7, color='#333', ha='center', arrowprops=dict(arrowstyle='->', lw=0.8))
axB.text(GFB.max() * 0.97, np.nanmax(FB['div']) * 0.05, f'MB at Gli1~{MB["Gli1_mRNA"]:.2f}  (off-scale $\\to$)',
         fontsize=7, color='#B2182B', ha='right')

fig.suptitle('v44: mitogen sensitivity in Hh-activity (Gli1) units, probing BELOW the division threshold (Hh titrated by GDC0449 from the floor up)',
             fontsize=10.5, y=1.0)
fig.tight_layout()
fig.savefig('simulations/fig_v44_mitogen_sensitivity.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v44_mitogen_sensitivity.pdf', bbox_inches='tight')
print('wrote simulations/fig_v44_mitogen_sensitivity.{png,pdf}')
print(f'  Hh threshold to cycle (Gli1): no-feedback {tNF:.4f}, feedback {tFB:.4f}; GNP full-Hh Gli1={GFB[-1]:.3f}, MB Gli1={MB["Gli1_mRNA"]:.3f}')
print(f'  floor CyclinD1: +fb {FB["Cd"][0]:.2f} (div {FB["div"][0]:.0f}), -fb {NF["Cd"][0]:.2f} (div {NF["div"][0]:.0f})')
