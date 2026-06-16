"""Proliferation-quiescence BIFURCATION via cyclin D1 / p27 heterogeneity (v44 + ensemble layer).

Spencer 2013 / Fan-Meyer 2021: at mitotic exit a population splits into immediate re-entry
(CDK2-inc, ~no G0) vs transient quiescence (CDK2-low, p27-high dwell), set by the mother's
cyclin D1 / p27 ratio crossing a sharp threshold. The deterministic v44 model gives every cell
one fate; here we add mother-set HETEROGENEITY in the cyclin D1 / p27 ratio and read out the
population bifurcation -- the transient-G0 fraction and dwell-time distribution -- and show that
SHH (mitogen) and EZH2 inhibition tune it.

Heterogeneity (paired across conditions: same draws, only SHH/EZH2i differ):
  - CyclinD1 setpoint  : k_Cd_translation = 0.801 * lognormal       [the D1 axis]
  - inherited p27      : P21_div          = lognormal(median 0.45, wide)  [the BIRTH-p27 axis]
  ratio_i = mean(Cd_i) / P21_div_i

Class per cell (settled limit cycle, p27 marker):
  arrest          : <2 divisions (deep quiescence)
  transient G0    : cycling, p27-high pre-S dwell >= 2 h   (CDK2-low branch)
  immediate re-entry : cycling, dwell < 2 h                (CDK2-inc branch)

No CDK2-reporter data yet -> distribution widths are plausible placeholders; the QUALITATIVE
bifurcation and its SHH/EZH2 tuning are the result. Run (background, ~10 min):
  ./venv/bin/python simulations/sim_g0_bifurcation.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

N = 140
T_END, N_PTS, SETTLE = 11000, 22000, 5500
DWELL_CUT = 2.0                       # h; transient-G0 vs immediate re-entry
rng = np.random.default_rng(7)

# mother-set heterogeneity, sampled ONCE (paired across conditions)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.68, N)), 0.25, 4.0)      # CyclinD1 setpoint multiplier (arrest axis)
# BIRTH p27 (P21_div): wide lognormal straddling the baked ultrasensitive threshold
# (probe_birthp27_dwell.py: p27<~0.2 -> immediate re-entry; 0.3-0.45 -> short G0 5-6h; >=0.6 -> ~12h plateau).
p21_div = np.clip(np.exp(rng.normal(np.log(0.45), 0.60, N)), 0.12, 2.5)  # immediate (low) <-> prolonged G0 (high)
KTL0 = 0.801                                                         # builder default k_Cd_translation (wide-search baked)

_RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)
SEL = ["time", "MPF", "P21", "aRc", "Dna", "Cd"]


def cell(ktl, p21d, shh, ezh2i, mycn, ptch1, p16=0.0, p18=0.464, ksyp21=0.002):
    _RR.reset()
    _RR['SHH'] = shh; _RR['EZH2i'] = ezh2i; _RR['MYCN_amplification'] = mycn
    _RR['Ptch1_copy_number'] = ptch1; _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl
    _RR['p16'] = p16; _RR['p18'] = p18; _RR['kSyP21'] = ksyp21   # CDK-inhibitor brake (GNP defaults / MB elevated)
    try:
        r = _RR.simulate(0, T_END, N_PTS, selections=SEL)
    except Exception:
        return dict(cls='arrest', dwell=np.nan, cd=np.nan)
    t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    cd = float(np.mean(r['Cd'][m]))
    if len(pk) < 2:
        return dict(cls='arrest', dwell=np.nan, cd=cd)
    per = float(np.mean(np.diff(tt[pk]))) / 60.0
    P21 = r['P21'][m]; aRc = r['aRc'][m]; Dna = r['Dna'][m]
    inS = (aRc > 0.05) & (Dna < 0.98); inG2 = Dna >= 0.98; preS = ~inS & ~inG2
    dwell = per * (preS & (P21 > 0.1)).mean() / ((preS | inS | inG2).mean() or 1)
    cls = 'transient_G0' if dwell >= DWELL_CUT else 'immediate'
    return dict(cls=cls, dwell=dwell, cd=cd)


CONDITIONS = {
    'GNP low-SHH':  dict(shh=0.35, ezh2i=0, mycn=1.0, ptch1=1.0),
    'GNP':          dict(shh=0.50, ezh2i=0, mycn=1.0, ptch1=1.0),
    'GNP high-SHH': dict(shh=1.00, ezh2i=0, mycn=1.0, ptch1=1.0),
    'GNP + EZH2i':  dict(shh=0.50, ezh2i=1, mycn=1.0, ptch1=1.0),
    'MB':           dict(shh=0.50, ezh2i=0, mycn=2.8, ptch1=0.1, p16=0.306, p18=1.553, ksyp21=0.002),
    'MB + EZH2i':   dict(shh=0.50, ezh2i=1, mycn=2.8, ptch1=0.1, p16=0.306, p18=1.553, ksyp21=0.002),
}

results = {}
for name, cond in CONDITIONS.items():
    cells = [cell(KTL0 * cd_scale[i], p21_div[i], **cond) for i in range(N)]
    results[name] = cells
    n_arr = sum(c['cls'] == 'arrest' for c in cells)
    n_tg0 = sum(c['cls'] == 'transient_G0' for c in cells)
    n_imm = sum(c['cls'] == 'immediate' for c in cells)
    dwells = [c['dwell'] for c in cells if c['cls'] == 'transient_G0']
    print(f"{name:14s}: arrest {100*n_arr/N:4.0f}%  transient-G0 {100*n_tg0/N:4.0f}%  "
          f"immediate {100*n_imm/N:4.0f}%   (mean G0 dwell {np.mean(dwells) if dwells else 0:.1f} h)")

# ---------------- figure ----------------
order = list(CONDITIONS)
fig, ax = plt.subplots(1, 3, figsize=(17, 5))

# (a) stacked fractions
imm = [100*np.mean([c['cls'] == 'immediate' for c in results[n]]) for n in order]
tg0 = [100*np.mean([c['cls'] == 'transient_G0' for c in results[n]]) for n in order]
arr = [100*np.mean([c['cls'] == 'arrest' for c in results[n]]) for n in order]
x = np.arange(len(order))
ax[0].bar(x, imm, color='#2ecc71', edgecolor='k', label='immediate re-entry (CDK2-inc)')
ax[0].bar(x, tg0, bottom=imm, color='#f39c12', edgecolor='k', label='transient G0 (CDK2-low)')
ax[0].bar(x, arr, bottom=np.array(imm)+np.array(tg0), color='#b03a2e', edgecolor='k', label='quiescent / arrest')
ax[0].set_xticks(x); ax[0].set_xticklabels(order, rotation=30, ha='right', fontsize=8)
ax[0].set_ylabel('% of cells'); ax[0].set_title('(a) Proliferation-quiescence bifurcation\n(SHH & EZH2i tune the split)', fontweight='bold', fontsize=11)
ax[0].legend(fontsize=8, loc='lower center')

# (b) dwell-time distributions (cycling cells), key conditions
for n, c in [('GNP low-SHH', '#1b9e77'), ('GNP', '#e7298a'), ('GNP + EZH2i', '#7570b3')]:
    d = [x['dwell'] for x in results[n] if x['cls'] in ('transient_G0', 'immediate')]
    ax[1].hist(d, bins=np.arange(0, 14, 1.0), alpha=0.55, color=c, label=n, edgecolor='k')
ax[1].axvline(DWELL_CUT, color='k', ls='--', alpha=0.6)
ax[1].set_xlabel('transient-G0 dwell (h)'); ax[1].set_ylabel('# cells')
ax[1].set_title('(b) G0 dwell-time distribution shifts\nwith mitogen / EZH2i', fontweight='bold', fontsize=11)
ax[1].legend(fontsize=8)

# (c) dwell vs cyclin D1 / p27 ratio (the single ultrasensitive threshold)
for n, c in [('GNP', '#e7298a'), ('GNP + EZH2i', '#7570b3'), ('MB', '#762A83')]:
    cells = results[n]
    ratio = [cells[i]['cd'] / p21_div[i] for i in range(N) if cells[i]['cls'] != 'arrest']
    dw = [cells[i]['dwell'] for i in range(N) if cells[i]['cls'] != 'arrest']
    ax[2].scatter(ratio, dw, s=22, color=c, alpha=0.6, label=n, edgecolor='none')
ax[2].set_xlabel('cyclin D1 / p27 ratio (Cd / inherited p27)'); ax[2].set_ylabel('transient-G0 dwell (h)')
ax[2].set_title('(c) G0 dwell set by the cyclin D1/p27 ratio\n(Fan-Meyer threshold)', fontweight='bold', fontsize=11)
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.2)

plt.tight_layout()
plt.savefig('simulations/fig_v44_g0_bifurcation.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_bifurcation.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_g0_bifurcation.png")
