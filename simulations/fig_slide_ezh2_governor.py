"""ONE PowerPoint-slide-sized (16:9) figure for a talk (JP 2026-09-21). Baked model, three beats:

  1  The model cycles, and cycling is Hedgehog-dependent  (GNP: CyclinB oscillates at SHH0.5, arrests at SHH0)
  2  Higher CDK-inhibitor expression gives MB a transient G0  (MB: 1x CDKI cycles; ~6x CDKI = long transient-G0 dwells)
  3  EZH2->H3K27me3 feedback on CyclinD1+CDK6 raises the mitogen threshold  (feedback ON vs OFF: divisions vs SHH)

(The S-phase-duration scaling beat is dropped for now -- S-phase is not the cycle-limiting step in the baked/EXCESS
regime, so the threshold does not scale with it. JP: revisit on the next sweep.)

Run:  ./venv/bin/python simulations/fig_slide_ezh2_governor.py    (~12 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

# ---- palette (talk-slide) ----
C_CYC = '#1f6f4a'   # cycling / Hh-on / feedback-off proliferates
C_ARR = '#9aa0a6'   # arrested / off
C_MB  = '#b0402f'   # MB transient G0
C_EZ  = '#6a2c8f'   # EZH2 feedback ON
plt.rcParams.update({'font.size': 12, 'axes.titlesize': 13, 'axes.labelsize': 12,
                     'axes.titleweight': 'bold', 'legend.fontsize': 10.5,
                     'xtick.labelsize': 10.5, 'ytick.labelsize': 10.5,
                     'font.family': 'DejaVu Sans'})


def new_rr(params=None):
    m = build_model_v44(params=params) if params else build_model_v44()
    rr = te.loada(m)
    rr.integrator.setValue('relative_tolerance', 1e-7)
    rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr


def set_gnp(rr, shh):
    rr['SHH'] = shh; rr['Ptch1_copy_number'] = 1.0
    rr['p18'] = 0.464; rr['p19'] = 0.36


def set_mb(rr, cdki=1.0):
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    ks = rr['kSyP21']
    rr['p18'] = 1.73 * cdki; rr['p19'] = 0.58 * cdki; rr['kSyP21'] = ks * cdki


def divisions(t, dna):
    return t[np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]]


fig = plt.figure(figsize=(13.33, 7.5))
gs = GridSpec(1, 3, figure=fig, wspace=0.30, left=0.055, right=0.985, top=0.80, bottom=0.135)

# ============================================================ (1) cycling is Hh-dependent
print('beat 1 ...', flush=True)
axA = fig.add_subplot(gs[0, 0])
rrA = new_rr()
for shh, col, lab in [(0.5, C_CYC, 'Hh ON (Shh 0.5)'), (0.0, C_ARR, 'Hh OFF (Shh 0)')]:
    set_gnp(rrA, shh); rrA.reset()
    r = rrA.simulate(0, 9000, 9000, selections=['time', 'Cb'])
    axA.plot(r['time'] / 60, r['Cb'], color=col, lw=1.9, label=lab)
axA.set_xlabel('time (h)'); axA.set_ylabel('CyclinB–CDK1 (mitotic driver)')
axA.set_xlim(0, 150); axA.legend(loc='upper right', framealpha=0.95)
axA.set_title('1  The model cycles —\nand cycling needs Hedgehog', loc='left')
axA.text(0.02, 0.94, 'GNP', transform=axA.transAxes, fontsize=11, fontweight='bold',
         color='#444', va='top')

# ============================================================ (2) MB transient G0 from higher CDKI
# Read out CyclinB-CDK1 (mitotic spikes): each spike = a division; the flat stretches between spikes are
# the transient-G0 dwells. CyclinD1 stays high throughout (the "massively expressed" point), so it does
# NOT show the dwells -- the mitotic driver does.
print('beat 2 ...', flush=True)
axB = fig.add_subplot(gs[0, 1])
rrB = new_rr()
Tb = 21000  # 350 h window
mb_div = {}
for cdki, col, lab in [(1.0, C_CYC, '1× CDK-inhibitor  (24 h cycle)'),
                       (6.0, C_MB, '6× CDK-inhibitor  (transient G0, ~71 h)')]:
    set_mb(rrB, cdki); rrB.reset()
    r = rrB.simulate(0, Tb, 16000, selections=['time', 'Cb', 'Dna'])
    axB.plot(r['time'] / 60, r['Cb'], color=col, lw=1.3, label=lab)
    mb_div[cdki] = divisions(r['time'], r['Dna']) / 60
axB.set_xlabel('time (h)'); axB.set_ylabel('CyclinB–CDK1  (each spike = a division)')
axB.set_xlim(0, 350); axB.set_ylim(0, 3.05); axB.legend(loc='upper right', framealpha=0.95)
axB.set_title('2  Higher CDK-inhibitors →\ntransient G0 in MB', loc='left')
axB.text(0.02, 0.94, 'MB', transform=axB.transAxes, fontsize=11, fontweight='bold',
         color='#444', va='top')
# double-arrow spanning one 6× G0 dwell (between the 1st and 2nd MB divisions)
d6 = mb_div[6.0]
if len(d6) >= 2:
    x0, x1 = d6[0], d6[1]; yb = 2.55
    axB.annotate('', xy=(x1, yb), xytext=(x0, yb),
                 arrowprops=dict(arrowstyle='<->', color=C_MB, lw=1.6))
    axB.text((x0 + x1) / 2, yb + 0.06, 'transient G0 dwell', ha='center', va='bottom',
             fontsize=9.5, color=C_MB, style='italic')

# ============================================================ (3) EZH2 feedback raises mitogen threshold
print('beat 3 ...', flush=True)
axC = fig.add_subplot(gs[0, 2])
shhs = np.array([0.0, 0.02, 0.035, 0.05, 0.08, 0.15, 0.30, 0.50])
T3 = 40000
for params, col, lab in [(None, C_EZ, 'EZH2 feedback ON'),
                         ({'f0_prc2': 1.0, 'f0_cdk6': 1.0}, C_ARR, 'EZH2 feedback OFF')]:
    rrC = new_rr(params)
    divs = []
    for shh in shhs:
        set_gnp(rrC, shh); rrC.reset()
        r = rrC.simulate(0, T3, T3 // 10, selections=['time', 'Dna'])
        divs.append(len(divisions(r['time'], r['Dna'])))
    axC.plot(shhs, divs, 'o-', color=col, lw=2.0, ms=6.5, label=lab)
    print('  ', lab, divs, flush=True)
axC.set_xlabel('Hedgehog input (Shh)'); axC.set_ylabel('divisions in ~28 days')
axC.legend(loc='center right', framealpha=0.95)
axC.set_title('3  EZH2 feedback raises\nthe mitogen threshold', loc='left')
axC.axvspan(0.0, 0.045, color=C_EZ, alpha=0.07)
axC.annotate('mark represses CyclinD1+CDK6\n→ Hh threshold to proliferate',
             xy=(0.045, 3), xytext=(0.22, 10), fontsize=9.5, color=C_EZ, ha='center',
             arrowprops=dict(arrowstyle='->', color=C_EZ))

fig.suptitle('An EZH2→H3K27me3 governor on CyclinD1/CDK6 sets a Hedgehog threshold for cell-cycle entry',
             fontsize=14.5, fontweight='bold', x=0.055, ha='left', y=0.955)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_slide_ezh2_governor')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight')
fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
