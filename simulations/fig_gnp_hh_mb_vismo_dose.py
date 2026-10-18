"""Two 3-pane dose-response strips from the baked model (JP 2026-09-21).
TOP row  — a GNP cell across Hedgehog doses (Shh): robust cycling -> near-threshold -> arrest.
BOTTOM row — an established MB cell treated with increasing vismodegib (HHi): vismo-RESISTANT cycling ->
             arrest only at high dose.
Each pane: CyclinB-CDK1 spikes (each spike = one division) over a filled CyclinD1 area (the mitogen drive).
Division RATE is ~flat while cycling; the dose sets the CyclinD1 LEVEL and the on/off threshold.

Run:  ./venv/bin/python simulations/fig_gnp_hh_mb_vismo_dose.py    (~6 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

C_DIV = '#155e46'   # CyclinB spikes (divisions)
C_CD  = '#e0902f'   # CyclinD1 drive (filled area)
C_ARR = '#b0392b'   # arrest annotation
plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12, 'axes.labelsize': 11,
                     'axes.titleweight': 'bold', 'xtick.labelsize': 10, 'ytick.labelsize': 10,
                     'font.family': 'DejaVu Sans'})

WIN = 9000   # 150 h shown per pane


def new_rr():
    rr = te.loada(build_model_v44())
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr


def run_gnp(shh):
    rr = new_rr()
    rr['SHH'] = shh; rr['Ptch1_copy_number'] = 1.0; rr['p18'] = 0.464; rr['p19'] = 0.36; rr['HHi'] = 0.0
    rr.reset()
    r = rr.simulate(0, WIN, 9000, selections=['time', 'Cb', 'Cd', 'Dna'])
    return dict(t=r['time'] / 60, Cb=r['Cb'], Cd=r['Cd'], Dna=r['Dna'])


def run_mb(hhi, T_est=8000):
    rr = new_rr()
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73; rr['p19'] = 0.58; rr['HHi'] = 0.0
    rr.reset()
    rr.simulate(0, T_est, 200, selections=['time'])   # establish the tumour untreated
    rr['HHi'] = hhi                                    # then add vismodegib
    r = rr.simulate(T_est, T_est + WIN, 9000, selections=['time', 'Cb', 'Cd', 'Dna'])
    return dict(t=(r['time'] - T_est) / 60, Cb=r['Cb'], Cd=r['Cd'], Dna=r['Dna'])


def n_div(d):
    return int(np.sum((d['Dna'][:-1] > 0.9) & (d['Dna'][1:] < 0.1)))


print('running GNP Hh doses ...', flush=True)
gnp = [(run_gnp(0.5), 'High Hh   (Shh 0.5)'),
       (run_gnp(0.05), 'Low Hh   (Shh 0.05)'),
       (run_gnp(0.0), 'No Hh   (Shh 0)')]
print('running MB vismo doses ...', flush=True)
mb = [(run_mb(0.0), 'No vismo   (HHi 0)'),
      (run_mb(0.5), 'Moderate vismo   (HHi 0.5)'),
      (run_mb(1.0), 'High vismo   (HHi 1.0)')]

fig = plt.figure(figsize=(14.5, 8.2))
gs = GridSpec(2, 3, figure=fig, hspace=0.52, wspace=0.34, left=0.06, right=0.945, top=0.86, bottom=0.085)


def draw(ax, d, title, cd_ymax):
    axr = ax.twinx()
    axr.fill_between(d['t'], 0, d['Cd'], color=C_CD, alpha=0.30, lw=0, zorder=1)
    axr.plot(d['t'], d['Cd'], color=C_CD, lw=1.3, zorder=2)
    axr.set_ylim(0, cd_ymax); axr.tick_params(axis='y', labelcolor=C_CD)
    ax.plot(d['t'], d['Cb'], color=C_DIV, lw=1.15, zorder=3)
    ax.set_ylim(0, 2.75); ax.set_xlim(0, 150)
    ax.set_zorder(axr.get_zorder() + 1); ax.patch.set_visible(False)   # spikes drawn over the CyclinD1 area
    ax.set_title(title, pad=7)
    bbox = dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.82)
    nd = n_div(d)
    if nd <= 3:
        ax.text(0.5, 0.55, 'arrested', transform=ax.transAxes, ha='center', color=C_ARR,
                fontsize=12, fontweight='bold', bbox=bbox)
    else:
        ax.text(0.03, 0.90, f'{nd} divisions', transform=ax.transAxes, ha='left', color=C_DIV,
                fontsize=9.5, fontweight='bold', bbox=bbox)
    ax.text(0.97, 0.90, f'CyclinD1 {d["Cd"][-3000:].mean():.1f}', transform=ax.transAxes, ha='right',
            color=C_CD, fontsize=9.5, fontweight='bold', bbox=bbox)
    return axr


CD_GNP, CD_MB = 3.6, 8.6
for j, (d, title) in enumerate(gnp):
    ax = fig.add_subplot(gs[0, j]); axr = draw(ax, d, title, CD_GNP)
    if j == 0: ax.set_ylabel('CyclinB–CDK1\n(each spike = division)')
    if j == 2: axr.set_ylabel('CyclinD1 (drive)', color=C_CD)
    ax.set_xlabel('time (h)')
for j, (d, title) in enumerate(mb):
    ax = fig.add_subplot(gs[1, j]); axr = draw(ax, d, title, CD_MB)
    if j == 0: ax.set_ylabel('CyclinB–CDK1\n(each spike = division)')
    if j == 2: axr.set_ylabel('CyclinD1 (drive)', color=C_CD)
    ax.set_xlabel('time (h)   (vismodegib from t = 0)')

fig.text(0.5, 0.965, 'GNP — increasing Hedgehog dose', ha='center', fontsize=14, fontweight='bold', color='#243b53')
fig.text(0.5, 0.475, 'Medulloblastoma — increasing vismodegib (established tumour, then treated)',
         ha='center', fontsize=14, fontweight='bold', color='#5b2333')

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_gnp_hh_mb_vismo_dose')
fig.savefig(out + '.png', dpi=200, bbox_inches='tight')
fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
print('GNP divisions:', [n_div(d) for d, _ in gnp], ' CyclinD1:', [round(float(d['Cd'][-3000:].mean()), 2) for d, _ in gnp])
print('MB  divisions:', [n_div(d) for d, _ in mb], ' CyclinD1:', [round(float(d['Cd'][-3000:].mean()), 2) for d, _ in mb])
