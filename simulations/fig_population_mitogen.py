"""Figure: many-cell POPULATION mitogen dose-response (calibrated per-cell distributions), EZH2->CyclinD1
repression ON vs OFF, GNP & MB. Reads sim_population_mitogen_cache.pkl. (JP 2026-07-20)
Run: ./venv/bin/python simulations/fig_population_mitogen.py"""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_population_mitogen_cache.pkl'), 'rb'))
SHH = R['shh']; c_on, c_off, c_cd, c_ez, c_mk = '#c0392b', '#2e86c1', '#0e6251', '#6c3483', '#1e8449'
DIV, CD, EZ, MK = 0, 1, 2, 3


def frac_prolif(A): return (A[:, :, DIV] * 24 > 0.3).mean(axis=0)

fig, ax = plt.subplots(2, 3, figsize=(16.6, 9.3))
fig.subplots_adjust(left=0.055, right=0.98, top=0.9, bottom=0.075, hspace=0.4, wspace=0.3)

# A: GNP population proliferating fraction vs mitogen, ON vs OFF (graded threshold)
a = ax[0, 0]
a.plot(SHH, frac_prolif(R['GNP_ON']) * 100, '-o', color=c_on, lw=2.2, ms=5, label='repression ON')
a.plot(SHH, frac_prolif(R['GNP_OFF']) * 100, '-o', color=c_off, lw=2.2, ms=5, label='repression OFF')
a.set_xscale('log'); a.set_xlabel('mitogen (Shh, log)'); a.set_ylabel('% of cells proliferating')
a.set_title('A  GNP population dose-response\n(graded threshold; repression shifts it right)', fontsize=10.5, fontweight='bold')
a.legend(fontsize=8.5, loc='upper left'); a.grid(alpha=0.2); a.set_ylim(-3, 103)

# B: MB population fraction
b = ax[0, 1]
b.plot(SHH, frac_prolif(R['MB_ON']) * 100, '-o', color=c_on, lw=2.2, ms=5, label='repression ON')
b.plot(SHH, frac_prolif(R['MB_OFF']) * 100, '-o', color=c_off, lw=2.2, ms=5, label='repression OFF')
b.set_xscale('log'); b.set_xlabel('mitogen (Shh, log)'); b.set_ylabel('% of cells proliferating')
b.set_title('B  MB population dose-response\n(constitutive Gli → mostly proliferating)', fontsize=10.5, fontweight='bold')
b.legend(fontsize=8.5, loc='lower right'); b.grid(alpha=0.2); b.set_ylim(-3, 103)

# C: GNP per-cell division-rate spread, ON (spaghetti + median)
c = ax[0, 2]; A = R['GNP_ON']
for i in range(min(A.shape[0], 60)):
    c.plot(SHH, A[i, :, DIV] * 24, color=c_on, lw=0.5, alpha=0.18)
c.plot(SHH, np.median(A[:, :, DIV], axis=0) * 24, color='k', lw=2.4, label='median')
c.set_xscale('log'); c.set_xlabel('mitogen (Shh, log)'); c.set_ylabel('division rate (÷/day)')
c.set_title('C  Per-cell heterogeneity (GNP, ON)\n(each line a cell; CyclinD1/p27 draws)', fontsize=10.5, fontweight='bold')
c.legend(fontsize=8.5); c.grid(alpha=0.2)

# D: readout distributions across cells at a mid mitogen (GNP ON) -- CyclinD1, EZH2, mark
d = ax[1, 0]; si = np.argmin(np.abs(SHH - 0.5)); A = R['GNP_ON']
for key, col, lab in [(CD, c_cd, 'CyclinD1'), (EZ, c_ez, 'EZH2')]:
    x = A[:, si, key]; x = x[np.isfinite(x)]
    d.hist(x / np.nanmedian(x), bins=18, alpha=0.5, color=col, label=lab, density=True)
d.set_xlabel('level / population median'); d.set_ylabel('cell density')
d.set_title(f'D  Cell-to-cell distributions at Shh={SHH[si]:.2f}\n(GNP, ON) — the sampled heterogeneity', fontsize=10.5, fontweight='bold')
d.legend(fontsize=8.5); d.grid(alpha=0.2)

# E: CyclinD1 setpoint (cd_scale) sets the cell's mitogen threshold (scatter)
e = ax[1, 1]; A = R['GNP_ON']; cds = R['cd_scale']
cell_thr = []
for i in range(A.shape[0]):
    pr = A[i, :, DIV] * 24 > 0.3
    cell_thr.append(SHH[np.argmax(pr)] if pr.any() else np.nan)
cell_thr = np.array(cell_thr)
e.scatter(cds, cell_thr, c=R['p21'], cmap='viridis', s=28, alpha=0.8)
e.set_yscale('log'); e.set_xlabel('CyclinD1 setpoint (cd_scale draw)'); e.set_ylabel("cell's mitogen threshold (log)")
e.set_title('E  Higher CyclinD1 setpoint → lower mitogen\nthreshold (color = birth-p27 draw)', fontsize=10.5, fontweight='bold')
cb = fig.colorbar(e.collections[0], ax=e, fraction=0.046); cb.set_label('birth-p27 (P21_div)', fontsize=8)
e.grid(alpha=0.2)

# F: population EZH2 & mark vs mitogen (median + IQR), ON
f = ax[1, 2]; A = R['GNP_ON']
for key, col, lab in [(EZ, c_ez, 'EZH2'), (MK, c_mk, 'H3K27me3')]:
    med = np.nanmedian(A[:, :, key], axis=0)
    q1 = np.nanpercentile(A[:, :, key], 25, axis=0); q3 = np.nanpercentile(A[:, :, key], 75, axis=0)
    norm = np.nanmax(med)
    f.plot(SHH, med / norm, '-', color=col, lw=2.2, label=lab)
    f.fill_between(SHH, q1 / norm, q3 / norm, color=col, alpha=0.18)
f.set_xscale('log'); f.set_xlabel('mitogen (Shh, log)'); f.set_ylabel('level (norm., median ± IQR)')
f.set_title('F  Population EZH2 & mark vs mitogen (GNP, ON)\n(median band = cell spread)', fontsize=10.5, fontweight='bold')
f.legend(fontsize=8.5); f.grid(alpha=0.2)

fig.suptitle('Many-cell population (calibrated CyclinD1 / EZH2 / birth-p27 distributions): the single-cell sharp threshold becomes a GRADED population dose-response; '
             'EZH2⊣CyclinD1 repression shifts the whole population threshold', fontsize=11.5, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_population_mitogen.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
