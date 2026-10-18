"""JP's SPEED-LIMITER framework, demonstrated in the mark equation (2026-07-21). 4 panels:
 A  the mark = writing/dilution balance -> a MAX SUSTAINABLE division rate (writer saturates)
 B  CyclinD1-as-target = a GOVERNOR: with it the mark is maintained; without it -> runaway -> differentiation
 C  the DIVISION decision is BISTABLE (go/no-go), the mark is graded
 D  regime map: locked-out / sustainable proliferation / too-fast -> differentiation, with GNP & MB
Imports the model from sim_speed_limiter.py. Run: ./venv/bin/python simulations/fig_speed_limiter.py"""
import os, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sim_speed_limiter as S

HERE = os.path.dirname(os.path.abspath(__file__))
TH = S.THETA_DIFF
C_SUS, C_DIFF, C_LOCK = '#1e8449', '#c0392b', '#7f8c8d'
C_GOV, C_NOG = '#2e86c1', '#e67e22'

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11.6))
fig.subplots_adjust(left=0.07, right=0.96, top=0.87, bottom=0.07, hspace=0.33, wspace=0.26)

# ================= A: mark maintenance / speed limit =================
a = ax[0, 0]
rr = np.linspace(0, 0.20, 400); PP = np.array([S.mark_of_rate(x) for x in rr])
Tc = np.where(rr > 0, 1 / np.maximum(rr, 1e-9), np.inf)
a.fill_between(rr, TH, 1.05, color=C_SUS, alpha=0.08)
a.fill_between(rr, -0.05, TH, color=C_DIFF, alpha=0.10)
a.plot(rr, PP, color='#5b2c6f', lw=3, zorder=4)
a.axhline(TH, color=C_DIFF, ls='--', lw=1.6); a.text(0.006, TH + 0.015, 'θ_diff — below this the marked DIFFERENTIATION genes de-repress', color=C_DIFF, fontsize=8.5)
rmax = rr[PP >= TH].max()
a.axvline(rmax, color='#333', ls=':', lw=1.4); a.text(rmax + 0.003, 0.1, f'max sustainable rate\n(Tc_min≈{1/rmax:.0f} h)', fontsize=8.5, color='#333')
for r0, lab, dx in [(0.0, 'quiescent\n(baseline EZH2)', 0.004), (0.062, 'GNP\n(Tc≈16 h, pushed hard)', -0.004), (0.13, 'too fast\n→ differentiation', 0.002)]:
    a.plot(r0, S.mark_of_rate(r0), 'o', ms=11, mfc='white', mec='#5b2c6f', mew=2.4, zorder=5)
    a.annotate(lab, (r0, S.mark_of_rate(r0)), xytext=(r0 + dx, S.mark_of_rate(r0) + (0.13 if r0 < 0.1 else -0.16)),
               fontsize=8.5, fontweight='bold', ha='center', color='#5b2c6f')
a.set_xlabel('division rate  r = 1/Tc  (1/h)'); a.set_ylabel('H3K27me3 mark  P = w(r) / (δ + ln2·r)')
a.set_ylim(-0.05, 1.05); a.set_xlim(0, 0.20)
a.set_title('A · The mark is a writing/dilution balance → a SPEED LIMIT\nwriter (EZH2) saturates, so dividing too fast collapses the mark', fontsize=11, fontweight='bold')
a.grid(alpha=0.2)

# ================= B: the CyclinD1 governor =================
b = ax[0, 1]; ms = np.linspace(0.05, 1.3, 60)
g = [S.solve(m, True, P0=1.0) for m in ms]; n = [S.solve(m, False) for m in ms]
rg = np.array([d['r'] for d in g]); Pg = np.array([d['P'] for d in g])
rn = np.array([d['r'] for d in n]); Pn = np.array([d['P'] for d in n])
b.plot(ms, rg, color=C_GOV, lw=2.8, label='rate — WITH governor (mark ⊣ CyclinD1)')
b.plot(ms, rn, color=C_NOG, lw=2.8, ls='-', label='rate — WITHOUT governor')
b.set_xlabel('mitogen  (SHH → CyclinD1 drive)'); b.set_ylabel('division rate  r  (1/h)', color='#333')
b.set_ylim(0, 0.20); b.set_title('B · CyclinD1-as-target = a GOVERNOR\nwith it: rate held sustainable, mark kept up. without it: runaway', fontsize=11, fontweight='bold')
b2 = b.twinx()
b2.plot(ms, Pg, color=C_GOV, lw=1.8, ls=':', alpha=0.9); b2.plot(ms, Pn, color=C_NOG, lw=1.8, ls=':', alpha=0.9)
b2.axhline(TH, color=C_DIFF, ls='--', lw=1.2); b2.fill_between(ms, -0.05, TH, color=C_DIFF, alpha=0.06)
b2.set_ylabel('mark P  (dotted)', color='#555'); b2.set_ylim(0, 1.05)
b2.text(0.07, TH - 0.06, 'P<θ_diff → differentiation', color=C_DIFF, fontsize=8)
b.legend(fontsize=8.5, loc='upper center'); b.grid(alpha=0.2)
b.text(0.98, 0.55, 'WITHOUT the brake,\nmitogen drives division\nstraight to max → mark\ncollapses → differentiation',
       transform=b.transAxes, ha='right', fontsize=8.5, color=C_NOG, style='italic',
       bbox=dict(boxstyle='round', fc='#fef5e7', ec=C_NOG))

# ================= C: bistable division decision =================
c = ax[1, 0]
rh = np.array([S.solve(m, True, P0=1.0)['r'] for m in ms])   # high-mark history
rl = np.array([S.solve(m, True, P0=0.0)['r'] for m in ms])   # low-mark history
bist = np.abs(rh - rl) > 0.01
c.plot(ms, rh, color=C_GOV, lw=2.8, label='from HIGH-mark history → governed (slow)')
c.plot(ms, rl, color=C_DIFF, lw=2.8, ls='--', label='from LOW-mark history → runaway (fast)')
if bist.any():
    c.axvspan(ms[bist].min(), ms[bist].max(), color='#f9e79f', alpha=0.7, zorder=0)
    c.text((ms[bist].min() + ms[bist].max()) / 2, 0.15, 'BISTABLE\ngo / no-go', ha='center', fontsize=10, fontweight='bold', color='#b9770e')
c.text(0.12, 0.008, 'LOCKED OUT\n(mark holds\nCyclinD1 down)', ha='center', fontsize=8.5, color=C_LOCK, fontweight='bold')
c.text(1.05, 0.15, 'committed', ha='center', fontsize=8.5, color='#333')
c.set_xlabel('mitogen  (SHH → CyclinD1 drive)'); c.set_ylabel('division rate  r  (1/h)'); c.set_ylim(-0.005, 0.20)
c.set_title('C · The DIVISION decision is BISTABLE, the mark is graded\nsame mitogen → divide or not, depending on mark history', fontsize=11, fontweight='bold')
c.legend(fontsize=8.5, loc='center right'); c.grid(alpha=0.2)

# ================= D: regime map =================
d = ax[1, 1]
mg = np.linspace(0.05, 1.2, 90); dg = np.linspace(0.6, 4.0, 80)   # mitogen x oncogenic CyclinD1 drive
Z = np.zeros((len(dg), len(mg)))
for iy, D in enumerate(dg):
    for ix, m in enumerate(mg):
        s = S.solve(m, True, P0=1.0, drive=D)     # governed (mark-maintaining) branch = the cell's protected state
        Z[iy, ix] = 0 if s['r'] < 0.005 else (1 if s['P'] >= TH else 2)   # locked / sustainable / differentiation
from matplotlib.colors import ListedColormap
d.pcolormesh(mg, dg, Z, cmap=ListedColormap([C_LOCK, C_SUS, C_DIFF]), shading='auto', alpha=0.85, vmin=0, vmax=2)
d.plot(0.5, 1.0, '*', ms=22, mfc='white', mec='#154360', mew=2); d.text(0.5, 1.18, 'GNP', ha='center', fontsize=11, fontweight='bold', color='#154360')
d.plot(0.5, 3.0, '*', ms=22, mfc='white', mec='#641e16', mew=2); d.text(0.5, 3.2, 'MB (constitutive Gli/MYCN)', ha='center', fontsize=10, fontweight='bold', color='#641e16')
d.set_xlabel('mitogen  (SHH)'); d.set_ylabel('oncogenic CyclinD1 drive  (Gli / MYCN)')
d.set_title('D · Regime map — GNP sits at the sustainable ceiling; MB overrides\n(pushed past the governor into the low-mark / too-fast regime)', fontsize=11, fontweight='bold')
for lab, col, xy in [('locked out', C_LOCK, (0.13, 3.5)), ('sustainable\nproliferation', '#145a32', (0.62, 1.9)), ('too fast →\ndifferentiation*', '#641e16', (0.95, 3.0))]:
    d.text(*xy, lab, fontsize=9, fontweight='bold', color=col, ha='center')
d.text(0.02, -0.13, '*in MB the "too-fast/low-mark" state is where ChIP shows MB<GNP H3K27me3; MB escapes differentiation by oncogenic means, not by the governor.',
       transform=d.transAxes, fontsize=8, style='italic', color='#555')

fig.suptitle("H3K27me3 as a division SPEED-LIMITER", fontsize=15, fontweight='bold', y=0.975)
fig.text(0.5, 0.925, "the mark is a graded writing/dilution rate-sensor · CyclinD1-as-a-target is the negative-feedback GOVERNOR · the bistable step is DIVIDE-or-NOT, not the mark",
         ha='center', fontsize=10.5, color='#333')
out = os.path.join(HERE, 'fig_speed_limiter.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
