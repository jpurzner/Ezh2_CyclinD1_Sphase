"""What happens to the FULL v44 model in the LATCH regime -- 4-panel consequences (JP 2026-07-21).
Reads sim_latch_consequences_cache.pkl. Run: ./venv/bin/python simulations/fig_latch_consequences.py"""
import os, pickle, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_latch_consequences_cache.pkl'), 'rb'))
C_BAK, C_LAT, C_ARR, C_DIV = '#2e86c1', '#c0392b', '#8e44ad', '#16a085'

# bistability test (from probe_latch_fullmodel.py; deterministic) -- Mk from HIGH vs LOW init
BIST = {  # cond: (baked_hi, baked_lo, latch_hi, latch_lo)
    'GNP+SHH\n(dividing)': (0.968, 0.968, 0.988, 0.001),
    'MB\n(dividing)':      (0.929, 0.929, 0.981, 0.003),
    'GNP+HHi\n(arrested)': (0.987, 0.987, 0.997, 0.001),
}

fig, ax = plt.subplots(2, 2, figsize=(15.5, 11.2))
fig.subplots_adjust(left=0.07, right=0.975, top=0.88, bottom=0.07, hspace=0.34, wspace=0.24)

# ---------- A: division-gating (the key result) ----------
a = ax[0, 0]; g = R['gating']
a.semilogx(g['a0'], g['arr_gap'], '-o', color=C_ARR, lw=2.6, ms=6, label='ARRESTED cell (GNP+HHi)')
a.semilogx(g['a0'], g['div_gap'], '-s', color=C_DIV, lw=2.6, ms=6, label='DIVIDING cell (GNP+SHH)')
a.axvline(0.00287, color='#ff2d55', ls='--', lw=1.6); a.text(0.00287, 0.5, ' v44 baked a0', color='#ff2d55', fontsize=8.5, rotation=90, va='center')
a.set_xlabel('nucleation  a0_prc2  (log)'); a.set_ylabel('bistability gap  |Mk(high init) − Mk(low init)|')
a.set_title('A · The latch is ARREST-GATED\nDividing cells can never latch — S-phase dilution always wins', fontsize=11, fontweight='bold')
a.legend(fontsize=9, loc='center left'); a.grid(alpha=0.25, which='both'); a.set_ylim(-0.05, 1.05)
a.text(0.97, 0.45, 'latch opens only in\nNON-dividing cells,\nbelow a0 ≈ 6×10⁻⁵', transform=a.transAxes,
       ha='right', fontsize=9, style='italic', color=C_ARR,
       bbox=dict(boxstyle='round', fc='#f5eef8', ec=C_ARR))

# ---------- B: bistability test (latch exists) ----------
b = ax[0, 1]; conds = list(BIST.keys()); x = np.arange(len(conds))
for i, c in enumerate(conds):
    bh, bl, lh, ll = BIST[c]
    b.plot([i - 0.16, i - 0.16], [bl, bh], '-', color=C_BAK, lw=2, alpha=0.4)
    b.plot(i - 0.16, bh, 'o', color=C_BAK, ms=9); b.plot(i - 0.16, bl, 'o', mfc='white', mec=C_BAK, mew=2, ms=9)
    b.plot([i + 0.16, i + 0.16], [ll, lh], '-', color=C_LAT, lw=2.5)
    b.plot(i + 0.16, lh, 'o', color=C_LAT, ms=11); b.plot(i + 0.16, ll, 'o', mfc='white', mec=C_LAT, mew=2.5, ms=11)
    if abs(lh - ll) > 0.1:
        b.annotate('', (i + 0.16, lh - 0.05), (i + 0.16, ll + 0.05), arrowprops=dict(arrowstyle='<->', color=C_LAT, lw=1.6))
b.set_xticks(x); b.set_xticklabels(conds, fontsize=9); b.set_ylabel('settled H3K27me3 mark  Mk')
b.set_title('B · The latch EXISTS (strong read-write config)\nsame condition, two heritable mark states', fontsize=11, fontweight='bold')
b.set_ylim(-0.05, 1.08); b.grid(alpha=0.25, axis='y')
b.plot([], [], 'o', color=C_BAK, label='BAKED: high & low init converge (monostable)')
b.plot([], [], 'o', color=C_LAT, label='LATCH: high vs low init SPLIT (bistable)')
b.legend(fontsize=8.5, loc='center left')

# ---------- C: validation damage ----------
c = ax[1, 0]; c.axis('off')
vb, vl = R['val_baked'], R['val_latch']
rows = [('CyclinD1 GNP (abs)', '≈1.4', vb['Cd_GNP'], vl['Cd_GNP'], vl['Cd_GNP'] > 2.5),
        ('CyclinD1 MB (abs)', '≈7', vb['Cd_MB'], vl['Cd_MB'], vl['Cd_MB'] > 12),
        ('CyclinD1 MB/GNP fold', '5.07', vb['fold'], vl['fold'], abs(vl['fold'] - 5.07) > 1.5),
        ('ChIP Mk MB/GNP (want <1)', '~0.5', vb['chip'], vl['chip'], vl['chip'] > 1.0),
        ('EZH2 MB/GNP', '2.05', vb['ez_fold'], vl['ez_fold'], abs(vl['ez_fold'] - 2.05) > 0.5),
        ('GNP divides (+SHH)', '>0', vb['div_GNP'], vl['div_GNP'], vl['div_GNP'] == 0),
        ('GNP arrest (−SHH)', '0', vb['div_GNPnoSHH'], vl['div_GNPnoSHH'], vl['div_GNPnoSHH'] > 0),
        ('GNP arrest (+HHi)', '0', vb['div_GNPHHi'], vl['div_GNPHHi'], vl['div_GNPHHi'] > 0)]
c.set_title('C · Validation damage — what the latch BREAKS', fontsize=11, fontweight='bold', loc='left')
yc = 0.90; c.text(0.02, 0.99, 'metric', fontsize=9, fontweight='bold', transform=c.transAxes)
c.text(0.52, 0.99, 'target', fontsize=9, fontweight='bold', transform=c.transAxes)
c.text(0.68, 0.99, 'BAKED', fontsize=9, fontweight='bold', color=C_BAK, transform=c.transAxes)
c.text(0.84, 0.99, 'LATCH', fontsize=9, fontweight='bold', color=C_LAT, transform=c.transAxes)
for name, tg, bkv, ltv, broken in rows:
    c.text(0.02, yc, name, fontsize=9, transform=c.transAxes, va='center')
    c.text(0.52, yc, tg, fontsize=9, transform=c.transAxes, va='center', color='#555')
    c.text(0.68, yc, f'{bkv:.2f}', fontsize=9, transform=c.transAxes, va='center', color=C_BAK)
    c.text(0.84, yc, f'{ltv:.2f}', fontsize=9.5, fontweight='bold' if broken else 'normal',
           color=(C_LAT if broken else '#2e7d32'), transform=c.transAxes, va='center')
    if broken:
        c.text(0.955, yc, '✗', fontsize=12, color=C_LAT, transform=c.transAxes, va='center', fontweight='bold')
    yc -= 0.115
c.text(0.02, -0.02, 'Biggest break: GNP proliferates WITHOUT mitogen (0→11 divisions) — Hedgehog-dependence LOST;\n'
                    'and the de-saturated low branch FLIPS the ChIP direction (MB now > GNP).',
       fontsize=9, style='italic', color=C_LAT, transform=c.transAxes, va='top')

# ---------- D: withdrawal-memory time course ----------
d = ax[1, 1]
wb, wl = R['wd_baked'], R['wd_latch']
b0 = np.cumsum([0] + list(wb['segbounds']))[:3] / 62.8
d.axvspan(b0[1], b0[2], color='#eaeaea', alpha=0.8, zorder=0)
d.text((b0[1] + b0[2]) / 2, 1.04, 'SHH withdrawn', ha='center', fontsize=9, color='#666')
d.text(b0[1] / 2, 1.04, 'proliferate', ha='center', fontsize=9, color='#666')
d.text((b0[2] + wb['t'][-1]) / 2, 1.04, 'SHH restored', ha='center', fontsize=9, color='#666')
d.plot(wb['t'], wb['Mk'], color=C_BAK, lw=2, label='BAKED  Mk')
d.plot(wl['t'], wl['Mk'], color=C_LAT, lw=2, label='LATCH  Mk')
for dt in wb['div_t']:
    d.plot(dt, -0.03, '|', color=C_BAK, ms=7, mew=1.4)
for dt in wl['div_t']:
    d.plot(dt, -0.07, '|', color=C_LAT, ms=7, mew=1.4)
d.annotate('LATCH: mark stays in the LOW basin\n(never latches; cell keeps dividing)', (wl['t'][-1] * 0.62, 0.03),
           xytext=(wl['t'][-1] * 0.36, 0.20), fontsize=8.5, color=C_LAT, fontweight='bold',
           arrowprops=dict(arrowstyle='->', color=C_LAT, lw=1.2))
d.set_xlabel('time (h)'); d.set_ylabel('H3K27me3 mark  Mk'); d.set_ylim(-0.11, 1.14)
d.set_title('D · Withdrawal test — NO spontaneous memory\nbaked arrests then re-proliferates; latch just keeps cycling', fontsize=11, fontweight='bold')
d.legend(fontsize=9, loc='center right'); d.grid(alpha=0.2)
d.text(0.02, 0.14, 'division ticks (│ blue=baked, red=latch) below the axis', transform=d.transAxes,
       ha='left', fontsize=8, style='italic', color='#666')

fig.suptitle('What happens to v44 in the LATCH regime (low a0 + cooperative read-write)?', fontsize=14, fontweight='bold', y=0.965)
fig.text(0.5, 0.925, 'A latch IS reachable — but only in ARRESTED cells, it never spontaneously engages, and it BREAKS the calibration (GNP loses Hedgehog-dependence).',
         ha='center', fontsize=10.5, color='#c0392b', fontweight='bold')
out = os.path.join(HERE, 'fig_latch_consequences.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
