"""Figure: cell-cycle duration sets H3K27me3 recovery from S-phase dilution (the missing buffering term). (JP 2026-07-20)
Reads sim_mark_cellcycle_recovery_cache.pkl. Run: ./venv/bin/python simulations/fig_mark_cellcycle_recovery.py"""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
R = pickle.load(open(os.path.join(HERE, 'sim_mark_cellcycle_recovery_cache.pkl'), 'rb'))
del_mk = R['del_mk']; c_s, c_l, c_mk = '#c0392b', '#2e86c1', '#1e8449'

fig, ax = plt.subplots(1, 4, figsize=(18, 4.6))
fig.subplots_adjust(left=0.045, right=0.99, top=0.82, bottom=0.16, wspace=0.34)

# A: within-cycle mark trajectory, short vs long Tc (S-dip + recovery)
a = ax[0]
for key, col, lab in (('long_Tc', c_l, f"long Tc ({R.get('long_Tc',{}).get('Tc',0):.0f} h)"),
                      ('short_Tc', c_s, f"short Tc ({R.get('short_Tc',{}).get('Tc',0):.0f} h)")):
    d = R.get(key)
    if d is None: continue
    a.plot(d['t'], d['mk'], color=col, lw=2, label=lab)
    dv = np.where((d['dna'][:-1] > 0.9) & (d['dna'][1:] < 0.1))[0]
    for i in dv: a.axvline(d['t'][i], color=col, ls=':', lw=1, alpha=0.5)
a.set_xlabel('time (h)'); a.set_ylabel('H3K27me3 (Mk)')
a.set_title('A  Within-cycle: S-phase dilution + recovery\n(÷ lines; short Tc recovers less each cycle)', fontsize=10.5, fontweight='bold')
a.legend(fontsize=8.5); a.grid(alpha=0.2)

# B: mark level vs Tc
b = ax[1]; sw = R['sweep']; Tc = sw[:, 1]; mk = sw[:, 4]
b.plot(Tc, mk, '-o', color=c_mk, lw=2.2, ms=6)
b.set_xlabel('cell-cycle duration Tc (h)'); b.set_ylabel('mean H3K27me3 (Mk)')
b.set_title('B  Longer Tc → more recovery → higher mark\n(the missing cell-cycle term)', fontsize=10.5, fontweight='bold')
b.grid(alpha=0.2)

# C: effective decay del_mk + ln2/Tc -- dilution dominates
c = ax[2]; dil = np.log(2) / (Tc * UPH)          # 1/model-unit
c.plot(Tc, dil, '-o', color='#c0392b', lw=2.2, ms=5, label='replication dilution  ln2/Tc')
c.axhline(del_mk, color='#7f8c8d', ls='--', lw=2, label=f'intrinsic turnover del_mk={del_mk:.4f}')
c.plot(Tc, del_mk + dil, '-s', color='k', lw=1.6, ms=4, label='effective decay  del_mk + ln2/Tc')
c.axhline(0.00044, color='#e67e22', ls=':', lw=1.6, label='default del_mk=0.00044')
c.plot(Tc, 0.00044 + dil, '-', color='#e67e22', lw=1.2, alpha=0.7)
c.set_yscale('log'); c.set_xlabel('cell-cycle duration Tc (h)'); c.set_ylabel('mark decay rate (1/model-time)')
c.set_title('C  Effective decay = del_mk + ln2/Tc\n(cell cycle adds ln2/Tc; ~ del_mk, dominates at fast Tc)', fontsize=10.5, fontweight='bold')
c.legend(fontsize=7.2); c.grid(alpha=0.2, which='both')

# D: fit -- Mk vs EZH2/(del_mk + ln2/Tc)
d = ax[3]; ez = sw[:, 5]
pred = ez / (del_mk + dil)
d.scatter(pred, mk, c=Tc, cmap='viridis', s=55, zorder=3)
z = np.polyfit(pred, mk, 1); xx = np.linspace(pred.min(), pred.max(), 20)
d.plot(xx, np.polyval(z, xx), 'k--', lw=1.4, alpha=0.7)
d.set_xlabel('writing / decay  =  EZH2 / (del_mk + ln2/Tc)'); d.set_ylabel('mean H3K27me3 (Mk)')
d.set_title(f'D  Mk ≈ w/(del_mk+ln2/Tc)\n(corr {np.corrcoef(mk,pred)[0,1]:.2f}; color = Tc)', fontsize=10.5, fontweight='bold')
cb = fig.colorbar(d.collections[0], ax=d, fraction=0.046); cb.set_label('Tc (h)', fontsize=8)
d.grid(alpha=0.2)

fig.suptitle('Cell-cycle duration sets H3K27me3 recovery from S-phase dilution — the MISSING buffering term: effective mark decay = del_mk + ln2/Tc, so the mark level & its buffering timescale depend on the CELL-CYCLE DURATION (dilution ~ del_mk; dominates at fast cycling)',
             fontsize=11, fontweight='bold', y=0.97)
out = os.path.join(HERE, 'fig_mark_cellcycle_recovery.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
