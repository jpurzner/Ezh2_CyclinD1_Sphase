"""Figure: the memory function on the current v44 default (Cdh1-EZH2 + G1-lengthening + Skp2) under mitogen RAMP
up/down (hysteresis) and FLUCTUATION (low-pass buffering). Reads sim_memory_ramp_fluctuation_cache.pkl. (JP 2026-07-20)
Run: ./venv/bin/python simulations/fig_memory_ramp_fluctuation.py"""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_memory_ramp_fluctuation_cache.pkl'), 'rb'))
c_cd, c_ez, c_mk, c_sk, c_up, c_dn, c_off = '#c0392b', '#6c3483', '#1e8449', '#d68910', '#2e86c1', '#e74c3c', '#95a5a6'

fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))
fig.subplots_adjust(left=0.055, right=0.985, top=0.9, bottom=0.075, hspace=0.36, wspace=0.30)


def norm(x): x = np.asarray(x); return x / (np.nanmax(np.abs(x)) + 1e-9)


def binned(shh, y, lo=0.06, hi=0.99, nb=28):
    """Cycle-average out the fast cell-cycle oscillation: mean y in Shh bins (reveals the hysteresis loop)."""
    shh = np.asarray(shh); y = np.asarray(y); edges = np.linspace(lo, hi, nb + 1)
    ctr = 0.5 * (edges[:-1] + edges[1:]); m = np.full(nb, np.nan)
    idx = np.digitize(shh, edges) - 1
    for i in range(nb):
        sel = idx == i
        if sel.any(): m[i] = np.mean(y[sel])
    ok = ~np.isnan(m)
    return ctr[ok], m[ok]

# ---- A: CyclinD1 hysteresis loop (fast ramp) -- full vs mark-off (Shh-binned = cycle-averaged) ----
a = ax[0, 0]; f = R['ramp_60_full']; o = R['ramp_60_markoff']
a.plot(*binned(f['sU'], f['U_Cd']), color=c_up, lw=2.4, label='UP-ramp (full)')
a.plot(*binned(f['sD'], f['D_Cd']), color=c_dn, lw=2.4, label='DOWN-ramp (full)')
a.plot(*binned(o['sU'], o['U_Cd']), color=c_off, lw=1.8, ls='--', alpha=0.9, label='mark OFF (up)')
a.plot(*binned(o['sD'], o['D_Cd']), color=c_off, lw=1.8, ls=':', alpha=0.9, label='mark OFF (down)')
a.set_xlabel('mitogen (Shh)'); a.set_ylabel('CyclinD1 (cycle-avg)'); a.set_title('A  CyclinD1 hysteresis (fast ramp 60h)\nloop needs the mark (+0.09) — mark-off has none', fontsize=10.5, fontweight='bold')
a.legend(fontsize=8, loc='upper left'); a.grid(alpha=0.2)

# ---- B: EZH2 + mark hysteresis (the slow memory variables lag) ----
b = ax[0, 1]
b.plot(*binned(f['sU'], f['U_EZH2']), color=c_ez, lw=2.2, label='EZH2 up')
b.plot(*binned(f['sD'], f['D_EZH2']), color=c_ez, lw=2.2, ls='--', label='EZH2 down')
b2 = b.twinx()
b2.plot(*binned(f['sU'], f['U_Mk']), color=c_mk, lw=2.2, label='H3K27me3 up')
b2.plot(*binned(f['sD'], f['D_Mk']), color=c_mk, lw=2.2, ls='--', label='H3K27me3 down')
b.set_xlabel('mitogen (Shh)'); b.set_ylabel('EZH2', color=c_ez); b.tick_params(axis='y', labelcolor=c_ez)
b2.set_ylabel('H3K27me3', color=c_mk); b2.tick_params(axis='y', labelcolor=c_mk)
b.set_title('B  Slow memory variables lag the ramp\n(EZH2 & mark up-branch ≠ down-branch)', fontsize=10.5, fontweight='bold')
b.legend(fontsize=7.5, loc='upper left'); b2.legend(fontsize=7.5, loc='lower right'); b.grid(alpha=0.2)

# ---- C: proliferation vs Shh, up vs down (does division persist on the way down?) ----
c = ax[0, 2]
c.plot(f['sU'], f['rU'], color=c_up, lw=2.2, label='UP-ramp')
c.plot(f['sD'], f['rD'], color=c_dn, lw=2.2, label='DOWN-ramp')
c.set_xlabel('mitogen (Shh)'); c.set_ylabel('division rate (1/h)')
c.set_title('C  Proliferation hysteresis\n(division persists to lower Shh on the way down)', fontsize=10.5, fontweight='bold')
c.legend(fontsize=8); c.grid(alpha=0.2)

# ---- D: fluctuation time series (period 33h) -- Cd tracks, EZH2 attenuated, mark flat ----
d = ax[1, 0]; fl = R['fluct_33']; t = fl['t'] / 62.8
d.plot(t, norm(fl['SHH']), color='#7f8c8d', lw=1.4, label='Shh (drive)')
d.plot(t, norm(fl['Cd']), color=c_cd, lw=2, label='CyclinD1 (tracks)')
d.plot(t, norm(fl['EZH2']), color=c_ez, lw=2, label='EZH2 (attenuated)')
d.plot(t, norm(fl['Mk']), color=c_mk, lw=2.2, label='H3K27me3 (nearly flat)')
d.set_xlabel('time (h)'); d.set_ylabel('normalized level'); d.set_xlim(t.min(), t.min() + 200)
d.set_title('D  Fluctuation (Shh period 33h): CyclinD1 follows,\nEZH2 buffered, mark is a stable memory', fontsize=10.5, fontweight='bold')
d.legend(fontsize=8, loc='lower right'); d.grid(alpha=0.2)

# ---- E: fluctuation low-pass -- CV of each variable vs drive period ----
e = ax[1, 1]
periods = [14, 33, 80]
cv = {k: [] for k in ('Cd', 'EZH2', 'Mk', 'Skp2')}
for T in periods:
    fl = R[f'fluct_{T}']; w = fl['t'] > T * 62.8
    for k in cv: cv[k].append(np.std(fl[k][w]) / (np.mean(fl[k][w]) + 1e-9))
for k, col, lab in [('Cd', c_cd, 'CyclinD1'), ('EZH2', c_ez, 'EZH2'), ('Skp2', c_sk, 'Skp2'), ('Mk', c_mk, 'H3K27me3')]:
    e.plot(periods, cv[k], '-o', color=col, lw=2, ms=6, label=lab)
e.set_xlabel('Shh fluctuation period (h)'); e.set_ylabel('response CV (variability)')
e.set_title('E  Low-pass filter: CyclinD1 varies most,\nthe mark is the deepest-buffered memory', fontsize=10.5, fontweight='bold')
e.legend(fontsize=8); e.grid(alpha=0.2)

# ---- F: Skp2 through the ramp (new proliferation-tracking variable, Shh-binned) ----
ff = ax[1, 2]
ff.plot(*binned(f['sU'], f['U_Skp2']), color=c_sk, lw=2.4, label='Skp2 up')
ff.plot(*binned(f['sD'], f['D_Skp2']), color=c_sk, lw=2.4, ls='--', label='Skp2 down')
cU_s, cU_v = binned(f['sU'], f['U_Cd'])
ff.plot(cU_s, cU_v / np.nanmax(cU_v) * np.nanmax(f['U_Skp2']), color=c_cd, lw=1.5, alpha=0.6, label='CyclinD1 (scaled)')
ff.set_xlabel('mitogen (Shh)'); ff.set_ylabel('Skp2')
ff.set_title('F  Skp2 through the ramp\n(mitogen-induced, tracks proliferation)', fontsize=10.5, fontweight='bold')
ff.legend(fontsize=8); ff.grid(alpha=0.2)

fig.suptitle('Memory function on the current model: a slow, stable H3K27me3 mark low-pass-filters mitogen fluctuations and creates CyclinD1/proliferation hysteresis under ramps '
             '— division persistence is a chromatin + cell-cycle memory', fontsize=12, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_memory_ramp_fluctuation.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
