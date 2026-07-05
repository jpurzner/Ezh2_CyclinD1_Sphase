"""Plot the high-resolution parameter sweep (sweep_param_exploration.py). Robust to partial/nan data
(reads whatever the checkpointed cache currently holds). Produces:
  fig_v44_paramsweep_entry.png       — entry mitogen-threshold heatmap (fine) + contours
  fig_v44_paramsweep_withdrawal.png  — withdrawal exit-Hh & divisions, small-multiples over decline rate
  fig_v44_paramsweep_cuts.png        — 1-D line cuts through the default operating point

Run:  ./venv/bin/python simulations/plot_param_exploration.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

CACHE = 'simulations/sweep_param_exploration_cache.npz'
z = np.load(CACHE, allow_pickle=True)
KEZ, F0, DEPTH, DRATES = z['KEZ'], z['F0'], z['DEPTH'], z['DRATES']
ZA, ZB_nd, ZB_ex = z['ZA'], z['ZB_nd'], z['ZB_ex']
KEZ0, F00, DEPTH0 = 0.022, 0.233, 1 - 0.233
nfrac = np.mean(np.isnan(ZA))
print(f'grids KEZ={len(KEZ)} F0={len(F0)} D={len(DRATES)};  entry {100*(1-nfrac):.0f}% filled;  '
      f'withdrawal {100*np.mean(~np.isnan(ZB_ex)):.0f}% filled')


def _mesh_edges(c):
    c = np.asarray(c, float); e = np.empty(len(c) + 1)
    e[1:-1] = 0.5 * (c[:-1] + c[1:]); e[0] = c[0] - (c[1] - c[0]) / 2; e[-1] = c[-1] + (c[-1] - c[-2]) / 2
    return e


KE, DE = _mesh_edges(KEZ), _mesh_edges(DEPTH)

# ---------- Fig 1: entry threshold ----------
fig, ax = plt.subplots(figsize=(8.5, 6.5))
pm = ax.pcolormesh(KE, DE, ZA, cmap='viridis', shading='flat')
try:
    cs = ax.contour(KEZ, DEPTH, np.ma.masked_invalid(ZA), levels=[0.2, 0.35, 0.5, 0.75, 1.0, 1.5], colors='w', linewidths=1.0)
    ax.clabel(cs, fmt='%.2f', fontsize=7)
except Exception:
    pass
ax.plot(KEZ0, DEPTH0, 'r*', ms=16, label='default')
ax.set_xlabel('EZH2 mitogen-responsiveness  (kEZE2f)'); ax.set_ylabel('CyclinD1 inhibition depth  (1 − f0)')
ax.set_title('Mitogen (SHH) threshold for first division\nboth EZH2 knobs raise the Hh gate', fontweight='bold')
fig.colorbar(pm, ax=ax, label='entry SHH threshold'); ax.legend(loc='lower right')
plt.tight_layout(); plt.savefig('simulations/fig_v44_paramsweep_entry.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_paramsweep_entry.pdf', bbox_inches='tight'); plt.close()

# ---------- Fig 2: withdrawal small-multiples ----------
nD = len(DRATES)
fig, axes = plt.subplots(2, nD, figsize=(3.0 * nD, 7.2))
exv = np.ma.masked_invalid(ZB_ex); ndv = np.ma.masked_invalid(ZB_nd)
exmin, exmax = np.nanmin(ZB_ex), np.nanmax(ZB_ex)
ndmax = np.nanmax(ZB_nd) if np.any(~np.isnan(ZB_nd)) else 1
for di, D in enumerate(DRATES):
    a = axes[0, di]
    pm = a.pcolormesh(KE, DE, ZB_ex[di], cmap='magma_r', shading='flat', vmin=exmin, vmax=exmax)
    a.plot(KEZ0, DEPTH0, 'c*', ms=9); a.set_title(f'{D/60:.0f} h', fontsize=10, fontweight='bold')
    if di == 0:
        a.set_ylabel('exit-Hh readout\ninhibition depth (1−f0)', fontsize=9)
    a.set_xticks([]);
    if di == nD - 1:
        fig.colorbar(pm, ax=a, label='exit SHH')
    b = axes[1, di]
    pm2 = b.pcolormesh(KE, DE, ZB_nd[di], cmap='cividis', shading='flat', vmin=0, vmax=ndmax)
    b.plot(KEZ0, DEPTH0, 'r*', ms=9); b.set_xlabel('kEZE2f', fontsize=8)
    if di == 0:
        b.set_ylabel('divisions readout\ninhibition depth (1−f0)', fontsize=9)
    if di == nD - 1:
        fig.colorbar(pm2, ax=b, label='withdrawal divisions')
fig.suptitle('Hh withdrawal vs (EZH2 responsiveness × inhibition depth), per decline rate (columns): higher/faster brake → exit at higher Hh, fewer divisions',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout(); plt.savefig('simulations/fig_v44_paramsweep_withdrawal.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_paramsweep_withdrawal.pdf', bbox_inches='tight'); plt.close()

# ---------- Fig 2b: CUMULATIVE PROLIFERATION (fold expansion 2^divisions) ----------
from scipy.ndimage import median_filter
CUM = np.full_like(ZB_nd, np.nan)
for di in range(nD):
    nds = median_filter(np.where(np.isnan(ZB_nd[di]), 0.0, ZB_nd[di]), size=3)   # de-speckle the discrete count
    CUM[di] = 2.0 ** nds
cmax = np.nanmax(CUM)
nc = 5; nr = int(np.ceil(nD / nc))
fig, axes = plt.subplots(nr, nc, figsize=(3.4 * nc, 3.1 * nr))
axes = np.atleast_2d(axes)
for di, D in enumerate(DRATES):
    a = axes[di // nc, di % nc]
    pm = a.pcolormesh(KE, DE, CUM[di], cmap='inferno', shading='flat', norm=LogNorm(vmin=1, vmax=cmax))
    a.plot(KEZ0, DEPTH0, 'c*', ms=10)
    a.set_title(f'Hh decline {D/60:.0f} h', fontsize=10, fontweight='bold')
    if di % nc == 0:
        a.set_ylabel('inhibition depth (1−f0)', fontsize=8)
    if di // nc == nr - 1:
        a.set_xlabel('kEZE2f', fontsize=8)
    fig.colorbar(pm, ax=a, label='fold expansion 2^div')
for k in range(nD, nr * nc):
    axes[k // nc, k % nc].axis('off')
fig.suptitle('Cumulative proliferation during Hh withdrawal (fold expansion = 2^divisions; symmetric GNP divisions): slow decline + weak/unresponsive EZH2 → exponentially more progeny',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout(); plt.savefig('simulations/fig_v44_paramsweep_cumprolif.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_paramsweep_cumprolif.pdf', bbox_inches='tight'); plt.close()

# ---------- Fig 3: line cuts through default ----------
ki = int(np.argmin(np.abs(KEZ - KEZ0))); fi = int(np.argmin(np.abs(F0 - F00)))
di_mid = int(np.argmin(np.abs(DRATES - 9000)))
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
ax[0].plot(DEPTH, ZA[:, ki], '-o', ms=3, label=f'kEZE2f={KEZ[ki]:.3f} (default col)')
ax[0].plot(DEPTH, ZA[:, int(0.85 * len(KEZ))], '-o', ms=3, label=f'kEZE2f={KEZ[int(0.85*len(KEZ))]:.3f} (high)')
ax[0].axvline(DEPTH0, color='k', ls=':', alpha=0.4)
ax[0].set_xlabel('inhibition depth (1−f0)'); ax[0].set_ylabel('entry SHH threshold'); ax[0].set_title('(A) Entry threshold vs inhibition', fontweight='bold'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2)
ax[1].plot(KEZ, ZA[fi, :], '-o', ms=3, color='#8e44ad', label=f'f0={F0[fi]:.2f} (default row)')
ax[1].axvline(KEZ0, color='k', ls=':', alpha=0.4)
ax[1].set_xlabel('EZH2 responsiveness (kEZE2f)'); ax[1].set_ylabel('entry SHH threshold'); ax[1].set_title('(B) Entry threshold vs EZH2 responsiveness', fontweight='bold'); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.2)
CUMc = 2.0 ** np.stack([median_filter(np.where(np.isnan(ZB_nd[d]), 0.0, ZB_nd[d]), size=3) for d in range(nD)])
ax[2].plot(DRATES / 60, CUMc[:, int(0.9 * len(F0)), ki], '-o', ms=4, color='#e67e22', label='weak inhib')
ax[2].plot(DRATES / 60, CUMc[:, fi, ki], '-o', ms=4, color='#117a65', label='default')
ax[2].plot(DRATES / 60, CUMc[:, int(0.15 * len(F0)), ki], '-o', ms=4, color='#c0392b', label='strong inhib')
ax[2].set_yscale('log'); ax[2].set_xlabel('Hh decline time (h)'); ax[2].set_ylabel('cumulative proliferation (2^div)'); ax[2].set_title('(C) Cumulative proliferation vs decline rate', fontweight='bold'); ax[2].legend(fontsize=8); ax[2].grid(alpha=0.2, which='both')
plt.tight_layout(); plt.savefig('simulations/fig_v44_paramsweep_cuts.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_paramsweep_cuts.pdf', bbox_inches='tight'); plt.close()
print('Saved: fig_v44_paramsweep_entry / _withdrawal / _cuts .png/.pdf')
