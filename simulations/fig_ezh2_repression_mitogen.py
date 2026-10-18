"""Figure: EZH2 -| CyclinD1 repression ON vs OFF -> the mitogen requirement (GNP & MB, multi-division), GNP ramp
hysteresis, and the 2D search for a regime with BOTH the mitogen-shift AND a ChIP-seq-correct (peaks-at-
proliferation) mark. Reads sim_ezh2_repression_mitogen_cache.pkl. (JP 2026-07-20)
Run: ./venv/bin/python simulations/fig_ezh2_repression_mitogen.py"""
import os, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_ezh2_repression_mitogen_cache.pkl'), 'rb'))
c_on, c_off, c_cd, c_ez, c_mk = '#c0392b', '#2e86c1', '#0e6251', '#6c3483', '#1e8449'

fig, ax = plt.subplots(2, 3, figsize=(16.8, 9.4))
fig.subplots_adjust(left=0.055, right=0.975, top=0.9, bottom=0.075, hspace=0.4, wspace=0.32)


def dose_panel(a, ct, title):
    on, off = R[f'dose_{ct}_ON'], R[f'dose_{ct}_OFF']
    a.plot(on['shh'], on['div_h'] * 24, '-o', color=c_on, lw=2.2, ms=5, label='repression ON')
    a.plot(off['shh'], off['div_h'] * 24, '-o', color=c_off, lw=2.2, ms=5, label='repression OFF')
    for d, col in ((on, c_on), (off, c_off)):
        if not np.isnan(d['thr']): a.axvline(d['thr'], color=col, ls=':', lw=1.4, alpha=0.8)
    sh = (on['thr'] / off['thr']) if off['thr'] else float('nan')
    a.set_xscale('log'); a.set_xlabel('mitogen (Shh, log)'); a.set_ylabel('division rate (÷ / day)')
    a.set_title(f'{title}\nrepression ON needs {sh:.2f}× more mitogen', fontsize=10.5, fontweight='bold')
    a.legend(fontsize=8.5, loc='upper left'); a.grid(alpha=0.2)


# A, B: dose-response GNP & MB
dose_panel(ax[0, 0], 'GNP', 'A  GNP mitogen dose-response')
dose_panel(ax[0, 1], 'MB', 'B  MB mitogen dose-response')

# C: GNP ramp hysteresis (fast rate), ON vs OFF -- division rate vs Shh
c = ax[0, 2]
key = 'ramp_60_ON' if 'ramp_60_ON' in R else ('ramp_120_ON' if 'ramp_120_ON' in R else None)
if key:
    dur = key.split('_')[1]
    for lab, col in (('ON', c_on), ('OFF', c_off)):
        rr = R[f'ramp_{dur}_{lab}']
        c.plot(rr['up']['shh'], rr['up']['div_h'] * 24, '-', color=col, lw=2, label=f'{lab} up')
        c.plot(rr['dn']['shh'], rr['dn']['div_h'] * 24, '--', color=col, lw=2, label=f'{lab} down')
c.set_xlabel('mitogen (Shh)'); c.set_ylabel('division rate (÷ / day)')
c.set_title('C  GNP ramp up/down hysteresis\n(memory; down-branch persists to lower Shh)', fontsize=10.5, fontweight='bold')
c.legend(fontsize=7.5, ncol=2); c.grid(alpha=0.2)

# D: multi-division trajectory (GNP ON) -- Cd/EZH2/mark + division ticks
d = ax[1, 0]
tr = R.get('traj_GNP_ON')
if tr is not None:
    t = tr['time'] / 62.8
    d.plot(t, tr['Cd'] / np.nanmax(tr['Cd']), color=c_cd, lw=1.5, label='CyclinD1')
    d.plot(t, tr['EZH2'] / np.nanmax(tr['EZH2']), color=c_ez, lw=1.5, label='EZH2')
    d.plot(t, tr['Mk'], color=c_mk, lw=1.8, label='H3K27me3')
    dv = np.where((tr['Dna'][:-1] > 0.9) & (tr['Dna'][1:] < 0.1))[0]
    for i in dv: d.axvline(t[i], color='#aaa', lw=0.8, alpha=0.6)
    d.set_xlim(0, min(300, t.max()))
d.set_xlabel('time (h)'); d.set_ylabel('level (norm.)')
d.set_title('D  GNP multi-division trajectory (ON)\n(grey = divisions; mark stays high)', fontsize=10.5, fontweight='bold')
d.legend(fontsize=8, loc='lower right'); d.grid(alpha=0.2)

# E: 2D sweep -- mitogen shift vs f0_prc2 x mark de-saturation (heatmap), 'BOTH' regime marked
e = ax[1, 1]
meta = R.get('sweep_meta'); off_ref = R.get('sweep_OFF_ref')
if meta and off_ref:
    KP, MARKS = meta['KP'], meta['MARKS']
    shift = np.full((len(KP), len(MARKS)), np.nan); slope = np.full_like(shift, np.nan)
    for ki in range(len(KP)):
        for mi in range(len(MARKS)):
            on = R.get(f'sweep_k{ki}_m{mi}')
            if on is None: continue
            shift[ki, mi] = (on['thr'] / off_ref['thr']) if off_ref['thr'] else np.nan
            slope[ki, mi] = on['Mk'][-3] - on['Mk'][2]
    im = e.imshow(np.clip(shift, 1, 8), aspect='auto', origin='lower', cmap='YlOrRd', vmin=1.0)
    for ki in range(len(KP)):
        for mi in range(len(MARKS)):
            if not np.isnan(slope[ki, mi]) and slope[ki, mi] > 0.005 and shift[ki, mi] > 1.5:
                e.plot(mi, ki, '*', color='#1e8449', ms=16, mec='k', mew=0.6)   # BOTH regime
            if not np.isnan(shift[ki, mi]):
                e.text(mi, ki, f"{shift[ki,mi]:.1f}", ha='center', va='center', fontsize=6.5,
                       color='k' if shift[ki, mi] < 4 else 'white')
    e.set_xticks(range(len(MARKS))); e.set_xticklabels([f'{i}' for i in range(len(MARKS))])
    e.set_yticks(range(len(KP))); e.set_yticklabels([f'{k:.4f}' for k in KP])
    e.set_xlabel('mark de-saturation idx (0=saturated → 4=strong)'); e.set_ylabel('K_prc2 (lower = repress at lower mark)')
    fig.colorbar(im, ax=e, fraction=0.046, label='mitogen shift (ON/OFF)')
e.set_title('E  2D search (f0=0.184): ★ = BOTH — mitogen\nshift>1.5× AND mark peaks at proliferation', fontsize=10.5, fontweight='bold')

# F: mark direction -- default (saturated, flat/down) vs a 'BOTH' regime (de-saturated + low K_prc2, tracks up)
f = ax[1, 2]
default = R.get('sweep_k0_m0')          # K=0.0074, saturated mark (current default-ish)
both_r = R.get('sweep_k1_m3') or R.get('sweep_k2_m3')   # K=0.004, de-saturated -> a BOTH candidate
if default is not None:
    f.plot(default['shh'], default['Mk'], '-', color='#95a5a6', lw=2.2, label='default (saturated, flat/down)')
if both_r is not None:
    f.plot(both_r['shh'], both_r['Mk'], '-', color=c_mk, lw=2.4, label="'both' regime (tracks proliferation UP)")
f.set_xscale('log'); f.set_xlabel('mitogen (Shh, log)'); f.set_ylabel('H3K27me3 @ CyclinD1 (Mk)')
f.set_title('F  Mark direction: default saturated-flat vs the\n"both" regime — mark UP with proliferation (ChIP-seq)', fontsize=10.5, fontweight='bold')
f.legend(fontsize=7.8); f.grid(alpha=0.2)

fig.suptitle('EZH2 ⊣ CyclinD1 repression raises the mitogen requirement (ON needs more mitogen); the search is for a regime with BOTH the requirement-shift AND a mark that peaks at max proliferation',
             fontsize=12, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_ezh2_repression_mitogen.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
