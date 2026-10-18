"""Transcript-MATCHED withdrawal (JP 2026-07-24): ON / OFF / EZH2i started at the SAME Ccnd1 transcript (~1.61,
each at its own Shh: ON 0.45, EZH2i 0.177, OFF 0.082), then Shh withdrawn at multiple decline rates over a
population of cells. Reads sim_matched_withdrawal_cache.npz. Run: ./venv/bin/python simulations/fig_matched_withdrawal.py"""
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
R = dict(np.load(os.path.join(HERE, 'sim_matched_withdrawal_cache.npz'), allow_pickle=True))
G = R['grid']; RATES = [int(x) for x in R['rates']]
MATCHED = R['matched'].item(); THRESH = R['thresh'].item()
CONDS = ['ON', 'OFF', 'EZH2i']; COL = {'ON': '#0e6251', 'OFF': '#c0392b', 'EZH2i': '#e67e22'}
RLAB = {RATES[0]: 'FAST', RATES[1]: 'MEDIUM', RATES[2]: 'SLOW'}
def band(a): return a.mean(0), a.mean(0) - a.std(0), a.mean(0) + a.std(0)

fig, ax = plt.subplots(2, 3, figsize=(21.5, 11.4))
fig.subplots_adjust(left=0.05, right=0.985, top=0.87, bottom=0.065, hspace=0.34, wspace=0.24)

# A/B/C: transcript decline for each decline rate
for j, D in enumerate(RATES):
    a = ax[0, j]
    for c in CONDS:
        mr = R[f'mr_{c}_{D}']; mu, lo, hi = band(mr)
        w = G >= -12
        a.plot(G[w], mu[w], color=COL[c], lw=2.6, label=c)
        a.fill_between(G[w], lo[w], hi[w], color=COL[c], alpha=0.13)
        last = float(np.median(R[f'last_{c}_{D}']))
        a.plot(last, np.interp(last, G, mu), 'v', color=COL[c], ms=9, mec='k', mew=0.5)
    a.axhline(1.61, color='#bbb', ls=':', lw=1.1); a.text(-11, 1.68, 'matched start', fontsize=7.5, color='#888')
    a.axvspan(0, D, color='#95a5a6', alpha=0.10); a.axvline(0, color='#333', ls=':', lw=1.1)
    a.text(D + 2, 0.15, 'Shh 0→0 over %dh' % D, fontsize=7.6, color='#555')
    a.set_xlabel('time from withdrawal onset (h)'); a.set_ylabel('CyclinD1 mRNA (Ccnd1)')
    a.set_xlim(-12, 175); a.set_ylim(0, 2.2)
    a.set_title('%s · %s withdrawal (Shh → 0 over %d h)\nmatched start; ▾ = median last division (arrest)' % ('ABC'[j], RLAB[D], D), fontsize=11.5, fontweight='bold')
    a.legend(fontsize=8.5, loc='upper right'); a.grid(alpha=0.2)

# D: the Shh trajectories (headroom) — medium rate
d = ax[1, 0]; Dm = RATES[1]
tt = np.linspace(-12, 175, 400)
for c in CONDS:
    m0 = MATCHED[c]; sh = np.where(tt < 0, m0, np.where(tt < Dm, m0 * np.clip(1 - tt / Dm, 0, 1), 0.0))
    d.plot(tt, sh, color=COL[c], lw=2.6, label='%s (start Shh %.2f)' % (c, m0))
    d.axhline(THRESH[c], color=COL[c], ls=':', lw=1.3)
    d.text(150, THRESH[c] + 0.005, '%s thr %.2f' % (c, THRESH[c]), fontsize=7.2, color=COL[c], ha='right')
d.axvspan(0, Dm, color='#95a5a6', alpha=0.10)
d.text(0.02, 0.55, 'MATCHED at the same transcript, the conditions sit at\nvery different Shh: ON has huge headroom above its\nthreshold, OFF almost none — yet all reach 0 together.', transform=d.transAxes, fontsize=8, color='#555')
d.set_xlabel('time from withdrawal onset (h)'); d.set_ylabel('Shh (each cell)'); d.set_xlim(-12, 175); d.set_ylim(0, 0.52)
d.set_title('D · The Shh they experience differs hugely (headroom)\neven at the SAME CyclinD1 — %s rate shown' % RLAB[Dm], fontsize=11.5, fontweight='bold')
d.legend(fontsize=8.3, loc='upper right'); d.grid(alpha=0.2)

# E: post-withdrawal divisions vs decline rate
e = ax[1, 1]; x = np.arange(len(RATES)); w = 0.25
for k, c in enumerate(CONDS):
    vals = [np.mean(R[f'post_{c}_{D}']) for D in RATES]; err = [np.std(R[f'post_{c}_{D}']) for D in RATES]
    e.bar(x + (k - 1) * w, vals, w, yerr=err, capsize=3, color=COL[c], alpha=0.9, label=c)
e.set_xticks(x); e.set_xticklabels(['%s\n(%dh)' % (RLAB[D], D) for D in RATES])
e.set_xlabel('Shh decline rate (ramp duration)'); e.set_ylabel('divisions after withdrawal (per cell)')
e.set_title('E · Divisions after withdrawal vs decline rate\n(matched start → similar counts; slower = a few more)', fontsize=11.5, fontweight='bold')
e.legend(fontsize=8.5, loc='upper left'); e.grid(alpha=0.2, axis='y')

# F: transcript floor after arrest vs decline rate (OFF stays elevated = cKO)
f = ax[1, 2]
i_end = np.argmin(np.abs(G - (G[-1] - 10)))
for c in CONDS:
    vals = [np.mean(R[f'mr_{c}_{D}'][:, i_end]) for D in RATES]; err = [np.std(R[f'mr_{c}_{D}'][:, i_end]) for D in RATES]
    f.plot(RATES, vals, '-o', color=COL[c], lw=2.4, ms=7, label=c)
    f.fill_between(RATES, np.array(vals) - np.array(err), np.array(vals) + np.array(err), color=COL[c], alpha=0.12)
f.set_xlabel('Shh decline rate (ramp duration, h)'); f.set_ylabel('CyclinD1 mRNA floor (arrested, end)')
f.text(0.03, 0.5, 'After arrest the de-repressed cell (OFF) keeps the\ntranscript ELEVATED (un-repressed basal) — ≈ the cKO\n(transcript high, but arrested/differentiating); the\ngoverned cell (ON) settles low.', transform=f.transAxes, fontsize=8, color='#555')
f.set_title('F · Post-arrest transcript FLOOR — de-repressed (OFF)\nstays high; governed (ON) settles low', fontsize=11.5, fontweight='bold')
f.legend(fontsize=8.5, loc='upper right'); f.grid(alpha=0.2)

fig.suptitle('Transcript-MATCHED withdrawal — ON/OFF/EZH2i started at the SAME Ccnd1 (~1.61; Shh 0.45/0.082/0.177): FAST decline → all arrest identically (memory-dominated); '
             'SLOW decline → de-repressed (OFF) keeps dividing far longer (threshold-dominated); OFF keeps the transcript elevated after arrest (≈ cKO)', fontsize=11.2, fontweight='bold', y=0.966)
out = os.path.join(HERE, 'fig_matched_withdrawal.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
for D in RATES:
    print(f"D={D}h: " + " | ".join(f"{c} post÷={np.mean(R[f'post_{c}_{D}']):.1f} last@={np.median(R[f'last_{c}_{D}']):.0f}h floor={np.mean(R[f'mr_{c}_{D}'][:,i_end]):.2f}" for c in CONDS))
