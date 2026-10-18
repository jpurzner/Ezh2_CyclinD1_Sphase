"""Poised bivalent differentiation gene -> PREMATURE GENE ACTIVATION (JP 2026-07-21). Serial me1->me2->me3
methylation; HIGH activation threshold; downstream the fired cell completes its final mitosis then exits to
G0 (differentiation itself is complex/G0-dependent, not modelled). 4 panels:
 A  normal GNP: serial chain builds heavy me3 + H3K4me3 poise -> gene OFF, cell keeps cycling
 B  divide TOO FAST -> me3 collapses below the poise -> PREMATURE ACTIVATION -> completes final mitosis -> G0
 C  EZH2i (writer loss) -> same premature activation + exit
 D  dose-response: high threshold -> gene stays OFF until division is pushed hard, then latches ON (cell exits)
Reads sim_diff_gene_cache.pkl. Run: ./venv/bin/python simulations/fig_diff_gene.py"""
import os, pickle, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_diff_gene_cache.pkl'), 'rb'))
C_ME3, C_DIFF, C_K4, C_DIV, C_ME2, C_ME1 = '#5b2c6f', '#c0392b', '#1e8449', '#2e86c1', '#a569bd', '#d2b4de'
H3K4 = 0.55

fig, ax = plt.subplots(2, 2, figsize=(15.4, 10.8))
fig.subplots_adjust(left=0.065, right=0.94, top=0.87, bottom=0.075, hspace=0.34, wspace=0.34)

def tcourse(a, d, title, onset=None, show_chain=False):
    a.axhline(H3K4, color=C_K4, ls='-', lw=1.6, alpha=0.85)
    a.text(d['t'][-1] * 0.99, H3K4 + 0.02, 'H3K4me3 poise', color=C_K4, fontsize=8, ha='right', fontweight='bold')
    if show_chain:
        a.plot(d['t'], d['me1'], color=C_ME1, lw=1.3, label='H3K27me1', alpha=0.9)
        a.plot(d['t'], d['me2'], color=C_ME2, lw=1.3, label='H3K27me2', alpha=0.9)
    a.plot(d['t'], d['Mk_diff'], color=C_ME3, lw=2.4, label='H3K27me3 (me3)')
    a.plot(d['t'], d['Diff'], color=C_DIFF, lw=2.4, label='gene activation (Diff)')
    if onset is not None:
        a.axvline(onset, color='#333', ls=':', lw=1.5)
    a.set_ylim(-0.03, 1.06); a.set_xlabel('time (h)'); a.set_ylabel('mark / activity')
    a.set_title(title, fontsize=10.6, fontweight='bold'); a.grid(alpha=0.2)
    a2 = a.twinx()
    a2.plot(d['t'], d['cumdiv'], color=C_DIV, lw=2.0, ls='--')
    a2.set_ylabel('cumulative divisions', color=C_DIV, fontsize=9); a2.tick_params(axis='y', labelcolor=C_DIV)
    a2.set_ylim(0, max(d['cumdiv'][-1], 5) * 1.12)
    a.legend(fontsize=7.6, loc='center left')
    return a2

tcourse(ax[0, 0], R['normal'], 'A · Normal GNP — bivalent & poised (serial me1→me2→me3)\nheavy me3 + H3K4me3 → gene OFF, cell keeps cycling', show_chain=True)
ax[0, 0].annotate('serial chain builds me3;\ngene repressed but "ready"', (R['normal']['t'][-1] * 0.52, 0.30), fontsize=8, ha='center',
                  color=C_ME3, fontweight='bold', bbox=dict(boxstyle='round', fc='#f4ecf7', ec=C_ME3))

tcourse(ax[0, 1], R['fast'], 'B · Divide TOO FAST → PREMATURE gene ACTIVATION\nme3 collapses; gene fires; completes final mitosis → G0', onset=R['onset_h'])
ax[0, 1].text(R['onset_h'] + 2, 0.42, 'push faster', fontsize=8.5, color='#333', rotation=90, va='center')
ax[0, 1].annotate('completes final\nmitoses → G0', (R['fast']['t'][-1] * 0.72, 0.72), fontsize=8.3, ha='center', color=C_DIV, fontweight='bold')

tcourse(ax[1, 0], R['ezh2i'], 'C · EZH2i (writer loss) → PREMATURE gene ACTIVATION\nme3 collapses; gene fires; completes final mitosis → G0', onset=R['onset_h'])
ax[1, 0].text(R['onset_h'] + 2, 0.42, 'EZH2i on', fontsize=8.5, color='#333', rotation=90, va='center')

# D: dose-response (high threshold)
d = ax[1, 1]; s = R['sweep']; order = np.argsort(s['mu']); mu = s['mu'][order] * 1e3
Df = s['Diff'][order]; dl = s['div_late'][order]
d.plot(mu, Df, '-s', color=C_DIFF, lw=2.6, ms=7, label='gene activation (final)')
d.set_xlabel('growth drive  μ ×10³   (→ harder push / faster division)'); d.set_ylabel('gene activation  Diff', color=C_DIFF)
d.tick_params(axis='y', labelcolor=C_DIFF); d.set_ylim(-0.03, 1.06)
d.set_title('D · HIGH threshold — gene stays OFF through normal & moderate\nspeed-up, then latches ON only when pushed hard (cell exits)', fontsize=10.6, fontweight='bold')
d.grid(alpha=0.2)
d2 = d.twinx(); d2.plot(mu, dl, '-o', color=C_DIV, lw=2.0, ms=6)
d2.set_ylabel('divisions after settling', color=C_DIV, fontsize=9); d2.tick_params(axis='y', labelcolor=C_DIV)
fired = mu[Df > 0.5]
if len(fired):
    d.axvspan(fired.min(), mu.max(), color=C_DIFF, alpha=0.06)
    d.text((fired.min() + mu.max()) / 2, 0.66, 'PREMATURE\nactivation →\nG0 exit', ha='center', fontsize=8.8, color=C_DIFF, fontweight='bold')
d.text(mu.min() * 1.04, 0.2, 'poised / cycling\n(robustly repressed)', fontsize=8.6, color=C_K4, fontweight='bold')
d.legend(fontsize=7.8, loc=(0.02, 0.5))

fig.suptitle('PREMATURE differentiation-gene ACTIVATION: a poised, heavily-me3 (serial me1→me2→me3) bivalent gene fires only when the writer can\'t keep up',
             fontsize=12, fontweight='bold', y=0.955)
fig.text(0.5, 0.915, 'HIGH threshold: fires only by dividing too fast (B — dilution outruns the slow serial writer) or losing the writer (C, EZH2i); '
                     'then the cell completes its final mitosis and exits to G0. Default OFF (28/29).',
         ha='center', fontsize=8.8, color='#444')
out = os.path.join(HERE, 'fig_diff_gene.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
