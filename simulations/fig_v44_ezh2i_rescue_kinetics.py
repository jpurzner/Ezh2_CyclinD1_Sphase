"""EZH2i fast-tracks CyclinD1 de-repression and rescues the vismodegib-arrested MB cell (calibrated
H3K27me3-MEMORY model). Vismo lowers Gli-driven CyclinD1 but the MYCN floor keeps it near threshold;
residual H3K27me3 at the CyclinD1 locus holds it just below, and EZH2i clears the mark to re-tip it.
EZH2i blocks H3K27me3 deposition -> the mark decays (active demethylation + replication-dilution) ->
CyclinD1 de-represses over ~24-48h -> cycle re-entry, AT WHATEVER TIME the inhibitor is given.

(A) CyclinD1 +/- EZH2i; (B) the H3K27me3 mark (the mechanism); (C) timed dosing @0/24/48h -> de-
represses whenever added (the fast-track); (D) rescue outcome (de-repression + divisions) vs dosing time.

NOTE (2026-07): the earlier "stable, EZH2-maintained lock-in" framing is SUPERSEDED. EZH2 is cell-cycle-
gated (Skp2/E2f), so it is not MYCN-sustained in a non-cycling cell; and the mark does NOT make G0
bistable -- a two-initial-condition test (fig_v44_g0_bistability_confirm) shows arrested G0 is MONOSTABLE
and reversible ("sticky"/metastable, not latched). The rescue result is unchanged: EZH2i de-represses the
MYCN-floored, mark-repressed borderline CyclinD1 and re-tips it into cycle.

(Mitogen-withdrawn GNPs are NOT shown: there EZH2 falls with cycle exit so the mark self-clears -- but by
the same cell-cycle-gating this now applies to the MB/vismo arrest too; EZH2i just accelerates it.)

Cache: fig_v44_ezh2i_rescue_kinetics_cache.npz  (--fresh). ~15 min fresh.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_ezh2i_rescue_kinetics_cache.npz'
rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_memory=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
SEL = ['time', 'Cd', 'H3K27_Cd', 'MPF']
T_POST = 9600                      # 160 h window
MB = dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002,
          k_Cd_translation=0.801 * 1.4, P21_div=0.5, EZH2i=0, HHi=0)
DOSE_MIN = {None: None, 0: 0, 24: 1440, 48: 2880, 72: 4320}   # dosing time label(h) -> minutes


def _set_tol(a, r):
    rr.integrator.setValue("absolute_tolerance", a); rr.integrator.setValue("relative_tolerance", r)


def _phase(dur):
    n = max(2, int(dur / 2.5))     # uniform 2.5-min output -> clean peak counting across phases
    try:
        return rr.simulate(0, dur, n, selections=SEL)
    except Exception:
        _set_tol(1e-6, 1e-4)
        try:
            return rr.simulate(0, dur, n, selections=SEL)
        finally:
            _set_tol(1e-9, 1e-6)


def trajectory(ezi_min):
    rr.reset()
    for k, v in MB.items():
        rr[k] = v
    rr.simulate(0, 6500, 13000)               # equilibrate cycling MB
    rr['HHi'] = 0.85                        # vismo at t=0
    if ezi_min == 0:
        rr['EZH2i'] = 1; segs = [(_phase(T_POST), 0.0)]
    elif ezi_min is None:
        segs = [(_phase(T_POST), 0.0)]
    else:
        segs = [(_phase(ezi_min), 0.0)]; rr['EZH2i'] = 1
        segs.append((_phase(T_POST - ezi_min), ezi_min))
    ts, cds, mks, mps = [], [], [], []
    for r, o in segs:
        sl = slice(1, None) if ts else slice(None)
        ts.append(r['time'][sl] + o); cds.append(r['Cd'][sl]); mks.append(r['H3K27_Cd'][sl]); mps.append(r['MPF'][sl])
    t = np.concatenate(ts) / 60.0; cd = np.concatenate(cds); mk = np.concatenate(mks); mpf = np.concatenate(mps)
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / 2.5))
    # divisions counted AFTER dosing (post-rescue re-entry)
    t0 = 0 if ezi_min in (0, None) else ezi_min / 60.0
    n_post = int(np.sum(t[pk] >= t0))
    return dict(t=t, cd=cd, mk=mk, pk=pk, n_post=n_post)


if FRESH or not os.path.exists(CACHE):
    print('computing MB vismo +/- EZH2i (timed dosing)...')
    save = {}
    for lab, mn in DOSE_MIN.items():
        tr = trajectory(mn)
        key = 'none' if lab is None else f'e{lab}'
        for f in ['t', 'cd', 'mk']:
            save[f'{key}_{f}'] = tr[f]
        save[f'{key}_pk'] = tr['pk']; save[f'{key}_npost'] = tr['n_post']
        print(f'  {key}: {len(tr["pk"])} total div, {tr["n_post"]} post-dosing')
    np.savez(CACHE, **save)
    print('cached ->', CACHE)
z = np.load(CACHE)

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(2, 2, figsize=(15, 11))
GREY, GRN = '#7f8c8d', '#1e8449'
TIMED = [('e0', 0, '#1e8449'), ('e24', 24, '#2980b9'), ('e48', 48, '#e67e22'), ('e72', 72, '#8e44ad')]

# (A) CyclinD1 +/- EZH2i (basic rescue + kinetics): vismo-only vs EZH2i@0h
a = ax[0, 0]
a.plot(z['none_t'], z['none_cd'], '-', color=GREY, lw=2.2, label='vismo only (EZH2-maintained arrest)')
a.plot(z['e0_t'], z['e0_cd'], '-', color=GRN, lw=2.2, label='vismo + EZH2i @0h')
a.axhline(z['none_cd'][z['none_t'] > 80].mean(), color=GREY, ls=':', alpha=0.5)
a.annotate('de-repression\n(~24–48 h)', xy=(36, 3.0), xytext=(60, 1.8), fontsize=9, color=GRN,
           arrowprops=dict(arrowstyle='->', color=GRN))
a.set_xlabel('time after vismo (h)', fontsize=11); a.set_ylabel('CyclinD1', fontsize=11)
a.set_title('(A) EZH2i de-represses CyclinD1 → rescues the\nvismo-arrested MB cell', fontweight='bold', fontsize=11)
a.legend(fontsize=9, loc='center right'); a.grid(alpha=0.15)

# (B) the mechanism: H3K27me3 mark
b = ax[0, 1]
b.plot(z['none_t'], z['none_mk'], '-', color=GREY, lw=2.2, label='vismo only (mark maintained)')
b.plot(z['e0_t'], z['e0_mk'], '-', color=GRN, lw=2.2, label='vismo + EZH2i (deposition blocked → decays)')
b.set_xlabel('time after vismo (h)', fontsize=11); b.set_ylabel('H3K27me3 at CyclinD1 locus', fontsize=11)
b.set_title('(B) Mechanism: EZH2i blocks H3K27me3 deposition\n→ mark cleared by demethylation + dilution', fontweight='bold', fontsize=11)
b.legend(fontsize=9, loc='upper right'); b.grid(alpha=0.15)

# (C) timed dosing -> de-represses whenever added (the fast-track)
c = ax[1, 0]
c.plot(z['none_t'], z['none_cd'], '-', color=GREY, lw=1.6, alpha=0.7, label='vismo only')
for key, h, col in TIMED:
    c.plot(z[f'{key}_t'], z[f'{key}_cd'], '-', color=col, lw=2, label=f'EZH2i @{h}h')
    c.axvline(h, color=col, ls=':', alpha=0.4)
c.set_xlabel('time after vismo (h)', fontsize=11); c.set_ylabel('CyclinD1', fontsize=11)
c.set_title('(C) Fast-track: EZH2i triggers CyclinD1 de-repression\nwhenever it is given (dotted = dosing time)', fontweight='bold', fontsize=11)
c.legend(fontsize=8.5, loc='center right'); c.grid(alpha=0.15)

# (D) rescue outcome vs dosing time
d = ax[1, 1]
hs = [h for _, h, _ in TIMED]
npost = [int(z[f'{key}_npost']) for key, _, _ in TIMED]
# final CyclinD1 (de-repression extent)
finalcd = [z[f'{key}_cd'][-1] for key, _, _ in TIMED]
d.bar(hs, npost, width=10, color=[col for _, _, col in TIMED], alpha=0.5, edgecolor='k')
for h, n in zip(hs, npost):
    d.text(h, n + 0.05, str(n), ha='center', fontweight='bold', fontsize=10)
d.axhline(0, color=GREY, lw=1)
d.text(72, 0.1, 'vismo only: 0 div (no rescue)', fontsize=8, color=GREY, ha='center')
d2 = d.twinx()
d2.plot(hs, finalcd, '-o', color='#c0392b', lw=2, ms=6, label='final CyclinD1 (de-repressed)')
d2.set_ylabel('final CyclinD1', color='#c0392b', fontsize=10); d2.set_ylim(0, 4)
d.set_xlabel('EZH2i dosing time after vismo (h)', fontsize=11)
d.set_ylabel('post-dosing divisions', fontsize=11)
d.set_title('(D) Rescue is available on-demand\n(de-repression completes whenever EZH2i is given)', fontweight='bold', fontsize=11)
d2.legend(fontsize=8, loc='upper right'); d.grid(alpha=0.15, axis='y')

plt.tight_layout()
plt.savefig('simulations/fig_v44_ezh2i_rescue_kinetics.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ezh2i_rescue_kinetics.pdf', bbox_inches='tight')
plt.close()
print('Saved: fig_v44_ezh2i_rescue_kinetics.png / .pdf')
