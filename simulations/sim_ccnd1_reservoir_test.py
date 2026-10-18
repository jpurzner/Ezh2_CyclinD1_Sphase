"""Does the MEASURED CyclinD1 mRNA half-life (3.4 h) still carry the Hh-withdrawal memory, if the
reservoir is large? (JP 2026-10-16). Memory is NOT in a slow half-life; it is a large CyclinD1 pool
draining past the commitment threshold. Coast time ~ tau * ln(pool/threshold). The model already has a
big reservoir (MB ~21x, GNP ~5x above K_CdRb). Three arms, SAME steady level (tx re-scaled to the half-
life so only the DRAIN TIME differs): slow 11 h (current) / measured 3.4 h / instantaneous (52 s). After
Hh-off (HHi=1) count post-withdrawal divisions = the coast. Run: venv/bin/python simulations/sim_ccnd1_reservoir_test.py
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
UPH = 62.8; LN2 = np.log(2); KCDRB = 0.4025; THR_M = KCDRB / 0.801   # Cd=0.801*Cd_mRNA -> mRNA threshold ~0.50
def divs(d): d = np.asarray(d); return np.where((d[:-1] > 0.9) & (d[1:] < 0.1))[0]
def smooth(x, t, w_h=8.0):
    dt = np.median(np.diff(t)); k = max(1, int(w_h / dt))
    return np.convolve(np.pad(x, k, mode='edge'), np.ones(k) / k, mode='same')[k:-k]

base = te.loada(build_model_v44())
TX0 = {p: base[p] for p in ['k_Cd_tx_basal', 'k_Cd_tx_Gli_max', 'k_Cd_tx_MYCN']}
KD0 = base['k_Cd_mRNA_deg']                        # 0.001 (t1/2 ~11 h)

def arm_over(t_half_h):
    k = LN2 / (t_half_h * UPH)                     # per model-time-unit
    s = k / KD0                                    # scale tx by the same ratio -> hold steady level fixed
    ov = {'k_Cd_mRNA_deg': k}
    ov.update({p: TX0[p] * s for p in TX0})
    return ov

ARMS = [('slow 11 h (current)', 11.0, '#0e6251'),
        ('measured 3.4 h', 3.4, '#c0392b'),
        ('instant (52 s)', 0.0138, '#e67e22')]
CELLS = {'GNP': dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0),
         'MB':  dict(Ptch1_copy_number=0.3)}
ONSET = 130.0; POST = 60.0

def run(cell_over, arm_thalf):
    rr = te.loada(build_model_v44(params=arm_over(arm_thalf)))
    rr.reset()
    for k, v in cell_over.items(): rr[k] = v
    if cell_over.get('Ptch1_copy_number', 1.0) < 1.0:        # MB: engage CyclinD2 expression
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
    rr.integrator.setValue('absolute_tolerance', 1e-7)
    T, CD, DNA = [], [], []; t0 = 0.0
    for ov, dur in [(dict(SHH=0.5), ONSET), (dict(HHi=1), POST)]:
        for k, v in ov.items(): rr[k] = v
        r = rr.simulate(t0, t0 + dur * UPH, max(2, int(dur * UPH / 4)), selections=['time', 'Cd_mRNA', 'Dna'])
        T += list(r['time']); CD += list(r['Cd_mRNA']); DNA += list(r['Dna']); t0 += dur * UPH
    t = np.array(T) / UPH; cd = np.array(CD); dna = np.array(DNA)
    cds = smooth(cd, t, 8.0)                                       # cycle-averaged reservoir envelope
    dvt = t[divs(dna)]; ndiv = int(np.sum(dvt > ONSET + 0.5))
    # reservoir drain: time after withdrawal for the smoothed Cd_mRNA to fall below the commitment threshold
    post = t > ONSET
    below = post & (cds < THR_M)
    tcross = (t[below][0] - ONSET) if below.any() else np.nan
    return t, cd, cds, dvt, ndiv, tcross

fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharex=True)
for ax, (cname, cover) in zip(axes, CELLS.items()):
    ratio = (8.34 if cname == 'MB' else 2.06) / KCDRB
    for (lab, th, col) in ARMS:
        try:
            t, cd, cds, dvt, ndiv, tcross = run(cover, th)
        except Exception as e:
            print(f"  {cname} {lab}: FAILED {str(e)[:60]}", flush=True); continue
        crs = f"{tcross:.0f} h" if tcross == tcross else ">60 h"
        print(f"  {cname:3s} {lab:20s}: Cd-below-threshold at {crs:>6}  (post-div={ndiv})", flush=True)
        m = t > ONSET - 45
        ax.plot(t[m], cds[m], color=col, lw=2.4, label=f'{lab}: crosses in {crs}')
        ax.plot(t[m], cd[m], color=col, lw=0.6, alpha=0.3)        # faint raw (cycle sawtooth)
        if tcross == tcross:
            ax.plot(ONSET + tcross, THR_M, 'o', color=col, ms=7, mec='k', mew=0.5, zorder=6)
    ax.axvline(ONSET, color='#333', ls=':', lw=1.3)
    ax.axhline(THR_M, color='#888', ls='--', lw=1.2); ax.text(ONSET - 43, THR_M * 1.15, 'commitment threshold', fontsize=8, color='#666')
    ax.text(ONSET + 1, ax.get_ylim()[1] * 0.95, 'Hh OFF', fontsize=8.5, rotation=90, va='top', color='#333')
    pred = 3.4 / LN2 * np.log(ratio)                             # tau*ln(ratio), tau = t_half/ln2, measured t_half 3.4h
    ax.set_title(f"{cname}: CyclinD1 reservoir {ratio:.0f}x above threshold\n"
                 f"measured-t1/2 drain predicted tau*ln({ratio:.0f}) = {pred:.0f} h", fontsize=10.5, fontweight='bold')
    ax.set_xlabel('time (h)'); ax.set_ylabel('CyclinD1 mRNA (reservoir)'); ax.legend(fontsize=8.6, loc='upper right'); ax.grid(alpha=0.2)
    ax.set_xlim(ONSET - 45, ONSET + POST); ax.set_ylim(0, None)
fig.suptitle("Accept the measured CyclinD1 mRNA half-life (3.4 h): the Hh-withdrawal memory survives as a LARGE-RESERVOIR drain\n"
             "(same reservoir level in all arms; only the drain time differs). Slow 11 h OVERSHOOTS; measured 3.4 h coasts ~1 cycle; instant arrests at once.",
             fontsize=11.5, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.9])
out = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/fig_ccnd1_reservoir_test.png"
plt.savefig(out, dpi=135); print('saved', out)
