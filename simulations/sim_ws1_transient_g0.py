"""WS1 (JP overnight 2026-07-30): detailed transient-G0 exploration on the BAKED dynamic-CDKI model.
Why MB shows transient G0 (p27-high / pRb-low pre-S dwell + re-entry) and GNP does not.

Mechanism (CDKI model): commitment is the CdP21 buffer -- p27 is sequestered on CyclinD1-CDK4/6 when Cd is
high (-> free p27 low -> CDK2 fires -> commit). MB has high INK4 (p16/p18/p19) that raises the CDK4/6 half-max
AND a birth-p27 (P21_div) so daughters are born p27-high; low Cd + high INK4 -> weak buffering -> free p27 stays
high -> pre-S pRb-low DWELL (transient G0) until Cd/CdP21 clears it -> re-entry. GNP: low INK4, no birth-p27 ->
free p27 cleared fast -> no G0.

Outputs: fig_ws1_transient_g0.png/.pdf + a console summary.
Run: ./venv/bin/python simulations/sim_ws1_transient_g0.py
"""
import os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

PRB_G0 = 1.5   # pre-S pRb(hyper) below this = hypophospho / G0 (validate convention)
M = build_model_v44()   # baked default = dynamic-CDKI
SEL = ['time', 'Dna', 'aRc', 'pRb', 'P21', 'CdP21', 'p18_prot', 'p19_prot', 'Cd', 'EZH2', 'MPF']
# cell-type operating points (baked tones)
GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, p19=0.36, kSyP21=0.00116, P21_div=0.0, MYCN_amplification=1.0)
MB  = dict(SHH=0.5, Ptch1_copy_number=0.1, p16=0.306, p18=1.73, p19=0.58, kSyP21=0.002, P21_div=0.6, MYCN_amplification=2.8)

def run(cond, T=40000, npts=16000):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            r = te.loada(M); r.reset()
            for k, v in cond.items():
                try: r[k] = v
                except Exception: pass
            r.integrator.setValue('absolute_tolerance', atol); r.integrator.setValue('relative_tolerance', 1e-6)
            r.integrator.setValue('maximum_num_steps', 4000000)
            res = r.simulate(0, T, npts, selections=SEL)
            return {k: np.asarray(res[k]) for k in SEL}
        except Exception: continue
    return None

def metrics(tr, settle=8000):
    """divisions, mean cycle length, per-birth G1/G0 dwell (birth->S entry), G0 time-fraction (pre-S & pRb<thr)."""
    if tr is None: return None
    t = tr['time']; d = tr['Dna']; aRc = tr['aRc']; pRb = tr['pRb']
    div_i = np.where((d[:-1] > 0.9) & (d[1:] < 0.1))[0]; div_t = t[div_i]
    div_t = div_t[div_t >= settle]
    # S-phase entry = aRc crosses up through 0.05 (aRc peaks ~0.1 in this model, not 0.5)
    s_i = np.where((aRc[:-1] < 0.05) & (aRc[1:] >= 0.05))[0]; s_t = t[s_i]; s_t = s_t[s_t >= settle]
    ncyc = len(div_t); cyc = float(np.median(np.diff(div_t))/60) if ncyc >= 2 else np.nan
    # G1/G0 dwell = time from each division to the next S entry
    dwell = []
    for dt in div_t[:-1]:
        nxt = s_t[s_t > dt]
        if len(nxt): dwell.append((nxt[0]-dt)/60)
    dwell = np.array(dwell)
    m = t >= settle
    preS = (aRc[m] < 0.05) & (d[m] < 0.9)   # pre-replication (not in S, not G2)
    g0frac = float(np.mean(preS & (pRb[m] < PRB_G0)))   # pre-S AND pRb-low = discrete G0 (NB pRb rarely dips this low in CDKI model)
    return dict(ncyc=ncyc, cyc=cyc, dwell=dwell, g0frac=g0frac,
                dwell_mean=float(np.mean(dwell)) if len(dwell) else np.nan,
                dwell_max=float(np.max(dwell)) if len(dwell) else np.nan)

print("=== WS1: transient G0, MB vs GNP (baked dynamic-CDKI) ===")
trG, trM = run(GNP), run(MB)
mG, mM = metrics(trG), metrics(trM)
for lab, mm in [('GNP', mG), ('MB', mM)]:
    if mm: print(f"  {lab}: {mm['ncyc']} div, cycle~{mm['cyc']:.1f}h, G1/G0 dwell mean={mm['dwell_mean']:.1f}h max={mm['dwell_max']:.1f}h, G0 time-fraction={mm['g0frac']*100:.1f}%")

# ---- birth-p27 sweep: how transient G0 emerges with birth-p27, GNP vs MB ----
p27_vals = [0.0, 0.3, 0.6, 1.0, 1.5, 2.2, 3.0]
sweep = {'GNP': [], 'MB': []}
for lab, base in [('GNP', GNP), ('MB', MB)]:
    for p in p27_vals:
        c = dict(base); c['P21_div'] = p
        mm = metrics(run(c, T=30000, npts=12000))
        sweep[lab].append(mm)
        if mm is None:
            print(f"  {lab} birth-p27={p}: FAIL")
        else:
            print(f"  {lab} birth-p27={p}: div={mm['ncyc']}, dwell_mean={mm['dwell_mean']:.1f}h, G0frac={mm['g0frac']*100:.0f}%")

# ---- figure ----
fig = plt.figure(figsize=(15.5, 10))
gs = GridSpec(3, 2, figure=fig, hspace=0.42, wspace=0.26, left=0.06, right=0.97, top=0.92, bottom=0.07)
C_G, C_M = '#2c7fb8', '#c0392b'
def hrs(t): return t/60.0
# A: GNP lineage (pRb + Dna) ; B: MB lineage (pRb + Dna) -- show MB pre-S pRb-low dwells
for ax, tr, lab, col in [(fig.add_subplot(gs[0,0]), trG, 'GNP', C_G), (fig.add_subplot(gs[0,1]), trM, 'MB', C_M)]:
    if tr is None: continue
    win = (tr['time'] >= 8000) & (tr['time'] <= 24000)
    ax.plot(hrs(tr['time'][win]), tr['pRb'][win], color=col, lw=1.3, label='pRb (hyper)')
    ax.axhline(PRB_G0, color='k', ls=':', lw=0.8); ax.fill_between(hrs(tr['time'][win]), 0, PRB_G0, color='#eee', zorder=0)
    ax2 = ax.twinx(); ax2.plot(hrs(tr['time'][win]), tr['Dna'][win], color='#999', lw=0.8, alpha=0.7)
    ax2.set_ylabel('Dna', color='#999', fontsize=8)
    ax.set_title(f'{lab}: pRb over time (grey band = G0, pRb<{PRB_G0})', fontsize=10, fontweight='bold', color=col)
    ax.set_xlabel('hours'); ax.set_ylabel('pRb (hyper)')
# C: mechanism during a window -- free p27, CdP21, INK4 (MB)
axC = fig.add_subplot(gs[1,0])
if trM is not None:
    win = (trM['time'] >= 8000) & (trM['time'] <= 24000)
    axC.plot(hrs(trM['time'][win]), trM['P21'][win], color='#c0392b', lw=1.2, label='free p27')
    axC.plot(hrs(trM['time'][win]), trM['CdP21'][win], color='#e67e22', lw=1.2, label='CdP21 (buffered p27)')
    axC.plot(hrs(trM['time'][win]), trM['p18_prot'][win]+trM['p19_prot'][win], color='#6c3483', lw=1.0, ls='--', label='INK4 (p18+p19)')
    axC.plot(hrs(trM['time'][win]), trM['Cd'][win], color='#16a085', lw=1.0, ls=':', label='CyclinD1 (Cd)')
axC.set_title('MB mechanism: free p27 high in the G0 dwell; CdP21 buffer + Cd clear it -> re-entry', fontsize=9.5, fontweight='bold')
axC.set_xlabel('hours'); axC.legend(fontsize=7); axC.grid(alpha=0.2)
# D: birth-p27 sweep -- G0 fraction GNP vs MB
axD = fig.add_subplot(gs[1,1])
for lab, col in [('GNP', C_G), ('MB', C_M)]:
    y = [ (m['g0frac']*100 if m else np.nan) for m in sweep[lab] ]
    axD.plot(p27_vals, y, 'o-', color=col, lw=1.6, label=lab)
axD.set_title('G0 time-fraction vs birth-p27: MB enters transient G0, GNP resists', fontsize=9.5, fontweight='bold')
axD.set_xlabel('birth-p27 (P21_div)'); axD.set_ylabel('G0 time-fraction (%)'); axD.legend(fontsize=8); axD.grid(alpha=0.2)
# E: birth-p27 sweep -- mean G1/G0 dwell
axE = fig.add_subplot(gs[2,0])
for lab, col in [('GNP', C_G), ('MB', C_M)]:
    y = [ (m['dwell_mean'] if m else np.nan) for m in sweep[lab] ]
    axE.plot(p27_vals, y, 'o-', color=col, lw=1.6, label=lab)
axE.set_title('mean G1/G0 dwell (birth->S) vs birth-p27', fontsize=9.5, fontweight='bold')
axE.set_xlabel('birth-p27 (P21_div)'); axE.set_ylabel('dwell (h)'); axE.legend(fontsize=8); axE.grid(alpha=0.2)
# F: dwell distributions (GNP vs MB at their operating points)
axF = fig.add_subplot(gs[2,1])
if mG and len(mG['dwell']): axF.hist(mG['dwell'], bins=15, color=C_G, alpha=0.6, label=f'GNP (n={len(mG["dwell"])})')
if mM and len(mM['dwell']): axF.hist(mM['dwell'], bins=15, color=C_M, alpha=0.6, label=f'MB (n={len(mM["dwell"])})')
axF.set_title('G1/G0 dwell distribution at operating point', fontsize=9.5, fontweight='bold')
axF.set_xlabel('dwell (h)'); axF.set_ylabel('count'); axF.legend(fontsize=8); axF.grid(alpha=0.2)
fig.suptitle('WS1 — Transient G0: MB (high INK4 + birth-p27) dwells pre-S / re-enters; GNP cycles continuously (baked dynamic-CDKI)',
             fontsize=12.5, fontweight='bold')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_ws1_transient_g0')
fig.savefig(out + '.png', dpi=130, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
