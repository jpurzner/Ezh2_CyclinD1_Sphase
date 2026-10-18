"""WS2 (JP overnight 2026-07-30): how mitogen (Shh) dose sets G1 duration in GNP vs MB, baked dynamic-CDKI.

Mechanism: Shh -> Gli -> CyclinD1(Cd). Higher Cd -> more p27 buffered on CyclinD1-CDK4/6 (CdP21) -> free p27
low -> CDK2 fires early -> SHORT G1. As Shh falls, Cd falls -> weaker buffering + the latent mu_eff G1-lengthening
(K_g1len) engages -> free p27 stays high longer -> G1 LENGTHENS graded -> eventually G0/arrest. MB (high INK4 +
Ptch-driven baseline Hh) sits at longer G1 for the same Shh and arrests at a higher Shh floor than GNP.

Outputs: fig_ws2_mitogen_g1.png/.pdf + console summary.  Run: ./venv/bin/python simulations/sim_ws2_mitogen_g1.py
"""
import os, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

M = build_model_v44()   # baked dynamic-CDKI default
SEL = ['time', 'Dna', 'aRc', 'pRb', 'P21', 'CdP21', 'Cd', 'Gli1']
GNP = dict(Ptch1_copy_number=1.0, p16=0.0, p18=0.464, p19=0.36, kSyP21=0.00116, P21_div=0.0, MYCN_amplification=1.0)
MB  = dict(Ptch1_copy_number=0.1, p16=0.306, p18=1.73, p19=0.58, kSyP21=0.002, P21_div=0.6, MYCN_amplification=2.8)
SHH_VALS = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0]

def run(cond, shh, T=30000, npts=12000):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            r = te.loada(M); r.reset()
            for k, v in cond.items():
                try: r[k] = v
                except Exception: pass
            r['SHH'] = shh
            r.integrator.setValue('absolute_tolerance', atol); r.integrator.setValue('relative_tolerance', 1e-6)
            r.integrator.setValue('maximum_num_steps', 4000000)
            res = r.simulate(0, T, npts, selections=SEL)
            return {k: np.asarray(res[k]) for k in SEL}
        except Exception: continue
    return None

def g1_metrics(tr, settle=8000):
    if tr is None: return None
    t = tr['time']; d = tr['Dna']; aRc = tr['aRc']
    div_t = t[np.where((d[:-1] > 0.9) & (d[1:] < 0.1))[0]]; div_t = div_t[div_t >= settle]
    s_t = t[np.where((aRc[:-1] < 0.05) & (aRc[1:] >= 0.05))[0]]; s_t = s_t[s_t >= settle]  # aRc peaks ~0.1, not 0.5
    ncyc = len(div_t); cyc = float(np.median(np.diff(div_t))/60) if ncyc >= 2 else np.nan
    dwell = [ (s_t[s_t > dt][0]-dt)/60 for dt in div_t[:-1] if len(s_t[s_t > dt]) ]
    m = t >= settle
    cd = float(np.mean(tr['Cd'][m])); p27 = float(np.mean(tr['P21'][m])); cdp = float(np.mean(tr['CdP21'][m]))
    return dict(ncyc=ncyc, cyc=cyc, g1=float(np.mean(dwell)) if dwell else np.nan, cd=cd, p27=p27, cdp=cdp)

print("=== WS2: mitogen(Shh) -> G1 duration, GNP vs MB (baked dynamic-CDKI) ===")
res = {'GNP': [], 'MB': []}
for lab, base in [('GNP', GNP), ('MB', MB)]:
    for shh in SHH_VALS:
        mm = g1_metrics(run(base, shh))
        res[lab].append(mm)
        s = 'FAIL' if mm is None else f"div={mm['ncyc']:2d} cyc={mm['cyc'] if not np.isnan(mm['cyc']) else 0:5.1f}h G1={mm['g1'] if not np.isnan(mm['g1']) else 0:5.1f}h Cd={mm['cd']:.2f} p27={mm['p27']:.3f}"
        print(f"  {lab} SHH={shh:.2f}: {s}")

fig = plt.figure(figsize=(15, 9)); gs = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.26, left=0.07, right=0.97, top=0.9, bottom=0.08)
C_G, C_M = '#2c7fb8', '#c0392b'
def series(lab, key): return [ (m[key] if (m and not (isinstance(m[key], float) and np.isnan(m[key]))) else np.nan) for m in res[lab] ]
# A: G1 duration vs Shh
axA = fig.add_subplot(gs[0,0])
for lab, col in [('GNP', C_G), ('MB', C_M)]:
    axA.plot(SHH_VALS, series(lab,'g1'), 'o-', color=col, lw=1.7, label=lab)
axA.set_title('A  G1 duration vs Shh dose — graded lengthening as mitogen falls', fontsize=10, fontweight='bold')
axA.set_xlabel('Shh dose'); axA.set_ylabel('G1 duration (birth→S, h)'); axA.legend(); axA.grid(alpha=0.2)
# B: divisions (cycling vs arrest) vs Shh
axB = fig.add_subplot(gs[0,1])
for lab, col in [('GNP', C_G), ('MB', C_M)]:
    axB.plot(SHH_VALS, series(lab,'ncyc'), 'o-', color=col, lw=1.7, label=lab)
axB.set_title('B  proliferation vs Shh — MB arrests at a higher Shh floor', fontsize=10, fontweight='bold')
axB.set_xlabel('Shh dose'); axB.set_ylabel('divisions in window'); axB.legend(); axB.grid(alpha=0.2)
# C: cycle length vs Shh
axC = fig.add_subplot(gs[1,0])
for lab, col in [('GNP', C_G), ('MB', C_M)]:
    axC.plot(SHH_VALS, series(lab,'cyc'), 'o-', color=col, lw=1.7, label=lab)
axC.set_title('C  cycle length vs Shh', fontsize=10, fontweight='bold')
axC.set_xlabel('Shh dose'); axC.set_ylabel('cycle length (h)'); axC.legend(); axC.grid(alpha=0.2)
# D: mechanism -- CyclinD1 (Cd) and free p27 vs Shh (GNP)
axD = fig.add_subplot(gs[1,1])
axD.plot(SHH_VALS, series('GNP','cd'), 'o-', color='#16a085', lw=1.6, label='CyclinD1 (Cd), GNP')
axD.plot(SHH_VALS, series('MB','cd'), 's--', color='#16a085', lw=1.2, alpha=0.7, label='Cd, MB')
ax2 = axD.twinx()
ax2.plot(SHH_VALS, series('GNP','p27'), 'o-', color='#c0392b', lw=1.6, label='free p27, GNP')
ax2.plot(SHH_VALS, series('MB','p27'), 's--', color='#c0392b', lw=1.2, alpha=0.7, label='free p27, MB')
ax2.set_ylabel('free p27', color='#c0392b')
axD.set_title('D  mechanism: Cd rises with Shh → free p27 falls → G1 shortens', fontsize=9.5, fontweight='bold')
axD.set_xlabel('Shh dose'); axD.set_ylabel('CyclinD1 (Cd)', color='#16a085')
l1,la1=axD.get_legend_handles_labels(); l2,la2=ax2.get_legend_handles_labels(); axD.legend(l1+l2, la1+la2, fontsize=7, loc='center right')
axD.grid(alpha=0.2)
fig.suptitle('WS2 — Mitogen(Shh) dose sets G1 duration: graded lengthening as Shh falls; MB longer + higher arrest floor (baked dynamic-CDKI)',
             fontsize=12, fontweight='bold')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_ws2_mitogen_g1')
fig.savefig(out + '.png', dpi=130, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
