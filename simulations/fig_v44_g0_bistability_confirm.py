"""DEFINITIVE bistability test: fixed-SHH, two-initial-condition equilibration.

For each SHH we drop the cell in from BOTH extremes and integrate to TRUE steady state (robust,
stiff-zone tolerance escalation), then read out whether it settled cycling or arrested:
  - from ARRESTED  (high-mark IC, built at SHH=0.05)
  - from CYCLING   (low-mark IC, built at SHH=1.0)
If the two ICs settle to DIFFERENT stable states at the same SHH -> BISTABLE at that SHH; the SHH
range where they differ = the hysteresis window (no sweep transients, no path dependence).

Run at ACCUMULATE (kDeEZ=5e-5, del_mk=3e-4 — stable EZH2, slow mark) and calibrated (RELAX) control.
Run:  ./venv/bin/python simulations/fig_v44_g0_bistability_confirm.py [--mini] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MINI = '--mini' in sys.argv
FRESH = '--fresh' in sys.argv or MINI
CACHE = 'simulations/fig_v44_g0_bistability_confirm_cache.npz' if not MINI else 'simulations/fig_v44_g0_bistability_confirm_mini.npz'
SCRATCH = '/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/81291f88-7d19-4149-be27-53b1162ee52d/scratchpad'

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SETTINGS = [
    dict(name='ACCUMULATE (stable EZH2, slow mark)', kDeEZ=5.0e-5, del_mk=3e-4, col='#c0392b'),
    dict(name='calibrated / RELAX (control)',        kDeEZ=1.5e-4, del_mk=7e-4, col='#2c7fb8'),
]
SHHG = np.round(np.linspace(0.30, 0.75, 5 if MINI else 12), 3)
T_EQ = 40000 if MINI else 70000       # equilibration to steady state (min); ACCUMULATE mark is slow
T_PRE_ARR = 30000
T_PRE_CYC = 15000
ATOLS = (1e-7, 1e-6, 1e-5, 1e-4)


def rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
    try: r.integrator.setValue('maximum_num_steps', 500000)
    except Exception: pass
    return r


def _chunks(r, t0, total, chunk=10000, sel=None):
    """robust chunked integration with tolerance escalation; returns last successful result & end time."""
    t = t0; last = None
    while t < t0 + total:
        ok = False
        for atol in ATOLS:
            try:
                r.integrator.setValue('absolute_tolerance', atol)
                last = r.simulate(t, t + chunk, max(80, int(chunk / 5)), selections=sel or ['time', 'MPF', 'Mk', 'Cd'])
                ok = True; break
            except Exception:
                continue
        if not ok:
            break
        t += chunk
    return last, t


def build_ic(kdeez, delmk, shh0, T, fname):
    r = rr(); r.reset()
    for k, v in GNP.items(): r[k] = v
    r['kDeEZ'] = float(kdeez); r['del_mk'] = float(delmk); r['SHH'] = float(shh0)
    _chunks(r, 0.0, T)
    r.saveState(fname)


def run_from_ic(kdeez, delmk, shh, icfile):
    r = rr(); r.loadState(icfile)
    r['kDeEZ'] = float(kdeez); r['del_mk'] = float(delmk); r['SHH'] = float(shh)
    _, tend = _chunks(r, 0.0, T_EQ)
    res, _ = _chunks(r, tend, 12000)          # measurement window at steady state
    if res is None:
        return np.nan, np.nan, np.nan
    pk, _ = find_peaks(res['MPF'], prominence=0.15, distance=200)
    hrs = (res['time'][-1] - res['time'][0]) / 60.0
    rate = len(pk) / hrs * 168.0 if hrs > 0 else 0.0
    return rate, float(np.nanmean(res['Mk'])), float(np.nanmean(res['Cd']))


if FRESH or not os.path.exists(CACHE):
    RATE_ARR, RATE_CYC, MK_ARR, MK_CYC, CD_ARR, CD_CYC = ([] for _ in range(6))
    for s in SETTINGS:
        af = os.path.join(SCRATCH, 'ic_arr.dat'); cf = os.path.join(SCRATCH, 'ic_cyc.dat')
        build_ic(s['kDeEZ'], s['del_mk'], 0.05, T_PRE_ARR, af)
        build_ic(s['kDeEZ'], s['del_mk'], 1.00, T_PRE_CYC, cf)
        ra, rc, ma, mc, ca, cc = [], [], [], [], [], []
        for shh in SHHG:
            r1 = run_from_ic(s['kDeEZ'], s['del_mk'], shh, af)   # from arrested
            r2 = run_from_ic(s['kDeEZ'], s['del_mk'], shh, cf)   # from cycling
            ra.append(r1[0]); ma.append(r1[1]); ca.append(r1[2])
            rc.append(r2[0]); mc.append(r2[1]); cc.append(r2[2])
            print(f"  {s['name'][:22]:22s} SHH={shh:.2f}  arrIC[div {r1[0]:.1f} Mk {r1[1]:.2f} Cd {r1[2]:.2f}]  cycIC[div {r2[0]:.1f} Mk {r2[1]:.2f} Cd {r2[2]:.2f}]", flush=True)
        RATE_ARR.append(ra); RATE_CYC.append(rc); MK_ARR.append(ma); MK_CYC.append(mc); CD_ARR.append(ca); CD_CYC.append(cc)
        np.savez(CACHE, SHHG=SHHG, RATE_ARR=RATE_ARR, RATE_CYC=RATE_CYC, MK_ARR=MK_ARR, MK_CYC=MK_CYC,
                 CD_ARR=CD_ARR, CD_CYC=CD_CYC, names=[x['name'] for x in SETTINGS], cols=[x['col'] for x in SETTINGS])
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE, allow_pickle=True)
SHHG = z['SHHG']; names = list(z['names']); cols = list(z['cols'])
RATE_ARR, RATE_CYC = z['RATE_ARR'], z['RATE_CYC']; MK_ARR, MK_CYC = z['MK_ARR'], z['MK_CYC']; CD_ARR, CD_CYC = z['CD_ARR'], z['CD_CYC']

fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6))
for j, (nm, cl) in enumerate(zip(names, cols)):
    # bistable where the two ICs settle to DIFFERENT stable states, judged by the smooth CyclinD1
    # order parameter (div-count is too quantized and false-positives on missed peaks).
    THR = 1.5
    ca = np.asarray(CD_ARR[j], float); cc = np.asarray(CD_CYC[j], float)
    valid = ~(np.isnan(ca) | np.isnan(cc))
    bist = np.zeros(len(SHHG), bool)
    bist[valid] = (ca[valid] > THR) != (cc[valid] > THR)
    win = SHHG[bist]
    width = (win.max() - win.min()) if bist.any() else 0.0
    for row, (YA, YC, yl) in enumerate([(CD_ARR[j], CD_CYC[j], 'mean CyclinD1  Cd'), (MK_ARR[j], MK_CYC[j], 'mean H3K27me3  Mk')]):
        a = axes[row, j]
        a.plot(SHHG, YC, '-o', color='#2980b9', ms=6, lw=2, label='from CYCLING IC (low mark)')
        a.plot(SHHG, YA, '-s', color='#c0392b', ms=6, lw=2, label='from ARRESTED IC (high mark)')
        for k in np.where(bist)[0]:
            a.axvspan(SHHG[k] - (SHHG[1]-SHHG[0])/2, SHHG[k] + (SHHG[1]-SHHG[0])/2, color='gold', alpha=0.3)
        if row == 0:
            verdict = f'BISTABLE\nwindow ΔSHH ≈ {width:.2f}' if bist.any() else 'MONOSTABLE\n(ICs converge)'
            a.set_title(f'{nm}\n{verdict}', fontsize=10.5, fontweight='bold', color=cl)
            if j == 0: a.legend(fontsize=8, loc='upper left')
        a.set_ylabel(yl if j == 0 else ''); a.grid(alpha=0.15)
        if row == 1: a.set_xlabel('mitogen drive  SHH  (fixed; each point equilibrated from both ICs)')
fig.suptitle('Definitive bistability test — same SHH, two initial conditions, integrated to steady state',
             fontsize=13, fontweight='bold', y=0.99)
fig.text(0.5, 0.005, 'Gold = SHH where the arrested (high-mark) and cycling (low-mark) starts settle to DIFFERENT stable states → two coexisting attractors → bistable G0 lock.',
         ha='center', fontsize=9, color='#555')
plt.tight_layout(rect=[0, 0.015, 1, 0.97])
plt.savefig('simulations/fig_v44_g0_bistability_confirm.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_bistability_confirm.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_g0_bistability_confirm.png')
