"""Probe: CyclinD1 de-repression KINETICS under vismo +/- EZH2i (does EZH2i fast-track it?).

MB cell equilibrated cycling (GDC=0), then vismo (GDC=0.85) applied at t=0. Arms:
  - vismo only          : feedback intact -> CyclinD1 held down by EZH2; de-repression only via slow EZH2 turnover
  - vismo + EZH2i @0h    : acute relief of repression -> fast CyclinD1 de-repression
  - vismo + EZH2i @40h   : delayed inhibitor -> CyclinD1 fast-tracks whenever added
Tracks CyclinD1(t), EZH2(t), divisions. Robust chunked integration (tol-retry).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MEM = '--instant' not in sys.argv     # default: H3K27me3-memory model (realistic kinetics)
rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, with_h3k27_memory=MEM))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
GDC_VISMO = 0.85
T_POST, CH = 7200, 100
_SEL = ['time', 'Cd', 'EZH2', 'MPF'] + (['H3K27_Cd'] if MEM else [])


def _set_tol(a, r):
    rr.integrator.setValue("absolute_tolerance", a); rr.integrator.setValue("relative_tolerance", r)


def _chunk(ch):
    try:
        return rr.simulate(0, ch, 20, selections=_SEL)
    except Exception:
        try:
            _set_tol(1e-6, 1e-4); out = rr.simulate(0, ch, 20, selections=_SEL); _set_tol(1e-9, 1e-6); return out
        except Exception:
            _set_tol(1e-9, 1e-6); return None


def _phase(dur):
    """one continuation phase (start=0 continues from current state); tol-retry on stiffness."""
    npts = max(2, int(dur / 5))
    try:
        return rr.simulate(0, dur, npts, selections=_SEL)
    except Exception:
        _set_tol(1e-6, 1e-4)
        try:
            return rr.simulate(0, dur, npts, selections=_SEL)
        finally:
            _set_tol(1e-9, 1e-6)


def run(ezi_time):
    rr.reset()
    rr['SHH'] = 0.5; rr['MYCN_amplification'] = 2.8; rr['Ptch1_copy_number'] = 0.1
    rr['p16'] = 0.306; rr['p18'] = 1.553; rr['kSyP21'] = 0.002
    rr['k_Cd_translation'] = 0.801 * 1.5; rr['P21_div'] = 0.5
    rr['EZH2i'] = 0; rr['GDC0449'] = 0
    rr.simulate(0, 6500, 13000)            # equilibrate cycling MB
    rr['GDC0449'] = GDC_VISMO              # apply vismo at t=0
    segs, off = [], 0.0
    if ezi_time == 0:
        rr['EZH2i'] = 1
        segs.append((_phase(T_POST), 0.0))
    elif ezi_time is None:
        segs.append((_phase(T_POST), 0.0))
    else:
        segs.append((_phase(ezi_time), 0.0))             # vismo-only phase
        rr['EZH2i'] = 1                                   # add inhibitor (persists for the next single simulate)
        segs.append((_phase(T_POST - ezi_time), ezi_time))  # EZH2i phase, continues from arrested state
    ts, cds, ezs, mps, mks = [], [], [], [], []
    for r, o in segs:
        sl = slice(1, None) if ts else slice(None)
        ts.append(r['time'][sl] + o); cds.append(r['Cd'][sl]); ezs.append(r['EZH2'][sl]); mps.append(r['MPF'][sl])
        mks.append(r['H3K27_Cd'][sl] if MEM else r['EZH2'][sl])
    t = np.concatenate(ts) / 60.0; cd = np.concatenate(cds); ez = np.concatenate(ezs); mp = np.concatenate(mps); mk = np.concatenate(mks)
    pk, _ = find_peaks(mp, prominence=0.15, distance=int(200 / (t[1] * 60 - t[0] * 60)))
    return t, cd, ez, mk, pk


ARMS = [('vismo only', None, '#8e44ad'), ('vismo + EZH2i @0h', 0.0, '#27ae60'), ('vismo + EZH2i @40h', 2400, '#e67e22')]
res = {}
mklab = 'H3K27me3 mark' if MEM else 'EZH2'
print(f"model: {'H3K27me3-MEMORY (kinetic de-repression)' if MEM else 'INSTANT (algebraic)'}")
print('arm                | CyclinD1 at t=  0h    5h   10h   20h   40h   60h  100h | divisions')
for name, ezt, col in ARMS:
    t, cd, ez, mk, pk = run(ezt)
    res[name] = (t, cd, ez, mk, pk, col)

    def at(h):
        return cd[np.argmin(np.abs(t - h))]
    print(f'{name:18s}| {at(0):14.2f} {at(5):5.2f} {at(10):5.2f} {at(20):5.2f} {at(40):5.2f} {at(60):5.2f} {at(100):5.2f} | {len(pk)}')

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for name, (t, cd, ez, mk, pk, col) in res.items():
    ax[0].plot(t, cd, '-', color=col, lw=1.8, label=name)
    if len(pk):
        ax[0].plot(t[pk], cd[pk], 'v', color=col, ms=6)
    ax[1].plot(t, mk, '-', color=col, lw=1.8, label=name)
ax[0].axvline(40, color='#e67e22', ls=':', alpha=0.5)
ax[0].set_xlabel('time after vismo (h)'); ax[0].set_ylabel('CyclinD1')
ax[0].set_title(f"CyclinD1 de-repression kinetics\n({'H3K27me3-memory: EZH2i de-represses over hours' if MEM else 'instant model'})", fontweight='bold'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.15)
ax[1].set_xlabel('time after vismo (h)'); ax[1].set_ylabel(mklab)
ax[1].set_title(f'{mklab} at CyclinD1 locus\n(EZH2i blocks deposition → mark decays)', fontweight='bold'); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.15)
plt.tight_layout(); plt.savefig('simulations/probe_derepression_kinetics.png', dpi=150, bbox_inches='tight'); plt.close()
print('Saved: probe_derepression_kinetics.png')
