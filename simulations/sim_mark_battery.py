"""The H3K27me3 "BATTERY" (JP mechanism): with a fast EZH2 (t1/2 10h, the 25/27 value) but a
self-reinforcing mark (moderate read-write k_w_mk), the H3K27me3 at the CyclinD1 promoter CHARGES UP
gradually over successive divisions during proliferation -> CyclinD1 gets progressively harder to push
(more repressed) -> the charged mark PERSISTS through G0 -> on withdrawal it slowly DISCHARGES, and the
discharge time (set by the mark turnover del_mk = the DURABILITY axis) is the withdrawal memory.

Single deterministic GNP cell: settle arrested -> cycle at SHH=0.55 (~300h, CHARGE) -> withdraw to 0.05
-> hold (~250h, DISCHARGE). Swept over del_mk (mark durability). kDeEZ fixed at the fast 10h value.

Run:  ./venv/bin/python simulations/sim_mark_battery.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
KDEEZ, KWMK = 0.0012, 0.0028           # fast EZH2 (t1/2 10h, 25/27) + weak-moderate read-write (gradual charge, sub-bistable)
DELMK = [0.0005, 0.0015, 0.004]        # mark turnover = DURABILITY axis (durable -> leaky)
COL = ['#08519c', '#3182bd', '#9ecae1']
SETTLE, CYCLE, HOLD, CHUNK, SHH_CYC = 3000.0, 15000.0, 18000.0, 30.0, 0.85   # high mitogen -> SUSTAIN cycling despite the charging mark
_R = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_R.integrator.setValue('relative_tolerance', 1e-6)
_R.integrator.setValue('maximum_num_steps', 1000000)
try: _R.integrator.setValue('maximum_time_step', 20.0)
except Exception: pass


def shh_of(t):
    if t < CYCLE: return SHH_CYC                     # proliferate (charge)
    return 0.15                                      # partial withdraw (coast/arrest -> discharge; avoids deep-arrest stiffness)


def run(delmk):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    _R['kDeEZ'] = KDEEZ; _R['k_w_mk'] = KWMK; _R['del_mk'] = delmk; _R['f0_mk'] = 0.145; _R['SHH'] = 0.05
    tm = 0.0
    while tm < SETTLE:
        try: _R.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: break
        tm += CHUNK
    T, SH, MK, CD, DN, EZ = [], [], [], [], [], []
    t = 0.0
    while t < CYCLE + HOLD:
        _R['SHH'] = float(shh_of(t)); r = None
        for atol, ms in [(1e-9, 20.0), (1e-8, 20.0), (1e-7, 10.0), (1e-6, 4.0), (1e-5, 1.5)]:
            try:
                _R.integrator.setValue('absolute_tolerance', atol); _R.integrator.setValue('maximum_time_step', ms)
                r = _R.simulate(tm, tm + CHUNK, 4, selections=['time', 'SHH', 'Mk', 'Cd_mRNA', 'Dna', 'EZH2']); break
            except Exception: r = None
        if r is None: break
        for A, s in [(T, 'time'), (SH, 'SHH'), (MK, 'Mk'), (CD, 'Cd_mRNA'), (DN, 'Dna'), (EZ, 'EZH2')]:
            A.append(r[s][1:])
        t += CHUNK; tm += CHUNK
    T = (np.concatenate(T) - SETTLE) / 60.0
    return dict(t=T, shh=np.concatenate(SH), mk=np.concatenate(MK), cd=np.concatenate(CD), dna=np.concatenate(DN), ez=np.concatenate(EZ))


runs = {dm: run(dm) for dm in DELMK}
ref = runs[DELMK[1]]                                 # middle del_mk for the Cd / per-division panels
div = np.where((ref['dna'][:-1] > 0.9) & (ref['dna'][1:] < 0.1))[0]
divt = ref['t'][div]
# per-division CyclinD1 peak (the "harder to push")
peaks = []
edges = np.concatenate([[0], div, [len(ref['t']) - 1]])
for a, b in zip(edges[:-1], edges[1:]):
    seg = ref['cd'][a:b];  peaks.append(np.max(seg) if len(seg) else np.nan)
peaks = np.array(peaks[1:])                           # per completed cycle
cyc_divt = divt[divt < CYCLE / 60.0]

fig, ax = plt.subplots(2, 2, figsize=(15, 9))
CYh = CYCLE / 60.0
# (A) SHH + mark charging/discharging
a = ax[0, 0]
a.axvspan(0, CYh, color='#2e8b57', alpha=0.06); a.axvspan(CYh, ref['t'][-1], color='#c0392b', alpha=0.06)
for dm, c in zip(DELMK, COL):
    a.plot(runs[dm]['t'], runs[dm]['mk'], color=c, lw=1.8, label=f'del_mk={dm:g} (t½ mark {np.log(2)/dm/60:.0f}h)')
a.set_ylabel('H3K27me3 mark (Mk) @ CyclinD1'); a.set_xlabel('time (h)')
a.set_title('(A) The BATTERY: mark CHARGES over divisions (green) then persists/DISCHARGES in G0 (red)', fontweight='bold', fontsize=10.5)
a.legend(fontsize=8, loc='center right'); a.grid(alpha=0.15); a.set_ylim(0, 1)
# (B) CyclinD1 over time (middle del_mk) with divisions
a = ax[0, 1]
a.axvspan(0, CYh, color='#2e8b57', alpha=0.06); a.axvspan(CYh, ref['t'][-1], color='#c0392b', alpha=0.06)
a.plot(ref['t'], ref['cd'], color='#8b1a1a', lw=1.1)
for x in cyc_divt: a.axvline(x, color='k', lw=0.4, alpha=0.18)
a.set_ylabel('CyclinD1 (Cd_mRNA)'); a.set_xlabel('time (h)')
a.set_title('(B) CyclinD1 sawtooth — envelope DECLINES as the mark charges (harder to push)', fontweight='bold', fontsize=10.5); a.grid(alpha=0.15)
# (C) per-division CyclinD1 peak (harder to push) + mark
a = ax[1, 0]
dn = np.arange(1, len(peaks) + 1)
ncyc = int(np.sum(cyc_divt < CYh)) or len(peaks)
a.plot(dn[:ncyc], peaks[:ncyc], '-o', color='#8b1a1a', ms=4, label='CyclinD1 peak / division')
a.set_xlabel('division number (during proliferation)'); a.set_ylabel('CyclinD1 peak (a.u.)', color='#8b1a1a')
ab = a.twinx()
mk_at_div = [ref['mk'][d] for d in div[:ncyc]]
ab.plot(dn[:ncyc], mk_at_div, '-s', color='#2e8b57', ms=4, label='mark at division')
ab.set_ylabel('mark Mk at division', color='#2e8b57'); ab.set_ylim(0, 1)
a.set_title('(C) "Harder to push": CyclinD1 peak falls as H3K27me3 charges over divisions', fontweight='bold', fontsize=10.5); a.grid(alpha=0.15)
# (D) withdrawal discharge (battery memory timescale) vs del_mk
a = ax[1, 1]
for dm, c in zip(DELMK, COL):
    d = runs[dm]; m = d['t'] >= CYh
    a.plot(d['t'][m] - CYh, d['mk'][m], color=c, lw=2.0, label=f'del_mk={dm:g}')
a.axhline(0.2, color='k', ls=':', lw=1, alpha=0.5); a.text(5, 0.21, 'repression threshold ~0.2', fontsize=8)
a.set_xlabel('time since withdrawal (h)'); a.set_ylabel('mark Mk (G0)')
a.set_title('(D) BATTERY DISCHARGE in G0 — durability (del_mk) sets the memory duration', fontweight='bold', fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.15); a.set_ylim(0, 1)
fig.suptitle('H3K27me3 "battery": fast EZH2 (t½10h, 25/27) + self-reinforcing mark — charge over divisions, persist through G0, discharge on withdrawal (GNP single cell)',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_mark_battery.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_mark_battery.pdf', bbox_inches='tight')
print('Saved sim_mark_battery.png')
print(f'divisions during proliferation: {ncyc}')
print('mark at successive divisions:', [f'{m:.2f}' for m in mk_at_div])
print('CyclinD1 peak per division:', [f'{p:.2f}' for p in peaks[:ncyc]])
for dm in DELMK:
    d = runs[dm]; m = d['t'] >= CYh
    mk_wd = d['mk'][m]; tt = d['t'][m] - CYh
    if len(mk_wd) == 0:
        print(f'del_mk={dm:g}: (no withdrawal window)'); continue
    below = np.where(mk_wd < 0.2)[0]
    print(f'del_mk={dm:g}: Mk at withdrawal {mk_wd[0]:.2f} -> discharge below 0.2 at {tt[below[0]]:.0f}h' if len(below) else f'del_mk={dm:g}: Mk {mk_wd[0]:.2f} -> still >0.2 at end ({tt[-1]:.0f}h)')
