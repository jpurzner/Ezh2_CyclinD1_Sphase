"""GNP mitogen ramp-ON (entry) -> sustain -> ramp-OFF (exit): EZH2, CyclinD1, and cell division,
WITH vs WITHOUT the H3K27me3-mediated repression of CyclinD1 (f0_prc2=1.0 removes it). Current 25/27
parameters. Single deterministic GNP cell; SHH staircased up from arrest, held, then down to arrest.

Shows: (i) entry -- the mark raises the Hh needed to start dividing; (ii) exit -- how EZH2/CyclinD1 fall
and when division collapses; the WITHOUT-mark cell coasts to a different (lower) exit mitogen.

Run:  ./venv/bin/python simulations/fig_v44_mitogen_ramp_onoff.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from src.build_model_v44_heldt import build_model_v44

SHH_LO, SHH_HI = 0.05, 0.60
SETTLE = 4000                       # settle arrested at SHH_LO
NUP, NDN, SUSTAIN = 30, 40, 6000    # ramp-up segs, ramp-down segs, sustain (min)
SEG = 300                           # min per ramp step
SEL = ["time", "SHH", "Cd_mRNA", "EZH2", "Dna", "P21"]


def ramp(model):
    rr = te.loada(model)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    try:
        rr.integrator.setValue("maximum_num_steps", 1000000); rr.integrator.setValue("maximum_time_step", 20.0)
    except Exception:
        pass
    rr['SHH'] = SHH_LO
    def seg(t0, t1, npts):
        for atol, ms in [(1e-9, 20.0), (1e-8, 20.0), (1e-7, 8.0), (1e-6, 3.0), (1e-5, 1.5)]:
            try:
                rr.integrator.setValue("absolute_tolerance", atol); rr.integrator.setValue("maximum_time_step", ms)
                return rr.simulate(t0, t1, npts, selections=SEL)
            except Exception:
                continue
        return None
    chunks = []; c = seg(0, SETTLE, 4000)
    if c is not None: chunks.append(c)
    t = SETTLE
    sched = ([('up', s) for s in np.linspace(SHH_LO, SHH_HI, NUP)]
             + [('hold', SHH_HI)] * int(SUSTAIN / SEG)
             + [('dn', s) for s in np.linspace(SHH_HI, SHH_LO, NDN)])
    for _, s in sched:
        rr['SHH'] = float(s); c = seg(t, t + SEG, 400)
        if c is None:
            break
        chunks.append(c); t += SEG
    return {k: np.concatenate([c[k] for c in chunks]) for k in SEL}


def divisions(d):
    Dna = d['Dna']; idx = np.where((Dna[:-1] > 0.9) & (Dna[1:] < 0.1))[0]
    return d['time'][idx] / 60.0


A = ramp(build_model_v44())                              # WITH mark (real model, 25/27)
B = ramp(build_model_v44(params={"f0_prc2": 1.0}))         # WITHOUT mark repression

fig, ax = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)
for row, (name, d, col_note) in enumerate([
        ('WITH H3K27me3 repression (real model, 25/27)', A, 'the mark raises the entry Hh and buffers CyclinD1'),
        ('WITHOUT H3K27me3 repression (f0_prc2=1.0)', B, 'no mark: CyclinD1 tracks mitogen directly')]):
    a = ax[row]
    m = d['time'] >= SETTLE                                   # drop the pre-ramp settle transient
    th = (d['time'][m] - SETTLE) / 60.0
    dv_all = divisions(d)
    if len(dv_all):
        ent, ex = dv_all.min(), dv_all.max()
        shh_ent = np.interp(ent * 60, d['time'], d['SHH']); shh_ex = np.interp(ex * 60, d['time'], d['SHH'])
    dv = (dv_all - SETTLE / 60.0); dv = dv[dv >= 0]           # shift division times to ramp-onset origin
    ent_s = ent - SETTLE / 60.0 if len(dv_all) else None; ex_s = ex - SETTLE / 60.0 if len(dv_all) else None
    a.plot(th, d['Cd_mRNA'][m], color='#117a65', lw=1.5, label='CyclinD1 (Cd_mRNA)')
    a.plot(th, d['EZH2'][m], color='#8e44ad', lw=1.7, label='EZH2 protein')
    a.plot(th, d['P21'][m], color='#b9770e', lw=1.1, alpha=0.7, label='p27')
    # division event ticks along the top
    ymax = max(np.nanmax(d['Cd_mRNA'][m]), 3) * 1.15
    a.set_ylim(0, ymax)
    for x in dv:
        a.plot([x, x], [ymax * 0.93, ymax], color='#e67e22', lw=1.1, alpha=0.85)
    a.plot([], [], color='#e67e22', lw=1.1, label=f'divisions (n={len(dv)})')
    # SHH on twin axis
    a2 = a.twinx(); a2.plot(th, d['SHH'][m], color='#888', lw=1.4, ls=':', label='SHH (mitogen)'); a2.set_ylim(0, 0.72)
    a2.set_ylabel('SHH (mitogen)', color='#666')
    if len(dv_all):
        a.axvline(ent_s, color='#2e7d32', ls='--', alpha=0.5); a.axvline(ex_s, color='#c0392b', ls='--', alpha=0.5)
        a.text(ent_s, ymax * 0.60, f' entry\n SHH~{shh_ent:.2f}', fontsize=8, color='#2e7d32')
        a.text(ex_s, ymax * 0.60, f' exit\n SHH~{shh_ex:.2f}', fontsize=8, color='#c0392b', ha='left')
    a.set_ylabel('level (a.u.)')
    a.set_title(f'{name}  —  {col_note}', fontsize=10.5, fontweight='bold')
    a.legend(fontsize=8, ncol=3, loc='upper left'); a.grid(alpha=0.18)
ax[1].set_xlabel('time (h)')
fig.suptitle('GNP mitogen ramp-ON → sustain → ramp-OFF: EZH2, CyclinD1, division — WITH vs WITHOUT H3K27me3 repression (25/27 params)',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_mitogen_ramp_onoff.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_mitogen_ramp_onoff.pdf', bbox_inches='tight')
plt.close()
for name, d in (('WITH mark', A), ('WITHOUT mark', B)):
    dv = divisions(d)
    if len(dv):
        se = np.interp(dv.min() * 60, d['time'], d['SHH']); sx = np.interp(dv.max() * 60, d['time'], d['SHH'])
        print(f'{name:13s}: {len(dv)} divisions | entry SHH~{se:.3f} | exit SHH~{sx:.3f} | '
              f'EZH2 cyc {d["EZH2"][(d["time"]>SETTLE)&(d["time"]<SETTLE+SUSTAIN)].mean():.2f}')
    else:
        print(f'{name:13s}: 0 divisions')
print('Saved fig_v44_mitogen_ramp_onoff.png')
