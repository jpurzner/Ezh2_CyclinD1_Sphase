"""Replicative dilution of H3K27me3 limits EZH2's repressive reach — mitogen sweep, dilution ON vs OFF.

The AUM module halves the Ccnd1 H3K27me3 mark (Mk) once per cycle at early S (Dna>0.05). Because the
cell-cycle period is set by the mitogen-driven CyclinD1, faster-cycling (higher-mitogen) cells dilute
the mark MORE often -> the mark cannot fully restore between divisions -> it ratchets below its
no-dilution steady state -> EZH2 represses Ccnd1 LESS. This is the endogenous-clock replicative-dilution
phenotype (Jadhav 2020): repression is proliferation-rate-dependent.

We compare the model WITH dilution (default) vs WITHOUT (the S-phase halving disabled: Mk=0.5*Mk ->
Mk=1.0*Mk) across a mitogen (SHH) sweep in the GNP context.

Panels:
  A  Mk(t) at a representative mitogen: ON = sawtooth (halved each division, re-methylates); OFF = smooth,
     higher. Division times marked -> visualizes the dilution events themselves.
  B  steady-state mark Mk vs mitogen, ON vs OFF -> dilution holds Mk below the no-dilution ceiling, and
     the gap WIDENS with mitogen (faster cycling = more-frequent dilution).
  C  steady-state CyclinD1 vs mitogen, ON vs OFF -> the mark relief translates into higher CyclinD1
     (EZH2's repressive reach is limited by proliferation).
  D  dilution frequency (divisions / 168 h) vs mitogen -> the driver: replication rate rises with mitogen,
     so the dilution the mark experiences rises with it.

Run:  ./venv/bin/python simulations/fig_v44_replicative_dilution.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

# dilution ON = default AUM; dilution OFF = disable the S-phase halving event
M_ON = build_model_v44(with_ezh2=True, with_hh=True)
M_OFF = M_ON.replace("Mk = 0.5*Mk", "Mk = 1.0*Mk")
assert M_OFF != M_ON, "failed to disable the dilution event"


def _rr(model):
    rr = te.loada(model)
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    return rr


RR_ON, RR_OFF = _rr(M_ON), _rr(M_OFF)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SEL = ['time', 'Mk', 'Cd', 'MPF', 'Dna', 'EZH2']


def run(rr, shh, T=13000, npts=26000):
    rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    rr['SHH'] = shh
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            return rr.simulate(0, T, npts, selections=SEL)
        except Exception:
            rr.reset()
            for k, v in GNP.items():
                rr[k] = v
            rr['SHH'] = shh
            rr.integrator.setValue("absolute_tolerance", atol)
    return None


def summarize(r, settle=7000):
    """steady mean Mk, mean Cd, and divisions/168h from a settled window."""
    if r is None:
        return np.nan, np.nan, 0.0
    m = r['time'] >= settle
    t = r['time'][m]; dt = t[1] - t[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    hours = (t[-1] - t[0]) / 60.0
    div_168 = len(pk) / hours * 168.0
    return float(r['Mk'][m].mean()), float(r['Cd'][m].mean()), div_168


SHH = np.array([0.2, 0.3, 0.45, 0.65, 0.9, 1.3, 1.8, 2.5])
mk_on, cd_on, dv_on = [], [], []
mk_off, cd_off, dv_off = [], [], []
for s in SHH:
    a = summarize(run(RR_ON, s)); mk_on.append(a[0]); cd_on.append(a[1]); dv_on.append(a[2])
    b = summarize(run(RR_OFF, s)); mk_off.append(b[0]); cd_off.append(b[1]); dv_off.append(b[2])
    print(f"SHH={s:.2f}  ON: Mk {a[0]:.3f} Cd {a[1]:5.2f} div {a[2]:.1f} | OFF: Mk {b[0]:.3f} Cd {b[1]:5.2f} div {b[2]:.1f}")
mk_on, cd_on, dv_on = map(np.array, (mk_on, cd_on, dv_on))
mk_off, cd_off, dv_off = map(np.array, (mk_off, cd_off, dv_off))

# time courses at a representative cycling mitogen
S_REP = 0.8
tc_on = run(RR_ON, S_REP); tc_off = run(RR_OFF, S_REP)

C_ON, C_OFF = '#2471a3', '#c0392b'
fig, ax = plt.subplots(2, 2, figsize=(13.5, 10))

# ---- A: Mk(t) sawtooth ----
a = ax[0, 0]
for r, col, lab in [(tc_off, C_OFF, 'dilution OFF'), (tc_on, C_ON, 'dilution ON')]:
    th = r['time'] / 60.0
    w = th >= (th[-1] - 90)     # last ~90 h window
    a.plot(th[w], r['Mk'][w], color=col, lw=1.8, label=lab)
# mark division times (MPF peaks) in the ON trace
th = tc_on['time'] / 60.0; w = th >= (th[-1] - 90)
pk, _ = find_peaks(tc_on['MPF'][w], prominence=0.15, distance=int(200 / (tc_on['time'][1] - tc_on['time'][0])))
for p in np.array(th[w])[pk]:
    a.axvline(p, color='#95a5a6', ls=':', lw=1.0, alpha=0.7)
a.plot([], [], color='#95a5a6', ls=':', label='division (S-phase halving)')
a.set_xlabel('time (h)'); a.set_ylabel('H3K27me3 mark  Mk')
a.set_title(f'(A) The dilution mechanism (GNP, SHH={S_REP})\nMk is halved each division then re-methylates (sawtooth)', fontweight='bold', fontsize=11)
a.legend(fontsize=8, loc='best'); a.grid(alpha=0.15)

# ---- B: steady Mk vs mitogen ----
b = ax[0, 1]
b.plot(SHH, mk_off, '-s', color=C_OFF, lw=2, label='dilution OFF (no-dilution ceiling)')
b.plot(SHH, mk_on, '-o', color=C_ON, lw=2, label='dilution ON')
b.fill_between(SHH, mk_on, mk_off, color='#f39c12', alpha=0.18, label='dilution relief')
b.set_xlabel('mitogen (SHH)'); b.set_ylabel('steady-state mark  Mk')
b.set_title('(B) Dilution holds the mark below its ceiling\n— the gap widens as mitogen speeds the cycle', fontweight='bold', fontsize=11)
b.legend(fontsize=8, loc='best'); b.grid(alpha=0.15)

# ---- C: steady CyclinD1 vs mitogen ----
c = ax[1, 0]
c.plot(SHH, cd_off, '-s', color=C_OFF, lw=2, label='dilution OFF')
c.plot(SHH, cd_on, '-o', color=C_ON, lw=2, label='dilution ON')
c.fill_between(SHH, cd_off, cd_on, color='#27ae60', alpha=0.15, label='extra CyclinD1 from dilution')
c.set_xlabel('mitogen (SHH)'); c.set_ylabel('steady-state CyclinD1 (Cd)')
c.set_title('(C) Relieved repression → higher CyclinD1\n(EZH2 reach limited by proliferation rate)', fontweight='bold', fontsize=11)
c.legend(fontsize=8, loc='best'); c.grid(alpha=0.15)

# ---- D: dilution frequency (divisions/168h) vs mitogen ----
d = ax[1, 1]
d.plot(SHH, dv_on, '-o', color=C_ON, lw=2, label='dilution ON')
d.plot(SHH, dv_off, '-s', color=C_OFF, lw=2, alpha=0.7, label='dilution OFF')
d.set_xlabel('mitogen (SHH)'); d.set_ylabel('divisions / 168 h  (= dilution frequency)')
d.set_title('(D) The driver: replication rate rises with mitogen\n→ more-frequent mark dilution at high mitogen', fontweight='bold', fontsize=11)
d.legend(fontsize=8, loc='best'); d.grid(alpha=0.15)

fig.suptitle('Replicative dilution of H3K27me3 limits EZH2 repression of CyclinD1 — proliferation-rate-dependent (AUM model)',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_replicative_dilution.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_replicative_dilution.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_replicative_dilution.png / .pdf")
