"""The H3K27me3 tension: EZH2-driven methylation (integrates S-phase duration + mitogen) vs replicative
dilution (frequency = division rate) — and how it sets the CyclinD1 transcript.

At steady cycling the Ccnd1 mark sits where DEPOSITION balances REMOVAL:
  deposition  =  EZH2·(a_rw_prc2·Mk + k0_mk)·(1−Mk)·(eviction)       [EZH2 integrates S-phase & mitogen dose]
  removal     =  del_mk·Mk   +   ½·Mk once per cycle (dilution, frequency = 1/T_cc, set by mitogen)
CyclinD1 transcript = drive · R(Mk),  R = f0 + (1−f0)/(1+(Mk/K)^n).  More mark → lower Cd_mRNA.

Controlling parameters:
  deposition ↑ (→ Cd_mRNA ↓): a_rw_prc2 (read-write), kEZE2f (EZH2 production), kDeEZ ↓ (EZH2 stability =
                              how much EZH2 integrates S-phase duration)
  removal   ↑ (→ Cd_mRNA ↑): del_mk (turnover), replicative dilution (½ per cycle)
  mediators:                 SHH (mitogen → faster cycle = more dilution, but also more EZH2) and
                              kSyDna (fork speed → S-phase duration → EZH2 integration + cycle length)

Panels:
  A  Cd_mRNA over (mitogen SHH × S-phase duration) — the physiological tension surface.
  B  S-phase-duration sweep: EZH2 ↑ and Cd_mRNA ↓ as S lengthens ("EZH2 integrates S-phase duration").
  C  control-parameter sweeps (× default): which knobs push Cd_mRNA up (removal) vs down (deposition).
  D  deposition vs dilution: a_rw_prc2 sweep, Cd_mRNA with replicative dilution ON vs OFF (gap = dilution relief).

Run:  ./venv/bin/python simulations/sim_ezh2_dilution_tension.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_dilution_tension_cache.npz'
M_ON = build_model_v44(with_ezh2=True, with_hh=True)
M_OFF = M_ON.replace("Mk = 0.5*Mk", "Mk = 1.0*Mk")


def _rr(m):
    rr = te.loada(m); rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    return rr


RR, RR_OFF = _rr(M_ON), _rr(M_OFF)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEF = {p: RR[p] for p in ['a_rw_prc2', 'del_mk', 'kDeEZ', 'kEZE2f', 'kSyDna']}
SEL = ['time', 'Mk', 'Cd_mRNA', 'EZH2', 'MPF', 'Dna']


def run(rr, shh=0.6, params=None, T=13000, npts=39000):
    rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    rr['SHH'] = shh
    if params:
        for k, v in params.items():
            rr[k] = v
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            return rr.simulate(0, T, npts, selections=SEL)
        except Exception:
            rr.reset()
            for k, v in GNP.items():
                rr[k] = v
            rr['SHH'] = shh
            if params:
                for k, v in params.items():
                    rr[k] = v
            rr.integrator.setValue("absolute_tolerance", atol)
    return None


def measure(r, settle=8000):
    if r is None:
        return dict(cdm=np.nan, mk=np.nan, ez=np.nan, per=np.nan, sdur=np.nan)
    m = r['time'] >= settle; t = r['time'][m]; dt = t[1] - t[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    per = float(np.mean(np.diff(t[pk])) / 60) if len(pk) > 1 else np.nan
    dna = r['Dna'][m]; sdur = float(np.mean((dna > 0.02) & (dna < 0.98)) * per) if per == per else np.nan
    return dict(cdm=float(r['Cd_mRNA'][m].mean()), mk=float(r['Mk'][m].mean()),
                ez=float(r['EZH2'][m].mean()), per=per, sdur=sdur)


if FRESH or not os.path.exists(CACHE):
    print("computing dilution-tension sweeps ...")
    # ---- A: (SHH x kSyDna) -> Cd_mRNA ----
    SHH = np.array([0.30, 0.45, 0.65, 0.90, 1.30, 1.90])
    KSY = np.array([0.066, 0.044, 0.030, 0.020, 0.014, 0.010])   # high->low fork speed = short->long S
    ZA = np.full((len(KSY), len(SHH)), np.nan)
    for i, ks in enumerate(KSY):
        for j, s in enumerate(SHH):
            ZA[i, j] = measure(run(RR, shh=s, params={'kSyDna': ks}))['cdm']
        print(f"  [A] kSyDna={ks:.3f} row done")
    # S-phase duration (h) for each kSyDna at reference mitogen (for the axis)
    SDUR = np.array([measure(run(RR, shh=0.65, params={'kSyDna': ks}))['sdur'] for ks in KSY])

    # ---- B: S-phase duration sweep -> EZH2, Mk, Cd_mRNA ----
    KSY_B = np.array([0.066, 0.044, 0.030, 0.022, 0.016, 0.012, 0.010])
    bEZ, bMK, bCD, bSD = [], [], [], []
    for ks in KSY_B:
        r = measure(run(RR, shh=0.65, params={'kSyDna': ks}))
        bEZ.append(r['ez']); bMK.append(r['mk']); bCD.append(r['cdm']); bSD.append(r['sdur'])
        print(f"  [B] kSyDna={ks:.3f} S~{r['sdur']:.1f}h done")

    # ---- C: control-parameter sweeps (x default) -> Cd_mRNA ----
    FOLD = np.array([0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    PARAMS_C = ['a_rw_prc2', 'kEZE2f', 'del_mk', 'kDeEZ']
    ZC = {}
    for p in PARAMS_C:
        row = [measure(run(RR, shh=0.65, params={p: DEF[p] * f}))['cdm'] for f in FOLD]
        ZC[p] = np.array(row); print(f"  [C] {p} sweep done")

    # ---- D: a_rw_prc2 sweep, dilution ON vs OFF ----
    KW = DEF['a_rw_prc2'] * FOLD
    dON = np.array([measure(run(RR, shh=0.65, params={'a_rw_prc2': k}))['cdm'] for k in KW])
    dOFF = np.array([measure(run(RR_OFF, shh=0.65, params={'a_rw_prc2': k}))['cdm'] for k in KW])
    print("  [D] a_rw_prc2 x dilution done")

    np.savez(CACHE, SHH=SHH, KSY=KSY, ZA=ZA, SDUR=SDUR, KSY_B=KSY_B, bEZ=bEZ, bMK=bMK, bCD=bCD, bSD=bSD,
             FOLD=FOLD, **{f'C_{p}': ZC[p] for p in PARAMS_C}, KW=KW, dON=dON, dOFF=dOFF)
    print("  cached ->", CACHE)

z = np.load(CACHE)
PARAMS_C = ['a_rw_prc2', 'kEZE2f', 'del_mk', 'kDeEZ']
LBL = {'a_rw_prc2': 'k_w (read-write methylation) ↓Cd', 'kEZE2f': 'kEZE2f (EZH2 production) ↓Cd',
       'del_mk': 'del_mk (mark turnover) ↑Cd', 'kDeEZ': 'kDeEZ (EZH2 degradation) ↑Cd'}
COL = {'a_rw_prc2': '#c0392b', 'kEZE2f': '#e67e22', 'del_mk': '#2471a3', 'kDeEZ': '#27ae60'}

fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))

# ---- A ----
a = ax[0, 0]
pm = a.pcolormesh(z['SHH'], np.arange(len(z['KSY'])), z['ZA'], shading='nearest', cmap='viridis')
a.set_yticks(np.arange(len(z['KSY'])))
a.set_yticklabels([f'{s:.0f}' for s in z['SDUR']])
a.set_xlabel('mitogen (SHH)'); a.set_ylabel('S-phase duration (h)')
a.set_title('(A) CyclinD1 transcript over the tension mediators\n(mitogen × S-phase duration)', fontweight='bold', fontsize=11)
fig.colorbar(pm, ax=a, label='CyclinD1 transcript (Cd_mRNA)')
a.annotate('long S + low mitogen\n→ mark wins → Cd LOW', xy=(0.4, 4.3), color='white', fontsize=7.5, ha='left')
a.annotate('short S + high mitogen\n→ dilution wins → Cd HIGH', xy=(1.0, 0.5), color='white', fontsize=7.5, ha='left')

# ---- B ----
b = ax[0, 1]; b2 = b.twinx()
order = np.argsort(z['bSD'])
sd = np.array(z['bSD'])[order]
b.plot(sd, np.array(z['bEZ'])[order], '-o', color='#8e44ad', lw=2, label='EZH2 (deposition driver)')
b.plot(sd, np.array(z['bMK'])[order], '-s', color='#c0392b', lw=2, label='H3K27me3 mark Mk')
b2.plot(sd, np.array(z['bCD'])[order], '-^', color='#117a65', lw=2.4, label='CyclinD1 transcript')
b.set_xlabel('S-phase duration (h)'); b.set_ylabel('EZH2 / mark', color='#333')
b2.set_ylabel('CyclinD1 transcript', color='#117a65')
b.set_title('(B) EZH2 integrates S-phase duration\nlonger S → more EZH2 → more mark → less CyclinD1', fontweight='bold', fontsize=11)
h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
b.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right'); b.grid(alpha=0.15)

# ---- C ----
c = ax[1, 0]
for p in PARAMS_C:
    c.plot(z['FOLD'], z[f'C_{p}'], '-o', color=COL[p], lw=2, ms=4, label=LBL[p])
c.axvline(1.0, color='k', ls=':', alpha=0.4); c.text(1.02, c.get_ylim()[0], ' default', fontsize=7)
c.set_xscale('log'); c.set_xticks([0.3, 0.5, 1, 2, 3]); c.set_xticklabels(['0.3', '0.5', '1', '2', '3'])
c.set_xlabel('parameter × default'); c.set_ylabel('CyclinD1 transcript (Cd_mRNA)')
c.set_title('(C) Control parameters of the tension\ndeposition knobs ↓ CyclinD1; removal knobs ↑ CyclinD1', fontweight='bold', fontsize=11)
c.legend(fontsize=7.5, loc='best'); c.grid(alpha=0.15)

# ---- D ----
d = ax[1, 1]
d.plot(z['KW'], z['dOFF'], '-s', color='#7f8c8d', lw=2, label='dilution OFF')
d.plot(z['KW'], z['dON'], '-o', color='#c0392b', lw=2, label='dilution ON (default)')
d.fill_between(z['KW'], z['dON'], z['dOFF'], color='#f39c12', alpha=0.18, label='dilution relief')
d.axvline(DEF['a_rw_prc2'], color='k', ls=':', alpha=0.4); d.text(DEF['a_rw_prc2'], d.get_ylim()[1], ' default', fontsize=7, va='top')
d.set_xlabel('a_rw_prc2 (read-write methylation strength)'); d.set_ylabel('CyclinD1 transcript (Cd_mRNA)')
d.set_title('(D) Deposition vs dilution\nstronger methylation ↓Cd, but dilution keeps clawing it back', fontweight='bold', fontsize=11)
d.legend(fontsize=8, loc='best'); d.grid(alpha=0.15)

fig.suptitle('H3K27me3 tension at Ccnd1: EZH2-driven methylation (integrates S-phase + mitogen) vs replicative dilution → CyclinD1 transcript',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_dilution_tension.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_dilution_tension.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_dilution_tension.png / .pdf")
