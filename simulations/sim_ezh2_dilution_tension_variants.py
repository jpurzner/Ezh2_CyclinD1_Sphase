"""Variants of the H3K27me3 deposition-vs-dilution tension surface (companion to
sim_ezh2_dilution_tension.py): (1) the HU-axis version — S-phase elongated by the EXPERIMENTAL lever
hydroxyurea (replication stress) instead of the internal fork-speed parameter; and (2) the MB/vismo
variant — the same tension in the tumour, where the mitogen axis is vismodegib titrating tonic Hh DOWN.

HU slows replication forks (vfork(HU)) -> S-phase elongates -> EZH2 (long-lived, S/G2-gated) integrates
the longer S -> more H3K27me3 -> lower CyclinD1 transcript. Period stays ~constant under HU, so HU
isolates the S-phase-integration arm of the tension from the dilution-frequency arm.

Panels:
  A  HU-axis GNP surface: CyclinD1 transcript over (mitogen SHH × HU dose).
  B  MB/vismo surface: CyclinD1 transcript over (vismodegib HHi × HU dose), MB context.
  C  HU decomposition (GNP): HU -> S-phase duration, EZH2, mark, CyclinD1 (the mechanism).
  D  vismo decomposition (MB): HHi -> EZH2, mark, CyclinD1 (the mitogen arm in the tumour).

Run:  ./venv/bin/python simulations/sim_ezh2_dilution_tension_variants.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_dilution_tension_variants_cache.npz'
rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)

GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002)
MB = dict(MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
SEL = ['time', 'Mk', 'Cd_mRNA', 'EZH2', 'Gli1', 'Gli_act', 'MPF', 'Dna']


def run(ctx, shh=0.6, hhi=0.0, hu=0.0, T=14000, npts=42000):
    base = dict(ctx);
    def _set():
        rr.reset()
        for k, v in base.items():
            rr[k] = v
        rr['SHH'] = shh; rr['HHi'] = hhi; rr['HU'] = hu; rr['EZH2i'] = 0
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        _set(); rr.integrator.setValue("absolute_tolerance", atol)
        try:
            return rr.simulate(0, T, npts, selections=SEL)
        except Exception:
            continue
    return None


def measure(r, settle=8500):
    if r is None:
        return dict(cdm=np.nan, mk=np.nan, ez=np.nan, gli1=np.nan, sdur=np.nan)
    m = r['time'] >= settle; t = r['time'][m]; dt = t[1] - t[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    per = float(np.mean(np.diff(t[pk])) / 60) if len(pk) > 1 else np.nan
    dna = r['Dna'][m]; sdur = float(np.mean((dna > 0.02) & (dna < 0.98)) * per) if per == per else np.nan
    return dict(cdm=float(r['Cd_mRNA'][m].mean()), mk=float(r['Mk'][m].mean()), ez=float(r['EZH2'][m].mean()),
                gli1=float((r['Gli1'][m] + r['Gli_act'][m]).mean()), sdur=sdur)


if FRESH or not os.path.exists(CACHE):
    print("computing tension variants ...")
    SHH = np.array([0.30, 0.45, 0.65, 0.90, 1.30, 1.90])
    HU = np.array([0.0, 0.15, 0.3, 0.45, 0.6])
    HHI = np.array([0.0, 0.3, 0.5, 0.7, 0.85, 0.92])
    # A: GNP (SHH x HU)
    ZA = np.full((len(HU), len(SHH)), np.nan)
    for i, hu in enumerate(HU):
        for j, s in enumerate(SHH):
            ZA[i, j] = measure(run(GNP, shh=s, hu=hu))['cdm']
        print(f"  [A] HU={hu:.2f} row done")
    # B: MB (HHi x HU)
    ZB = np.full((len(HU), len(HHI)), np.nan)
    for i, hu in enumerate(HU):
        for j, h in enumerate(HHI):
            ZB[i, j] = measure(run(MB, shh=0.5, hhi=h, hu=hu))['cdm']
        print(f"  [B] HU={hu:.2f} row done")
    # C: GNP HU decomposition
    HU_C = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65])
    cSD, cEZ, cMK, cCD = [], [], [], []
    for hu in HU_C:
        r = measure(run(GNP, shh=0.6, hu=hu)); cSD.append(r['sdur']); cEZ.append(r['ez']); cMK.append(r['mk']); cCD.append(r['cdm'])
        print(f"  [C] HU={hu:.2f} S~{r['sdur']:.1f}h done")
    # D: MB vismo decomposition
    HHI_D = np.array([0.0, 0.2, 0.4, 0.6, 0.75, 0.85, 0.92])
    dEZ, dMK, dCD, dG1 = [], [], [], []
    for h in HHI_D:
        r = measure(run(MB, shh=0.5, hhi=h)); dEZ.append(r['ez']); dMK.append(r['mk']); dCD.append(r['cdm']); dG1.append(r['gli1'])
        print(f"  [D] HHi={h:.2f} done")
    np.savez(CACHE, SHH=SHH, HU=HU, HHI=HHI, ZA=ZA, ZB=ZB,
             HU_C=HU_C, cSD=cSD, cEZ=cEZ, cMK=cMK, cCD=cCD,
             HHI_D=HHI_D, dEZ=dEZ, dMK=dMK, dCD=dCD, dG1=dG1)
    print("  cached ->", CACHE)

z = np.load(CACHE)
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))

# ---- A: GNP SHH x HU ----
a = ax[0, 0]
pm = a.pcolormesh(z['SHH'], z['HU'], z['ZA'], shading='nearest', cmap='viridis')
a.set_xlabel('mitogen (SHH)'); a.set_ylabel('HU dose  (replication stress → longer S)')
a.set_title('(A) HU-axis tension surface (GNP)\nCyclinD1 transcript over mitogen × HU', fontweight='bold', fontsize=11)
fig.colorbar(pm, ax=a, label='CyclinD1 transcript (Cd_mRNA)')
a.annotate('↑HU → longer S → more EZH2\n→ more mark → less CyclinD1', xy=(0.4, 0.5), color='white', fontsize=7.5)

# ---- B: MB HHi x HU ----
b = ax[0, 1]
pm2 = b.pcolormesh(z['HHI'], z['HU'], z['ZB'], shading='nearest', cmap='viridis')
b.set_xlabel('vismodegib (HHi) → less mitogen'); b.set_ylabel('HU dose  (replication stress → longer S)')
b.set_title('(B) MB/vismo tension surface\nCyclinD1 transcript over vismo × HU', fontweight='bold', fontsize=11)
fig.colorbar(pm2, ax=b, label='CyclinD1 transcript (Cd_mRNA)')
b.annotate('vismo + HU both raise the mark\n→ deepest CyclinD1 repression', xy=(0.35, 0.5), color='white', fontsize=7.5)

# ---- C: GNP HU decomposition ----
c = ax[1, 0]; c2 = c.twinx()
sd = np.array(z['cSD'])
c.plot(sd, z['cEZ'], '-o', color='#8e44ad', lw=2, label='EZH2 (deposition driver)')
c.plot(sd, z['cMK'], '-s', color='#c0392b', lw=2, label='H3K27me3 mark Mk')
c2.plot(sd, z['cCD'], '-^', color='#117a65', lw=2.4, label='CyclinD1 transcript')
c.set_xlabel('S-phase duration (h)  [set by HU]'); c.set_ylabel('EZH2 / mark')
c2.set_ylabel('CyclinD1 transcript', color='#117a65')
c.set_title('(C) HU elongates S → EZH2 integrates it\n(GNP; the experimental replication-stress lever)', fontweight='bold', fontsize=11)
h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels()
c.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right'); c.grid(alpha=0.15)

# ---- D: MB vismo decomposition ----
d = ax[1, 1]; d2 = d.twinx()
hh = z['HHI_D']
d.plot(hh, z['dEZ'], '-o', color='#8e44ad', lw=2, label='EZH2')
d.plot(hh, z['dMK'], '-s', color='#c0392b', lw=2, label='H3K27me3 mark Mk')
d2.plot(hh, z['dCD'], '-^', color='#117a65', lw=2.4, label='CyclinD1 transcript')
d.set_xlabel('vismodegib (HHi)'); d.set_ylabel('EZH2 / mark')
d2.set_ylabel('CyclinD1 transcript', color='#117a65')
d.set_title('(D) Vismo (less mitogen) in MB\nGli-drive ↓ → CyclinD1 ↓ toward the MYCN floor', fontweight='bold', fontsize=11)
h1, l1 = d.get_legend_handles_labels(); h2, l2 = d2.get_legend_handles_labels()
d.legend(h1 + h2, l1 + l2, fontsize=8, loc='center right'); d.grid(alpha=0.15)

fig.suptitle('H3K27me3 tension — HU-axis (replication stress → S-phase) and MB/vismo variants → CyclinD1 transcript',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_dilution_tension_variants.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_dilution_tension_variants.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_dilution_tension_variants.png / .pdf")
