"""CyclinD1 -> transient-G0 mechanism (CORRECTED, p27 marker).

(a) Transient G0 (p27-high, pre-S) duration vs mean CyclinD1 protein -- SHH / k_Cd_translation /
    EZH2-repression knobs collapse onto one curve; below a threshold the cell is permanently G0.
(b) Phase partition (G0/G1/S/G2) vs CyclinD1 -- G0 shrinks into G1/S, period ~constant.
(c) One-cycle anatomy with the transient G0 (p27-high) shaded: phospho-Rb saturates early, but
    E2f / CyclinE+A / p27 show the real ~10 h growth-timed G0 before commitment.
(d) EZH2 -| CyclinD1 -| G0: EZH2 inhibitor de-represses CyclinD1 and collapses the transient G0.

Run:  ./venv/bin/python simulations/fig_v44_g0_mechanism.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "Cd", "P21", "E2f", "Ce", "Ca", "aRc", "Dna", "MPF", "pRb", "mass"]
T_END, N_PTS, SETTLE, P27 = 16000, 64000, 7000, 0.1


def sim(params=None, shh=0.5, mycn=1.0, ptch1=1.0, ezh2i=0.0):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=(params or None)))
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = shh; rr['MYCN_amplification'] = mycn; rr['Ptch1_copy_number'] = ptch1; rr['EZH2i'] = ezh2i
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def metrics(res):
    t = res['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1]-tt[0]
    pk, _ = find_peaks(res['MPF'][m], prominence=0.15, distance=int(200/dt))
    cd = float(np.mean(res['Cd'][m]))
    P21 = res['P21'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = Dna >= 0.98; preS = ~in_S & ~in_G2
    G0 = preS & (P21 > P27); G1 = preS & (P21 <= P27)
    if len(pk) < 2:
        return dict(cd=cd, arrest=True, period=np.nan, G0=np.nan, G1=np.nan, S=np.nan, G2=np.nan)
    period = float(np.mean(np.diff(tt[pk])))/60.0
    fr = {k: v.mean() for k, v in dict(G0=G0, G1=G1, S=in_S, G2=in_G2).items()}
    s = sum(fr.values()) or 1.0
    return dict(cd=cd, arrest=False, period=period, **{k: period*v/s for k, v in fr.items()})


# sweeps
shh_s = [(metrics(sim(shh=s))) for s in [0.35, 0.45, 0.6, 0.8, 1.0]]
tl_s = [(metrics(sim({"k_Cd_translation": k}))) for k in [0.35, 0.45, 0.55, 0.65, 0.75, 0.95]]
ke_s = [(metrics(sim({"f0_prc2": k}))) for k in [0.05, 0.2, 0.4, 0.7]]   # AUM: f0_prc2 leaky floor (lower = stronger repression)
mb = metrics(sim(shh=0.5, ptch1=0.1, mycn=2.8))
mb_e = metrics(sim(shh=0.5, ptch1=0.1, mycn=2.8, ezh2i=1.0))
gnp = metrics(sim(shh=0.5)); gnp_e = metrics(sim(shh=0.5, ezh2i=1.0))

fig, ax = plt.subplots(2, 2, figsize=(13.5, 10.5))

# (a) G0 vs CyclinD1 protein
for s, lab, c, mk in [(shh_s, 'SHH dose', '#1b9e77', 'o'), (tl_s, 'translation', '#d95f02', 's'),
                      (ke_s, 'EZH2 repression', '#7570b3', '^')]:
    pts = [(a['cd'], a['G0']) for a in s if not a['arrest']]
    if pts:
        x, y = np.array(pts).T; ax[0, 0].scatter(x, y, s=70, color=c, marker=mk, edgecolor='k', label=lab, zorder=3)
arr_cd = max([a['cd'] for a in tl_s + ke_s if a['arrest']] + [0])
ax[0, 0].axvspan(0, arr_cd, color='#fdecea', alpha=0.7, zorder=0)
ax[0, 0].text(arr_cd*0.5, ax[0, 0].get_ylim()[1]*0.5 if False else 6, 'permanent\nG0 (arrest)', ha='center',
              fontsize=9, color='#c0392b', fontweight='bold')
ax[0, 0].scatter([mb['cd']], [mb['G0']], s=130, color='#762A83', marker='*', edgecolor='k', zorder=4, label='MB')
ax[0, 0].set_xlabel('CyclinD1 protein (mean, a.u.)'); ax[0, 0].set_ylabel('transient G0 (p27-high), h/cycle')
ax[0, 0].set_title('(a) CyclinD1 sets transient-G0 duration\n(all knobs collapse; arrest below threshold)',
                   fontweight='bold', fontsize=11)
ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.2)

# (b) phase partition vs CyclinD1 (translation sweep, cycling only)
cyc = [a for a in tl_s if not a['arrest']]
xs = [a['cd'] for a in cyc]
bot = np.zeros(len(cyc))
for key, col in [('G0', '#fbeee6'), ('G1', '#aed6f1'), ('S', '#a3e4d7'), ('G2', '#fdebd0')]:
    vals = [a[key] for a in cyc]
    ax[0, 1].bar(range(len(cyc)), vals, bottom=bot, label=key, color=col, edgecolor='k')
    bot += np.array(vals)
ax[0, 1].set_xticks(range(len(cyc))); ax[0, 1].set_xticklabels([f'{c:.2f}' for c in xs])
ax[0, 1].set_xlabel('CyclinD1 protein (a.u.)'); ax[0, 1].set_ylabel('phase duration (h)')
ax[0, 1].set_title('(b) Phase partition vs CyclinD1\n(G0 shrinks into G1/S; period ~constant)',
                   fontweight='bold', fontsize=11)
ax[0, 1].legend(fontsize=8, ncol=4, loc='upper center')

# (c) one-cycle anatomy (p27-marked G0)
r = sim(shh=0.5); t = r['time']; mm = r['mass']
drops = np.where((mm[:-1]-mm[1:]) > 0.2*mm[:-1])[0]; drops = drops[t[drops] >= SETTLE]
a0, a1 = t[drops[0]], t[drops[1]]
seg = (t >= a0) & (t <= a1)
th = (t[seg]-a0)/60.0
def nz(x): x = x[seg]; return x/np.nanmax(x) if np.nanmax(x) > 0 else x
ax[1, 0].plot(th, nz(r['Cd']), color='#117a65', lw=1.8, label='CyclinD1 (norm)')
ax[1, 0].plot(th, nz(r['pRb']), color='#c0392b', lw=1.5, ls='--', label='phospho-Rb (norm; saturates early!)')
ax[1, 0].plot(th, r['P21'][seg], color='#b9770e', lw=2.0, label='p27 (=P21)')
ax[1, 0].plot(th, nz(r['E2f']), color='#2471a3', lw=1.6, label='E2f (norm)')
ax[1, 0].plot(th, nz(r['aRc']), color='#17a2b8', lw=1.4, label='replication aRc (norm)')
# shade p27-high pre-S window
p27seg = r['P21'][seg]; aRcseg = r['aRc'][seg]; Dnaseg = r['Dna'][seg]
preS = ~((aRcseg > 0.05) & (Dnaseg < 0.98)) & (Dnaseg < 0.98)
g0mask = preS & (p27seg > P27)
if g0mask.any():
    g0_end = th[g0mask][-1]
    ax[1, 0].axvspan(0, g0_end, color='#fbeee6', alpha=0.7, zorder=0)
    ax[1, 0].text(g0_end/2, 0.92, f'transient G0\n(p27-high, ~{g0_end:.0f} h)', ha='center',
                  fontsize=9, color='#b9770e', fontweight='bold')
ax[1, 0].set_xlabel('time since division (h)'); ax[1, 0].set_ylabel('level (norm / a.u.)')
ax[1, 0].set_title('(c) One cycle (GNP): the transient G0 is p27-high / E2f-low,\nNOT marked by phospho-Rb',
                   fontweight='bold', fontsize=10.5)
ax[1, 0].legend(fontsize=7.5, loc='center right'); ax[1, 0].grid(alpha=0.2)

# (d) EZH2 -| CyclinD1 -| G0
groups = ['GNP', 'GNP\n+EZH2i', 'MB', 'MB\n+EZH2i']
g0v = [gnp['G0'], gnp_e['G0'], mb['G0'], mb_e['G0']]
cdv = [gnp['cd'], gnp_e['cd'], mb['cd'], mb_e['cd']]
xb = np.arange(4)
b = ax[1, 1].bar(xb, g0v, color=['#1b9e77', '#66c2a5', '#762A83', '#c2a5cf'], edgecolor='k')
for i, (g, c) in enumerate(zip(g0v, cdv)):
    ax[1, 1].text(i, g+0.15, f'{g:.1f} h\n(Cd={c:.2f})', ha='center', fontsize=8, fontweight='bold')
ax[1, 1].set_xticks(xb); ax[1, 1].set_xticklabels(groups)
ax[1, 1].set_ylabel('transient G0 (p27-high), h/cycle')
ax[1, 1].set_title('(d) EZH2 inhibition de-represses CyclinD1\nand collapses the transient G0',
                   fontweight='bold', fontsize=11)
ax[1, 1].grid(alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig('simulations/fig_v44_g0_mechanism.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_mechanism.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_g0_mechanism.png")
print(f"GNP G0={gnp['G0']:.1f}h ({gnp['G0']/gnp['period']*100:.0f}%)  +EZH2i={gnp_e['G0']:.1f}h  | "
      f"MB G0={mb['G0']:.1f}h ({mb['G0']/mb['period']*100:.0f}%)  +EZH2i={mb_e['G0']:.1f}h")
