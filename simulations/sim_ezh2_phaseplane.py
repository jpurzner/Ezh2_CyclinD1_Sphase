"""Ferrell-style phase-plane / bifurcation analysis of the EZH2 -| CyclinD1 feedback (v44 baked).

The feedback is a two-node loop: CyclinD1 -> (Rb-E2f / cell cycle, + mitogen-dose) -> EZH2 -| CyclinD1.
CyclinD1 is FAST (~1 min) with no cycle-dependent production except EZH2 repression -> its nullcline
is exact & analytic. EZH2 is SLOW (~4.6 cycles) -> its nullcline is traced numerically by breaking the
loop (EZH2i=1) and sweeping the CyclinD1 drive, so EZH2 becomes a pure readout of Cd.

  Cd-nullcline (dCd/dt=0):   Cd = A * Krep/(Krep+EZH2),  A = ktl*D/(k_Cd_deg*k_Cd_mRNA_deg)
                             D = mitogen drive (Gli/MYCN, constant per condition)   [decreasing]
  EZH2-nullcline (dEZH2/dt=0): feedback OFF, sweep ktl -> (mean Cd, mean EZH2)        [increasing]
  intersection = operating point; the real limit cycle (feedback ON) orbits it.

Panels:
  A  GNP nullcline portrait: both nullclines, fast/slow direction field, limit-cycle orbit, operating
     point, and the (un-repressed) no-feedback point -> the negative feedback is monostable & buffering.
  B  condition overlay: GNP / MB / +EZH2i -- how the mitogen drive (D) and feedback removal move the
     operating point along / off the nullclines.
  C  bifurcation vs CyclinD1 drive (ktl), +/- feedback, up & down continuation -> saddle-node /
     hysteresis: where the proliferation-quiescence switch sits and how the feedback shifts it.
  D  bifurcation vs mitogen (SHH, in Gli1 units), +/- feedback -> the feedback raises the mitogen
     threshold for proliferation (the mitogen-gated quiescence boundary).

Cache: fig_v44_ezh2_phaseplane_cache.npz  (recompute with --fresh).  ~30-40 min fresh.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_ezh2_phaseplane_cache.npz'

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)

C = {k: rr[k] for k in ['k_Cd_translation', 'k_Cd_deg', 'k_Cd_mRNA_deg', 'K_EZH2_repression', 'kDeEZ',
                         'k_Cd_tx_basal', 'k_Cd_tx_Gli_max', 'K_Gli_act_CycD', 'K_Gli_rep_CycD',
                         'n_Gli_act', 'n_Gli_rep', 'k_Cd_tx_MYCN', 'K_MYCN_Cd', 'n_MYCN_Cd']}
KREP = C['K_EZH2_repression']

GNP = dict(shh=0.5, ezh2i=0, mycn_amp=1.0, ptch1=1.0, p16=0.0, p18=0.464, ksyp21=0.002)
MB = dict(shh=0.5, ezh2i=0, mycn_amp=2.8, ptch1=0.1, p16=0.306, p18=1.553, ksyp21=0.002)


def drive_D(ga, g1, grep, mycn):
    gli = ga + g1
    gli_term = C['k_Cd_tx_Gli_max'] * gli**C['n_Gli_act'] / (C['K_Gli_act_CycD']**C['n_Gli_act'] + gli**C['n_Gli_act']) \
        * (C['K_Gli_rep_CycD']**C['n_Gli_rep'] / (C['K_Gli_rep_CycD']**C['n_Gli_rep'] + grep**C['n_Gli_rep']))
    mycn_term = C['k_Cd_tx_MYCN'] * mycn**C['n_MYCN_Cd'] / (C['K_MYCN_Cd']**C['n_MYCN_Cd'] + mycn**C['n_MYCN_Cd'])
    return C['k_Cd_tx_basal'] + gli_term + mycn_term


SEL = ['time', 'Cd', 'EZH2', 'Gli_act', 'Gli1', 'Gli_rep', 'MYCN', 'MPF']


def run(shh, ezh2i, mycn_amp, ptch1, p16, p18, ksyp21, ktl=None, T=14000, settle=8000, reset=True):
    if reset:
        rr.reset()
    rr['SHH'] = shh; rr['EZH2i'] = ezh2i; rr['MYCN_amplification'] = mycn_amp
    rr['Ptch1_copy_number'] = ptch1; rr['p16'] = p16; rr['p18'] = p18; rr['kSyP21'] = ksyp21
    if ktl is not None:
        rr['k_Cd_translation'] = ktl
    r = rr.simulate(0, T, int(T / 0.5), selections=SEL)
    m = r['time'] >= settle
    return {k: r[k][m] for k in SEL}


SEL_BIF = ['time', 'MPF', 'E2f']


def _meas(seg):
    """(proliferation rate div/100h, mean E2f) over a settled segment."""
    t = seg['time']; dt = t[1] - t[0]
    pk, _ = find_peaks(seg['MPF'], prominence=0.15, distance=int(200 / dt))
    hours = (t[-1] - t[0]) / 60.0
    rate = len(pk) / hours * 100.0
    return rate, float(np.mean(seg['E2f']))


def cond_pack(cond, label):
    """operating point + analytic Cd-nullcline coeff A for a condition (feedback ON)."""
    o = run(**cond)
    D = drive_D(o['Gli_act'].mean(), o['Gli1'].mean(), o['Gli_rep'].mean(), o['MYCN'].mean())
    A = C['k_Cd_translation'] * D / (C['k_Cd_deg'] * C['k_Cd_mRNA_deg'])
    return dict(label=label, cd=o['Cd'], ez=o['EZH2'], cd_op=float(o['Cd'].mean()),
                ez_op=float(o['EZH2'].mean()), A=float(A), D=float(D))


def ezh2_nullcline(base, ktls, T=14000, settle=8000):
    """feedback OFF (ezh2i=1), sweep ktl -> (meanCd, meanEZH2)."""
    cds, ezs = [], []
    for k in ktls:
        c = dict(base); c['ezh2i'] = 1
        o = run(**c, ktl=k, T=T, settle=settle)
        cds.append(float(o['Cd'].mean())); ezs.append(float(o['EZH2'].mean()))
    return np.array(cds), np.array(ezs)


_PMAP = {'shh': 'SHH', 'ezh2i': 'EZH2i', 'mycn_amp': 'MYCN_amplification',
         'ptch1': 'Ptch1_copy_number', 'p16': 'p16', 'p18': 'p18', 'ksyp21': 'kSyP21'}


def _apply(base):
    for k, key in _PMAP.items():
        rr[key] = base[k]


def _last_seg(r, frac):
    m = r['time'] >= r['time'][-1] - (r['time'][-1] - r['time'][0]) * frac
    return {'time': r['time'][m], 'MPF': r['MPF'][m], 'E2f': r['E2f'][m]}


def bifurc(base, param, values, T=8000, settle_frac=0.4):
    """reset-per-point sweep in `param` (robust; each point from the standard IC, no continuation
    stiffness). Returns (proliferation-rate array, mean-E2f array). mean E2f = smooth commitment
    order parameter (Rb-E2f), the clean Ferrell response readout."""
    rates, e2fs = [], []
    for v in values:
        rr.reset(); _apply(base); rr[param] = v
        try:
            r = rr.simulate(0, T, int(T / 0.5), selections=SEL_BIF)
            seg = _last_seg(r, settle_frac)
        except Exception:
            try:
                rr.reset(); _apply(base); rr[param] = v
                r = rr.simulate(0, 12000, 24000, selections=SEL_BIF)
                seg = _last_seg(r, 0.35)
            except Exception:
                rates.append(np.nan); e2fs.append(np.nan); continue
        rt, ef = _meas(seg)
        rates.append(rt); e2fs.append(ef)
    return np.array(rates), np.array(e2fs)


# ------------------------------------------------------------------ compute / cache
NCACHE = 'simulations/fig_v44_ezh2_phaseplane_stage1_cache.npz'

# stage 1 (expensive ~16 min): operating points + nullclines -> cached separately so a later
# bifurcation failure can't discard them.
if FRESH or not os.path.exists(NCACHE):
    print("stage 1: operating points + nullclines...")
    gnp = cond_pack(GNP, 'GNP')
    mb = cond_pack(MB, 'MB')
    gnp_ezi = cond_pack(dict(GNP, ezh2i=1), 'GNP+EZH2i')
    ktls_null = np.geomspace(0.08, 4.5, 16)
    nc_gnp_cd, nc_gnp_ez = ezh2_nullcline(GNP, ktls_null)
    nc_mb_cd, nc_mb_ez = ezh2_nullcline(MB, ktls_null)
    np.savez(NCACHE,
             gnp_cd=gnp['cd'], gnp_ez=gnp['ez'], gnp_op=[gnp['cd_op'], gnp['ez_op'], gnp['A']],
             mb_cd=mb['cd'], mb_ez=mb['ez'], mb_op=[mb['cd_op'], mb['ez_op'], mb['A']],
             gnpezi_cd=gnp_ezi['cd'], gnpezi_ez=gnp_ezi['ez'], gnpezi_op=[gnp_ezi['cd_op'], gnp_ezi['ez_op'], gnp_ezi['A']],
             ktls_null=ktls_null, nc_gnp_cd=nc_gnp_cd, nc_gnp_ez=nc_gnp_ez, nc_mb_cd=nc_mb_cd, nc_mb_ez=nc_mb_ez)
    print("  stage-1 cached ->", NCACHE)
nz = np.load(NCACHE)

# stage 2: bifurcations (robust reset-per-point; mean-E2f response curves)
if FRESH or not os.path.exists(CACHE):
    print("stage 2: bifurcation sweeps...")
    ktl_bif = np.linspace(0.08, 1.3, 16)
    cF_rate, cF_e2f = bifurc(GNP, 'k_Cd_translation', ktl_bif)                 # feedback ON
    cN_rate, cN_e2f = bifurc(dict(GNP, ezh2i=1), 'k_Cd_translation', ktl_bif)  # feedback OFF
    shh_bif = np.linspace(0.0, 3.0, 16)
    dF_rate, dF_e2f = bifurc(GNP, 'SHH', shh_bif)
    dN_rate, dN_e2f = bifurc(dict(GNP, ezh2i=1), 'SHH', shh_bif)
    gli1_at = []
    for s in shh_bif:
        o = run(**dict(GNP, shh=s), T=6000, settle=4000)
        gli1_at.append(float((o['Gli_act'] + o['Gli1']).mean()))
    gli1_at = np.array(gli1_at)
    np.savez(CACHE, **{k: nz[k] for k in nz.files},
             ktl_bif=ktl_bif, cF_rate=cF_rate, cF_e2f=cF_e2f, cN_rate=cN_rate, cN_e2f=cN_e2f,
             shh_bif=shh_bif, gli1_at=gli1_at, dF_rate=dF_rate, dF_e2f=dF_e2f, dN_rate=dN_rate, dN_e2f=dN_e2f,
             KREP=KREP, kDeEZ=C['kDeEZ'], k_Cd_deg=C['k_Cd_deg'])
    print("  cached ->", CACHE)

z = np.load(CACHE)
KREP = float(z['KREP'])

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(2, 2, figsize=(15, 12))


def cd_nullcline(A, ez):
    return A * KREP / (KREP + ez)


# ---- Panel A: GNP nullcline portrait ----
a = ax[0, 0]
cd_op, ez_op, A = z['gnp_op']
ez_axis = np.linspace(0.001, 2.2, 200)
a.plot(cd_nullcline(A, ez_axis), ez_axis, '-', color='#c0392b', lw=2.5, label='Cd-nullcline (EZH2 represses Cd)')
order = np.argsort(z['nc_gnp_cd'])
a.plot(z['nc_gnp_cd'][order], z['nc_gnp_ez'][order], '-', color='#2471a3', lw=2.5, label='EZH2-nullcline (Cd induces EZH2)')
a.plot(z['nc_gnp_cd'][order], z['nc_gnp_ez'][order], 'o', color='#2471a3', ms=4)
# fast/slow direction field (normalized)
from scipy.interpolate import interp1d
ez_of_cd = interp1d(z['nc_gnp_cd'][order], z['nc_gnp_ez'][order], bounds_error=False, fill_value=(z['nc_gnp_ez'][order][0], z['nc_gnp_ez'][order][-1]))
gx, gy = np.meshgrid(np.linspace(0.3, 6.0, 22), np.linspace(0.05, 2.1, 18))
dC = float(z['k_Cd_deg']) * (cd_nullcline(A, gy) - gx)          # fast
dE = float(z['kDeEZ']) * (ez_of_cd(gx) - gy)                    # slow
# scale to comparable visual magnitude per-component then normalize direction
u = dC / (np.abs(dC).max() + 1e-9); v = dE / (np.abs(dE).max() + 1e-9)
nrm = np.hypot(u, v) + 1e-9
a.quiver(gx, gy, u / nrm, v / nrm, color='#95a5a6', alpha=0.55, width=0.0035, scale=28, pivot='mid')
a.plot(z['gnp_cd'], z['gnp_ez'], '-', color='#27ae60', lw=1.3, alpha=0.9, label='limit-cycle orbit (feedback ON)')
a.plot(cd_op, ez_op, '*', color='k', ms=20, zorder=6, label=f'operating point ({cd_op:.2f}, {ez_op:.2f})')
# no-feedback point: where Cd sits with repression removed (on EZH2-nullcline at full drive A)
a.plot(A, float(ez_of_cd(A)), 'D', color='#e67e22', ms=10, zorder=6, label=f'no-feedback Cd ({A:.2f})')
a.annotate('', xy=(cd_op, ez_op), xytext=(A, float(ez_of_cd(A))),
           arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2, ls='--'))
a.text((A + cd_op) / 2, float(ez_of_cd(A)) + 0.08, 'EZH2 buffers\nCyclinD1', color='#e67e22', fontsize=9, ha='center')
a.set_xlim(0, 6.2); a.set_ylim(0, 2.2)
a.set_xlabel('CyclinD1 (Cd)', fontsize=11); a.set_ylabel('EZH2', fontsize=11)
a.set_title('(A) EZH2–CyclinD1 phase plane (GNP)\nnegative feedback → single stable node (monostable, buffering)', fontweight='bold', fontsize=11)
a.legend(fontsize=8, loc='upper right'); a.grid(alpha=0.15)

# ---- Panel B: condition overlay ----
b = ax[0, 1]
for pack_cd, pack_op, col, lab, ncd, ncz in [
    ('gnp_cd', 'gnp_op', '#27ae60', 'GNP', 'nc_gnp_cd', 'nc_gnp_ez'),
    ('mb_cd', 'mb_op', '#8e44ad', 'MB', 'nc_mb_cd', 'nc_mb_ez'),
    ('gnpezi_cd', 'gnpezi_op', '#e67e22', 'GNP+EZH2i', 'nc_gnp_cd', 'nc_gnp_ez')]:
    cdv, ezv, Av = z[pack_op]
    b.plot(cd_nullcline(Av, ez_axis), ez_axis, '--', color=col, lw=1.6, alpha=0.8)
    b.plot(z[pack_cd], z['%s' % pack_cd.replace('_cd', '_ez')], '-', color=col, lw=0.6, alpha=0.18)
    b.plot(cdv, ezv, '*', color=col, ms=18, zorder=6, label=f'{lab}  ({cdv:.2f}, {ezv:.2f})')
ordm = np.argsort(z['nc_mb_cd'])
b.plot(z['nc_mb_cd'][ordm], z['nc_mb_ez'][ordm], '-', color='#2471a3', lw=2.0, alpha=0.9, label='EZH2-nullcline (MB CKIs)')
ordg = np.argsort(z['nc_gnp_cd'])
b.plot(z['nc_gnp_cd'][ordg], z['nc_gnp_ez'][ordg], '-', color='#5dade2', lw=2.0, alpha=0.9, label='EZH2-nullcline (GNP CKIs)')
b.set_xlim(0, 12); b.set_ylim(0, 2.2)
b.set_xlabel('CyclinD1 (Cd)', fontsize=11); b.set_ylabel('EZH2', fontsize=11)
b.set_title('(B) Mitogen drive & feedback removal move the operating point\n(dashed = Cd-nullcline per condition; EZH2i → un-repressed)', fontweight='bold', fontsize=11)
b.legend(fontsize=8, loc='lower right'); b.grid(alpha=0.15)

def bifpoint(x, rate):
    """first x where proliferation turns on (rate > 0.5 div/100h) = the bifurcation point."""
    idx = np.where(np.nan_to_num(rate) > 0.5)[0]
    return float(x[idx[0]]) if len(idx) else None


# ---- Panel C: bifurcation vs CyclinD1 drive, +/- feedback ----
c = ax[1, 0]
ktl = z['ktl_bif']
c.plot(ktl, z['cN_e2f'], '-o', color='#7f8c8d', lw=2, ms=4, label='no feedback (EZH2i)')
c.plot(ktl, z['cF_e2f'], '-o', color='#c0392b', lw=2, ms=4, label='EZH2 feedback ON')
bN, bF = bifpoint(ktl, z['cN_rate']), bifpoint(ktl, z['cF_rate'])
if bN:
    c.axvline(bN, color='#7f8c8d', ls='--', alpha=0.8)
if bF:
    c.axvline(bF, color='#c0392b', ls='--', alpha=0.8)
if bN and bF:
    yy = np.nanmax(z['cN_e2f']) * 0.5
    c.annotate('', xy=(bF, yy), xytext=(bN, yy), arrowprops=dict(arrowstyle='->', color='k', lw=1.6))
    c.text((bN + bF) / 2, yy * 1.06, f'Δ threshold\n+{bF - bN:.2f}', ha='center', fontsize=8.5, fontweight='bold')
c.axvline(C['k_Cd_translation'], color='k', ls=':', alpha=0.4)
c.text(C['k_Cd_translation'], np.nanmax(z['cN_e2f']) * 0.03, ' baked', fontsize=7.5)
c.set_xlabel('CyclinD1 drive  (k_Cd_translation)', fontsize=11)
c.set_ylabel('commitment readout  (mean E2F)', fontsize=11)
c.set_title('(C) Proliferation–quiescence bifurcation vs CyclinD1 drive\nEZH2 feedback shifts the threshold to higher drive (dashed = onset)', fontweight='bold', fontsize=11)
c.legend(fontsize=8, loc='lower right'); c.grid(alpha=0.15)

# ---- Panel D: bifurcation vs mitogen (SHH / Gli1), +/- feedback ----
d = ax[1, 1]
gli1 = z['gli1_at']
d.plot(gli1, z['dN_e2f'], '-o', color='#7f8c8d', lw=2, ms=4, label='no feedback (EZH2i)')
d.plot(gli1, z['dF_e2f'], '-o', color='#c0392b', lw=2, ms=4, label='EZH2 feedback ON')
gN, gF = bifpoint(gli1, z['dN_rate']), bifpoint(gli1, z['dF_rate'])
if gN:
    d.axvline(gN, color='#7f8c8d', ls='--', alpha=0.8)
if gF:
    d.axvline(gF, color='#c0392b', ls='--', alpha=0.8)
if gN and gF:
    d.axvspan(gN, gF, color='#c0392b', alpha=0.08)
    d.text((gN + gF) / 2, np.nanmax(z['dN_e2f']) * 0.5, 'feedback-\ngated\nquiescence', ha='center', fontsize=8, color='#c0392b')
d.set_xlabel('mitogen  (Gli1 transcript, Hh activity)', fontsize=11)
d.set_ylabel('commitment readout  (mean E2F)', fontsize=11)
d.set_title('(D) Mitogen bifurcation: EZH2 feedback raises the\nmitogen threshold for proliferation (gated quiescence)', fontweight='bold', fontsize=11)
d.legend(fontsize=8, loc='lower right'); d.grid(alpha=0.15)

plt.tight_layout()
plt.savefig('simulations/fig_v44_ezh2_phaseplane.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ezh2_phaseplane.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_ezh2_phaseplane.png / .pdf")
