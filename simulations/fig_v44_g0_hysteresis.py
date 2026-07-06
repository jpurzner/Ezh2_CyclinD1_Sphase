"""Does the H3K27me3 mark make G0 arrest BISTABLE?  Hysteresis test.

*** SUPERSEDED (2026-07) by fig_v44_g0_bistability_confirm.py ***
This quasi-static SHH sweep gave a FALSE-POSITIVE hysteresis in the ACCUMULATE regime: the mark there is
so slow (EZH2 t1/2 ~230 h) that a ~200 h/step sweep cannot equilibrate the arrested branch, so it LOOKS
locked but is only a long transient. The rigorous two-initial-condition test (equilibrate each SHH from
both a high-mark and a low-mark start to true steady state) shows G0 is MONOSTABLE / reversible
(metastable "sticky" arrest, not a bistable lock). Kept only as a methodological cautionary record --
do NOT cite as bistability evidence.

For each parameter setting we sweep the mitogen drive (SHH) through the arrest threshold along TWO
branches that carry their state (no reset between steps):
  - ASCENDING  : start ARRESTED with a fully built-up mark (low SHH, long pre-equilibration), raise SHH.
  - DESCENDING : start CYCLING with a diluted (low) mark (high SHH), lower SHH.
If the two branches diverge over a range of SHH -> two coexisting stable states -> BISTABLE (a real
epigenetic G0 lock, hysteresis width = the gap). If they coincide -> MONOSTABLE (reversible G0).

Settings:
  (1) calibrated (RELAX)   kDeEZ=1.5e-4, del_mk=7e-4, mark ON  -> expect ~monostable
  (2) mark OFF (control)   same, f0_mk=1.0 (mark cannot repress) -> the INTRINSIC commitment hysteresis
  (3) ACCUMULATE           kDeEZ=5e-5 (stable EZH2), del_mk=3e-4 (slow turnover) -> expect BISTABLE
The mark's contribution = hysteresis in (1)/(3) BEYOND the commitment baseline (2).

Run:  ./venv/bin/python simulations/fig_v44_g0_hysteresis.py [--mini] [--fresh]
Cache: simulations/fig_v44_g0_hysteresis_cache.npz
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
CACHE = 'simulations/fig_v44_g0_hysteresis_cache.npz' if not MINI else 'simulations/fig_v44_g0_hysteresis_mini.npz'

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SETTINGS = [
    dict(name='calibrated (RELAX)',  kDeEZ=1.5e-4, del_mk=7e-4, f0=0.233, col='#2c7fb8'),
    dict(name='mark OFF (control)',  kDeEZ=1.5e-4, del_mk=7e-4, f0=1.00, col='#7f7f7f'),
    dict(name='ACCUMULATE regime',   kDeEZ=5.0e-5, del_mk=3e-4, f0=0.233, col='#c0392b'),
]
SHH = np.round(np.linspace(0.05, 1.0, 7 if MINI else 16), 3)
SETTLE = 6000 if MINI else 12000        # min per SHH step (carry state)
PREEQ_ARR = 12000 if MINI else 26000    # arrested pre-equilibration (build the mark)
PREEQ_CYC = 8000 if MINI else 15000     # cycling pre-equilibration (dilute the mark)


def rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
    try: r.integrator.setValue('maximum_num_steps', 200000)
    except Exception: pass
    return r


def _sim(r, t, dt):
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            return r.simulate(t, t + dt, max(60, int(dt / 3)), selections=['time', 'MPF', 'Mk', 'Cd'])
        except Exception:
            continue
    return None


def branch(kdeez, delmk, f0, ascending):
    r = rr(); r.reset()
    for k, v in GNP.items(): r[k] = v
    r['kDeEZ'] = float(kdeez); r['del_mk'] = float(delmk); r['f0_mk'] = float(f0)
    shh_seq = SHH if ascending else SHH[::-1]
    # pre-equilibrate on the appropriate branch
    r['SHH'] = float(shh_seq[0])
    t = 0.0
    _sim(r, t, PREEQ_ARR if ascending else PREEQ_CYC); t += (PREEQ_ARR if ascending else PREEQ_CYC)
    rate = {}; mk = {}; cd = {}
    for shh in shh_seq:
        r['SHH'] = float(shh)
        res = _sim(r, t, SETTLE); t += SETTLE
        if res is None:
            rate[shh] = np.nan; mk[shh] = np.nan; cd[shh] = np.nan; continue
        w = res['time'] >= res['time'][0] + SETTLE * 0.4      # score the settled tail
        pk, _ = find_peaks(res['MPF'][w], prominence=0.15, distance=200)
        hrs = (res['time'][w][-1] - res['time'][w][0]) / 60.0
        rate[shh] = len(pk) / hrs * 168.0 if hrs > 0 else 0.0
        mk[shh] = float(np.nanmean(res['Mk'][w])); cd[shh] = float(np.nanmean(res['Cd'][w]))
    return (np.array([rate[s] for s in SHH]), np.array([mk[s] for s in SHH]), np.array([cd[s] for s in SHH]))


if FRESH or not os.path.exists(CACHE):
    RATE_A, RATE_D, MK_A, MK_D, CD_A, CD_D = ([] for _ in range(6))
    for s in SETTINGS:
        ra, ma, ca = branch(s['kDeEZ'], s['del_mk'], s['f0'], ascending=True)
        rd, md, cd = branch(s['kDeEZ'], s['del_mk'], s['f0'], ascending=False)
        RATE_A.append(ra); RATE_D.append(rd); MK_A.append(ma); MK_D.append(md); CD_A.append(ca); CD_D.append(cd)
        print(f"  {s['name']:22s} ascending div-rate {np.round(ra,1)}", flush=True)
        print(f"  {' '*22} descending div-rate {np.round(rd,1)}", flush=True)
    np.savez(CACHE, SHH=SHH, RATE_A=RATE_A, RATE_D=RATE_D, MK_A=MK_A, MK_D=MK_D, CD_A=CD_A, CD_D=CD_D,
             names=[s['name'] for s in SETTINGS], cols=[s['col'] for s in SETTINGS])
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE, allow_pickle=True)
SHH = z['SHH']; RATE_A, RATE_D = z['RATE_A'], z['RATE_D']; MK_A, MK_D = z['MK_A'], z['MK_D']
names = list(z['names']); cols = list(z['cols'])

CD_A, CD_D = z['CD_A'], z['CD_D']
THR_CD = 1.8      # CyclinD1 order-parameter midpoint: > = cycling, < = arrested
fig, axes = plt.subplots(2, 3, figsize=(15, 8.4))
for j, (nm, cl) in enumerate(zip(names, cols)):
    # hysteresis width from the continuous Cd order parameter (smoother than the quantized div-count)
    def cross(cd_branch):
        idx = np.where(cd_branch >= THR_CD)[0]
        return SHH[idx[0]] if len(idx) else np.nan
    s_up = cross(CD_A[j])              # ascending (from arrested): SHH where it ESCAPES to high-Cd cycling
    s_dn = cross(CD_D[j])             # descending (from cycling): SHH where it drops below (COLLAPSES)
    width = (s_up - s_dn) if (s_up == s_up and s_dn == s_dn) else np.nan
    bist = (width == width) and (width > 0.02)

    a = axes[0, j]
    a.plot(SHH, CD_D[j], '-o', color='#2980b9', ms=5, lw=2, label='descending (from cycling)')
    a.plot(SHH, CD_A[j], '-s', color='#c0392b', ms=5, lw=2, label='ascending (from arrested)')
    a.axhline(THR_CD, color='k', ls=':', lw=0.9, alpha=0.5)
    if bist:
        a.axvspan(s_dn, s_up, color='gold', alpha=0.3)
        a.text(0.5 * (s_dn + s_up), THR_CD * 1.35, f'BISTABLE\nΔSHH={width:.2f}', ha='center', fontsize=9, color='#8a6d00', fontweight='bold')
    else:
        a.text(0.5, 0.9, 'monostable', transform=a.transAxes, ha='center', fontsize=9, color='#555', style='italic')
    a.set_title(f'{nm}', fontsize=11, fontweight='bold', color=cl)
    a.set_ylabel('mean CyclinD1  Cd' if j == 0 else ''); a.grid(alpha=0.15)
    if j == 0: a.legend(fontsize=8, loc='upper left')

    b = axes[1, j]
    b.plot(SHH, MK_D[j], '-o', color='#2980b9', ms=5, lw=2)
    b.plot(SHH, MK_A[j], '-s', color='#c0392b', ms=5, lw=2)
    if bist:
        b.axvspan(s_dn, s_up, color='gold', alpha=0.3)
    b.set_xlabel('mitogen drive  SHH'); b.set_ylabel('mean H3K27me3  Mk' if j == 0 else ''); b.grid(alpha=0.15)
fig.suptitle('Does the H3K27me3 mark make G0 arrest bistable?  —  hysteresis of the arrested (↑) vs cycling (↓) branch',
             fontsize=13, fontweight='bold', y=0.98)
fig.text(0.5, 0.005, 'Two branches carry their state through the SHH sweep. Overlap (gold) = a range of drive where BOTH a cycling and an arrested state are stable = bistable G0 lock.',
         ha='center', fontsize=9, color='#555')
plt.tight_layout(rect=[0, 0.01, 1, 0.97])
plt.savefig('simulations/fig_v44_g0_hysteresis.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_hysteresis.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_g0_hysteresis.png')
