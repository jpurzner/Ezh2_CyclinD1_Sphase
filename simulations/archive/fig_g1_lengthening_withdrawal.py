"""Withdrawal memory via GRADED G1-LENGTHENING (JP hypothesis 2026-07-24): as Shh falls, cells keep dividing with
progressively LENGTHENING G1 (graded slowdown), NOT abrupt exit. The lever is the latent mu_eff growth-slowdown
(mu_min_frac / K_g1len): growth slows as CyclinD1 falls -> G1 (growth-limited) lengthens. Baked-latent (K_g1len=0.6
below operating Cd); tuning K_g1len=0.8, mu_min_frac=0.1 engages it. VALIDATION-COMPATIBLE (28/29).
JP's EZH2-reduction driver (faster kDeEZ -> EZH2 tracks falling mitogen -> de-represses CyclinD1 -> pool sustains
more divisions) also works but BREAKS validation (23/29, kDeEZ was pinned) -> shown as an option.
3 panels: A period-vs-Shh dose-response (current abrupt vs graded); B withdrawal transient (division rasters +
CyclinD1); C EZH2-reduction option (default-turnover flat vs fast-turnover reduces + more divisions).
Compute caches to _cache.npz. Run: ./venv/bin/python simulations/fig_g1_lengthening_withdrawal.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
CACHE = os.path.join(HERE, 'fig_g1_lengthening_withdrawal_cache.npz'); FORCE = os.environ.get('FORCE', '0') == '1'
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
CURRENT = {}; GRADED = {'mu_min_frac': 0.1, 'K_g1len': 0.8}
kDeEZ0 = 0.0002766221763090641; kTlEZ0 = 0.015098862630302894
EZRED = {'mu_min_frac': 0.1, 'K_g1len': 0.8, 'kDeEZ': kDeEZ0 * 3, 'kTlEZ': kTlEZ0 * 3}
def divs(dna): dna = np.asarray(dna); return np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
def sim_robust(rr, t0, t1, n, sel):
    for atol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        try: rr.integrator.setValue('absolute_tolerance', atol); return rr.simulate(t0, t1, n, selections=sel)
        except Exception: rr.integrator.setValue('maximum_time_step', 1.0)
    return None
def newrr(ov):
    rr = te.loada(build_model_v44(params=ov)); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('absolute_tolerance', 1e-8); rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('maximum_num_steps', 5000000)
    return rr
def period(ov, shh, T=900):
    rr = newrr(ov); rr['SHH'] = shh
    r = sim_robust(rr, 0, T * UPH, int(T * 3), ['time', 'Dna'])
    if r is None: return np.inf
    t = r['time'] / UPH; dvt = t[divs(r['Dna'])]
    if len(dvt) < 3: return np.inf
    iv = np.diff(dvt)
    # steady only if the cell is still dividing at the end (last div within the last cycle-length of T)
    if (T - dvt[-1]) > 1.6 * iv[-1]: return np.inf          # arrested (no recent division)
    return float(np.median(iv[-3:]))                         # settled period = median of last 3 intervals
def transient(ov, ramp=320.0, Twarm=130, Tpost=140):
    rr = newrr(ov); rr['SHH'] = 0.5; sim_robust(rr, 0, Twarm * UPH, 300, ['time', 'Dna'])
    T_, DNA, CD, EZ, SH = [], [], [], [], []; t = Twarm * UPH; chunk = 8.0; Tend = (Twarm + ramp + Tpost) * UPH
    while t < Tend - 1e-6:
        dt = t / UPH - Twarm; shh = max(0.0, 0.5 * (1 - dt / ramp)) if dt < ramp else 0.0; rr['SHH'] = shh
        r = sim_robust(rr, t, t + chunk, 3, ['time', 'Dna', 'Cd', 'EZH2'])
        if r is not None:
            T_ += list(np.asarray(r['time'])[1:]); DNA += list(np.asarray(r['Dna'])[1:]); CD += list(np.asarray(r['Cd'])[1:]); EZ += list(np.asarray(r['EZH2'])[1:]); SH += [shh] * (len(r['time']) - 1)
        t += chunk
    T_ = (np.array(T_) / UPH - Twarm)
    return T_, np.array(CD), np.array(EZ), np.array(SH), T_[divs(np.array(DNA))]

if os.path.exists(CACHE) and not FORCE:
    Z = dict(np.load(CACHE, allow_pickle=True)); print('loaded', CACHE)
else:
    Z = {}; SHHS = np.array([0.5, 0.35, 0.25, 0.18, 0.14, 0.11, 0.09, 0.07])
    Z['shhs'] = SHHS
    for tag, ov in (('cur', CURRENT), ('grad', GRADED)):
        Z['per_' + tag] = np.array([period(ov, s) for s in SHHS]); print('period', tag, 'done', flush=True)
    for tag, ov in (('cur', CURRENT), ('grad', GRADED), ('ezred', EZRED)):
        t, cd, ez, sh, dv = transient(ov)
        Z['t_' + tag] = t; Z['cd_' + tag] = cd; Z['ez_' + tag] = ez; Z['sh_' + tag] = sh; Z['dv_' + tag] = dv
        print('transient', tag, 'done  #div=', len(dv), flush=True)
    np.savez(CACHE, **Z)

C_CUR, C_GRAD, C_EZ = '#c0392b', '#0e6251', '#8e44ad'
def smooth(t, x, w_h=11.0):
    if len(t) < 3: return x
    dt = np.median(np.diff(t)); k = max(1, int(w_h / dt))
    return np.convolve(np.pad(x, k, mode='edge'), np.ones(k) / k, mode='same')[k:-k]
fig, ax = plt.subplots(1, 3, figsize=(20.5, 6.4))
fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.14, wspace=0.26)

# A: interval SEQUENCE — the direct graded-lengthening demonstration (robust)
a = ax[0]
for tag, col, lab in (('cur', C_CUR, 'CURRENT — even intervals then stop'), ('grad', C_GRAD, 'GRADED — intervals RAMP UP (G1 lengthens)')):
    dv = Z['dv_' + tag]; dv = dv[dv > -20]; iv = np.diff(dv)
    a.plot(np.arange(1, len(iv) + 1), iv, '-o', color=col, lw=2.6, ms=7, label=lab)
a.set_xlabel('division number (during withdrawal)'); a.set_ylabel('inter-division interval Tc (h)')
a.text(0.04, 0.80, 'as mitogen falls, the GRADED model lengthens each\ncycle (growth slows with CyclinD1) — "lengthen G1\nbut not exit"; CURRENT keeps ~even cycles then stops', transform=a.transAxes, fontsize=8.2, color='#555')
a.set_title('A · Division intervals RAMP UP during withdrawal\n(graded G1-lengthening) vs even-then-stop', fontsize=11, fontweight='bold')
a.legend(fontsize=8.2, loc='lower right'); a.grid(alpha=0.2)

# B: CyclinD1 (cycle-smoothed) + Shh during withdrawal
b = ax[1]
for tag, col, lab in (('cur', C_CUR, 'CyclinD1 — current'), ('grad', C_GRAD, 'CyclinD1 — graded')):
    t = Z['t_' + tag]; w = t > -8
    b.plot(t[w], smooth(t, Z['cd_' + tag])[w], color=col, lw=2.4, label=lab)
    for x in Z['dv_' + tag]:
        if x > -8: b.plot(x, -0.35 - 0.35 * ('cur grad'.split().index(tag)), '|', color=col, ms=12, mew=2.0)
b.plot(Z['t_grad'], Z['sh_grad'] * 3, color='#95a5a6', lw=1.1, ls='--', label='Shh × 3 (withdrawal ramp)')
b.axhline(0, color='#ccc', lw=0.8)
b.text(0.03, 0.10, 'ticks = divisions (cycle-smoothed CyclinD1). GRADED keeps CyclinD1\ndeclining smoothly, dividing slower; both stay Shh-DEPENDENT (arrest)', transform=b.transAxes, fontsize=7.8, color='#555')
b.set_xlabel('time from start of withdrawal (h)'); b.set_ylabel('CyclinD1 (smoothed) / Shh×3'); b.set_ylim(-0.95, None); b.set_xlim(-8, 300)
b.set_title('B · CyclinD1 + Shh during the gradual withdrawal', fontsize=11, fontweight='bold')
b.legend(fontsize=7.9, loc='upper right'); b.grid(alpha=0.2, axis='x')

# C: the EZH2-reduction TENSION — EZH2 sampled at each division (consistent phase → clean trend)
c = ax[2]
t = Z['t_grad']; dv = Z['dv_grad']; dv = dv[dv > -8]
ezd = np.array([Z['ez_grad'][np.argmin(np.abs(t - x))] for x in dv])
c.plot(dv, ezd / ezd[0], '-o', color=C_GRAD, lw=2.6, ms=7, label='EZH2 per cycle (graded, validation-OK)')
c.axhline(1.0, color='#c0392b', lw=1.4, ls='--'); c.text(150, 1.005, 'flat — EZH2 does NOT reduce', color='#c0392b', fontsize=8, ha='center')
c.text(0.04, 0.10, "JP's mechanism NEEDS EZH2 to REDUCE during withdrawal\n(to de-repress the CyclinD1 pool). But the model BANKS EZH2\n(slow τ~14 h + it tracks the cycle, not the mitogen) → it stays\nFLAT while the cell keeps dividing. Making it reduce needs faster\nturnover, which breaks validation (23/29). So the phenotype comes\nfrom GROWTH-slowdown (A/B), NOT EZH2 reduction — JP's specific\ndriver is UNRESOLVED (would need a joint re-opt).", transform=c.transAxes, fontsize=7.6, color='#555')
c.set_xlabel('time from start of withdrawal (h)'); c.set_ylabel('EZH2 at each division (norm. to 1st)'); c.set_xlim(-8, 175); c.set_ylim(0.9, 1.06)
c.set_title('C · The EZH2-reduction TENSION — the calibrated model\nBANKS EZH2 (flat), contra JP\'s required driver', fontsize=11, fontweight='bold')
c.legend(fontsize=8.2, loc='lower left'); c.grid(alpha=0.2)

fig.suptitle('Withdrawal memory via GRADED G1-LENGTHENING: the latent mu_eff growth-slowdown (mu_min_frac 0.4→0.1, K_g1len 0.6→0.8) makes division intervals RAMP UP as mitogen falls, then arrest '
             '— "lengthen G1 but not exit", VALIDATION-COMPATIBLE (28/29). JP\'s EZH2-reduction DRIVER is not supported as-is (model banks EZH2; forcing it breaks validation)', fontsize=10.6, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_g1_lengthening_withdrawal.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
print('intervals cur:', np.round(np.diff(Z['dv_cur'][Z['dv_cur']>-20]),0))
print('intervals grad:', np.round(np.diff(Z['dv_grad'][Z['dv_grad']>-20]),0))
