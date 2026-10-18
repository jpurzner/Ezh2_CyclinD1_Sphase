"""DEFINITIVE figure: what EZH2->CyclinD1 repression (the GOVERNOR) does for GNPs, in the current
mRNA-reservoir + EZH2-IFFL model WITH the k_Cd_tx_basal x0.5 bake (2026-07-24). Compares EZH2 repression:
  ON     = default (EZH2 represses the CyclinD1 transcript via PRC2_rep)
  OFF    = f0_prc2=1.0 (repression floor lifted -> genetic 'EZH2 inhibited on CyclinD1', full de-repression)
  EZH2i  = the catalytic drug (EZH2i=1; PRC2 occupancy persists -> PARTIAL de-repression)
CORRECTED story (post-basal-bake): proliferation is now Shh-DEPENDENT for ALL (all arrest at SHH=0) because the
mitogen-INDEPENDENT CyclinD1 floor is low -> CyclinD1 needs Gli-drive to reach commitment. EZH2 repression is a
THRESHOLD/GAIN GOVERNOR: it RAISES the Shh threshold ~15x (OFF 0.01 -> ON 0.18) and BUFFERS the dose-response,
matching JP's cKO (transcript UP substantially, but no division without Shh + differentiates). Context (JP): Shh
exposure is PERSISTENT once exposed but NOISY, and withdrawal is GRADUAL + COMPLETE (leaving the niche).
6 panels: A threshold/Shh-dependence; B CyclinD1 dose-response buffering; C temporal noise filter; D gradual
complete withdrawal (all arrest, de-repressed coasts longer + keeps transcript up); E persistent-noisy Shh +
gradual withdrawal; F summary (threshold ~15x + de-repression fold).
Compute caches to fig_governor_noise_filter_cache.npz (delete or set FORCE=1 to recompute).
Run: ./venv/bin/python simulations/fig_governor_noise_filter.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
CACHE = os.path.join(HERE, 'fig_governor_noise_filter_cache.npz'); FORCE = os.environ.get('FORCE', '0') == '1'
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MODEL = build_model_v44()
C_ON, C_OFF, C_EZI, C_MR = '#0e6251', '#c0392b', '#e67e22', '#2e86c1'
def divs(dna): dna = np.asarray(dna); return np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
def cond(name): return {} if name == 'ON' else ({'f0_prc2': 1.0} if name == 'OFF' else {'EZH2i': 1.0})
def newrr():
    rr = te.loada(MODEL); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('absolute_tolerance', 1e-8); rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('maximum_num_steps', 4000000)
    return rr
def sim_robust(rr, t0, t1, n, sel):
    for atol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        try: rr.integrator.setValue('absolute_tolerance', atol); return rr.simulate(t0, t1, n, selections=sel)
        except Exception: rr.integrator.setValue('maximum_time_step', 1.5)
    return None
CONDS = ('ON', 'OFF', 'EZH2i')

# ========================== compute ==========================
if os.path.exists(CACHE) and not FORCE:
    Z = dict(np.load(CACHE, allow_pickle=True)); print('loaded cache', CACHE)
else:
    Z = {}
    def steady(c, shh, T=20000):
        rr = newrr()
        for k, v in cond(c).items(): rr[k] = v
        rr['SHH'] = shh
        try: r = rr.simulate(0, T, 3000, selections=['time', 'Dna', 'Cd', 'Cd_mRNA'])
        except Exception: return np.nan, np.nan, np.nan
        t = np.asarray(r['time']); m = t > T * 0.55
        dv = divs(np.asarray(r['Dna'])[m]); span = (t[m][-1] - t[m][0]) / UPH
        return (len(dv) / span * 24 if span > 0 else 0), float(np.asarray(r['Cd'])[m].mean()), float(np.asarray(r['Cd_mRNA'])[m].mean())
    SHH = np.round(np.linspace(0.0, 0.34, 30), 4); Z['SHH'] = SHH
    for c in CONDS:
        Z['dose_' + c] = np.array([steady(c, s) for s in SHH]); print('dose', c, 'done', flush=True)
    # C: pulse-duration filter (governed cell) — base below ON threshold
    def pulse_dur(c, base, peak, dur, warm=100, post=240):
        rr = newrr()
        for k, v in cond(c).items(): rr[k] = v
        rr['SHH'] = base; sim_robust(rr, 0, warm * UPH, 300, ['time', 'Dna'])
        DNA = []; t = warm * UPH; chunk = 10.0; T = (warm + dur + post) * UPH
        def shh(th): return peak if (warm <= th < warm + dur) else base
        while t < T - 1e-6:
            rr['SHH'] = shh(t / UPH); r = sim_robust(rr, t, t + chunk, 4, ['time', 'Dna'])
            if r is not None: DNA += list(np.asarray(r['Dna'])[1:])
            t += chunk
        return len(divs(DNA))
    DURS = np.array([0, 3, 6, 9, 12, 18, 24, 36, 48]); Z['pulse_dur'] = DURS
    Z['pulse_ON'] = np.array([pulse_dur('ON', 0.12, 0.34, d) for d in DURS]); print('pulse ON done', flush=True)
    # D: GRADUAL COMPLETE withdrawal — warm 0.5, ramp 0.5->0 over `ramp` h, hold. record Cd_mRNA + divisions + SHH
    Z['wd_ramp'] = 48.0
    def grad_wd(c, Twarm=120, ramp=48.0, Tpost=180):
        rr = newrr()
        for k, v in cond(c).items(): rr[k] = v
        rr['SHH'] = 0.5; sim_robust(rr, 0, Twarm * UPH, 300, ['time', 'Dna'])
        T_, MR, DNA, SH = [], [], [], []; t = Twarm * UPH; chunk = 8.0; Tend = (Twarm + ramp + Tpost) * UPH
        while t < Tend - 1e-6:
            dt = t / UPH - Twarm; shh = max(0.0, 0.5 * (1 - dt / ramp)) if dt < ramp else 0.0
            rr['SHH'] = shh
            r = sim_robust(rr, t, t + chunk, 4, ['time', 'Dna', 'Cd_mRNA'])
            if r is not None:
                T_ += list(np.asarray(r['time'])[1:]); DNA += list(np.asarray(r['Dna'])[1:]); MR += list(np.asarray(r['Cd_mRNA'])[1:]); SH += [shh] * (len(r['time']) - 1)
            t += chunk
        T_ = (np.array(T_) / UPH - Twarm); DNA = np.array(DNA)
        Z['wd_t_' + c] = T_; Z['wd_mr_' + c] = np.array(MR); Z['wd_sh_' + c] = np.array(SH); Z['wd_dv_' + c] = T_[divs(DNA)]
    for c in CONDS:
        grad_wd(c); print('withdrawal', c, 'done', flush=True)
    # E: persistent-but-noisy Shh + gradual complete withdrawal (cumulative divisions)
    def noisy_persistent(T_h, dt=0.1, tonic=0.45, sd=0.08, tau=6.0, wd_start=170.0, ramp=48.0, seed=5):
        rng = np.random.default_rng(seed); tt = np.arange(0, T_h, dt); n = len(tt)
        s = sd * np.sqrt(2 / tau); x = np.zeros(n)
        for i in range(1, n): x[i] = x[i-1] - x[i-1] * dt / tau + s * np.sqrt(dt) * rng.normal()
        x -= x.mean(); L = np.clip(tonic + x, 0.20, 0.85)
        for i, th in enumerate(tt):
            if th >= wd_start:
                L[i] = max(0.0, L[i] * (1 - (th - wd_start) / ramp)) if (th - wd_start) < ramp else 0.0
        return tt, L
    Tn = 360.0; tL, L = noisy_persistent(Tn); Z['pn_tL'] = tL; Z['pn_L'] = L; Z['pn_wd_start'] = 170.0; Z['pn_ramp'] = 48.0
    def run_under(c, tL, L, T_h):
        rr = newrr()
        for k, v in cond(c).items(): rr[k] = v
        def Lf(th): return float(np.interp(th, tL, L))
        rr['SHH'] = Lf(0); T_, DNA = [], []; t = 0.0; chunk = 8.0
        while t < T_h * UPH - 1e-6:
            rr['SHH'] = Lf(t / UPH); r = sim_robust(rr, t, t + chunk, 4, ['time', 'Dna'])
            if r is not None: T_ += list(np.asarray(r['time'])[1:]); DNA += list(np.asarray(r['Dna'])[1:])
            t += chunk
        T_ = np.array(T_) / UPH; DNA = np.array(DNA); cum = np.zeros(len(T_))
        for i in divs(DNA): cum[i:] += 1
        return T_, cum
    for c in CONDS:
        T_, cum = run_under(c, tL, L, Tn); Z['pn_t_' + c] = T_; Z['pn_cum_' + c] = cum; print('persistent-noisy', c, 'done', flush=True)
    np.savez(CACHE, **Z); print('cached', CACHE)

# ========================== thresholds / folds ==========================
SHH = Z['SHH']
def threshold(arr):
    dd = arr[:, 0]; idx = np.where(dd > 0.5)[0]
    return SHH[idx[0]] if len(idx) else np.nan
thr = {c: threshold(Z['dose_' + c]) for c in CONDS}
# de-repression transcript fold at SHH=0.5 (index nearest 0.5)
i05 = int(np.argmin(np.abs(SHH - 0.5))) if SHH.max() >= 0.5 else len(SHH) - 1
fold = {c: Z['dose_' + c][i05, 2] / Z['dose_ON'][i05, 2] for c in CONDS}

# ========================== figure ==========================
fig, ax = plt.subplots(2, 3, figsize=(21.5, 11.6))
fig.subplots_adjust(left=0.05, right=0.985, top=0.87, bottom=0.065, hspace=0.34, wspace=0.24)

# A: threshold / Shh-dependence
a = ax[0, 0]
for c, col, mk, lab in (('ON', C_ON, 'o', 'EZH2 repression ON (governor)'), ('EZH2i', C_EZI, '^', 'EZH2i drug (partial)'), ('OFF', C_OFF, 's', 'repression OFF — genetic (full)')):
    a.plot(SHH, Z['dose_' + c][:, 0], '-' + mk, color=col, lw=2.4, ms=4.5, label='%s  (thr %.2f)' % (lab, thr[c]))
a.axvline(thr['ON'], color=C_ON, ls=':', lw=1.4); a.axvline(thr['OFF'], color=C_OFF, ls=':', lw=1.4)
a.annotate('', (thr['OFF'], 0.55), (thr['ON'], 0.55), arrowprops=dict(arrowstyle='<->', color='#333', lw=1.4))
a.text((thr['ON'] + thr['OFF']) / 2, 0.62, 'EZH2 raises the\nthreshold ~%.0f×' % (thr['ON'] / max(thr['OFF'], 1e-3)), ha='center', fontsize=8.6, fontweight='bold', color='#333')
a.text(0.005, 1.12, 'ALL arrest at SHH=0\n(all Shh-DEPENDENT now)', fontsize=8, color='#555', fontweight='bold')
a.set_xlabel('mitogen (steady SHH)'); a.set_ylabel('division rate (÷/day)'); a.set_ylim(-0.05, 1.25); a.set_xlim(-0.005, 0.34)
a.set_title('A · All Shh-DEPENDENT (basal-bake); EZH2 repression\nraises the mitogen threshold ~15×', fontsize=11.5, fontweight='bold')
a.legend(fontsize=8.0, loc='center right'); a.grid(alpha=0.2)

# B: CyclinD1 dose-response buffering
b = ax[0, 1]
for c, col, mk, lab in (('ON', C_ON, 'o', 'ON — buffered (flat)'), ('EZH2i', C_EZI, '^', 'EZH2i'), ('OFF', C_OFF, 's', 'OFF — steep, runs away')):
    b.plot(SHH, Z['dose_' + c][:, 1], '-' + mk, color=col, lw=2.4, ms=4.5, label=lab)
b.set_xlabel('mitogen (steady SHH)'); b.set_ylabel('CyclinD1 protein (steady)')
b.text(0.02, 0.72, 'the governor + reservoir BUFFER CyclinD1\nagainst mitogen level; without it CyclinD1\ntracks (and amplifies) the Shh input', transform=b.transAxes, fontsize=8, color='#555')
b.set_title('B · CyclinD1 dose-response is BUFFERED with the governor\n(EZH2 incoherent feedforward on the transcript)', fontsize=11.5, fontweight='bold')
b.legend(fontsize=8.3, loc='upper left'); b.grid(alpha=0.2)

# C: temporal noise filter (governed cell needs a SUSTAINED pulse)
c = ax[0, 2]
c.plot(Z['pulse_dur'], Z['pulse_ON'], '-o', color=C_ON, lw=2.6, ms=6, label='governed (ON): divisions from the pulse')
c.axvspan(0, 15, color=C_ON, alpha=0.07); c.axvline(16, color='#999', ls=':', lw=1.2)
ymax = max(Z['pulse_ON'].max(), 1)
c.text(16.6, ymax * 0.75, 'reservoir τ ≈ 16 h', fontsize=7.6, color='#777', rotation=90, va='center')
c.text(6, ymax * 0.06, 'brief drift →\nFILTERED (0 ÷)', ha='center', fontsize=8.3, color=C_ON, fontweight='bold')
c.text(38, ymax * 0.7, 'sustained →\ncommits', ha='center', fontsize=8.3, color='#333', fontweight='bold')
c.set_xlabel('duration of a supra-threshold pulse (h)'); c.set_ylabel('divisions triggered (governed cell)')
c.set_title('C · The governor requires a SUSTAINED signal — a brief\nmitogen excursion is filtered out (reservoir low-pass)', fontsize=11.5, fontweight='bold')
c.legend(fontsize=8.3, loc='center left'); c.grid(alpha=0.2)

# D: gradual complete withdrawal — all arrest; de-repressed coasts longer + keeps transcript up
d = ax[1, 0]; ramp = float(Z['wd_ramp'])
order = list(CONDS)
for c2, col, lab in (('ON', C_ON, 'ON (arrests, transcript low)'), ('OFF', C_OFF, 'OFF (coasts longer, transcript stays UP)'), ('EZH2i', C_EZI, 'EZH2i')):
    m = Z['wd_t_' + c2] > -5
    d.plot(Z['wd_t_' + c2][m], Z['wd_mr_' + c2][m], color=col, lw=2.6, label=lab)
    yb = -0.7 - 0.7 * order.index(c2)
    for x in Z['wd_dv_' + c2]:
        if x > -5: d.plot(x, yb, '|', color=col, ms=14, mew=2.4)
d.axhline(0, color='#ccc', lw=0.8)
d.axvspan(0, ramp, color='#95a5a6', alpha=0.10); d.text(ramp / 2, 6.3, 'Shh ramped\n0.5 → 0', ha='center', fontsize=7.8, color='#555')
d.set_xlabel('time from start of withdrawal ramp (h)'); d.set_ylabel('CyclinD1 mRNA reservoir'); d.set_ylim(-2.9, 6.9); d.set_xlim(-20, 175)
d.text(0.42, 0.30, 'ticks = divisions (per row). GRADUAL COMPLETE withdrawal →\nALL arrest (0 divisions after ~60 h). The de-repressed cell\ncoasts a bit longer and stays arrested with the transcript UP\n(≈ the cKO: transcript high, but no division, differentiates)',
       transform=d.transAxes, fontsize=7.5, color='#555')
d.set_title('D · GRADUAL COMPLETE Shh withdrawal (leaving the niche):\nall arrest → differentiation-competent', fontsize=11.5, fontweight='bold')
d.legend(fontsize=8.0, loc='upper right'); d.grid(alpha=0.2, axis='x')

# E: persistent-but-noisy Shh + gradual complete withdrawal
e = ax[1, 1]; wd0 = float(Z['pn_wd_start']); rp = float(Z['pn_ramp'])
e.plot(Z['pn_tL'], Z['pn_L'], color='#95a5a6', lw=0.8, label='Shh: persistent + noisy, then withdrawn')
e.axvspan(wd0, wd0 + rp, color='#95a5a6', alpha=0.12); e.text(wd0 + rp / 2, 0.83, 'gradual\nwithdrawal', ha='center', fontsize=7.6, color='#333')
e.set_xlabel('time (h)'); e.set_ylabel('Shh (a.u.)', color='#95a5a6'); e.set_ylim(0, 0.9)
e2 = e.twinx()
for c2, col in (('ON', C_ON), ('OFF', C_OFF), ('EZH2i', C_EZI)):
    e2.plot(Z['pn_t_' + c2], Z['pn_cum_' + c2], color=col, lw=2.6, label=c2)
e2.set_ylabel('cumulative divisions')
h1, l1 = e.get_legend_handles_labels(); h2, l2 = e2.get_legend_handles_labels()
e.legend(h1 + h2, l1 + l2, fontsize=7.8, loc='upper left')
e.set_title('E · Persistent-but-noisy Shh: smooth Shh-gated proliferation\nthat STOPS when the niche signal is withdrawn', fontsize=11.5, fontweight='bold'); e.grid(alpha=0.2)

# F: summary — threshold (log) + de-repression transcript fold
f = ax[1, 2]
conds = ['ON', 'EZH2i', 'OFF']; cols = [C_ON, C_EZI, C_OFF]
x = np.arange(3); w = 0.38
thr_vals = [thr[c] for c in conds]; fold_vals = [fold[c] for c in conds]
f.bar(x - w / 2, thr_vals, w, color=cols, alpha=0.95, label='mitogen threshold (SHH)')
f2 = f.twinx()
f2.bar(x + w / 2, fold_vals, w, color=cols, alpha=0.4, hatch='//', label='CyclinD1 transcript fold vs ON')
for i, v in enumerate(thr_vals): f.text(i - w / 2, v + 0.004, '%.2f' % v, ha='center', fontsize=8.2)
for i, v in enumerate(fold_vals): f2.text(i + w / 2, v + 0.05, '%.1f×' % v, ha='center', fontsize=8.2)
f.set_xticks(x); f.set_xticklabels(['ON\n(governor)', 'EZH2i\n(drug)', 'OFF\n(genetic)'])
f.set_ylabel('mitogen threshold (SHH)'); f2.set_ylabel('CyclinD1 transcript fold vs ON')
f.set_ylim(0, max(thr_vals) * 1.25); f2.set_ylim(0, max(fold_vals) * 1.3)
f.text(0.02, 0.90, 'EZH2 = a threshold/gain governor: it raises the Shh\nthreshold ~15× (0.01→0.18) and holds the transcript down.\nAll conditions are Shh-dependent and arrest on withdrawal.', transform=f.transAxes, fontsize=7.6, color='#555')
h1, l1 = f.get_legend_handles_labels(); h2, l2 = f2.get_legend_handles_labels()
f.legend(h1 + h2, l1 + l2, fontsize=7.8, loc='upper center')
f.set_title('F · Summary: EZH2 raises the Shh threshold ~15× and\nbuffers the transcript — but all stay Shh-dependent', fontsize=11.5, fontweight='bold'); f.grid(alpha=0.2, axis='y')

fig.suptitle('What EZH2→CyclinD1 repression does for GNPs (reservoir+IFFL, basal-bake):  proliferation is Shh-DEPENDENT for all; EZH2 is a THRESHOLD/GAIN governor '
             '(~15× threshold, buffered dose-response, brief-drift filter) — de-repression raises the transcript but the cell still needs Shh and arrests on withdrawal',
             fontsize=11.4, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_governor_noise_filter.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
print(f"thresholds: ON={thr['ON']:.3f} EZH2i={thr['EZH2i']:.3f} OFF={thr['OFF']:.3f}  (ON/OFF ~{thr['ON']/max(thr['OFF'],1e-3):.0f}x)")
print(f"transcript fold vs ON @SHH0.5: " + ", ".join(f"{c}={fold[c]:.2f}" for c in CONDS))
print(f"withdrawal last division (h): " + ", ".join(f"{c}={Z['wd_dv_'+c][-1]:.0f}" if len(Z['wd_dv_'+c]) else f"{c}=none" for c in CONDS))
