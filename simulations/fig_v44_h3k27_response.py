"""Intuition plots for the H3K27me3 -> CyclinD1 response characteristics of v44.

TWO FIGURES:
  fig_v44_h3k27_response.png  -- STATIC transfer functions + parameter sensitivity (how the
    variables shape transcription): the leaky-Hill readout R(Mk), the methylation dose-response
    Mk_ss(EZH2) and how cell-cycle DILUTION shifts it, the CyclinD1 "throttle" R(Mk_ss(EZH2)),
    and per-parameter small-multiples (K_mk, n_mk, f0_mk, g_mk, del_mk, k_w_mk).
  fig_v44_h3k27_timing.png    -- TIMING & cell-cycle impact (full model runs): de-repression
    kinetics after EZH2i (t1/2 set by del_mk + dilution), the replicative-dilution sawtooth,
    mean mark vs cell-cycle period (dilution-vs-restoration race), and entry-timing brake on/off.

Run:  ./venv/bin/python simulations/fig_v44_h3k27_response.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

# ---- default H3K27me3 params (from build_model_v44_heldt.py) ----
P = dict(k_w_mk=0.00160, k0_mk=0.000258, del_mk=0.00070, g_mk=0.30,
         K_tx_mk=5.0, p_tx_mk=2.0, K_mk=0.305, n_mk=4.15, f0_mk=0.233)
TX_REP = 2.6**2 / (5.0**2 + 2.6**2)   # representative transcription Hill (Cd_mRNA~2.6) for the eviction arm
LN2 = np.log(2.0)

C_MARK = '#7a5aa6'; C_TX = '#b0402f'; C_DRIVE = '#8a94a1'
C_ARR = '#c0392b'; C_SLOW = '#e67e22'; C_FAST = '#117a65'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
                     'axes.titleweight': 'bold', 'figure.dpi': 110})


def R_of_Mk(Mk, K=P['K_mk'], n=P['n_mk'], f0=P['f0_mk']):
    """leaky repression readout: fraction of transcriptional drive that survives the mark."""
    return f0 + (1 - f0) / (1 + (Mk / K) ** n)


def mk_ss(EZH2, T=1 - P['g_mk'] * TX_REP, ddil=0.0, kw=P['k_w_mk'], k0=P['k0_mk'], dl=P['del_mk']):
    """steady-state H3K27me3 occupancy for a given EZH2, transcription-antagonism T, and extra
    linear loss ddil (= ln2/period, the mean-field replicative dilution). Root in [0,1] of
    EZH2*T*(kw*Mk+k0)*(1-Mk) = (dl+ddil)*Mk."""
    EZH2 = np.atleast_1d(EZH2).astype(float)
    a = EZH2 * T
    L = dl + ddil
    out = np.empty_like(a)
    for i, ai in enumerate(a):
        if kw <= 0:                     # de-novo only -> linear
            out[i] = ai * k0 / (L + ai * k0)
        else:
            A = ai * kw; B = L - ai * (kw - k0); C = -ai * k0
            disc = B * B - 4 * A * C
            out[i] = (-B + np.sqrt(max(disc, 0.0))) / (2 * A)
    return np.clip(out, 0, 1)


# =====================================================================================
# FIGURE 1 -- static response characteristics + sensitivity
# =====================================================================================
fig = plt.figure(figsize=(15, 9))
gs = GridSpec(3, 6, figure=fig, height_ratios=[1.25, 1.0, 1.0], hspace=0.55, wspace=0.75)
EZ = np.linspace(0, 3, 400)

# --- A: readout R(Mk) ---
axA = fig.add_subplot(gs[0, 0:2])
Mk = np.linspace(0, 1, 400)
for n, col in [(1, '#c9b8e0'), (2, '#a487c9'), (P['n_mk'], C_MARK), (8, '#4a2f6b')]:
    axA.plot(Mk, R_of_Mk(Mk, n=n), color=col, lw=2.3 if n == P['n_mk'] else 1.5,
             label=f"n={n:g}" + (" (default)" if n == P['n_mk'] else ""))
axA.axhline(P['f0_mk'], ls=':', color='#555', lw=1); axA.text(0.02, P['f0_mk'] + 0.02, 'floor f0 = 0.233', fontsize=8, color='#555')
axA.axvline(P['K_mk'], ls=':', color='#555', lw=1); axA.text(P['K_mk'] + 0.02, 0.9, 'K_mk', fontsize=8, color='#555')
axA.set_xlabel('H3K27me3 occupancy  Mk'); axA.set_ylabel('R = surviving transcription\n(fraction of drive)')
axA.set_title('(A) Readout: the mark throttles,\nnever silences (leaky Hill)'); axA.set_ylim(0, 1.02)
axA.legend(fontsize=8, loc='upper right'); axA.grid(alpha=0.15)

# --- B: methylation dose-response Mk_ss(EZH2), dilution shifts it ---
axB = fig.add_subplot(gs[0, 2:4])
axB.plot(EZ, mk_ss(EZ, ddil=0.0), color=C_ARR, lw=2.4, label='arrested (no dilution)')
axB.plot(EZ, mk_ss(EZ, ddil=LN2 / (34 * 60)), color=C_SLOW, lw=2.0, label='cycling, 34 h period')
axB.plot(EZ, mk_ss(EZ, ddil=LN2 / (18 * 60)), color=C_FAST, lw=2.0, label='cycling, 18 h period')
axB.plot(EZ, mk_ss(EZ, kw=0.0), color='#888', lw=1.6, ls='--', label='no read-write (k_w=0)')
axB.set_xlabel('EZH2 level'); axB.set_ylabel('steady-state mark  Mk_ss')
axB.set_title('(B) Writing vs erasing:\nfaster cycling dilutes the mark'); axB.set_ylim(0, 1.02)
axB.legend(fontsize=7.5, loc='lower right'); axB.grid(alpha=0.15)

# --- C: the throttle -- transcript vs EZH2 (buffering) ---
axC = fig.add_subplot(gs[0, 4:6])
for ddil, col, lab in [(0.0, C_ARR, 'arrested'), (LN2 / (34 * 60), C_SLOW, '34 h'), (LN2 / (18 * 60), C_FAST, '18 h')]:
    axC.plot(EZ, R_of_Mk(mk_ss(EZ, ddil=ddil)), color=col, lw=2.2, label=f'brake, {lab}')
axC.axhline(1.0, color=C_DRIVE, lw=2.0, ls='--', label='no brake (drive)')
axC.fill_between(EZ, R_of_Mk(mk_ss(EZ)), 1.0, color=C_MARK, alpha=0.10)
axC.annotate('repression\n= buffering', xy=(2.2, 0.6), fontsize=8.5, color=C_MARK, ha='center')
axC.set_xlabel('EZH2 level  (rises with mitogen)'); axC.set_ylabel('relative CyclinD1 transcript')
axC.set_title('(C) Throttle: more EZH2 →\nmore mark → less transcript'); axC.set_ylim(0, 1.05)
axC.legend(fontsize=7.5, loc='lower left'); axC.grid(alpha=0.15)

# --- D-I: per-parameter sensitivity small-multiples (impact on transcription) ---
SENS = [('K_mk', 'K_mk', 'mark→repression midpoint'),
        ('n_mk', 'n_mk', 'switch steepness'),
        ('f0_mk', 'f0_mk', 'leak floor'),
        ('g_mk', 'g_mk', 'transcription eviction'),
        ('del_mk', 'del_mk', 'demethylation / turnover'),
        ('k_w_mk', 'k_w_mk', 'read-write autocatalysis')]
for k, (pname, sym, desc) in enumerate(SENS):
    ax = fig.add_subplot(gs[1 + k // 3, 2 * (k % 3):2 * (k % 3) + 2])
    base = P[pname]
    for scale, col, lw in [(0.5, '#4a90c2', 1.6), (1.0, '#333', 2.3), (2.0, '#c25a4a', 1.6)]:
        kw = dict(kw=P['k_w_mk'], k0=P['k0_mk'], dl=P['del_mk'], T=1 - P['g_mk'] * TX_REP)
        n, K, f0 = P['n_mk'], P['K_mk'], P['f0_mk']
        val = base * scale
        if pname == 'del_mk': kw['dl'] = val
        elif pname == 'k_w_mk': kw['kw'] = val
        elif pname == 'g_mk': kw['T'] = 1 - val * TX_REP
        elif pname == 'n_mk': n = val
        elif pname == 'K_mk': K = val
        elif pname == 'f0_mk': f0 = val
        y = R_of_Mk(mk_ss(EZ, **kw), K=K, n=n, f0=f0)
        ax.plot(EZ, y, color=col, lw=lw, label=f'{scale:g}×')
    ax.set_title(f'{sym}  —  {desc}', fontsize=9.0, fontweight='bold')
    ax.set_ylim(0, 1.05); ax.grid(alpha=0.15)
    if k % 3 == 0: ax.set_ylabel('rel. transcript')
    if k // 3 == 1: ax.set_xlabel('EZH2 level')
    ax.legend(fontsize=7, loc='lower left', ncol=3, columnspacing=0.8, handlelength=1.0)

fig.suptitle('H3K27me3 → CyclinD1 response characteristics — static transfer functions & parameter sensitivity',
             fontsize=13, fontweight='bold', y=0.98)
fig.text(0.5, 0.005, 'Panels D–I: relative CyclinD1 transcript R(Mk_ss(EZH2)) as each H3K27me3 parameter is scaled 0.5× / 1× / 2× '
         '(transcription-eviction held at a representative Cd_mRNA).', ha='center', fontsize=8, color='#666')
plt.savefig('simulations/fig_v44_h3k27_response.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_h3k27_response.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_h3k27_response.png')


# =====================================================================================
# FIGURE 2 -- timing & cell-cycle impact (full model)
# =====================================================================================
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)


def rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
    r.integrator.setValue('absolute_tolerance', 1e-8)
    try: r.integrator.setValue('maximum_num_steps', 60000)
    except Exception: pass
    return r


def setgnp(r, shh=1.0):
    r.reset()
    for k, v in GNP.items(): r[k] = v
    r['SHH'] = shh


def movavg(t, x, win=1320):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


fig2 = plt.figure(figsize=(15, 9))
gs2 = GridSpec(2, 2, figure=fig2, hspace=0.34, wspace=0.26)

# --- A: de-repression kinetics after EZH2i (timing set by del_mk + dilution) ---
axA = fig2.add_subplot(gs2[0, 0])
t0 = 4000.0; span = 9000.0
for scale, col in [(0.5, C_SLOW), (1.0, C_MARK), (2.0, C_FAST)]:
    r = rr(); setgnp(r, shh=1.0); r['del_mk'] = P['del_mk'] * scale
    a = r.simulate(0, t0, 4000, selections=['time', 'Mk'])
    mk0 = a['Mk'][-1]
    r['EZH2i'] = 1
    b = r.simulate(t0, t0 + span, 9000, selections=['time', 'Mk'])
    tt = (b['time'] - t0) / 60.0
    frac = b['Mk'] / mk0
    axA.plot(tt, frac, color=col, lw=2.2, label=f"del_mk × {scale:g}")
    # t1/2
    below = np.where(frac <= 0.5)[0]
    if len(below):
        thalf = tt[below[0]]
        axA.plot(thalf, 0.5, 'o', color=col, ms=6)
        if scale == 1.0:
            axA.annotate(f't½ ≈ {thalf:.0f} h', xy=(thalf, 0.5), xytext=(thalf + 6, 0.62),
                         fontsize=9, color=col, arrowprops=dict(arrowstyle='->', color=col))
axA.axhline(0.5, ls=':', color='#888', lw=1)
axA.set_xlabel('time after EZH2i (h)'); axA.set_ylabel('mark remaining  Mk / Mk₀')
axA.set_title('(A) De-repression timing: how fast the mark clears\nwhen the writer is inhibited')
axA.legend(fontsize=8.5, title='turnover rate', title_fontsize=8); axA.grid(alpha=0.15); axA.set_ylim(0, 1.02)

# --- B: replicative-dilution sawtooth (mark vs cycle) ---
axB = fig2.add_subplot(gs2[0, 1])
r = rr(); setgnp(r, shh=1.0)
c = r.simulate(0, 9000, 18000, selections=['time', 'Mk', 'Cd', 'MPF', 'Dna'])
tt = c['time'] / 60.0
axB.plot(tt, c['Mk'], color=C_MARK, lw=1.6, label='Mk (H3K27me3)')
axB.set_xlabel('time (h)'); axB.set_ylabel('Mk', color=C_MARK); axB.tick_params(axis='y', colors=C_MARK)
axB.set_title('(B) Replicative dilution: the mark is halved every S phase,\nthen rebuilds — a sawtooth clocked by the cell cycle')
# mark division times (MPF peaks)
pk, _ = find_peaks(c['MPF'], prominence=0.15, distance=200)
for p in pk:
    axB.axvline(tt[p], color='#bbb', lw=0.8, ls='--')
axBt = axB.twinx()
axBt.plot(tt, c['Cd'], color=C_TX, lw=1.3, alpha=0.7, label='CyclinD1')
axBt.set_ylabel('CyclinD1 (Cd)', color=C_TX); axBt.tick_params(axis='y', colors=C_TX)
axB.grid(alpha=0.12)
axB.plot([], [], color='#bbb', ls='--', label='division'); axB.legend(fontsize=8, loc='upper right')

# --- C: mean mark & transcript vs cell-cycle period (dilution-restoration race) ---
axC = fig2.add_subplot(gs2[1, 0])
MUS = np.linspace(0.00030, 0.00095, 12)
periods, meanMk, meanTx = [], [], []
for mu in MUS:
    r = rr(); setgnp(r, shh=1.0); r['mu'] = float(mu)
    try:
        d = r.simulate(0, 14000, 28000, selections=['time', 'Mk', 'Cd_mRNA', 'MPF'])
    except Exception:
        periods.append(np.nan); meanMk.append(np.nan); meanTx.append(np.nan); continue
    t = d['time']; m = t >= 5000
    pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
    tp = t[m][pk]
    per = np.median(np.diff(tp)) / 60.0 if len(tp) >= 2 else np.nan
    periods.append(per); meanMk.append(np.nanmean(d['Mk'][m])); meanTx.append(np.nanmean(d['Cd_mRNA'][m]))
periods = np.array(periods); meanMk = np.array(meanMk); meanTx = np.array(meanTx)
ok = ~np.isnan(periods)
axC.plot(periods[ok], meanMk[ok], '-o', color=C_MARK, lw=2.2, ms=5, label='mean Mk')
axC.set_xlabel('cell-cycle period (h)'); axC.set_ylabel('mean H3K27me3  Mk', color=C_MARK)
axC.tick_params(axis='y', colors=C_MARK)
axC.set_title('(C) Cell-cycle speed sets the mark:\nslower cycling → less dilution → stronger brake')
axCt = axC.twinx()
axCt.plot(periods[ok], meanTx[ok], '-s', color=C_TX, lw=2.0, ms=5, label='mean transcript')
axCt.set_ylabel('mean CyclinD1 transcript', color=C_TX); axCt.tick_params(axis='y', colors=C_TX)
axC.grid(alpha=0.15)

# --- D: entry-timing, brake on vs off ---
axD = fig2.add_subplot(gs2[1, 1])
t0r, D = 500.0, 9000.0     # SHH ramp 0.05 -> 1.0 over 150 h
def shh_ramp(tt): return 0.05 if tt < t0r else (1.0 if tt > t0r + D else 0.05 + (1.0 - 0.05) * (tt - t0r) / D)
for ezi, col, lab in [(0, C_MARK, 'brake ON (default)'), (1, C_FAST, 'brake OFF (EZH2i)')]:
    r = rr(); setgnp(r, shh=0.05); r['EZH2i'] = ezi
    T, MP, SH = [], [], []; t = 0.0; chunk = 120.0
    while t < t0r + D + 1500:
        r['SHH'] = float(shh_ramp(t))
        try:
            s = r.simulate(t, t + chunk, 60, selections=['time', 'MPF', 'SHH'])
            T.append(s['time'][1:]); MP.append(s['MPF'][1:]); SH.append(s['SHH'][1:])
        except Exception:
            break
        t += chunk
    T = np.concatenate(T); MP = np.concatenate(MP); SH = np.concatenate(SH)
    axD.plot((T - t0r) / 60.0, SH, color=col, lw=1.4, alpha=0.5)
    pk, _ = find_peaks(MP, prominence=0.15, distance=200)
    if len(pk):
        te_h = (T[pk[0]] - t0r) / 60.0; shh_e = SH[pk[0]]
        axD.axvline(te_h, color=col, lw=2.0)
        axD.plot(te_h, shh_e, 'v', color=col, ms=11)
        axD.annotate(f'{lab}\n1st division: {te_h:.0f} h, SHH {shh_e:.2f}',
                     xy=(te_h, shh_e), xytext=(te_h + 3, shh_e + 0.12 - 0.24 * ezi),
                     fontsize=8.5, color=col, arrowprops=dict(arrowstyle='->', color=col))
axD.set_xlabel('time along rising-Hh ramp (h)'); axD.set_ylabel('SHH input')
axD.set_title('(D) Entry timing: the brake delays\ncell-cycle onset to higher Hh / later time')
axD.grid(alpha=0.15); axD.set_xlim(-5, 165)

fig2.suptitle('H3K27me3 → CyclinD1 timing & cell-cycle impact (full v44 model)', fontsize=13, fontweight='bold', y=0.98)
plt.savefig('simulations/fig_v44_h3k27_timing.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_h3k27_timing.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_h3k27_timing.png')
print('periods(h):', np.round(periods, 1))
print('meanMk:', np.round(meanMk, 3))
