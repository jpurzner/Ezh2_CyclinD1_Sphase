"""Biologically realistic Shh mitogen input, DECOMPOSITION of the Gli signal, and the demonstration that the
CyclinD1 mRNA RESERVOIR low-pass-filters it (JP 2026-07-21; grounded in docs/hh-signalling-timing.md).
Reads sim_realistic_mitogen_cache.pkl. Run: ./venv/bin/python simulations/fig_realistic_mitogen.py"""
import os, pickle, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_realistic_mitogen_cache.pkl'), 'rb'))
C_LIG, C_SLOW = '#95a5a6', '#2c3e50'
C1, C2, C3, C4 = '#bdc3c7', '#2e86c1', '#16a085', '#c0392b'
C_GLI, C_MRNA, C_PROT = '#8e44ad', '#0e6251', '#e67e22'
def divmarks(t, dna): d = np.asarray(dna); return t[np.where((d[:-1] > 0.9) & (d[1:] < 0.1))[0]]

fig, ax = plt.subplots(2, 3, figsize=(21.5, 11.4))
fig.subplots_adjust(left=0.045, right=0.99, top=0.87, bottom=0.065, hspace=0.34, wspace=0.22)

# ---------------------------------------------------------------------------------------------------------------
# A: the evidence-driven ligand input
a = ax[0, 0]; s = R['lig_stats']
a.plot(R['t_lig'], R['L_full'], color=C_LIG, lw=0.8, label='Shh ligand L(t): tonic + slow drift + fast noise')
a.plot(R['t_lig'], R['L_slow'], color=C_SLOW, lw=2.4, label='slow component (paracrine drift — the "real" signal)')
a.set_xlabel('time (h)'); a.set_ylabel('Shh ligand (a.u.)'); a.set_ylim(0, max(R['L_full']) * 1.12)
a.text(0.02, 0.06, 'CV = %.2f  =  slow %.2f (τ≈40 h, passband)\n              +  fast %.2f (τ≈0.5 h, stopband)\nno fast pulse train (that is ERK, not Hh)'
       % (s['cv'], s['cv_slow'], s['cv_fast']), transform=a.transAxes, fontsize=8, color='#555')
a.set_title('A · Input: paracrine Shh — tonic + slow drift + fast noise\n(evidence: slowly-varying, low-pass-able; NOT a pulse train)', fontsize=11.5, fontweight='bold')
a.legend(fontsize=8.2, loc='upper right'); a.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# B: Gli-signal decomposition — Gli is quasi-static, so it TRACKS the ligand (carries the noise)
b = ax[0, 1]
t1, o1 = R['d1_tonic']; t2, o2 = R['d2_noise']; t3, o3 = R['d3_cilium']; t4, o4 = R['d4_full']
b.plot(t1, o1['Gli_act'], color=C1, lw=1.6, ls='--', label='1. tonic ligand → flat Gli')
b.plot(t2, o2['Gli_act'], color=C2, lw=1.5, label='2. + ligand noise → Gli TRACKS it (quasi-static, noisy)')
b.plot(t3, o3['Gli_act'], color=C3, lw=1.5, alpha=0.9, label='3. + cilium duty cycle → G2/M dips')
b.plot(t4, o4['Gli_act'], color=C4, lw=2.3, label='4. + ~20 h adaptation → slow setpoint drop (= full)')
b.set_xlabel('time (h)'); b.set_ylabel('Gli activity (Gli_act)')
b.text(0.02, 0.05, 'v44 Gli is fast (protein τ ~40 min–2 h « the ~22 h cycle):\nit does NOT filter — the noise reaches Gli (filtered downstream, panel C/D)', transform=b.transAxes, fontsize=7.6, color='#555')
b.set_title('B · Gli signal DECOMPOSED — what each component contributes\n(Gli is quasi-static → it carries the ligand noise)', fontsize=11.5, fontweight='bold')
b.legend(fontsize=7.6, loc='upper right'); b.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# C: the reservoir removes the FAST ligand noise (high-frequency component; isolates noise from the cell-cycle
# oscillation that dominates Cd_mRNA's TOTAL variance). Uses d2_noise (cilium off: Gli = pure ligand readout).
c = ax[1, 0]
tcn, ocn = R['d2_noise']
def highpass(t, x, w_h=3.0):                                # subtract a ~3 h rolling mean -> keep only fast fluctuations
    dt = np.median(np.diff(t)); k = max(1, int(w_h / dt))
    base = np.convolve(np.pad(x, k, mode='edge'), np.ones(k) / k, mode='same')[k:-k]
    return (x - base) / x.mean()                            # relative deviation
gli_hf = highpass(tcn, ocn['Gli_act']); mrna_hf = highpass(tcn, ocn['Cd_mRNA'])
g_amp, m_amp = gli_hf.std(), mrna_hf.std()
c.plot(tcn, gli_hf, color=C_GLI, lw=0.9, alpha=0.9, label='Gli fast noise (amp %.1f%%) — carried through' % (100 * g_amp))
c.plot(tcn, mrna_hf, color=C_MRNA, lw=2.4, label='CyclinD1 mRNA fast noise (amp %.1f%%) — filtered out' % (100 * m_amp))
c.axhline(0, color='#ccc', lw=0.8, ls=':')
c.text(0.02, 0.08, 'high-frequency (period < 3 h) component only, so the slow\ncell-cycle oscillation is removed: the slow (~16 h) transcript\nRESERVOIR attenuates the fast ligand noise %.0f× → the\nproliferation decision is buffered against ligand noise'
       % (g_amp / max(m_amp, 1e-9)), transform=c.transAxes, fontsize=7.9, color=C_MRNA, fontweight='bold')
c.set_xlabel('time (h)'); c.set_ylabel('fast fluctuation (Δ from 3 h baseline, /mean)')
c.set_title('C · The CyclinD1 mRNA RESERVOIR removes the fast ligand noise\n(temporal averaging by slow mRNA turnover — the baked memory)', fontsize=11.5, fontweight='bold')
c.legend(fontsize=8.3, loc='upper right'); c.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# D: Bode transfer function ligand -> {Gli, Cd_mRNA} — the quantitative low-pass, with the analytic first-order curve
d = ax[1, 1]; B = R['bode']; per = np.asarray(B['period']); tau = R['reservoir_tau_h']
gli_n = B['Gli_act'] / B['Gli_act'][-1]; mrna_n = B['Cd_mRNA'] / B['Cd_mRNA'][-1]
Tcell = 22.8                                                # ~cycle period: lock-in is contaminated here
d.axvspan(0.3, 3, color=C2, alpha=0.07); d.axvspan(40, 200, color=C_MRNA, alpha=0.07)
d.axvspan(Tcell / 1.25, Tcell * 1.25, color='#999', alpha=0.12)
d.text(0.85, 1.62, 'fast noise\n(stopband)', fontsize=7.3, color=C2, ha='center')
d.text(95, 1.62, 'slow drift\n(passband)', fontsize=7.3, color=C_MRNA, ha='center')
d.text(Tcell, 0.06, 'cell cycle\n(lock-in\ncontaminated)', fontsize=6.6, color='#777', ha='center')
# analytic first-order low-pass |H| = 1/sqrt(1+(2*pi*tau/period)^2), normalised to the longest measured period
pp = np.logspace(np.log10(0.5), np.log10(200), 300)
Hn = lambda P: 1.0 / np.sqrt(1 + (2 * np.pi * tau / P) ** 2)
d.semilogx(pp, Hn(pp) / Hn(per[-1]), color=C_MRNA, lw=1.4, ls='--', alpha=0.75, label='first-order low-pass, τ=%.0f h (theory)' % tau)
d.semilogx(per, gli_n, '-o', color=C_GLI, lw=2.2, ms=5, label='ligand → Gli (flat = quasi-static, no filtering)')
d.semilogx(per, mrna_n, 's', color=C_MRNA, ms=6, label='ligand → CyclinD1 mRNA (measured)')
for ta, lab in ((40 / 60, 'Gli1 t½ 40′'), (2.0, 'Smo t½ 2 h')):
    d.axvline(ta, color='#bbb', lw=0.8, ls=':'); d.text(ta, 0.02, lab, rotation=90, fontsize=6.6, color='#999', va='bottom', ha='right')
d.set_xlabel('ligand fluctuation period (h)  — fast ← → slow'); d.set_ylabel('transfer gain (norm. to slow limit)')
d.set_ylim(0, 1.85)
d.set_title('D · Transfer function: Gli passes all frequencies; the reservoir\nis a first-order LOW-PASS (fast fluctuations strongly attenuated)', fontsize=11.5, fontweight='bold')
d.legend(fontsize=7.7, loc='center left'); d.grid(alpha=0.2, which='both')

# ---------------------------------------------------------------------------------------------------------------
# E: cilium duty cycle
e = ax[0, 2]
tcc, occ = R['d3_cilium']; m = (tcc > 40) & (tcc < 95)
e.plot(tcc[m], occ['cilium_gate'][m], color='#16a085', lw=2.3, label='cilium gate (Shh reception)')
e.plot(tcc[m], occ['Cb'][m] / (max(occ['Cb'][m]) + 1e-9), color='#7f8c8d', lw=1.2, ls=':', label='CyclinB (G2/M marker, scaled)')
for dt in divmarks(tcc, occ['Dna']):
    if 40 < dt < 95: e.axvline(dt, color='#c0392b', lw=0.9, alpha=0.5)
e.text(0.02, 0.10, 'red = divisions; gate drops through G2/M (cilium resorbed →\nShh reception OFF), restored after ÷. ON-fraction = %.2f (~94%% ciliated)'
       % R['cilium_on_fraction'], transform=e.transAxes, fontsize=8, color='#555')
e.set_xlabel('time (h)'); e.set_ylabel('cilium gate / CyclinB'); e.set_ylim(-0.03, 1.08)
e.set_title('E · Cilium DUTY CYCLE — reception phase-gated OFF at G2/M\n(~94%% ciliated progenitors, pre-mitotic resorption)', fontsize=11.5, fontweight='bold')
e.legend(fontsize=8.2, loc='center right'); e.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# F: cell-to-cell heterogeneity — noisy per-cell Gli, smoother per-cell reservoir
f = ax[1, 2]
for i, (tc2, g, mrna, Lc) in enumerate(R['cells']):
    w = tc2 > 12                                            # trim the initialization transient
    gm = g[w].mean(); mm = mrna[w].mean()
    f.plot(tc2[w], g[w] / gm, lw=0.8, alpha=0.5, color=C_GLI)
    f.plot(tc2[w], mrna[w] / mm, lw=1.6, alpha=0.8, color=C_MRNA)
f.plot([], [], color=C_GLI, lw=1.2, label='per-cell Gli (noisy)')
f.plot([], [], color=C_MRNA, lw=2.0, label='per-cell CyclinD1 mRNA (reservoir-smoothed)')
f.set_xlabel('time (h)'); f.set_ylabel('level (per cell, norm. to mean)')
f.set_title('F · Cell-to-cell heterogeneity — noisy single-cell Hh response,\nreservoir-smoothed downstream (precision emerges downstream)', fontsize=11.5, fontweight='bold')
f.legend(fontsize=8.3, loc='upper right'); f.grid(alpha=0.2)

fig.suptitle('Realistic Shh input & Gli-signal decomposition:  noisy paracrine ligand → quasi-static Gli (carries the noise) → the slow CyclinD1 mRNA RESERVOIR low-pass-filters it '
             '(cilium duty cycle + ~20 h adaptation layered on)', fontsize=12, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_realistic_mitogen.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
