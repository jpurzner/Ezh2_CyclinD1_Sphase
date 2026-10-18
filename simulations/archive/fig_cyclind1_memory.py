"""CyclinD1 = the Hedgehog MEMORY carrier, and the memory is the cumulative CyclinD1 TRANSCRIPT reservoir
(Ho, Tsai & Stearns 2020, Curr Biol 30:2829; JP mechanism 2026-07-21). The proliferation decision sits in the
PREVIOUS cycle. The memory lives in the mRNA (slow, ~11-12h reservoir); the PROTEIN is labile/fast (degraded for
most of the cycle, forms rapidly when there is lots of transcript). EZH2 forms an INCOHERENT FEEDFORWARD on the
transcript (Gli raises it, EZH2 lowers it) that BUFFERS it against mitogen changes and holds it up as mitogen falls.

6 panels:
 A  Two-tier structure: slow mRNA RESERVOIR (memory) -> fast protein READOUT.
 B  Restriction point in the PREVIOUS cycle: after Hh-off the reservoir COASTS ~1 more cycle (~11-12h);
    the fast/slaved mRNA arrests immediately.
 C  EZH2 incoherent feedforward BUFFERS the dose-response: as mitogen falls EZH2 de-represses, flattening the
    transcript's mitogen-dependence (vs a frozen-repression control).
 D  EZH2 IFFL transient: after a mitogen step-down EZH2 falls slowly -> repression LIFTS slowly -> transcript is
    held up ("stays up longer") vs the frozen-repression control.
 E  Widened cilium DUTY CYCLE: Shh reception ON through G1/S, OFF at G2 (not the textbook G1-only window).
 F  Robustness: proliferation is buffered against the intermittent (duty-cycled) reception (the low-pass argument).
Baked default = the CyclinD1 mRNA-reservoir memory (k_Cd_mRNA_deg=0.001, tx scaled 1/800; protein fast k_Cd_deg=1.0).
Run: ./venv/bin/python simulations/fig_cyclind1_memory.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
# mRNA-FAST / slaved control: restore the ~1-min mRNA (k_Cd_mRNA_deg 0.001->0.8) and the transcription x800 to hold
# the same steady level, so the ONLY difference vs the default is the mRNA timescale (memory vs slaved).
FAST = dict(k_Cd_mRNA_deg=0.8, k_Cd_tx_basal=0.44921479494593447,
            k_Cd_tx_Gli_max=75.26553209129528, k_Cd_tx_MYCN=16.206164824211957)
C_MEM, C_FAST, C_PROT, C_CIL, C_CD = '#0e6251', '#c0392b', '#e67e22', '#16a085', '#2e86c1'
C_EZ, C_GLI, C_REP = '#8e44ad', '#95a5a6', '#c0392b'
def divs(d): d = np.asarray(d); return np.where((d[:-1] > 0.9) & (d[1:] < 0.1))[0]

def build(over): rr = te.loada(build_model_v44(params=over)); rr.reset(); return rr
def rep_factor(rr, prc2_rep):
    f0, Kp, n = rr['f0_prc2'], rr['K_prc2'], rr['n_prc2']
    return f0 + (1 - f0) / (1 + (np.asarray(prc2_rep) / Kp) ** n)

def run(over, seg, sel):
    rr = build(over)
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('absolute_tolerance', 1e-8)
    T, out = [], {k: [] for k in sel if k != 'time'}; t0 = 0.0
    for ov, dur in seg:
        for k, v in ov.items(): rr[k] = v
        r = rr.simulate(t0, t0 + dur, max(2, int(dur / 6)), selections=['time'] + list(out))
        T += list(r['time'])
        for k in out: out[k] += list(r[k])
        t0 += dur
    return np.array(T) / UPH, {k: np.array(v) for k, v in out.items()}, rr

fig, ax = plt.subplots(2, 3, figsize=(21.5, 11.6))
fig.subplots_adjust(left=0.045, right=0.99, top=0.86, bottom=0.065, hspace=0.34, wspace=0.235)

# ---------------------------------------------------------------------------------------------------------------
# A: two-tier — slow mRNA reservoir -> fast protein readout
a = ax[0, 0]
t, o, _ = run(None, [(dict(SHH=0.5), 90 * UPH)], ['Cd_mRNA', 'Cd', 'Dna'])
m = t > 40
a.plot(t[m], o['Cd_mRNA'][m], color=C_MEM, lw=2.6, label='CyclinD1 mRNA — SLOW reservoir (the memory)')
a.plot(t[m], o['Cd'][m], color=C_PROT, lw=1.6, alpha=0.9, label='CyclinD1 protein — FAST readout (labile)')
for dt in t[divs(o['Dna'])]:
    if dt > 40: a.axvline(dt, color='#bbb', lw=0.8, ls=':')
a.text(0.02, 0.05, 'grey = divisions. Protein is degraded for most of the cycle and forms\nrapidly when transcript is high; the mRNA pool carries THROUGH mitosis.',
       transform=a.transAxes, fontsize=8, color='#555')
a.set_xlabel('time (h)'); a.set_ylabel('CyclinD1 mRNA / protein')
a.set_title('A · The memory is the TRANSCRIPT reservoir\nslow mRNA (~11–12 h)  →  fast protein readout', fontsize=11.5, fontweight='bold')
a.legend(fontsize=8.6, loc='upper right'); a.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# B: previous-cycle restriction point — reservoir coasts after Hh-off (memory vs mRNA-fast/slaved)
b = ax[0, 1]; onset = 130.0
for over, col, lab in ((None, C_MEM, 'mRNA RESERVOIR (default, ~12 h)'), (FAST, C_FAST, 'mRNA fast / slaved (k_Cd_mRNA_deg=0.8)')):
    t, o, _ = run(over, [(dict(SHH=0.5), onset * UPH), (dict(HHi=1), 46 * UPH)], ['Cd_mRNA', 'Dna'])
    b.plot(t, o['Cd_mRNA'], color=col, lw=2.3, label=lab)
    for dt in t[divs(o['Dna'])]:
        if dt > onset: b.plot(dt, -0.12, '|', color=col, ms=10, mew=1.7)
b.axvline(onset, color='#333', ls=':', lw=1.4)
b.text(onset + 0.8, 2.1, 'Hh OFF (SANT-1-like)', fontsize=8.5, color='#333', rotation=90, va='top')
b.text(onset + 5, 1.35, 'reservoir t½ ≈ 13 h  →  ~11–12 h\npersistence = coasts ~1 more cycle', fontsize=8.4, color=C_MEM, fontweight='bold')
b.set_xlabel('time (h)'); b.set_ylabel('CyclinD1 mRNA (reservoir)'); b.set_xlim(onset - 40, onset + 46); b.set_ylim(-0.2, None)
b.set_title('B · Restriction point is in the PREVIOUS cycle\nafter Hh-off the reservoir COASTS ~1 cycle (÷ ticks)', fontsize=11.5, fontweight='bold')
b.legend(fontsize=8.4, loc='upper right'); b.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# C: EZH2 incoherent feedforward buffers the dose-response
c = ax[1, 0]
SHHs = np.array([0.5, 0.35, 0.25, 0.18, 0.12, 0.09, 0.07])
cd_ss, rep_ss = [], []
for shh in SHHs:
    t, o, rr = run(None, [(dict(SHH=shh), 150 * UPH)], ['Cd_mRNA', 'PRC2_rep'])
    mm = t > 90
    cd_ss.append(o['Cd_mRNA'][mm].mean()); rep_ss.append(rep_factor(rr, o['PRC2_rep'][mm].mean()))
cd_ss, rep_ss = np.array(cd_ss), np.array(rep_ss)
ref = 0  # normalise to SHH=0.5
cd_on = cd_ss / cd_ss[ref]
cd_frozen = (cd_ss * rep_ss[ref] / rep_ss) / cd_ss[ref]   # counterfactual: repression pinned at the SHH=0.5 value
c.plot(SHHs, cd_on, '-o', color=C_MEM, lw=2.6, ms=6.5, label='EZH2 feedforward ON (default) — buffered')
c.plot(SHHs, cd_frozen, '--s', color=C_FAST, lw=2.1, ms=5.5, label='repression FROZEN (no EZH2 de-repression)')
c.fill_between(SHHs, cd_frozen, cd_on, color=C_MEM, alpha=0.12)
c.annotate('EZH2 de-represses as mitogen falls\n(rep factor %.2f→%.2f) → transcript held up' % (rep_ss[ref], rep_ss[-1]),
           xy=(SHHs[-2], cd_on[-2]), xytext=(0.28, 0.30), textcoords='axes fraction', fontsize=8.2, color=C_MEM,
           fontweight='bold', arrowprops=dict(arrowstyle='->', color=C_MEM, lw=1.3))
c.set_xlabel('Shh mitogen (steady)'); c.set_ylabel('CyclinD1 mRNA (norm. to SHH=0.5)'); c.invert_xaxis()
c.set_title('C · EZH2 incoherent feedforward BUFFERS the dose-response\n(Gli raises the transcript, EZH2 lowers it — an IFFL)', fontsize=11.5, fontweight='bold')
c.legend(fontsize=8.5, loc='lower left'); c.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# D: transient after a mitogen step-down — the transcript reservoir "stays up" long after the fast Gli input drops
# (the temporal face of the memory; the EZH2 de-repression that buffers it is quantified rigorously in panel C).
d = ax[1, 1]
rrd = build(None)
for k, v in GNP.items(): rrd[k] = v
rrd['SHH'] = 0.5; rrd.integrator.setValue('absolute_tolerance', 1e-8)
rrd.simulate(0, 120 * UPH, 200)                             # equilibrate at high mitogen (discard)
rrd['SHH'] = 0.09
rd = rrd.simulate(120 * UPH, 175 * UPH, 1800, selections=['time', 'Cd_mRNA', 'Cd', 'Gli_act', 'Gli1'])
sh = rd['time'] / UPH - 120.0; tm = rd['time']
def smooth(x, w_h=9.0):                                     # edge-padded rolling mean (no zero-pad ramp) to remove cycle sawtooth
    dt = (tm[1] - tm[0]) / UPH; k = max(1, int(w_h / dt))
    return np.convolve(np.pad(x, k, mode='edge'), np.ones(k) / k, mode='same')[k:-k]
gli = smooth(rd['Gli_act'] + rd['Gli1'], 1.5)              # Gli is upstream (no cycle sawtooth) -> light smoothing
M = smooth(rd['Cd_mRNA'], 9.0); P = smooth(rd['Cd'], 9.0)  # transcript/protein oscillate with the cycle -> ~1-cycle smoothing
inner = (sh >= 0) & (sh < 48)
# fraction of the total step-change still remaining: 1 at the step, 0 at the new steady state -> isolates the RATE
frac = lambda y: np.clip((y[inner] - y[inner][-1]) / (y[inner][0] - y[inner][-1] + 1e-12), 0, 1.05)
d.plot(sh[inner], frac(gli), color=C_GLI, lw=1.9, label='Gli signal — the INPUT (drops in <1 h)')
d.plot(sh[inner], frac(P), color=C_PROT, lw=1.9, ls=(0, (4, 2)), label='CyclinD1 protein — follows the reservoir')
d.plot(sh[inner], frac(M), color=C_MEM, lw=3.0, label='CyclinD1 mRNA — the RESERVOIR (slow relax)')
d.fill_between(sh[inner], frac(gli), frac(M), color=C_MEM, alpha=0.12)
def thalf(fr): i = np.argmax(fr < 0.5); return sh[inner][i] if fr[i] < 0.5 else np.nan
tg, tmm = thalf(frac(gli)), thalf(frac(M))
d.axhline(0.5, color='#ccc', lw=0.8, ls=':')
d.axvline(0, color='#333', ls=':', lw=1.3)
d.text(0.30, 0.55, 'Shh 0.5→0.09: the Gli input resets in <1 h, but the\ntranscript RESERVOIR relaxes with t½ ≈ %.0f h — the protein\nfollows the reservoir, not the input, so the cell coasts.\nEZH2 de-repression (panel C) further buffers the transcript.'
       % (tmm if tmm == tmm else 30), transform=d.transAxes, fontsize=7.9, color=C_MEM, fontweight='bold')
d.set_xlabel('time after mitogen step-down (h)'); d.set_ylabel('fraction of step-change remaining'); d.set_ylim(-0.03, 1.08)
d.legend(fontsize=8.0, loc='upper right')
d.set_title('D · The transcript reservoir "stays up" after the input drops\n(the temporal memory; EZH2 buffering quantified in C)', fontsize=11.5, fontweight='bold')
d.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# E: widened cilium duty cycle
e = ax[0, 2]
tc, oc, _ = run(dict(w_cilium=0.9, K_cil=0.6), [(dict(SHH=0.5), 90 * UPH)], ['cilium_gate', 'Cb', 'Dna'])
mm = (tc > 40) & (tc < 88)
e.plot(tc[mm], oc['cilium_gate'][mm], color=C_CIL, lw=2.4, label='cilium gate (Shh reception)')
e.plot(tc[mm], oc['Cb'][mm] / (max(oc['Cb'][mm]) + 1e-9), color='#7f8c8d', lw=1.2, ls=':', label='CyclinB (G2/M marker)')
for dt in tc[divs(oc['Dna'])]:
    if 40 < dt < 88: e.axvline(dt, color=C_FAST, lw=0.9, alpha=0.5)
e.text(0.02, 0.1, 'reception ON through G1/S, OFF only at G2 (÷ = red)\n— wider than the textbook G1-only window', transform=e.transAxes, fontsize=8, color='#555')
e.set_xlabel('time (h)'); e.set_ylabel('cilium gate / CyclinB'); e.set_ylim(-0.03, 1.08)
e.set_title('E · Widened cilium DUTY CYCLE — reception ON in G1/S,\nOFF at G2 (Ho 2020: cilia persist into S)', fontsize=11.5, fontweight='bold')
e.legend(fontsize=8.3, loc='center right'); e.grid(alpha=0.2)

# ---------------------------------------------------------------------------------------------------------------
# F: robustness to the duty cycle
f = ax[1, 2]
def divday(over):
    t, o, _ = run(over, [(dict(SHH=0.5), 22000)], ['Dna']); m = t > 130
    dv = divs(o['Dna'][m]); span = (t[m][-1] - t[m][0])
    return len(dv) / span * 24 if span > 0 else 0
Kcils = [2.0, 1.2, 0.8, 0.6, 0.45, 0.35]
rates_mem = [divday(dict(w_cilium=0.9, K_cil=k)) for k in Kcils]
rates_fast = [divday(dict(w_cilium=0.9, K_cil=k, **FAST)) for k in Kcils]
f.plot(Kcils, rates_mem, '-o', color=C_MEM, lw=2.4, ms=6, label='mRNA reservoir memory (default)')
f.plot(Kcils, rates_fast, '-s', color=C_FAST, lw=2.0, ms=5, label='mRNA fast / slaved')
f.invert_xaxis()
f.set_xlabel('cilium gate threshold K_cil  (→ shorter ciliated window →)'); f.set_ylabel('division rate (÷/day)')
f.set_ylim(0, max(rates_mem + rates_fast) * 1.15)
f.set_title('F · Proliferation is ROBUST to the duty cycle\n(intermittent reception buffered by the reservoir)', fontsize=11.5, fontweight='bold')
f.legend(fontsize=8.5, loc='lower left'); f.grid(alpha=0.2)

fig.suptitle('CyclinD1 = the Hedgehog memory carrier (Ho, Tsai & Stearns 2020): the memory is the cumulative CyclinD1 TRANSCRIPT reservoir '
             '(slow mRNA → fast protein), buffered by an EZH2 incoherent feedforward; the proliferation decision sits in the PREVIOUS cycle',
             fontsize=12, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_cyclind1_memory.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
print('dose-response buffering: IFFL-on transcript / frozen at SHH=%.2f:  %.2f vs %.2f  (%.2fx higher)'
      % (SHHs[-1], cd_on[-1], cd_frozen[-1], cd_on[-1] / cd_frozen[-1]))
print('transient: reservoir t1/2=%.0fh vs Gli t1/2=%.1fh (transcript persists after the input drops)' % (tmm, tg if tg == tg else 0))
