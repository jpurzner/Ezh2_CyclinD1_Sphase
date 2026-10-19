"""How does protein respond to a transcriptional burst, in our model? (JP 2026-10-19)

Controlled pulse experiment on the baked two-step model (with_cdki_mrna): hold a cycling cell at steady
state, then apply a defined square BURST to a CDKI's transcription knob (*_tx: 1 -> AMP for DUR, then back
to 1), and record the transcription input, the mRNA species, and the protein species.

We use p18 as the clean demonstrator: its expression module (transcription -> p18_m -> p18_prot) is fully
DECOUPLED from the cell cycle (no E2f / Skp2 / Cdh1 / division terms), so the traces are the PURE two-stage
low-pass with no confounds. Panel A: the canonical response (input -> mRNA -> protein). Panel B: the SAME
burst with the protein half-life set to 0.6h (p21-like), 3.9h (p18 model-tuned) and 23.9h (p18 measured) --
mRNA stage held fixed -- to show that a stable protein filters the burst out. p21/p27 map onto the fast/slow
ends of this axis; p21's short protein is why it (not p18/p27) transmits bursts into transient-G0.
"""
import os, sys, numpy as np, tellurium as te
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

LN2 = 0.6931471805599453
WARM = 9000.0          # long warmup so even a 24h-protein settles to steady state
REC = 4800.0           # 80 h record window
DTR = 15.0             # fine recording step (min)
PULSE_T0 = 600.0       # burst starts 10 h into the record
PULSE_DUR = 480.0      # 8 h burst
AMP = 4.0              # transcription x4 during the burst
SEL = ['p18_m', 'p18_prot']


def run_pulse(kDe_p18_override=None):
    """Warm up a cycling MB cell, then apply the p18 transcription burst; return t, tx(t), mRNA(t), prot(t)."""
    rr = te.loada(build_model_v44())                 # baked default = two-step
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 20000000)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    if kDe_p18_override is not None:
        rr['kDe_p18'] = kDe_p18_override             # protein decay (also scales synthesis; self-normalizing -> ss = p18)
    rr['p18_tx'] = 1.0
    rr.simulate(0, WARM, 50)                         # settle
    t, tx, mr, pr = [], [], [], []
    clock = 0.0
    while clock < REC:
        in_pulse = PULSE_T0 <= clock < PULSE_T0 + PULSE_DUR
        rr['p18_tx'] = AMP if in_pulse else 1.0
        r = rr.simulate(clock, clock + DTR, 2, selections=SEL)
        t.append(clock + DTR); tx.append(AMP if in_pulse else 1.0)
        mr.append(float(r['p18_m'][-1])); pr.append(float(r['p18_prot'][-1]))
        clock += DTR
    t = np.array(t) / 60.0                            # hours
    mr = np.array(mr); pr = np.array(pr)
    base = slice(0, int(PULSE_T0 / DTR) - 1)          # pre-pulse baseline for fold-normalisation
    return t, np.array(tx), mr / np.mean(mr[base]), pr / np.mean(pr[base])


# Panel A: canonical response (model-tuned p18, t1/2 3.9h)
tA, txA, mrA, prA = run_pulse()
# Panel B: protein-stability sweep (same mRNA stage; protein t1/2 = 0.6 / 3.9 / 23.9 h)
cases = [('p21-like 0.6 h', LN2 / (0.6 * 60), '#c0392b'),
         ('p18 tuned 3.9 h', LN2 / (3.9 * 60), '#e8a33d'),
         ('p18 measured 23.9 h', LN2 / (23.9 * 60), '#2c6fa6')]
sweep = [(lab, col, run_pulse(k)) for lab, k, col in cases]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(14, 5.4))
fig.suptitle("How protein responds to a transcriptional burst (two-step model): transcription → mRNA → protein, each filtered by its half-life",
             fontsize=12.5, fontweight='bold')

# A
axA.axvspan(PULSE_T0/60, (PULSE_T0+PULSE_DUR)/60, color='#ddd', alpha=0.6, lw=0, label='burst (tx ×4, 8 h)')
axA.plot(tA, txA, color='#777', lw=1.4, ls=':', label='transcription rate (input)')
axA.plot(tA, mrA, color='#c0392b', lw=2.0, label='p18 mRNA (t½ 3.7 h)')
axA.plot(tA, prA, color='#1f6f54', lw=2.2, label='p18 protein (t½ 3.9 h)')
axA.set_xlabel('time (h)'); axA.set_ylabel('fold vs baseline'); axA.legend(fontsize=9, loc='upper right')
axA.set_title("(A) A step in transcription → mRNA rises/decays at its t½,\nprotein lags & smooths further — a two-stage low-pass", fontsize=10.5)
axA.axhline(1, color='#ccc', lw=0.8)

# B
axB.axvspan(PULSE_T0/60, (PULSE_T0+PULSE_DUR)/60, color='#ddd', alpha=0.6, lw=0, label='burst (tx ×4, 8 h)')
for lab, col, (t, tx, mr, pr) in sweep:
    axB.plot(t, pr, color=col, lw=2.2, label=f'protein {lab}')
axB.plot(tA, mrA, color='#555', lw=1.3, ls='--', label='mRNA (shared, t½ 3.7 h)')
axB.axhline(1, color='#ccc', lw=0.8)
axB.set_xlabel('time (h)'); axB.set_ylabel('protein fold vs baseline'); axB.legend(fontsize=9, loc='upper right')
axB.set_title("(B) Same burst, protein t½ varied (mRNA stage fixed):\nstable protein filters the burst out — why p21 (fast) transmits it, p18/p27 (stable) buffer", fontsize=10.5)

plt.tight_layout(rect=[0, 0, 1, 0.93])
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'overnight', 'fig_pulse_response.png')
plt.savefig(out, dpi=140); print('saved', out)
# quick numeric summary
peak_m = mrA.max(); peak_p = prA.max()
print(f"p18 (tuned): mRNA peak {peak_m:.2f}x, protein peak {peak_p:.2f}x fold")
for lab, col, (t, tx, mr, pr) in sweep:
    print(f"  protein peak [{lab}]: {pr.max():.2f}x")
