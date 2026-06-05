"""Is the model's ~3x CyclinD1 de-repression supported by the Ezh2 cKO data?

Raw experiment: cKO 2.46x (Fig 3A, log2FC 1.301), Tazemetostat 2.2x (Fig 3C). Model: ~3x.
Two reasons the measured bulk fold UNDER-estimates the true per-cell feedback (JP's caveat):
  (1) Incomplete Math1-Cre recombination -> bulk RNA-seq mixes recombined + WT-escaper cells, diluting
      the de-repression. Per-cell fold X from a bulk fold B at recombination efficiency f:
          B = f*X + (1-f)*1   =>   X = (B - (1-f)) / f
  (2) The cKO deletes the SET domain (methyltransferase) but leaves EZH2 PROTEIN -> possible residual
      (non-catalytic) repression; and H3K27me3 dilution needs cell divisions, so 24 h Tazemetostat is
      also incomplete. Both make even the per-cell fold a LOWER BOUND.

Panel (a): per-cell fold vs Math1-Cre efficiency (model 3x marked).
Panel (b): model feedback strength (K_EZH2_repression) -> measured de-repression fold + baseline
           CyclinD1, so we can see what feedback strength the data imply.

Run:  ./venv/bin/python simulations/fig_v44_ckoo_feedback_strength.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import mean_settled, SEL

# ---- (1) recombination correction ----
f = np.linspace(0.5, 1.0, 26)
def percell(B): return (B - (1 - f)) / f
X_cko = percell(2.46)      # cKO bulk
X_taz = percell(2.20)      # Tazemetostat bulk
print("Per-cell CyclinD1 de-repression implied by the bulk cKO (2.46x):")
for ff in (1.0, 0.9, 0.8, 0.7, 0.6):
    print(f"  Math1-Cre efficiency {ff:.0%}: per-cell fold = {(2.46-(1-ff))/ff:.2f}x")

# ---- (2) model feedback strength sweep ----
def run(K, shh=0.5, ezh2i=0):
    for atol in (1e-9, 1e-8, 1e-7):
        rr = te.loada(build_model_v44(params={"K_EZH2_repression": K}))
        rr.integrator.setValue("absolute_tolerance", atol); rr.integrator.setValue("relative_tolerance", 1e-6)
        rr['SHH'] = shh; rr['MYCN_amplification'] = 1.0; rr['Ptch1_copy_number'] = 1.0; rr['EZH2i'] = ezh2i
        try:
            return rr.simulate(0, 10080, 20160, selections=SEL)
        except Exception:
            continue
    return None

Ks = [0.20, 0.30, 0.40, 0.50, 0.70, 1.00]
fold_m, base_m = [], []
for K in Ks:
    on = run(K, ezh2i=0); off = run(K, ezh2i=1)     # ezh2i=1 removes repression -> the cKO/inhibitor readout
    cd_on = mean_settled(on, 'Cd_mRNA'); cd_off = mean_settled(off, 'Cd_mRNA')
    fold_m.append(cd_off / cd_on); base_m.append(cd_on)
    print(f"  K_EZH2_repression={K:.2f}: CyclinD1 baseline={cd_on:.2f}, de-repression fold={cd_off/cd_on:.2f}x")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))

ax[0].plot(f, X_cko, '#762A83', lw=2.4, label='from cKO bulk 2.46x (Fig 3A)')
ax[0].plot(f, X_taz, '#E08214', lw=2.0, ls='--', label='from Tazemetostat 2.2x (Fig 3C)')
ax[0].axhline(3.0, color='#2c3e50', ls=':', lw=1.8, label='model (~3x)')
ax[0].axvspan(0.8, 0.95, color='gray', alpha=0.12)
ax[0].text(0.875, ax[0].get_ylim()[1]*0.2 if False else 4.6, 'typical Math1-Cre\nefficiency', ha='center', fontsize=8, color='#555')
ax[0].set_xlabel('Math1-Cre recombination efficiency'); ax[0].set_ylabel('implied per-cell CyclinD1 fold')
ax[0].set_title('(a) The bulk cKO fold UNDER-estimates per-cell de-repression\n(escaper cells dilute it)', fontweight='bold', fontsize=11)
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.2)
ax[0].invert_xaxis()

ax2 = ax[1]; ax2b = ax2.twinx()
ax2.plot(Ks, fold_m, '#2c3e50', lw=2.4, marker='o', label='de-repression fold (cKO/EZH2i readout)')
ax2b.plot(Ks, base_m, '#1b9e77', lw=2.0, marker='s', ls='--', label='baseline CyclinD1')
ax2.axhspan(2.2, 3.1, color='#fdebd0', alpha=0.5)
ax2.set_xlabel('K_EZH2_repression  (smaller = STRONGER feedback)'); ax2.set_ylabel('de-repression fold', color='#2c3e50')
ax2b.set_ylabel('baseline CyclinD1', color='#1b9e77')
ax2.set_title('(b) Model feedback strength sets the measurable fold\n(stronger feedback = bigger fold, lower baseline)', fontweight='bold', fontsize=11)
ax2.invert_xaxis(); ax2.legend(fontsize=8, loc='upper left'); ax2b.legend(fontsize=8, loc='lower left'); ax2.grid(alpha=0.2)

fig.suptitle('Is the model ~3x CyclinD1 supported by the Ezh2 cKO? (accounting for incomplete cKO)',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_ckoo_feedback_strength.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ckoo_feedback_strength.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_ckoo_feedback_strength.png")
