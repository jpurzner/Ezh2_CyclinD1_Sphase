"""CORRECTED transient-G0 analysis (p27 marker, alignment-free).

Earlier pRb-based metric was wrong: phospho-Rb jumps to ~max within ~20 min of division (fast
CyclinD-driven phosphorylation; the single Rb species conflates mono- and hyper-phosphorylation),
but the FUNCTIONAL commitment (E2f release, CyclinE/A activation, p27 clearance via the Skp2-p27
feedforward) happens ~10 h later. So the transient G0 = the p27-HIGH / E2f-LOW / pre-S window.

Define (manuscript: G0 = p27-positive):
  G0 = fraction of cycle with (P21 > 0.1  AND  pre-replication)   x period   [alignment-free]
  G1 = pre-replication AND p27 cleared (committed, pre-S)
  S, G2 as before. Arrest = stays p27-high forever (never ignites the feedforward).

Run:  ./venv/bin/python simulations/probe_transient_g0_v2.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "Cd", "P21", "E2f", "aRc", "Dna", "MPF", "mass"]
T_END, N_PTS, SETTLE, P27_HI = 16000, 64000, 7000, 0.1


def sim(params=None, shh=0.5, mycn=1.0, ptch1=1.0, ezh2i=0.0):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=(params or None)))
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = shh; rr['MYCN_amplification'] = mycn; rr['Ptch1_copy_number'] = ptch1; rr['EZH2i'] = ezh2i
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def metrics(res):
    t = res['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    mpf = res['MPF'][m]; P21 = res['P21'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / dt))
    cd = float(np.mean(res['Cd'][m]))
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = Dna >= 0.98; preS = ~in_S & ~in_G2
    G0 = preS & (P21 > P27_HI); G1 = preS & (P21 <= P27_HI)
    if len(pk) < 2:
        # non-cycling: report p27-high fraction (≈1 if permanently arrested in G0)
        return dict(cd=cd, period=np.nan, arrest=True, g0frac=float(G0.mean()))
    period = float(np.mean(np.diff(tt[pk]))) / 60.0
    fr = {k: v.mean() for k, v in dict(G0=G0, G1=G1, S=in_S, G2=in_G2).items()}
    s = sum(fr.values()) or 1.0
    dur = {k: period * v / s for k, v in fr.items()}
    return dict(cd=cd, period=period, arrest=False, g0frac=float(G0.mean()), **dur)


def line(tag, a):
    if a['arrest']:
        return f"{tag}  ARREST (p27-high frac {a['g0frac']*100:.0f}% = permanent G0)"
    return (f"{tag}  Cd={a['cd']:.3f}  period={a['period']:.1f}h  "
            f"G0(p27)={a['G0']:.1f}h  G1={a['G1']:.1f}h  S={a['S']:.1f}h  G2={a['G2']:.1f}h "
            f"(G0={a['G0']/a['period']*100:.0f}%)")


print("=" * 86)
print("CORRECTED transient G0 (p27-high, pre-S) vs CyclinD1 drive — GNP")
print("-" * 86)
print("k_Cd_translation sweep (protein scaled, transcript ~fixed):")
for kt in [0.18, 0.22, 0.26, 0.32, 0.40]:
    print("  " + line(f"kTl={kt:.2f}", metrics(sim({"k_Cd_translation": kt}))))
print("\nK_EZH2_repression sweep (higher = weaker repression = more CyclinD1):")
for ke in [0.3, 0.5, 0.8, 1.5]:
    print("  " + line(f"K_EZH2={ke:.1f}", metrics(sim({"K_EZH2_repression": ke}))))
print("\nSHH dose sweep (mitogen):")
for s in [0.35, 0.5, 0.75, 1.0]:
    print("  " + line(f"SHH={s:.2f}", metrics(sim(shh=s))))

print("\n" + "=" * 86)
print("MB context (Ptch1=0.1, MYCN_amp=2.8) and the EZH2-inhibitor (de-represses CyclinD1):")
print("  " + line("MB         ", metrics(sim(shh=0.5, ptch1=0.1, mycn=2.8))))
print("  " + line("MB + EZH2i ", metrics(sim(shh=0.5, ptch1=0.1, mycn=2.8, ezh2i=1.0))))
print("  " + line("GNP        ", metrics(sim(shh=0.5))))
print("  " + line("GNP + EZH2i", metrics(sim(shh=0.5, ezh2i=1.0))))
