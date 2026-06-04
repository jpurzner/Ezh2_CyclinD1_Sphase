"""Why is the transient G0 ~0, and what controls its duration?

Hypothesis: the daughter is born at mass ~1.5 but the commitment size gate M_commit=1.3, so the
gate is ALREADY OPEN at birth -> the Skp2-p27-Rb feedforward fires immediately -> no growth-delayed
p27-high G0. Test: raise M_commit above birth mass and measure the resulting transient G0 (p27-high
window), and show how CyclinD1 protein modulates it once the gate opens.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "Cd", "pRb", "P21", "Skp2", "aRc", "Dna", "MPF", "mass"]
T_END, N_PTS, SETTLE, P21_LO = 16000, 64000, 7000, 0.1


def sim(params, shh=0.5, mycn=1.0, ptch1=1.0):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=params))
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = shh; rr['MYCN_amplification'] = mycn; rr['Ptch1_copy_number'] = ptch1
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def divtimes(res):
    t = res['time']; mass = res['mass']
    d = np.where((mass[:-1] - mass[1:]) > 0.2 * mass[:-1])[0]
    if not len(d): return np.array([])
    keep = [d[0]]
    for x in d[1:]:
        if x - keep[-1] > 50: keep.append(x)
    return t[keep]


def metrics(res):
    t = res['time']; m = t >= SETTLE; tt = t[m]
    P21 = res['P21'][m]; pRb = res['pRb'][m]; mass = res['mass'][m]
    mpf = res['MPF'][m]; dt = tt[1]-tt[0]
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200/dt))
    if len(pk) < 2:
        return None
    period = float(np.mean(np.diff(tt[pk])))/60.0
    divs = divtimes(res); divs = divs[divs >= SETTLE]
    g0, birth = [], []
    for i in range(len(divs)-1):
        a, b = divs[i], divs[i+1]
        seg = (tt >= a) & (tt < b)
        ts, p27 = tt[seg], P21[seg]
        birth.append(mass[seg][0] if seg.any() else np.nan)
        below = np.where(p27 < P21_LO)[0]
        g0.append((ts[below[0]]-a)/60.0 if len(below) else (b-a)/60.0)
    # phase fractions for preS/S/G2
    aRc = res['aRc'][m]; Dna = res['Dna'][m]
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = Dna >= 0.98; preS = ~in_S & ~in_G2
    return dict(period=period, g0=float(np.mean(g0)), birth=float(np.mean(birth)),
                cd=float(np.mean(res['Cd'][m])), preS=period*preS.mean(),
                S=period*in_S.mean(), G2=period*in_G2.mean())


print("=" * 84)
print("M_commit sweep (GNP, SHH=0.5): transient G0 appears once M_commit > birth mass")
print(f"{'M_commit':>9} {'birth_mass':>11} {'meanCd':>7} {'period':>7} {'G0(p27)':>8} {'preS':>6} {'S':>6} {'G2':>6}")
for mc in [1.3, 1.6, 1.8, 2.0, 2.2, 2.4]:
    r = sim({"M_commit": mc})
    a = metrics(r)
    if a is None:
        print(f"{mc:>9.1f}   ARREST / non-cycling"); continue
    print(f"{mc:>9.1f} {a['birth']:>11.2f} {a['cd']:>7.3f} {a['period']:>6.1f}h "
          f"{a['g0']:>7.2f}h {a['preS']:>5.1f}h {a['S']:>5.1f}h {a['G2']:>5.1f}h")

print("\n" + "=" * 84)
print("At M_commit=2.0 (a real growth-gated G0), vary CyclinD1 via k_Cd_translation:")
print(f"{'kTl':>6} {'meanCd':>7} {'period':>7} {'G0(p27)':>8} {'preS':>6}")
for kt in [0.20, 0.26, 0.32, 0.40]:
    r = sim({"M_commit": 2.0, "k_Cd_translation": kt})
    a = metrics(r)
    if a is None:
        print(f"{kt:>6.2f}   ARREST"); continue
    print(f"{kt:>6.2f} {a['cd']:>7.3f} {a['period']:>6.1f}h {a['g0']:>7.2f}h {a['preS']:>5.1f}h")
