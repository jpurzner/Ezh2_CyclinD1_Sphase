"""Calibrate the GLI1-autoregulation slow epigenetic memory (Gli1_epi) into the main v44 model.

Goal: with the memory ON, a vismo-treated (previously cycling) MB cell holds Gli1 up then declines
gradually over several cell cycles -> COASTS through extra divisions, while:
  - GNP is unchanged (memory gated to broken feedback),
  - MB BASELINE targets are preserved (Gli1 MB/GNP ~6.9x, Cd MB/GNP ~7.58x) via re-attribution
    (g_smo_Gli1_broken reduces MB's Smo-driven Gli1 so the memory supplies the rest),
  - the MB+HHi vismo ENDPOINT is preserved (in that condition vismo is on from t=0 so the memory never
    charges -> Cd -> the MYCN floor, 0.144xMB).

Search (g_smo_Gli1_broken, k_Gli1_auto) at fixed memory kinetics; report baseline + coasting + checks.
Run:  ./venv/bin/python simulations/v44_calibrate_gli_memory.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

GNP = dict(f=1.0, mycn=1.0,  p16=0.0,  p18=0.4, ksy=0.0020)
MB  = dict(f=0.1, mycn=2.86, p16=0.88, p18=1.2, ksy=0.0040)
# fixed memory kinetics (engagement + slow discharge); search the two re-attribution/strength knobs
MEM = dict(K_Gli1_auto=0.3, n=4, k_epi_on=0.02, k_epi_off=0.0005)


def _mk(P):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    rr.integrator.setValue("relative_tolerance", 1e-7)
    rr['K_Gli1_auto'] = P['K_Gli1_auto']; rr['n_Gli1_auto'] = P['n']
    rr['k_epi_on'] = P['k_epi_on']; rr['k_epi_off'] = P['k_epi_off']
    rr['k_Gli1_auto'] = P['k_Gli1_auto']; rr['g_smo_Gli1_broken'] = P['g_smo']
    return rr


def _seg(rr, t0, t1, n):
    for a in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        try:
            rr.integrator.setValue("absolute_tolerance", a)
            rr.integrator.setValue("maximum_num_steps", 2000000)
            return rr.simulate(t0, t1, n, selections=["time", "MPF", "Cd", "Gli1", "Gli1_epi"])
        except Exception:
            continue
    return None


def _set(rr, c, gdc):
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = c['f']; rr['MYCN_amplification'] = c['mycn']
    rr['p16'] = c['p16']; rr['p18'] = c['p18']; rr['kSyP21'] = c['ksy']; rr['EZH2i'] = 0; rr['GDC0449'] = gdc


def baseline(P, c, T=9000):
    rr = _mk(P); rr.reset(); _set(rr, c, 0.0)
    r = _seg(rr, 0, T, T * 3)
    if r is None: return np.nan, np.nan
    m = r['time'] >= 5000
    return float(np.mean(r['Gli1'][m])), float(np.mean(r['Cd'][m]))


def endpoint_hhi(P, c, T=12000):
    rr = _mk(P); rr.reset(); _set(rr, c, 1.0)         # vismo from t=0: memory never charges
    r = _seg(rr, 0, T, T * 3)
    if r is None: return np.nan
    return float(np.mean(r['Cd'][r['time'] >= T - 3000]))


def withdrawal(P, c, T1=8000, TEND=46000):
    rr = _mk(P); rr.reset(); _set(rr, c, 0.0)
    s1 = _seg(rr, 0, T1, T1 * 3)
    if s1 is None: return None
    rr['GDC0449'] = 1.0
    s2 = _seg(rr, T1, TEND, (TEND - T1) * 3)
    if s2 is None: return None
    t = np.concatenate([s1['time'], s2['time']]); mpf = np.concatenate([s1['MPF'], s2['MPF']])
    gli = np.concatenate([s1['Gli1'], s2['Gli1']])
    dt = np.median(np.diff(t)); pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / dt))
    g0 = np.mean(gli[(t >= 5000) & (t <= T1)])
    samp = [float(gli[np.argmin(np.abs(t - (T1 + d)))]) / g0 for d in (1400, 4200, 8400, 16000)]
    return int(np.sum(t[pk] >= T1)), samp


# GNP reference (memory gated off -> unchanged)
gnp_g, gnp_cd = baseline(dict(MEM, k_Gli1_auto=0.0, g_smo=1.0), GNP)
print(f"GNP baseline: Gli1={gnp_g:.4f}  Cd={gnp_cd:.3f}  (targets: MB/GNP Gli1 6.9x -> MB Gli1~{6.9*gnp_g:.3f})")
gnp_wd = withdrawal(dict(MEM, k_Gli1_auto=0.0, g_smo=1.0), GNP)
print(f"GNP vismo withdrawal: divisions_to_arrest={gnp_wd[0]}\n")

MB_GLI_TGT = 6.9 * gnp_g
print(f"Search (k_Gli1_auto, g_smo) to hit MB baseline Gli1 ~ {MB_GLI_TGT:.3f}:")
print(f"{'k_auto':>7}{'g_smo':>7}{'MBGli1':>8}{'MBCd':>7}{'MB/GNP':>7}{'HHiCd/MB':>9}{'div_aft':>8}  Gli1 +1/+3/+6/+12 cyc")
for kauto in [0.06, 0.10, 0.14, 0.18]:
    # bisect g_smo in [0,1] to hit MB baseline Gli1 target
    lo, hi = 0.0, 1.0
    for _ in range(7):
        g = 0.5 * (lo + hi)
        mbg, _ = baseline(dict(MEM, k_Gli1_auto=kauto, g_smo=g), MB)
        if np.isnan(mbg): break
        if mbg > MB_GLI_TGT: hi = g
        else: lo = g
    g = 0.5 * (lo + hi)
    P = dict(MEM, k_Gli1_auto=kauto, g_smo=g)
    mbg, mbcd = baseline(P, MB)
    hhi = endpoint_hhi(P, MB)
    wd = withdrawal(P, MB)
    if wd is None or np.isnan(mbg):
        print(f"{kauto:7.2f}{g:7.3f}  (sim failed)"); continue
    samp = "/".join(f"{x:.2f}" for x in wd[1])
    print(f"{kauto:7.2f}{g:7.3f}{mbg:8.3f}{mbcd:7.2f}{mbg/gnp_g:7.2f}{hhi/mbcd:9.3f}{wd[0]:8d}  {samp}")
