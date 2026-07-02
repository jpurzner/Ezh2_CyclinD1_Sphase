"""Calibrate the mean-field AUM H3K27me3 module (transcription->PRC2 eviction arm + EZH2-scaled
methylation + lit-review-paced turnover). Co-varies {del_mk, k_w_mk} (turnover vs read-write, which
together set the mark t1/2 and M_ss), holding the other params, and scores the four validate
CyclinD1 targets at the SAME HHi=1.0 / stiff-retry protocol validate_v44.py uses, plus the EZH2i
de-repression half-time (lit-review target ~24 h). Goal: slowest turnover (most dilution-dependent /
lit-review-paced) that keeps all four CyclinD1 targets in band.

Usage: ./venv/bin/python simulations/calibrate_h3k27_aum.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

MB = dict(ptch1_cn=0.1, mycn_amp=2.8, p16=0.306, p18=1.553, ksyp21=0.002)
GNP = dict(ptch1_cn=1.0, mycn_amp=1.0, p16=0.0, p18=0.464, ksyp21=0.002)
# validate target bands (target, tol)
BANDS = {'GNP+HHi/GNP': (0.157, 0.30), 'MB+HHi/MB': (0.144, 0.30),
         'MB/GNP': (5.07, 0.35), 'EZH2i fold': (2.2, 0.30)}


def _run(rr, hhi=0.0, ezh2i=0.0, cond=GNP, t_end=12000, n_pts=48000):
    last = None
    for te_end, te_pts in ((t_end, n_pts), (8000, 32000), (6000, 24000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr.reset()
            try:
                rr.integrator.setValue("maximum_num_steps", 300000)
            except Exception:
                pass
            rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = cond['ptch1_cn']; rr['HHi'] = hhi
            rr['EZH2i'] = ezh2i; rr['MYCN_amplification'] = cond['mycn_amp']
            rr['p16'] = cond['p16']; rr['p18'] = cond['p18']; rr['kSyP21'] = cond['ksyp21']
            rr.integrator.setValue("absolute_tolerance", atol)
            try:
                return rr.simulate(0, te_end, te_pts, selections=['time', 'Cd', 'Mk'])
            except Exception as e:
                last = e
    raise last


def cd(rr, **kw):
    r = _run(rr, **kw); return float(np.mean(r['Cd'][r['time'] >= 4000]))


def derep_thalf(rr):
    """MB+HHi+EZH2i: equilibrate vismo (HHi=0.95) -> add EZH2i -> time for Mk to halve (h)."""
    rr.reset()
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.1; rr['MYCN_amplification'] = 2.8
    rr['p16'] = 0.306; rr['p18'] = 1.553; rr['kSyP21'] = 0.002; rr['HHi'] = 0.95; rr['EZH2i'] = 0
    rr.integrator.setValue("absolute_tolerance", 1e-9)
    rr.simulate(0, 8000, 16000); mk0 = rr['Mk']
    rr['EZH2i'] = 1
    r = rr.simulate(0, 5000, 10000, selections=['time', 'Mk'])
    th = r['time'] / 60.0; hit = np.where(r['Mk'] <= mk0 / 2)[0]
    return mk0, (th[hit[0]] if len(hit) else float('nan'))


def score(delm, kw, g=0.30, ktx=5.0, k0=0.000258, K=0.305, n=4.15, f0=0.233):
    p = dict(del_mk=delm, k_w_mk=kw, g_mk=g, K_tx_mk=ktx, k0_mk=k0, K_mk=K, n_mk=n, f0_mk=f0)
    rr = te.loada(build_model_v44(with_hh=True, with_h3k27_dilution=True, params=p))
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    g_cd = cd(rr, cond=GNP); mb_cd = cd(rr, cond=MB)
    vals = {'GNP+HHi/GNP': cd(rr, hhi=1.0, cond=GNP) / g_cd,
            'MB+HHi/MB': cd(rr, hhi=1.0, cond=MB) / mb_cd,
            'MB/GNP': mb_cd / g_cd,
            'EZH2i fold': cd(rr, ezh2i=1.0, cond=GNP) / g_cd}
    mk0, th = derep_thalf(rr)
    return vals, th


if __name__ == "__main__":
    print(f"{'del_mk':>7} {'k_w':>7} | " + " ".join(f"{k:>12}" for k in BANDS) + f" | {'t1/2(h)':>7}  pass")
    for delm, kw in [(0.00117, 0.00261), (0.00090, 0.00200), (0.00070, 0.00160),
                     (0.00060, 0.00150), (0.00050, 0.00130), (0.00080, 0.00180)]:
        vals, th = score(delm, kw)
        oks = []
        cells = []
        for k, (tg, tol) in BANDS.items():
            ok = abs(vals[k] - tg) / tg <= tol
            oks.append(ok)
            cells.append(f"{vals[k]:7.3f}{'*' if ok else ' '}    ")
        print(f"{delm:.5f} {kw:.5f} | " + " ".join(cells) + f" | {th:7.1f}  {sum(oks)}/4")
