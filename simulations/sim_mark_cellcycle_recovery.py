"""The cell-cycle duration sets how much H3K27me3 recovers from S-phase dilution -> the MISSING term in the
buffering equation (JP 2026-07-20). The mark is HALVED at S (replication) and re-methylated over the rest of the
cycle; so the steady-state mark depends on the recovery time = the cell-cycle length Tc. The dilution acts as an
effective decay of ln2/Tc ON TOP of the intrinsic turnover del_mk. In matched units the two are COMPARABLE
(del_mk~0.028/h vs ln2/Tc~0.032/h at Tc=22h); dilution DOMINATES at fast cycling (short Tc), turnover at slow,
crossover Tc~ln2/del_mk. Either way the cell cycle is a first-order contributor the old del_mk-only equation missed.

  effective mark decay:   delta_eff = del_mk + ln2/Tc      (comparable terms; dilution wins at fast Tc)
  steady mark:            P* ~ w / (del_mk + ln2/Tc)       (w = writing rate, EZH2*PRC2)
  since Tc is mitogen-dependent (graded G1-lengthening), the mark & its buffering timescale are cell-cycle-coupled.

We (A) show the within-cycle mark trajectory (S-dip + recovery) for SHORT vs LONG Tc, in a de-saturated mark regime
where it is visible; (B) sweep Tc (via growth rate mu) and check P* ~ w/(del_mk + ln2/Tc). Run:
  ./venv/bin/python simulations/sim_mark_cellcycle_recovery.py
"""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
# de-saturated mark so the S-dip + recovery are VISIBLE (default mark is saturated ~0.98, masks it)
DESAT = dict(a_rw_prc2=0.0005, del_mk=0.0024, K_prc2=0.004)
M = build_model_v44()


def run(mu, over, shh=0.5, T=26000):
    rr = te.loada(M); rr.reset(); rr['SHH'] = shh; rr['mu_min_frac'] = 1.0; rr['mu'] = mu
    for k, v in over.items(): rr[k] = v
    for atol in (1e-8, 1e-7, 1e-6, 1e-5):
        try: rr.integrator.setValue('absolute_tolerance', atol); r = rr.simulate(0, T, int(T / 2), selections=['time', 'Dna', 'Mk', 'EZH2']); break
        except Exception: r = None
    return r


if __name__ == '__main__':
    res = {}
    print("Cell-cycle duration -> H3K27me3 recovery from S-phase dilution.\n")

    # ---- A. within-cycle trajectories (de-saturated mark), short vs long Tc ----
    print("=== A. within-cycle mark trajectory (de-saturated), short vs long Tc ===")
    for lab, mu in (('short_Tc', 0.0011), ('long_Tc', 0.00035)):
        r = run(mu, DESAT)
        t = np.asarray(r['time']); m = t > 12000
        tt = t[m]; dna = np.asarray(r['Dna'])[m]; mk = np.asarray(r['Mk'])[m]
        dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
        Tc = float(np.mean(np.diff(tt[dv])) / UPH) if len(dv) >= 2 else np.nan
        # extract ~2 cycles from the last divisions for the trace
        if len(dv) >= 3:
            i0, i1 = dv[-3], min(dv[-1] + int((dv[-1]-dv[-2])*0.3), len(tt)-1)
            res[lab] = dict(t=(tt[i0:i1] - tt[i0]) / UPH, mk=mk[i0:i1], dna=dna[i0:i1], Tc=Tc, mk_mean=float(mk.mean()))
        print(f"  {lab}: Tc={Tc:.1f}h  mean Mk={mk.mean():.3f}  (dip-to-recover each cycle)")

    # ---- B. sweep Tc (via mu): mark level vs Tc; check P* ~ w/(del_mk + ln2/Tc) ----
    print("\n=== B. mark level vs cell-cycle duration Tc (de-saturated) ===")
    print(f"  {'mu':>8} {'Tc(h)':>7} {'ln2/Tc':>7} {'del_mk+ln2/Tc':>13} {'Mk':>7} {'EZH2(~w)':>9}")
    B = []
    for mu in (0.0012, 0.0009, 0.0007, 0.0005, 0.00035, 0.00025):
        r = run(mu, DESAT); t = np.asarray(r['time']); m = t > 12000
        tt = t[m]; dna = np.asarray(r['Dna'])[m]
        dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
        Tc = float(np.mean(np.diff(tt[dv])) / UPH) if len(dv) >= 2 else np.nan
        mk = np.asarray(r['Mk'])[m].mean(); ez = np.asarray(r['EZH2'])[m].mean()
        dil = np.log(2) / Tc if Tc == Tc else np.nan
        B.append((mu, Tc, dil, DESAT['del_mk'] + dil / UPH, mk, ez))   # ln2/Tc in 1/h; del_mk in 1/model-unit -> convert
        print(f"  {mu:>8.5f} {Tc:>7.1f} {dil:>7.3f} {DESAT['del_mk'] + dil/UPH:>13.4f} {mk:>7.3f} {ez:>9.3f}")
    res['sweep'] = np.array(B); res['del_mk'] = DESAT['del_mk']

    with open(os.path.join(HERE, 'sim_mark_cellcycle_recovery_cache.pkl'), 'wb') as fh:
        pickle.dump(res, fh)
    # quick fit check: Mk vs w/(del_mk + ln2/Tc)  (w ~ EZH2, both in model units; ln2/Tc converted to 1/model-unit)
    Bm = np.array(B); Tc = Bm[:, 1]; mk = Bm[:, 4]; ez = Bm[:, 5]
    dil_mu = np.log(2) / (Tc * UPH)                                    # dilution rate in 1/model-unit
    pred = ez / (DESAT['del_mk'] + dil_mu); pred = pred / pred.max() * mk.max()
    print(f"\n  fit P* ~ EZH2/(del_mk + ln2/Tc): corr(Mk, pred) = {np.corrcoef(mk, pred)[0,1]:.3f}")
    print("  -> figure: fig_mark_cellcycle_recovery.py")
