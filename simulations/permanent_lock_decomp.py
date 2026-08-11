"""What holds the PERMANENT-LOCK arrest, and does weakening H3K27me3->CyclinD1 give a TRANSIENT G0?
(JP 2026-08-09). Uses the two-compartment module in the permanent-lock config (a_P=0.5, f_P=0.3).

Protocol: arrest at 6x CDKI, then drop CDKI to 1x; measure re-entry (dwell, None = permanent lock) and the
LOCKED-STATE readout (what is repressed). Then independently weaken the CyclinD1 arm (w_cd_prox) and the
CDK6 arm (w_cdk6_prox) of the proximal repression (1 = full, 0 = that arm un-repressed).

Run:  PYTHONPATH=. ./venv/bin/python simulations/permanent_lock_decomp.py
"""
import os, numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)
PERM = dict(a_P=0.50, f_P=0.30)   # permanent-lock config


def run(extra, mult_hi=6.0, t_switch=35000, T=115000):
    P = dict(BP); P.update(PERM); P.update(extra)
    rr = te.loada(build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_proximal_distal=True, params=P))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 15000000)
    def setp(mult):
        for k, v in P.items():
            try: rr[k] = v
            except Exception: pass
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
        rr['p18'] = 1.73 * mult; rr['p19'] = 0.58 * mult; rr['kSyP21'] = BP['kSyP21'] * mult
    sel = ['time', 'Dna', 'P_prox', 'prox_amp', 'Cd', 'cdk6', 'E2f', 'Rb', 'Rbm', 'RbmE2f', 'pRb']
    setp(mult_hi); rr.reset()
    r1 = rr.simulate(0, t_switch, t_switch // 10, selections=sel)
    setp(1.0)   # drop CDKI to baseline; species state (loaded P_prox) carries over
    r2 = rr.simulate(t_switch, T, (T - t_switch) // 10, selections=sel)
    dna = r2['Dna']; dd = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    dwell = float((r2['time'][dd[0]] - t_switch) / 60) if len(dd) else None
    lock = {k: float(r2[k][-1]) for k in ['P_prox', 'prox_amp', 'Cd', 'cdk6', 'E2f', 'Rb', 'Rbm', 'RbmE2f', 'pRb']}
    return dwell, lock


# baseline cycling reference (a normally cycling MB cell) for comparison
def cycling_ref():
    P = dict(BP); P.update(PERM)
    rr = te.loada(build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_proximal_distal=True, params=P))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 15000000)
    for k, v in P.items():
        try: rr[k] = v
        except Exception: pass
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73; rr['p19'] = 0.58; rr['kSyP21'] = BP['kSyP21']
    rr.reset(); r = rr.simulate(0, 40000, 20000, selections=['time', 'Cd', 'cdk6', 'E2f', 'Rbm', 'RbmE2f', 'pRb'])
    mm = r['time'] >= 30000
    return {k: float(r[k][mm].mean()) for k in ['Cd', 'cdk6', 'E2f', 'Rbm', 'RbmE2f', 'pRb']}

ref = cycling_ref()
print("cycling ref (mean): " + "  ".join(f"{k}={v:.3f}" for k, v in ref.items()))
print()
d, lk = run({})
print("PERMANENT LOCK (w_cd=1, w_cdk6=1):  re-entry =", d)
print("  locked state: " + "  ".join(f"{k}={lk[k]:.3f}" for k in lk))
print()
print("--- weaken H3K27me3->CyclinD1 repression only (w_cd_prox: 1 -> 0) ---")
for w in (1.0, 0.75, 0.5, 0.25, 0.0):
    d, lk = run(dict(w_cd_prox=w))
    v = 'TRANSIENT' if d is not None else 'permanent'
    print(f"  w_cd={w:<4}: {v:9s} re-entry={d}   locked Cd={lk['Cd']:.2f} cdk6={lk['cdk6']:.2f} E2f={lk['E2f']:.3f}")
print()
print("--- weaken H3K27me3->CDK6 repression only (w_cdk6_prox: 1 -> 0) ---")
for w in (1.0, 0.75, 0.5, 0.25, 0.0):
    d, lk = run(dict(w_cdk6_prox=w))
    v = 'TRANSIENT' if d is not None else 'permanent'
    print(f"  w_cdk6={w:<4}: {v:9s} re-entry={d}   locked Cd={lk['Cd']:.2f} cdk6={lk['cdk6']:.2f} E2f={lk['E2f']:.3f}")
