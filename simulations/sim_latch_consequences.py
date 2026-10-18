"""What happens to the FULL v44 model in the LATCH regime -- the consequences (JP 2026-07-21).
Three experiments, BAKED (monostable) vs LATCH (a0_prc2=1e-5, n_rw=3, a_rw_prc2=0.02):
  (1) WITHDRAWAL-MEMORY time course: proliferate (SHH on) -> withdraw SHH -> RESTORE SHH.
      Does the cell re-proliferate (no memory) or stay arrested (latched memory)?
  (2) MITOGEN hysteresis: SHH swept up then down (continuation) -> loop = memory.
  (3) key VALIDATION metrics baked vs latch (what breaks).
Caches sim_latch_consequences_cache.pkl. Run: ./venv/bin/python simulations/sim_latch_consequences.py"""
import os, sys, pickle, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8

LATCH = dict(a0_prc2=1e-5, n_rw=3, a_rw_prc2=0.02)
GNP = dict(HHi=0, MYCN_amplification=1.0, Ptch1_copy_number=1.0)
MB = dict(SHH=0.5, HHi=0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, P21_div=0.6)

def divs_from(dna):
    dna = np.asarray(dna); return (dna[:-1] > 0.9) & (dna[1:] < 0.1)

# ---------- (1) withdrawal protocol (segmented, state carried over) ----------
def withdrawal(over, seg=(('SHH', 0.5, 11000), ('SHH', 0.03, 13000), ('SHH', 0.5, 14000)), mk0=0.02):
    rr = te.loada(build_model_v44(params=over)); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr['Mk'] = mk0; rr['m2_me'] = 0.05; rr['m1_me'] = 0.05
    rr.integrator.setValue('absolute_tolerance', 1e-7)
    T, Mk, Cd, Dna, EZ, SH = [], [], [], [], [], []
    t0 = 0.0
    for key, val, dur in seg:
        rr[key] = val
        r = rr.simulate(t0, t0 + dur, int(dur / 6), selections=['time', 'Mk', 'Cd', 'Dna', 'EZH2'])
        T += list(r['time']); Mk += list(r['Mk']); Cd += list(r['Cd']); Dna += list(r['Dna']); EZ += list(r['EZH2'])
        SH += [val] * len(r['time']); t0 += dur
    T = np.array(T) / UPH
    dv = np.where(divs_from(Dna))[0]
    return dict(t=T, Mk=np.array(Mk), Cd=np.array(Cd), Dna=np.array(Dna), EZ=np.array(EZ), SH=np.array(SH),
                div_t=T[dv], segbounds=[s[2] for s in seg])

# ---------- (2) SHH continuation (up then down) ----------
def continuation(over, shhs, settle=9000):
    rr = te.loada(build_model_v44(params=over)); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr['Mk'] = 0.02; rr['m2_me'] = 0.05; rr['m1_me'] = 0.05
    rr.integrator.setValue('absolute_tolerance', 1e-7)
    def sweep(seq):
        mk, cd = [], []
        for s in seq:
            rr['SHH'] = s
            try:
                r = rr.simulate(0, settle, 1500, selections=['time', 'Mk', 'Cd'])
            except Exception:
                r = rr.simulate(0, settle, 1500, selections=['time', 'Mk', 'Cd'])
            m = np.asarray(r['time']) > settle * 0.6
            mk.append(float(np.asarray(r['Mk'])[m].mean())); cd.append(float(np.asarray(r['Cd'])[m].mean()))
        return np.array(mk), np.array(cd)
    up_mk, up_cd = sweep(shhs)          # low -> high (carry state)
    dn_mk, dn_cd = sweep(shhs[::-1])    # high -> low (carry state)
    return dict(shh=shhs, up_mk=up_mk, up_cd=up_cd, dn_mk=dn_mk[::-1], dn_cd=dn_cd[::-1])

# ---------- (3) key validation metrics ----------
def settle_state(over, cond, mk0=0.15, T=18000):
    rr = te.loada(build_model_v44(params=over)); rr.reset()
    for k, v in cond.items(): rr[k] = v
    rr['Mk'] = mk0
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, T, 3500, selections=['time', 'Mk', 'Cd', 'EZH2', 'Dna'])
            m = np.asarray(r['time']) > T * 0.55
            return dict(Mk=float(np.asarray(r['Mk'])[m].mean()), Cd=float(np.asarray(r['Cd'])[m].mean()),
                        EZ=float(np.asarray(r['EZH2'])[m].mean()), div=int(np.sum(divs_from(r['Dna']))))
        except Exception:
            rr.reset()
            for k, v in cond.items(): rr[k] = v
            rr['Mk'] = mk0
    return None

def val_metrics(over):
    g = settle_state(over, dict(SHH=0.5, **GNP)); mb = settle_state(over, MB)
    gh = settle_state(over, dict(SHH=0.5, HHi=1, MYCN_amplification=1.0, Ptch1_copy_number=1.0))
    gm = settle_state(over, dict(SHH=0.0, **GNP))
    return dict(Cd_GNP=g['Cd'], Cd_MB=mb['Cd'], fold=mb['Cd'] / g['Cd'], hhi_ratio=gh['Cd'] / g['Cd'],
                chip=mb['Mk'] / g['Mk'], ez_fold=mb['EZ'] / g['EZ'],
                div_GNP=g['div'], div_MB=mb['div'], div_GNPHHi=gh['div'], div_GNPnoSHH=gm['div'])

# ---------- (0) division-gating: at which a0 does the latch open, dividing vs arrested? ----------
def bistable_gap(over, cond, mk0):
    hi = settle_state(over, cond, mk0=0.90); lo = settle_state(over, cond, mk0=0.02)
    if hi is None or lo is None: return np.nan, np.nan
    return abs(hi['Mk'] - lo['Mk']), hi['div'] if hi['div'] is not None else np.nan

def gating_sweep(a0s, n_rw=3):
    div_gap, arr_gap = [], []
    for a0 in a0s:
        ov = dict(a0_prc2=a0, n_rw=n_rw)
        g, _ = bistable_gap(ov, dict(SHH=0.5, **GNP), 0.02)          # DIVIDING (GNP+SHH)
        a, _ = bistable_gap(ov, dict(SHH=0.5, HHi=1, MYCN_amplification=1.0, Ptch1_copy_number=1.0), 0.02)  # ARRESTED
        div_gap.append(g); arr_gap.append(a)
    return dict(a0=a0s, div_gap=np.array(div_gap), arr_gap=np.array(arr_gap))

if __name__ == '__main__':
    res = {}
    print("=== (0) division-gating: latch (hi/lo Mk gap) vs a0, dividing vs arrested ===")
    a0s = np.logspace(-6, -2.2, 16)
    res['gating'] = gating_sweep(a0s)
    for a0, dg, ag in zip(a0s, res['gating']['div_gap'], res['gating']['arr_gap']):
        print(f"  a0={a0:.2e}  dividing gap={dg:.3f}  arrested gap={ag:.3f}")

    print("\n=== (1) withdrawal-memory protocol ===")
    res['wd_baked'] = withdrawal({}); res['wd_latch'] = withdrawal(LATCH)
    for nm, d in (('BAKED', res['wd_baked']), ('LATCH', res['wd_latch'])):
        b = np.cumsum([0] + list(d['segbounds']))[:3] / UPH
        p1 = np.sum(d['div_t'] < b[1]); p2 = np.sum((d['div_t'] >= b[1]) & (d['div_t'] < b[2])); p3 = np.sum(d['div_t'] >= b[2])
        print(f"  {nm}: divisions  proliferate={p1}  withdraw={p2}  RESTORE={p3}   (memory if RESTORE stays 0)")

    print("\n=== (2) SHH hysteresis (continuation) ===")
    shhs = np.linspace(0.02, 0.6, 22)
    res['cont_baked'] = continuation({}, shhs); res['cont_latch'] = continuation(LATCH, shhs)
    for nm, d in (('BAKED', res['cont_baked']), ('LATCH', res['cont_latch'])):
        gap = np.max(np.abs(d['up_mk'] - d['dn_mk']))
        print(f"  {nm}: max |up-down| Mk gap = {gap:.3f}  ({'HYSTERESIS/memory' if gap > 0.1 else 'reversible/no memory'})")

    print("\n=== (3) key validation metrics ===")
    res['val_baked'] = val_metrics({}); res['val_latch'] = val_metrics(LATCH)
    hdr = ['Cd_GNP', 'Cd_MB', 'fold', 'hhi_ratio', 'chip', 'ez_fold', 'div_GNP', 'div_MB', 'div_GNPHHi', 'div_GNPnoSHH']
    tgt = dict(Cd_GNP='~1.4', Cd_MB='~7', fold='5.07', hhi_ratio='0.16', chip='~0.5(<1)', ez_fold='2.05',
               div_GNP='>0', div_MB='>0', div_GNPHHi='0', div_GNPnoSHH='0')
    print(f"  {'metric':13s} {'target':>9} {'BAKED':>9} {'LATCH':>9}")
    for k in hdr:
        print(f"  {k:13s} {tgt[k]:>9} {res['val_baked'][k]:>9.3f} {res['val_latch'][k]:>9.3f}")

    pickle.dump(res, open(os.path.join(HERE, 'sim_latch_consequences_cache.pkl'), 'wb'))
    print("\n-> cached. figure: fig_latch_consequences.py")
