"""Kinetics of EZH2 / PRC2 inhibition on CyclinD1, the role of H3K27me3, and how HISTORICAL model variants differ.

Four experiments, all in an MB cell, tracking the CyclinD1 fold (Cd relative to pre-perturbation baseline, smoothed
over one cycle to remove the cell-cycle oscillation) and the H3K27me3 mark over 72h after a perturbation at t=0:

  A. EZH2i (catalytic: SET-domain block, complex stays) vs EED KO (complex removed) in the CURRENT chain model.
     EZH2i -> gradual + PARTIAL de-repression (mark decays but PRC2 occupancy floor a0 remains); EED KO -> fuller.
  B. How much H3K27me3 is there: sweep the read-write strength a_rw_prc2 (the mark's self-sustaining contribution
     to PRC2), then apply EZH2i -> more mark -> more buffered / slower de-repression.
  C. HISTORICAL vs CURRENT: EZH2i de-repression kinetics under each modeling generation --
     instant (direct-EZH2, with_h3k27_dilution=False), explicit-mark memory (with_h3k27_memory),
     PRC2-occupancy lumped (with_h3k27_chain=False), serial me-chain (current default).

Time -> hours via 62.8 units/hour. Run:  ./venv/bin/python simulations/sim_ezh2_prc2_kinetics.py
"""
import sys, os, json, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

@contextlib.contextmanager
def _quiet():
    fd = sys.stderr.fileno(); s = os.dup(fd); dn = os.open(os.devnull, os.O_WRONLY)
    try: os.dup2(dn, fd); yield
    finally: os.dup2(s, fd); os.close(dn); os.close(s)

UPH = 62.8; SETTLE = 8000.0; TRACK_H = 72.0; TRACK = TRACK_H * UPH
CYCLE_U = 1385.0                         # for cycle-smoothing the CyclinD1 oscillation
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553)

def _smooth_cycle(t, y):
    """moving average over ~one cycle to get the CyclinD1 envelope (de-repression signal, not the oscillation)."""
    dt = np.median(np.diff(t)); w = max(3, int(CYCLE_U / dt) | 1)
    k = np.ones(w) / w
    return np.convolve(np.pad(y, w // 2, mode='edge'), k, mode='valid')[:len(y)]

def run_case(build_kwargs, perturb, mark_species='Mk', overrides=None):
    """Settle MB, record baseline CyclinD1, apply `perturb` dict at t=0, track Cd fold + mark over TRACK."""
    os.environ.setdefault('TWO_STEP_RB', '1')
    m = build_model_v44(**build_kwargs)
    rr = te.loada(m); rr['SHH'] = 0.5
    for k, v in MB.items(): rr[k] = v
    if overrides:
        for k, v in overrides.items():
            try: rr[k] = v
            except Exception: pass
    rr.reset(); rr['SHH'] = 0.5
    for k, v in MB.items(): rr[k] = v
    if overrides:
        for k, v in overrides.items():
            try: rr[k] = v
            except Exception: pass
    has_mark = mark_species in rr.getFloatingSpeciesIds()
    sel = ['time', 'Cd'] + ([mark_species] if has_mark else [])
    with _quiet():
        r0 = rr.simulate(0, SETTLE, 3200, selections=sel)
        base = float(np.mean(np.asarray(r0['Cd'])[-800:]))         # baseline CyclinD1 (last ~1.5 cycles)
        for k, v in perturb.items():
            try: rr[k] = v
            except Exception: pass
        r1 = rr.simulate(SETTLE, SETTLE + TRACK, 2400, selections=sel)
    t = (np.asarray(r1['time']) - SETTLE) / UPH
    cd = _smooth_cycle(np.asarray(r1['time']), np.asarray(r1['Cd'])) / base
    mk = np.asarray(r1[mark_species]) if has_mark else np.zeros_like(t)
    return t.tolist(), cd.tolist(), mk.tolist(), base

CHAIN = dict()                          # current default
out = {'track_h': TRACK_H, 'A': {}, 'B': {}, 'C': {}}

print("A. EZH2i (catalytic) vs EED KO  [current chain model]")
t, cd, mk, b = run_case(CHAIN, perturb=dict(EZH2i=1.0)); out['A']['EZH2i (catalytic)'] = dict(t=t, cd=cd, mk=mk)
print(f"   EZH2i:  CyclinD1 fold @12h {np.interp(12,t,cd):.2f}, @24h {np.interp(24,t,cd):.2f}, @72h {cd[-1]:.2f}")
t, cd, mk, b = run_case(CHAIN, perturb=dict(a0_prc2=0.0, a_rw_prc2=0.0)); out['A']['EED KO (complex removed)'] = dict(t=t, cd=cd, mk=mk)
print(f"   EED KO: CyclinD1 fold @12h {np.interp(12,t,cd):.2f}, @24h {np.interp(24,t,cd):.2f}, @72h {cd[-1]:.2f}")

print("B. mark amount sweep (a_rw_prc2 read-write) then EZH2i")
base_arw = 0.003731831639056823
for scale in (0.3, 1.0, 2.0, 3.0):
    t, cd, mk, b = run_case(CHAIN, perturb=dict(EZH2i=1.0), overrides=dict(a_rw_prc2=base_arw * scale))
    out['B'][f'a_rw x{scale}'] = dict(t=t, cd=cd, mk=mk)
    print(f"   a_rw x{scale}: mark_ss {mk[0]:.3f}, EZH2i fold @24h {np.interp(24,t,cd):.2f}, @72h {cd[-1]:.2f}")

print("C. historical vs current model, EZH2i de-repression kinetics")
VARIANTS = {
    'instant (direct EZH2)':   dict(with_h3k27_dilution=False, with_prc2=False, with_h3k27_memory=False),
    'explicit-mark memory':    dict(with_h3k27_memory=True),
    'PRC2-occupancy (lumped)': dict(with_h3k27_chain=False),
    'serial me-chain (current)': dict(),
}
for name, kw in VARIANTS.items():
    ms = 'H3K27_Cd' if 'memory' in name else 'Mk'
    t, cd, mk, b = run_case(kw, perturb=dict(EZH2i=1.0), mark_species=ms)
    out['C'][name] = dict(t=t, cd=cd, mk=mk)
    print(f"   {name:26}: EZH2i fold @1h {np.interp(1,t,cd):.2f}, @12h {np.interp(12,t,cd):.2f}, @72h {cd[-1]:.2f}")

with open(os.path.join(os.path.dirname(__file__), 'ezh2_prc2_kinetics_results.json'), 'w') as f:
    json.dump(out, f)
print("  -> simulations/ezh2_prc2_kinetics_results.json")
