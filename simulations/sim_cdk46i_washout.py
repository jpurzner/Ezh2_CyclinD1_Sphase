"""CDK4/6i washout kinetics in MB: 72h palbociclib -> washout, tracking the pRb-Ser807/811⁺ cycling fraction over
time. Reproduces JP's experiment (1 µM reversible, rebounds; 5 µM saturating, slower to re-enter). Palbo = graded
residual CDK4/6 activity (kPhRbCd=KP*resid): 1 µM -> resid 0.60 (sub-saturating), 5 µM -> resid 0.23 (saturating,
from the K~1.5 µM fit). Time is mapped to hours via the model's calibration (62.8 units/hour; cycle ~22h).

Run:  ./venv/bin/python simulations/sim_cdk46i_washout.py [--n=40] [--workers=8]
"""
import sys, os, argparse, json, contextlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, multiprocessing as mp, tellurium as te
from src.build_model_v44_heldt import build_model_v44

@contextlib.contextmanager
def _quiet():
    fd = sys.stderr.fileno(); s = os.dup(fd); dn = os.open(os.devnull, os.O_WRONLY)
    try: os.dup2(dn, fd); yield
    finally: os.dup2(s, fd); os.close(dn); os.close(s)

UPH = 62.8                              # model units per hour (measured, G2+M~3.3h)
SETTLE = 4000.0; PALBO_H = 72.0; WASH_H = 96.0
PALBO_DUR = PALBO_H * UPH; WASH_DUR = WASH_H * UPH
CHUNK = 40.0; CYCLE_U = 1385.0                # cycling = divided within the last cycle
KTL0, CD_SDLOG, INK4_SDLOG, CIP_SDLOG = 0.801, 0.70, 0.40, 0.40
RESERVE_KSYP21, RESERVE_FRAC = 0.11, 0.20
MB = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002)
DOSES = {'1 µM (sub-saturating)': 0.61, '5 µM (saturating)': 0.23}   # resid mapping (K~1.5 µM)

_RR = None; _KP = None
def _init(_):
    global _RR, _KP
    os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
    _RR = te.loada(build_model_v44()); _KP = float(_RR['kPhRbCd'])

def _seg(rr, t0, t1, npts):
    """Robust single-segment integration; returns (time, Dna) arrays (empty on failure)."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(t0, t1, npts, selections=['time', 'Dna'])
            return np.asarray(r['time']), np.asarray(r['Dna'])
        except Exception:
            continue
    return np.array([]), np.array([])

def _cell(args):
    resid, ktl, p16, p18, kcip, reserve = args
    rr = _RR
    rr.reset(); rr['SHH'] = 0.5
    for k, v in MB.items(): rr[k] = v
    rr['k_Cd_translation'] = ktl; rr['p16'] = p16; rr['p18'] = p18
    rr['kSyP21'] = RESERVE_KSYP21 if reserve else kcip
    t1 = SETTLE; t2 = SETTLE + PALBO_DUR; t3 = SETTLE + PALBO_DUR + WASH_DUR
    npts = lambda a, b: max(2, int((b - a) / CHUNK))
    with _quiet():
        rr['kPhRbCd'] = _KP;         ts, ps = _seg(rr, 0.0, t1, npts(0, t1))      # settle
        rr['kPhRbCd'] = _KP * resid; tp, pp = _seg(rr, t1, t2, npts(t1, t2))      # palbo (single integration)
        rr['kPhRbCd'] = _KP;         tw, pw = _seg(rr, t2, t3, npts(t2, t3))      # washout
    T = np.concatenate([a for a in (ts, tp, tw) if a.size])
    P = np.concatenate([a for a in (ps, pp, pw) if a.size])
    if T.size == 0: return np.array([SETTLE, t3]), np.array([0.0, 0.0])
    return T, P

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=40); ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    rng = np.random.default_rng(20260716); n = a.n
    ktl = KTL0 * np.exp(rng.normal(0, CD_SDLOG, n))
    p16 = MB['p16'] * np.exp(rng.normal(0, INK4_SDLOG, n))
    p18 = MB['p18'] * np.exp(rng.normal(0, INK4_SDLOG, n))
    kcip = 0.002 * np.exp(rng.normal(0, CIP_SDLOG, n))
    is_res = np.arange(n) < int(round(RESERVE_FRAC * n))
    # common time grid in HOURS, zeroed at palbo start
    T_end = SETTLE + PALBO_DUR + WASH_DUR
    grid_u = np.arange(0, T_end, CHUNK)
    tgrid_h = (grid_u - SETTLE) / UPH
    def smooth(y, w=9):
        k = np.ones(w) / w
        return np.convolve(np.pad(y, w // 2, mode='edge'), k, mode='valid')[:len(y)]
    out = {'palbo_h': PALBO_H, 'wash_h': WASH_H, 'settle_h': SETTLE / UPH, 'time_h': tgrid_h.tolist(), 'doses': {}}
    with mp.Pool(a.workers, initializer=_init, initargs=(None,)) as pool:
        for label, resid in DOSES.items():
            args = [(resid, float(ktl[i]), float(p16[i]), float(p18[i]), float(kcip[i]), bool(is_res[i])) for i in range(n)]
            traces = pool.map(_cell, args)
            frac = np.zeros(len(grid_u))
            for T, D in traces:                              # D = Dna(t); cycling = divided within the last cycle
                if T.size < 2: continue
                div_t = T[:-1][(D[:-1] > 0.9) & (D[1:] < 0.1)]   # division times
                if div_t.size == 0: continue
                cyc = np.array([np.any((div_t > gt - CYCLE_U) & (div_t <= gt)) for gt in grid_u])
                frac += cyc.astype(float)
            frac = smooth(100.0 * frac / n)
            out['doses'][label] = frac.tolist()
            base = np.mean(frac[(tgrid_h > -30) & (tgrid_h < 0)])
            plateau = np.mean(frac[(tgrid_h >= PALBO_H * 0.5) & (tgrid_h <= PALBO_H)])   # palbo plateau (2nd half)
            rec_24 = np.interp(PALBO_H + 24, tgrid_h, frac); rec_end = frac[-1]
            print(f"  {label:24}: baseline {base:.0f}% -> palbo plateau {plateau:.0f}% -> washout +24h {rec_24:.0f}% -> +{WASH_H:.0f}h {rec_end:.0f}%")
    with open(os.path.join(os.path.dirname(__file__), 'cdk46i_washout_results.json'), 'w') as f:
        json.dump(out, f)
    print("  -> simulations/cdk46i_washout_results.json")

if __name__ == '__main__':
    main()
