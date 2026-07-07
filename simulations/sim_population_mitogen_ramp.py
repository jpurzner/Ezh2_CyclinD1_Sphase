"""Population mitogen ramp-ON (entry) and ramp-OFF (exit) — WITH vs WITHOUT the H3K27me3 mark repression
of CyclinD1, on the data-grounded heterogeneous population.

Each GNP cell carries the three measured-variation axes (CyclinD1 abundance CV 0.70, birth p27 CV 0.33,
EZH2 abundance CV ~0.55; scale-free). We ramp mitogen (SHH) UP from arrest (entry) and DOWN from cycling
(exit), and record the SHH at which each cell first divides (entry threshold) / last divides (exit
threshold). We do this WITH the mark repression (default, f0_mk=0.233) and WITHOUT it (f0_mk=1.0, the
H3K27me3 mark cannot repress CyclinD1). The mark should GATE the population — raising the Hh needed to
enter and the Hh at which cells exit (buffering: protection against premature entry AND premature exit).

Run:  ./venv/bin/python simulations/sim_population_mitogen_ramp.py [--full] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FULL = '--full' in sys.argv
FRESH = '--fresh' in sys.argv or FULL
N = 120 if FULL else 40
CACHE = 'simulations/sim_population_mitogen_ramp_cache.npz'
KTL0, KTLEZ0 = 0.801, 0.004
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.633, 0.32, 0.51    # DATA CVs (0.70 / 0.33 / ~0.55)
D, PREEQ = 9000.0, 3000.0                             # ramp duration (150 h) + pre-equilibration
CHUNK = 150.0

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass


def setcell(cd, p21d, ez, f0):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    _R['k_Cd_translation'] = float(KTL0 * cd); _R['P21_div'] = float(p21d)
    _R['kTlEZ'] = float(KTLEZ0 * ez); _R['f0_mk'] = float(f0)


def _chunk_sim(t, dt):
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            return _R.simulate(t, t + dt, max(40, int(dt / 4)), selections=['time', 'MPF', 'SHH'])
        except Exception:
            continue
    return None


def ramp(cd, p21d, ez, f0, entry):
    setcell(cd, p21d, ez, f0)
    shh0, shh1 = (0.05, 1.10) if entry else (1.10, 0.03)
    _R['SHH'] = shh0
    t = 0.0
    _chunk_sim(t, PREEQ); t = PREEQ                      # pre-equilibrate (arrested / cycling)
    shh_fn = lambda tt: shh0 + (shh1 - shh0) * min(max((tt - PREEQ) / D, 0.0), 1.0)
    T, MP, SH = [], [], []; T_end = PREEQ + D + 1500
    while t < T_end:
        _R['SHH'] = float(shh_fn(t))
        res = _chunk_sim(t, CHUNK)
        if res is None:
            break
        T.append(res['time'][1:]); MP.append(res['MPF'][1:]); SH.append(res['SHH'][1:]); t += CHUNK
    if not T:
        return np.nan
    T = np.concatenate(T); MP = np.concatenate(MP); SH = np.concatenate(SH)
    pk, _ = find_peaks(MP, prominence=0.15, distance=200)
    pk = pk[T[pk] >= PREEQ]                               # divisions during the ramp only
    if len(pk) == 0:
        return shh1 if entry else shh0                   # never entered / never exited (edge)
    return float(SH[pk[0]]) if entry else float(SH[pk[-1]])   # entry = 1st div SHH; exit = last div SHH


if FRESH or not os.path.exists(CACHE):
    rng = np.random.default_rng(7)
    cd = np.clip(np.exp(rng.normal(0, CD_SDLOG, N)), 0.2, 6.0)
    p21 = np.clip(np.exp(rng.normal(np.log(0.85), P27_SDLOG, N)), 0.4, 2.5)
    ez = np.clip(np.exp(rng.normal(0, EZ_SDLOG, N)), 0.25, 4.0)
    ent_on = np.full(N, np.nan); ent_off = np.full(N, np.nan)
    exi_on = np.full(N, np.nan); exi_off = np.full(N, np.nan)
    for i in range(N):
        ent_on[i] = ramp(cd[i], p21[i], ez[i], 0.233, entry=True)
        ent_off[i] = ramp(cd[i], p21[i], ez[i], 1.000, entry=True)
        exi_on[i] = ramp(cd[i], p21[i], ez[i], 0.233, entry=False)
        exi_off[i] = ramp(cd[i], p21[i], ez[i], 1.000, entry=False)
        print(f"  cell {i+1:3d}/{N}  entry SHH  mark {ent_on[i]:.2f} / no-mark {ent_off[i]:.2f}   exit SHH  mark {exi_on[i]:.2f} / no-mark {exi_off[i]:.2f}", flush=True)
    np.savez(CACHE, ent_on=ent_on, ent_off=ent_off, exi_on=exi_on, exi_off=exi_off, N=N)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
ent_on, ent_off, exi_on, exi_off = z['ent_on'], z['ent_off'], z['exi_on'], z['exi_off']
C_ON, C_OFF = '#7a5aa6', '#e67e22'


def cdf(v, xs):
    v = v[np.isfinite(v)]
    return np.array([np.mean(v <= x) for x in xs])


xs = np.linspace(0.0, 1.1, 200)
fig, ax = plt.subplots(2, 2, figsize=(14, 9))

# (A) ENTRY: cumulative fraction of the population that has entered (÷ rising SHH)
ax[0, 0].plot(xs, cdf(ent_on, xs), color=C_ON, lw=2.6, label='with H3K27me3 repression')
ax[0, 0].plot(xs, cdf(ent_off, xs), color=C_OFF, lw=2.6, label='without (mark cannot repress)')
ax[0, 0].axhline(0.5, color='k', ls=':', lw=0.8, alpha=0.5)
ax[0, 0].set_xlabel('mitogen SHH (rising →)'); ax[0, 0].set_ylabel('fraction of population entered (cycling)')
ax[0, 0].set_title('(A) Ramp-ON entry: the mark raises the Hh\nthreshold to start proliferating', fontweight='bold')
ax[0, 0].legend(fontsize=9, loc='lower right'); ax[0, 0].grid(alpha=0.15)

# (B) EXIT: fraction still cycling at SHH x = fraction whose exit point is below x (CDF of exit-SHH)
ax[0, 1].plot(xs, cdf(exi_on, xs), color=C_ON, lw=2.6, label='with H3K27me3 repression')
ax[0, 1].plot(xs, cdf(exi_off, xs), color=C_OFF, lw=2.6, label='without (mark cannot repress)')
ax[0, 1].axhline(0.5, color='k', ls=':', lw=0.8, alpha=0.5)
ax[0, 1].set_xlabel('mitogen SHH (falling ←)'); ax[0, 1].set_ylabel('fraction of population still cycling')
ax[0, 1].invert_xaxis()
ax[0, 1].set_title('(B) Ramp-OFF exit: the mark makes cells exit\nat higher Hh (earlier on withdrawal)', fontweight='bold')
ax[0, 1].legend(fontsize=9, loc='lower left'); ax[0, 1].grid(alpha=0.15)

# (C) entry-SHH distributions
parts = ax[1, 0].violinplot([ent_on[np.isfinite(ent_on)], ent_off[np.isfinite(ent_off)]], showmedians=True)
for pc, c in zip(parts['bodies'], [C_ON, C_OFF]): pc.set_facecolor(c); pc.set_alpha(0.6)
ax[1, 0].set_xticks([1, 2]); ax[1, 0].set_xticklabels(['with mark', 'without'])
ax[1, 0].set_ylabel('entry SHH threshold'); ax[1, 0].set_title(f'(C) Entry threshold distribution\nmedian {np.nanmedian(ent_on):.2f} vs {np.nanmedian(ent_off):.2f}', fontweight='bold'); ax[1, 0].grid(alpha=0.15, axis='y')

# (D) exit-SHH distributions
parts = ax[1, 1].violinplot([exi_on[np.isfinite(exi_on)], exi_off[np.isfinite(exi_off)]], showmedians=True)
for pc, c in zip(parts['bodies'], [C_ON, C_OFF]): pc.set_facecolor(c); pc.set_alpha(0.6)
ax[1, 1].set_xticks([1, 2]); ax[1, 1].set_xticklabels(['with mark', 'without'])
ax[1, 1].set_ylabel('exit SHH threshold'); ax[1, 1].set_title(f'(D) Exit threshold distribution\nmedian {np.nanmedian(exi_on):.2f} vs {np.nanmedian(exi_off):.2f}', fontweight='bold'); ax[1, 1].grid(alpha=0.15, axis='y')

fig.suptitle(f'Population mitogen ramp — H3K27me3 repression gates cell-cycle ENTRY (↑Hh) and EXIT (↑Hh) (N={len(ent_on)}; 3 measured-variation axes)',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_population_mitogen_ramp.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_population_mitogen_ramp.pdf', bbox_inches='tight')
plt.close()
print(f'\nENTRY median SHH: mark {np.nanmedian(ent_on):.2f}  no-mark {np.nanmedian(ent_off):.2f}')
print(f'EXIT  median SHH: mark {np.nanmedian(exi_on):.2f}  no-mark {np.nanmedian(exi_off):.2f}')
print('Saved sim_population_mitogen_ramp.png')
