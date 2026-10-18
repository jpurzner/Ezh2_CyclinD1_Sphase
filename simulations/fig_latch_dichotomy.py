"""What IS the latch? A didactic MONOSTABLE-vs-BISTABLE figure for the mark equation (JP 2026-07-21).
Same reduced equation as fig_memory_parameter_pairs.py:

    dM/dt = writing(M) - erasure(M),   writing = λ·W0·(a0 + a_rw·M^n_rw)·(1-M),   erasure = δ·M

Left column  = MONOSTABLE (v44 baked: high nucleation a0)         -> no memory.
Right column = BISTABLE   (low a0 + cooperative read-write n_rw=3) -> a LATCH (memory).
Rows: (1) rate balance & fixed points, (2) potential landscape, (3) the mitogen-sweep MEMORY test.
Run: ./venv/bin/python simulations/fig_latch_dichotomy.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
DEL = 0.00386          # baked δ_eff (del_mk + eraser@Gli0.5 + ln2/Tc + h33)
W0 = 18.674            # calibrated writing scale (M*≈0.96 at baked)

MONO = dict(a0=0.0028749753, a_rw=0.0022657, n_rw=2.0)   # v44 baked -> monostable
BIST = dict(a0=2e-6,        a_rw=0.010,     n_rw=3.0)     # low a0 + cooperative -> latch

def writing(M, p, lam=1.0):
    return lam * W0 * (p['a0'] + p['a_rw'] * M ** p['n_rw']) * (1 - M)

def rate(M, p, lam=1.0):
    return writing(M, p, lam) - DEL * M

def stable_fps(p, lam):
    M = np.linspace(1e-5, 1 - 1e-5, 4000)
    f = rate(M, p, lam)
    idx = np.where(np.diff(np.sign(f)) != 0)[0]
    fps = []
    for i in idx:
        m = M[i] - f[i] * (M[i + 1] - M[i]) / (f[i + 1] - f[i])
        fps.append((m, f[i] > 0 and f[i + 1] < 0))     # (M*, stable?)
    return fps

def potential(p, lam=1.0):
    M = np.linspace(0, 1, 800)
    V = -np.concatenate([[0], np.cumsum((rate(M[1:], p, lam) + rate(M[:-1], p, lam)) / 2 * np.diff(M))])
    return M, V

# --- colours ---
C_W, C_E, C_STA, C_UNS = '#1e8449', '#c0392b', '#154360', '#c0392b'
C_UP, C_DN, C_MEM = '#2e86c1', '#e67e22', '#f9e79f'

fig, ax = plt.subplots(3, 2, figsize=(13.6, 14.2))
fig.subplots_adjust(left=0.08, right=0.965, top=0.83, bottom=0.05, hspace=0.42, wspace=0.24)
COLS = [('MONOSTABLE — no memory\n(v44 baked: high nucleation a0)', MONO, ax[:, 0]),
        ('BISTABLE — a LATCH (memory)\n(low a0 + cooperative read-write n_rw=3)', BIST, ax[:, 1])]

Mgrid = np.linspace(0, 1, 500)
for title, p, col in COLS:
    # ---------- ROW 1: rate balance ----------
    a = col[0]
    a.plot(Mgrid, writing(Mgrid, p), color=C_W, lw=2.6, label='writing  λ·W0·(a0 + a_rw·M^n_rw)·(1−M)')
    a.plot(Mgrid, DEL * Mgrid, color=C_E, lw=2.6, label='erasure  δ_eff·M')
    for m, stab in stable_fps(p, 1.0):
        a.plot(m, DEL * m, 'o', ms=13, mfc=(C_STA if stab else 'white'), mec=(C_STA if stab else C_UNS),
               mew=2.2, zorder=5)
        a.annotate('stable' if stab else 'unstable\n(tipping pt)', (m, DEL * m),
                   textcoords='offset points', xytext=(6, 12 if stab else -26), fontsize=8.5,
                   color=(C_STA if stab else C_UNS), fontweight='bold')
    a.set_title(title, fontsize=12, fontweight='bold', pad=10)
    a.set_xlabel('H3K27me3 mark  M'); a.set_ylabel('rate')
    a.set_xlim(0, 1); a.set_ylim(0, DEL * 1.05); a.legend(fontsize=7.8, loc='upper center')
    a.text(0.55, 0.10, 'fixed points = where\nwriting crosses erasure', transform=a.transAxes,
           fontsize=8.5, style='italic', va='bottom', ha='center', color='#555')

    # ---------- ROW 2: potential landscape ----------
    b = col[1]; M, V = potential(p, 1.0)
    b.plot(M, V, color='#5b2c6f', lw=2.8)
    for m, stab in stable_fps(p, 1.0):
        vv = np.interp(m, M, V)
        b.plot(m, vv, 'o', ms=15, mfc=(C_STA if stab else 'white'), mec=(C_STA if stab else C_UNS), mew=2.2, zorder=5)
    b.set_xlabel('H3K27me3 mark  M'); b.set_ylabel('potential  V(M)  =  −∫ dM/dt')
    b.set_xlim(0, 1)
    nwell = sum(1 for _, s in stable_fps(p, 1.0) if s)
    b.text(0.5, 0.93, ('ONE valley — the mark always rolls to the\nsame state (no memory of the past)'
                       if nwell < 2 else 'TWO valleys — the mark STAYS in whichever\nvalley it is in = it REMEMBERS (latch)'),
           transform=b.transAxes, ha='center', va='top', fontsize=9.5, fontweight='bold',
           color=('#7f8c8d' if nwell < 2 else '#5b2c6f'),
           bbox=dict(boxstyle='round', fc='#fdf2e9' if nwell >= 2 else '#f4f6f6', ec='none'))

    # ---------- ROW 3: mitogen sweep (memory test) ----------
    c = col[2]
    lams = np.linspace(0.1, 4.0, 240)
    up, dn = [], []
    m_prev = 0.001
    for lam in lams:                      # UP-sweep (mitogen increasing)
        fps = [m for m, s in stable_fps(p, lam) if s] or [m_prev]
        m_prev = min(fps, key=lambda z: abs(z - m_prev)); up.append(m_prev)
    m_prev = 0.999
    for lam in lams[::-1]:                 # DOWN-sweep (mitogen decreasing)
        fps = [m for m, s in stable_fps(p, lam) if s] or [m_prev]
        m_prev = min(fps, key=lambda z: abs(z - m_prev)); dn.append(m_prev)
    dn = dn[::-1]
    up, dn = np.array(up), np.array(dn)
    c.fill_between(lams, up, dn, where=np.abs(up - dn) > 0.02, color=C_MEM, alpha=0.9, zorder=1, label='memory (hysteresis)')
    c.plot(lams, up, color=C_UP, lw=2.6, label='mitogen ↑ (turning ON)')
    c.plot(lams, dn, color=C_DN, lw=2.6, ls='--', label='mitogen ↓ (turning OFF)')
    c.axvline(1.0, color='#999', lw=1, ls=':')
    c.set_xlabel('mitogen  λ  (∝ SHH · EZH2 writing)'); c.set_ylabel('steady H3K27me3 mark  M')
    c.set_ylim(-0.03, 1.03); c.set_xlim(lams[0], lams[-1]); c.legend(fontsize=8, loc='center right')
    if np.max(np.abs(up - dn)) > 0.05:
        li = int(np.argmax(np.abs(up - dn)))
        c.annotate('same mitogen,\nTWO outcomes\n= MEMORY', (lams[li], (up[li] + dn[li]) / 2),
                   fontsize=9.5, fontweight='bold', ha='center', color='#b9770e',
                   bbox=dict(boxstyle='round', fc='white', ec='#e67e22'))
        c.set_title('mark can be HIGH or LOW at the same mitogen → LATCH', fontsize=9.5, color='#b9770e', fontweight='bold')
    else:
        c.set_title('mark just tracks mitogen (up = down) → no memory', fontsize=9.5, color='#7f8c8d', fontweight='bold')

# ----- explanatory banner -----
fig.text(0.5, 0.985, 'What makes the H3K27me3 mark a MEMORY LATCH?', ha='center', fontsize=14, fontweight='bold')
fig.text(0.5, 0.963,
         'The mark writes MORE of itself where it already exists — read-write  a_rw·M^n_rw  is POSITIVE FEEDBACK.',
         ha='center', fontsize=10.5, fontweight='bold', color='#1e8449')
fig.text(0.5, 0.943,
         'a0 (nucleation): seeds the mark REGARDLESS of state — high a0 always fills the empty state, so no OFF state can exist → monostable.',
         ha='center', fontsize=9.6, color='#444')
fig.text(0.5, 0.927,
         'n_rw (cooperativity): n≥2 makes the feedback switch-like.   The LATCH needs all three: LOW a0  +  strong a_rw  +  n_rw ≥ 2.',
         ha='center', fontsize=9.6, color='#444')
fig.text(0.5, 0.909,
         'v44 baked has a0 ~300× too high → MONOSTABLE (a passive buffer, left column), not a latch (right column).',
         ha='center', fontsize=9.6, style='italic', color='#7f8c8d')

# row labels
for r, lab in enumerate(['1 · RATE BALANCE', '2 · LANDSCAPE', '3 · MITOGEN SWEEP']):
    ax[r, 0].text(-0.16, 0.5, lab, transform=ax[r, 0].transAxes, rotation=90, va='center', ha='center',
                  fontsize=10.5, fontweight='bold', color='#34495e')
    ax[r, 0].grid(alpha=0.15); ax[r, 1].grid(alpha=0.15)

out = os.path.join(HERE, 'fig_latch_dichotomy.png')
fig.savefig(out, dpi=145); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
for title, p, _ in COLS:
    fps = stable_fps(p, 1.0)
    desc = [(round(m, 3), 'S' if s else 'U') for m, s in fps]
    print('  {:24s}: {} stable state(s)  FPs={}'.format(title.splitlines()[0], sum(1 for _, s in fps if s), desc))
