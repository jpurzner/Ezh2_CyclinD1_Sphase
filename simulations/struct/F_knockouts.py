"""STRUCTURE ANALYSIS — step F: causal knockout map.

Zero each loop/module knob (via VALIDATE_PARAMS in a clean subprocess -> no reset() param-leak) and record which
of the 28 validated behaviors flip pass->fail, plus the CyclinD1 fold, MB G0, and division counts. Builds a
KO x behavior breakage matrix -> the load-bearing ranking + which modules are separable. Also perturbs (2x, not 0)
K_EZH2_repression to confirm it is DYNAMICALLY DEAD in the chain default.

Run: ./venv/bin/python simulations/struct/F_knockouts.py
"""
import os, sys, json, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VALIDATE = os.path.join(ROOT, 'simulations', 'validate_v44.py')
OUT = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix='ko_')

# KO knob -> (param, value, loop/module it tests)
KOS = [
    ('baseline',        None,               None,  'baseline'),
    ('E2f synth (L1)',  'kSyE2f',           0.0,   'Rb-E2f'),
    ('E2f autoact(L1)', 'kSyE2fE2f',        0.0,   'Rb-E2f'),
    ('CyclinE (L1)',    'kSyCe',            0.0,   'CDK2'),
    ('CyclinA',         'kSyCa',            0.0,   'CDK2'),
    ('Skp2 (L2)',       'kSySkp2',          0.0,   'Skp2-commit'),
    ('p27 clear (L2)',  'kDeP21Cd',         0.0,   'Skp2-commit'),
    ('p27 synth',       'kSyP21',           0.0,   'CDK2-p27'),
    ('CyclinD1 (mito)', 'k_Cd_translation', 0.0,   'CyclinD1'),
    ('CycD-CDK4/6 Rb',  'kPhRbCd',          0.0,   'commitment'),
    ('EZH2 fb (N2)',    'kTlEZ',            0.0,   'EZH2'),
    ('EZH2 cyc-synth',  'kEZE2f',           0.0,   'EZH2'),
    ('Gli eraser',      'k_jmjd3_gli',      0.0,   'H3K27'),
    ('mark rd-wr (L5)', 'a_rw_prc2',        0.0,   'H3K27'),
    ('mark nucleate',   'a0_prc2',          0.0,   'H3K27'),
    ('S-phase (Dna)',   'kSyDna',           0.0,   'Replication'),
    ('Cdc20 (N1)',      'kaCdc20',          0.0,   'MPF'),
    ('growth (N5)',     'mu',               0.0,   'Mass'),
    ('EZH2-direct*2',   'K_EZH2_repression', 0.9694, 'DEAD-test'),  # 2x; chain default -> expect NO change
]


def run(name, param, val):
    rp = os.path.join(TMP, f'{name.replace(" ","_").replace("/","")}.json')
    env = dict(os.environ, TWO_STEP_RB='1', H3K27_CHAIN='1', VALIDATE_RESULTS=rp)
    if param is not None:
        env['VALIDATE_PARAMS'] = json.dumps({param: val})
    try:
        subprocess.run([PY, VALIDATE], env=env, cwd=ROOT, capture_output=True, timeout=600)
        return json.load(open(rp))
    except Exception as e:
        print(f"  {name}: ERR {repr(e)[:80]}")
        return None


base = run('baseline', None, None)
base_checks = {c['name']: c['pass'] for c in base['checks']}
base_ck = {c['name']: c for c in base['checks']}
names = list(base_checks.keys())
print(f"baseline: {sum(base_checks.values())}/{len(names)} pass\n")

rows = []
for name, param, val, mod in KOS:
    if name == 'baseline':
        continue
    out = run(name, param, val)
    if out is None:
        continue
    ck = {c['name']: c for c in out['checks']}
    flipped = [n for n in names if base_checks.get(n) and not ck.get(n, {}).get('pass', True)]
    fold = float(ck.get('CyclinD1 MB/GNP', {}).get('actual', float('nan')))
    g0 = float(out.get('phase_dmso_duration', {}).get('G0', float('nan')))
    passed = int(out.get('passed', 0))
    ndiv = out.get('divisions', {})
    rows.append(dict(name=name, param=param, module=mod, passed=passed, n_broken=len(flipped),
                     broken=flipped, fold=fold, mb_g0=g0))
    print(f"  {name:18} [{param}={val}] -> {passed}/28  broke {len(flipped):2d}  fold {fold:.2f}  MBG0 {g0:.1f}%"
          + ("  <- BREAKS: " + ", ".join(b[:22] for b in flipped[:4]) if flipped else ""))

rows.sort(key=lambda r: r['n_broken'], reverse=True)
print("\n=== LOAD-BEARING RANKING (most behaviors broken first) ===")
for r in rows:
    print(f"  {r['n_broken']:2d} broken  {r['name']:18} ({r['module']})")

# dead-wire verdict
dead = [r for r in rows if r['name'] == 'EZH2-direct*2']
if dead:
    print(f"\nDEAD-WIRE TEST: doubling K_EZH2_repression broke {dead[0]['n_broken']} behaviors "
          + ("(CONFIRMED DEAD in chain default)" if dead[0]['n_broken'] == 0 else "(has effect)"))

# ---- heatmap: KO x behavior breakage ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
kos = [r for r in rows if r['name'] != 'EZH2-direct*2']
M = np.zeros((len(kos), len(names)))
for i, r in enumerate(kos):
    for n in r['broken']:
        M[i, names.index(n)] = 1
fig, ax = plt.subplots(figsize=(16, 7))
ax.imshow(M, cmap='Reds', aspect='auto', vmin=0, vmax=1)
ax.set_yticks(range(len(kos))); ax.set_yticklabels([f"{r['name']} ({r['module']})" for r in kos], fontsize=8)
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=90, fontsize=6)
ax.set_title('Knockout x behavior breakage — red = this KO breaks that validated behavior', fontsize=12)
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'F_knockouts.png'), dpi=140, bbox_inches='tight')
print("wrote F_knockouts.png")

json.dump(dict(baseline_pass=sum(base_checks.values()), n_checks=len(names), rows=rows),
          open(os.path.join(OUT, 'F_knockouts.json'), 'w'), indent=2)
print("wrote F_knockouts.json")
