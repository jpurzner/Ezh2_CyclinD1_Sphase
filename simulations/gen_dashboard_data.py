"""Generate trajectory data for the interactive single-cell dashboard.

Simulates a handful of GNP and MB lineages (heterogeneous via the calibrated abundance draws), records the key
biomolecules over several divisions, computes cell-cycle phase per timepoint, and exports a compact JSON that the
dashboard HTML embeds. Deterministic lineages (one draw per cell) so each cell is a clean, interpretable trajectory.

Run:  ./venv/bin/python simulations/gen_dashboard_data.py
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

KTL0 = 0.801
T_END, N_PTS = 8400, 1400            # ~140 h, ~6 min/point
SEL = ['time', 'Cd', 'EZH2', 'P21', 'Skp2', 'pRb', 'E2f', 'Ce', 'Ca', 'MPF', 'Dna', 'mass', 'aRc']

# display metadata: model name -> (label, color, group, short description)
SPECIES = [
    ('Cd',   'CyclinD1',       '#16a34a', 'mitogen',    'Mitogen/growth signal (Hh+MYCN); drives commitment.'),
    ('EZH2', 'EZH2',           '#9333ea', 'mitogen',    'PRC2 methyltransferase; represses CyclinD1.'),
    ('P21',  'p27',            '#d97706', 'commit',     'CDK inhibitor. HIGH = quiescent/G0; the fate marker.'),
    ('Skp2', 'Skp2',           '#dc2626', 'commit',     'Degrades p27. Rises at commitment (the R-point switch).'),
    ('pRb',  'phospho-Rb',     '#be123c', 'commit',     'Hyper-phosphorylated Rb = past the restriction point.'),
    ('E2f',  'E2F',            '#1d4ed8', 'commit',     'Proliferation transcription factor; released by pRb.'),
    ('Ce',   'CyclinE',        '#0d9488', 'cyclin',     'CDK2 partner; fires S-phase entry.'),
    ('Ca',   'CyclinA',        '#2563eb', 'cyclin',     'CDK2 partner through S/G2.'),
    ('MPF',  'CyclinB-CDK1',   '#ea580c', 'cyclin',     'Mitotic trigger (MPF). Spikes at mitosis -> division.'),
    ('Dna',  'DNA content',    '#64748b', 'cell',       'Replication 0->1; resets to 0 at division.'),
    ('mass', 'cell mass',      '#92400e', 'cell',       'Cell size; grows, halves at division.'),
    ('aRc',  'replication forks','#0891b2','cell',      'Active replication (S-phase).'),
]
GROUPS = {'mitogen': 'Mitogen / epigenetic', 'commit': 'Commitment (R-point)', 'cyclin': 'Cyclins / CDK', 'cell': 'Cell state'}

# hand-picked cells spanning the behaviour (cd multiplier, birth-p27, label)
GNP_CELLS = [(1.4, 0.5, 'fast'), (1.0, 0.6, 'typical'), (0.7, 0.9, 'hesitant')]
MB_CELLS  = [(1.4, 1.4, 'driven'), (1.0, 1.8, 'typical'), (0.7, 2.4, 'quiescent-prone')]
MB_PARAMS = dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002)


def simulate(cd_mult, p21div, params):
    r = te.loada(build_model_v44())
    r.integrator.setValue('absolute_tolerance', 1e-8); r.integrator.setValue('relative_tolerance', 1e-6)
    r.integrator.setValue('maximum_num_steps', 2000000); r.integrator.setValue('maximum_time_step', 5.0)
    r['SHH'] = 0.5
    for k, v in params.items(): r[k] = v
    r['k_Cd_translation'] = KTL0 * cd_mult
    r['P21_div'] = p21div
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            return r.simulate(0, T_END, N_PTS, selections=SEL)
        except Exception:
            continue
    return r.simulate(0, T_END, N_PTS, selections=SEL)


def phase_of(dna, arc, p21, mpf):
    ph = np.empty(len(dna), dtype=int)   # 0 G0, 1 G1, 2 S, 3 G2, 4 M
    for i in range(len(dna)):
        if mpf[i] > 0.3:                    ph[i] = 4
        elif dna[i] >= 0.9:                 ph[i] = 3
        elif arc[i] > 0.05 and dna[i] < 0.9:ph[i] = 2
        elif p21[i] > 0.1:                  ph[i] = 0
        else:                               ph[i] = 1
    return ph


def build():
    cells = []
    smax = {s[0]: 1e-9 for s in SPECIES}
    for typ, spec in (('GNP', GNP_CELLS), ('MB', MB_CELLS)):
        for cd_mult, p21div, tag in spec:
            d = simulate(cd_mult, p21div, {} if typ == 'GNP' else MB_PARAMS)
            t = d['time'] / 60.0
            dna = d['Dna']
            divs = list(np.round(t[np.where((dna[:-1] > 0.9) & (dna[1:] < 0.5))[0] + 1], 2))
            ph = phase_of(d['Dna'], d['aRc'], d['P21'], d['MPF'])
            traj = {s[0]: [round(float(x), 4) for x in d[s[0]]] for s in SPECIES}
            for s in SPECIES: smax[s[0]] = max(smax[s[0]], max(traj[s[0]]))
            cells.append(dict(id=f"{typ}-{tag}", type=typ, tag=tag, cd_mult=cd_mult, p21div=p21div,
                              traj=traj, phase=[int(x) for x in ph], divisions=divs))
    out = dict(
        time=[round(float(x), 2) for x in (np.arange(N_PTS) * (T_END / (N_PTS - 1)) / 60.0)],
        species=[dict(key=s[0], label=s[1], color=s[2], group=s[3], desc=s[4]) for s in SPECIES],
        groups=GROUPS, smax={k: round(v, 4) for k, v in smax.items()},
        phases=['G0', 'G1', 'S', 'G2', 'M'],
        phase_colors=['#f59e0b', '#93c5fd', '#22c55e', '#fb923c', '#ef4444'],
        cells=cells)
    path = os.path.join(os.path.dirname(__file__), 'dashboard_data.json')
    with open(path, 'w') as f:
        json.dump(out, f, separators=(',', ':'))
    kb = os.path.getsize(path) / 1024
    print(f"wrote {path}  ({kb:.0f} KB, {len(cells)} cells, {N_PTS} pts)")
    for c in cells:
        print(f"  {c['id']:22s} cd_mult={c['cd_mult']} birth-p27={c['p21div']}  divisions={len(c['divisions'])}")


if __name__ == '__main__':
    build()
