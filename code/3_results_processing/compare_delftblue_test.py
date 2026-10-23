# -*- coding: utf-8 -*-
"""
SOLVER CHECK: DelftBlue cluster (Abaqus 2024) against the workstation
(Abaqus 2025), one design point.

The test case is the A350 Baseline condition at a fill level of 0.10, the
smallest production model (126,816 elements at the sized wall of 5.674 mm).
It exists in the workstation results, so the cluster run can be compared
number by number before the cluster is used for the remaining stages. The
recorded outcome is in results/solver_check_delftblue/.

ON DELFTBLUE (after the upload of the thesis folder):
    mkdir -p /scratch/$USER/dbtest/code/1_abaqus_model
    cp /scratch/$USER/thesis/code/1_abaqus_model/* /scratch/$USER/dbtest/code/1_abaqus_model/
    cp /scratch/$USER/thesis/sized_walls.csv /scratch/$USER/dbtest/
    cd /scratch/$USER/dbtest
    sbatch --job-name=dbtest --time=01:00:00 \
      --export=ALL,WORKDIR=/scratch/$USER/dbtest,AC=A350,STAGE=production,E190_ONLY_LC=Baseline,E190_ONLY_FILL=0.1 \
      code/1_abaqus_model/run_final_delftblue.sh

ON THE WORKSTATION, when it has finished (new Command Prompt):
    cd <local working folder>
    mkdir delftblue_test
    scp <netid>@login.delftblue.tudelft.nl:/scratch/<netid>/dbtest/A350_production*_results.csv delftblue_test\\
    abaqus python code/3_results_processing/compare_delftblue_test.py

Standard library only.
"""
import csv
import glob
import os

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
KEYS = [('n_el', 'elements'), ('T_min', 'T min [K]'), ('T_max', 'T max [K]'),
        ('dT', 'dT [K]'), ('hfl_max_Wm2', 'peak heat flux [W/m2]'),
        ('mises_MPa', 'far-field von Mises [MPa]'), ('fos', 'yield fos'),
        ('u_max_mm', 'max displacement [mm]'), ('sumRF_N', 'total reaction [N]'),
        ('mean_S11_mem_MPa', 'mean axial membrane [MPa]'),
        ('mean_S22_mem_MPa', 'mean hoop membrane [MPa]'),
        ('probe_S22_sneg_MPa', 'crown probe hoop [MPa]'),
        ('RF_aft_kN', 'aft ring [kN]'), ('RF_mid_kN', 'mid ring [kN]'),
        ('RF_fwd_kN', 'fwd ring [kN]'), ('ovality_mm', 'ovality [mm]')]


def pick(paths):
    for p in paths:
        for r in csv.DictReader(open(p)):
            try:
                if (r['load_case'] == 'Baseline' and abs(float(r['fill']) - 0.1) < 1e-6
                        and abs(float(r['t_wall_mm']) - 5.674) < 1e-6
                        and abs(float(r['p_gauge_bar']) - 1.5) < 1e-4):
                    return r
            except (KeyError, ValueError):
                pass
    return None


lap = pick(sorted(glob.glob(os.path.join(RP.PRODUCTION, 'A350_production*_results.csv'))))
dbl = pick(sorted(glob.glob(os.path.join(RP.DELFTBLUE, 'A350_production*_results.csv'))))
if lap is None or dbl is None:
    raise SystemExit('missing: %s' % ('laptop row' if lap is None else 'DelftBlue row'))

print('%-28s %14s %14s %9s' % ('quantity', 'laptop 2025', 'DelftBlue 2024', 'diff %'))
worst = 0.0
for k, name in KEYS:
    a, b = float(lap[k]), float(dbl[k])
    d = 100.0 * (b - a) / a if a else 0.0
    if k != 'ovality_mm':
        worst = max(worst, abs(d))
    print('%-28s %14.4f %14.4f %+9.3f' % (name, a, b, d))
print('\nlargest difference (excluding ovality, a small local quantity): %.3f %%' % worst)
if dbl['n_el'].split('.')[0] == lap['n_el'].split('.')[0] and worst < 0.5:
    print('SAME MESH AND SAME RESULTS: DelftBlue can run the remaining stages.')
elif worst < 2.0:
    print('Mesh or results differ slightly: acceptable for the check stages, state it in the report.')
else:
    print('RESULTS DIFFER: do not use DelftBlue for the thesis stages.')
