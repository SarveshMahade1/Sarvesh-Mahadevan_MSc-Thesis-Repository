# -*- coding: utf-8 -*-
"""
================================================================================
OUTER VACUUM JACKET: REQUIRED THICKNESS UNDER EXTERNAL PRESSURE
TU Delft MSc thesis - Sarvesh Mahadevan
================================================================================

Backs the jacket table of Chapter 3 (Section "Outer Vacuum Jacket
Configuration") and entry 35 of the calculation appendix. The jacket holds the
vacuum, so it carries one atmosphere of external pressure; it is not part of
the finite-element model.

Unstiffened jacket: long thin cylinder, elastic collapse (Bresse)
    p_cr = E / (4 (1 - nu^2)) * (t / R)^3
Ring-stiffened jacket: Windenburg and Trilling (1934), ring spacing L
    p_cr = 2.42 E (t/D)^(5/2) / ( (1 - nu^2)^(3/4) * (L/D - 0.45 (t/D)^(1/2)) )
Both solved for t with p_cr / FOS = 1 atm, FOS 1.5, aluminium E = 73.1 GPa.

Jacket radius = inner-vessel mid-surface radius + half the inner wall + the
evacuated MLI gap of that aircraft. The sized inner walls are read from
results/sized_walls.csv and the MLI thicknesses from
results/processed/insulation.json (written by insulation_sizing.py), so the
jacket follows the sized vessel and its own blanket.

RUN (from anywhere):   python code/2_requirements_and_sizing/outer_jacket_closed_form.py
================================================================================
"""
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402

E, NU, P_EXT, FOS = 73.1e9, 0.33, 101325.0, 1.5
R_MID = {'E190': 1.3015, 'A320': 1.7080, 'A350': 2.5771}
T_WALL = {}
for r in csv.DictReader(open(RP.SIZED_WALLS)):
    T_WALL[r['aircraft']] = float(r['t_wall_mm']) / 1000.0
_INS = json.load(open(os.path.join(RP.PROCESSED, 'insulation.json')))
MLI = {ac: _INS[ac]['mli_mm'] / 1000.0 for ac in R_MID}
SPACINGS = (0.25, 0.5, 1.0, 2.0)


def r_out(ac):
    return R_MID[ac] + T_WALL[ac] / 2.0 + MLI[ac]


def t_unstiffened(ac):
    return r_out(ac) * (4.0 * (1.0 - NU * NU) * P_EXT * FOS / E) ** (1.0 / 3.0)


def t_stiffened(ac, spacing):
    D = 2.0 * r_out(ac)
    lo, hi = 1e-5, 0.1
    for _ in range(200):
        t = 0.5 * (lo + hi)
        pc = 2.42 * E * (t / D) ** 2.5 / ((1.0 - NU * NU) ** 0.75 *
                                          (spacing / D - 0.45 * (t / D) ** 0.5))
        if pc / FOS >= P_EXT:
            hi = t
        else:
            lo = t
    return hi


if __name__ == '__main__':
    print('External pressure %.0f Pa, FOS %.1f, E %.1f GPa\n' % (P_EXT, FOS, E / 1e9))
    print('%-5s %9s %13s  %s' % ('', 'R_out m', 'unstiffened', 'ring-stiffened at spacing ' +
                                  ', '.join('%.2f' % s for s in SPACINGS) + ' m'))
    for ac in ('E190', 'A320', 'A350'):
        print('%-5s %9.4f %10.1f mm  %s' % (ac, r_out(ac), t_unstiffened(ac) * 1e3,
              ', '.join('%.2f mm' % (t_stiffened(ac, s) * 1e3) for s in SPACINGS)))
