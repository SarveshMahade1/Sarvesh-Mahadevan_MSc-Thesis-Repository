# -*- coding: utf-8 -*-
"""
Reads code/design_basis.json, the single source of every design input, and
derives the quantities that follow from it. Used by the plain-Python scripts
and by the Abaqus model, so it runs on Python 2.7 and 3.

    import design_basis as DB
    DB.AIRCRAFT                 ['E190', 'A320', 'A350']
    DB.geometry('E190')         R, L, h_cap, V, A, slenderness, ...
    DB.heat_budget('E190')      702.6 W
    DB.MAT['E_Pa'], DB.TH['h_dry_W_per_m2K'], DB.LOAD_CASES, DB.FILLS

To run a study with other inputs without editing the thesis file, copy it and
point the environment variable LH2_DESIGN_BASIS at the copy.
"""
import json
import math
import os
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
PATH = os.environ.get('LH2_DESIGN_BASIS', '').strip() or os.path.join(HERE, 'design_basis.json')
with open(PATH) as _f:
    DATA = json.load(_f, object_pairs_hook=OrderedDict)


def _public(d):
    """The records of a section, without its '_doc' entries."""
    return OrderedDict((k, v) for k, v in d.items() if not k.startswith('_'))


MISSION = _public(DATA['mission'])
MAT = _public(DATA['material'])
CRIT = _public(DATA['criteria'])
PRESS = _public(DATA['pressures'])
TH = _public(DATA['thermal'])
BUDGET = _public(DATA['heat_budget'])
ANALYSIS = _public(DATA['analysis'])
FILLS = list(DATA['fill_levels'])
_AC = _public(DATA['aircraft'])
AIRCRAFT = [k for k in DATA['aircraft'] if not k.startswith('_')]      # file order
LOAD_CASES = dict((k, (v['nx'], v['ny'], v['nz'])) for k, v in _public(DATA['load_cases']).items())
LOAD_CASE_ORDER = [k for k in DATA['load_cases'] if not k.startswith('_')]
LOAD_CASE_LABEL = dict((k, v.get('label', k)) for k, v in _public(DATA['load_cases']).items())

P_DESIGN = PRESS['p_design_gauge_Pa']                              # Pa gauge, strength
P_STABILITY = PRESS['p_fill_abs_Pa'] - PRESS['p_ambient_Pa']       # Pa gauge, stability
SIGMA_ALLOW = MAT['Fty_Pa'] / CRIT['yield_margin']                 # Pa
G = 9.81


def record(ac):
    return _AC[ac]


def label(ac):
    return _AC[ac].get('label', ac)


def mission_geometry(ac):
    """Chapter 3.1 rule: energy equivalence, radius from the fuselage, length from the volume."""
    r = _AC[ac]
    m_lh2 = r['jet_a_mass_kg'] * MISSION['lhv_jet_a_MJ_per_kg'] / MISSION['lhv_lh2_MJ_per_kg']
    V = m_lh2 / (MISSION['rho_lh2_kg_m3'] * MISSION['mission_fill_level'])
    R = round(MISSION['tank_to_fuselage_diameter_ratio'] * r['fuselage_diameter_m'] / 2.0, 4)
    h = R / MISSION['cap_aspect_ratio_a_over_b']
    L = round((V - 4.0 / 3.0 * math.pi * R * R * h) / (math.pi * R * R), 4)
    return dict(m_lh2_kg=m_lh2, V_required_m3=V, R=R, L=L)


def cap_area(R, h):
    """Exact area of the two semi-ellipsoidal caps (one oblate spheroid)."""
    e = math.sqrt(1.0 - (h / R) ** 2)
    return 2.0 * math.pi * R * R * (1.0 + (1.0 - e * e) / e * math.atanh(e))


def geometry(ac):
    """Model geometry: the stored R_mid_m and L_barrel_m, or derived from the mission when null."""
    r = _AC[ac]
    if r.get('R_mid_m') is not None and r.get('L_barrel_m') is not None:
        R, L, src = r['R_mid_m'], r['L_barrel_m'], 'design_basis.json'
    else:
        mg = mission_geometry(ac)
        R, L, src = mg['R'], mg['L'], 'derived from the mission inputs'
    a_b = MISSION['cap_aspect_ratio_a_over_b']
    h = R / a_b
    V = math.pi * R * R * L + 4.0 / 3.0 * math.pi * R * R * h
    A = 2.0 * math.pi * R * L + cap_area(R, h)
    out = dict(R=R, L=L, h_cap=h, cap_aspect=a_b, L_total=L + 2.0 * h, V=V, A=A,
               slenderness=L / (2.0 * R), source=src)
    if r.get('fuselage_length_m'):
        out['share_of_fuselage'] = (L + 2.0 * h) / r['fuselage_length_m']
    return out


def heat_budget(ac):
    """Heat leak budget [W]: the stored Q_budget_W, or the boil-off rule of design_basis.json."""
    q = _AC[ac].get('Q_budget_W')
    if q is not None:
        return q
    V = geometry(ac)['V']
    return round(BUDGET['boiloff_fraction_per_h'] * V * MISSION['rho_lh2_kg_m3']
                 * BUDGET['h_fg_J_per_kg'] / 3600.0, 1)


def sigma_cr(R, t):
    """NASA SP-8007 axial buckling stress of an unpressurised cylinder [Pa]: (sigma_cr, gamma, sigma_cl)."""
    E, nu = MAT['E_Pa'], MAT['nu']
    s_cl = E * t / (R * math.sqrt(3.0 * (1.0 - nu * nu)))
    phi = math.sqrt(R / t) / 16.0
    gam = 1.0 - 0.901 * (1.0 - math.exp(-phi))
    return gam * s_cl, gam, s_cl


if __name__ == '__main__':
    print('Design basis: %s' % PATH)
    for ac in AIRCRAFT:
        g = geometry(ac)
        print('%-9s R %.4f m  L %.4f m  V %.2f m3  A %.1f m2  L/2R %.2f  Q %.1f W  (%s)'
              % (label(ac), g['R'], g['L'], g['V'], g['A'], g['slenderness'], heat_budget(ac), g['source']))
    print('Load cases: %s' % ', '.join('%s %s' % (k, LOAD_CASES[k]) for k in LOAD_CASE_ORDER))
    print('Fill levels: %s' % FILLS)
    print('Allowable: %.1f MPa; stability pressure difference %.5f bar' % (SIGMA_ALLOW / 1e6, P_STABILITY / 1e5))
