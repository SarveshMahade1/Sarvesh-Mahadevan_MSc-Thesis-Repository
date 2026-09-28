# -*- coding: utf-8 -*-
"""
FINAL NUMBERS FOR THE THESIS REPORT (28 Sep 2026)
=================================================

Every number quoted in the report text, tables and graphs is taken from the
file this script writes, results/processed/final_numbers.json, and from the
tables beside it. Nothing is typed by hand.

Inputs (all written by E190_MASTER_thermo_structural.py):
    sized_walls.csv                                  FE sizing record
    results/wall_sizing_iterations/<AC>/<AC>_sizing_results.csv           sizing iterations
    <AC>_production*_results.csv                     1.5 bar, 20 cases
    <AC>_stability*_results.csv                      0.187 bar, 20 cases
    <AC>_convergence_{pdesign,pmin}_results.csv      3 mesh levels each

Only rows at the recorded sized wall and at the stage pressure are accepted.
Geometry and the free-surface depth are recomputed here from the same
AIRCRAFT_DATA and the same formulas as the MASTER script.

Run with any Python 3 (standard library only):
    python python_scripts/final_numbers.py
"""
import csv
import glob
import json
import math
import os

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
OUT = RP.PROCESSED
if not os.path.isdir(OUT):
    os.makedirs(OUT)

AC = ['E190', 'A320', 'A350']
LCS = ['Baseline', 'LC1_Maneuver', 'LC3_Combined', 'LC7_Emergency']
FILLS = [0.1, 0.3, 0.5, 0.75, 0.9]
P_DES, P_MIN = 1.5, 0.18675

# ---- identical to the MASTER script -------------------------------------
GEO = {'E190': dict(R=1.3015, L=13.9675, Q=702.6),
       'A320': dict(R=1.7080, L=10.0120, Q=919.3),
       'A350': dict(R=2.5771, L=33.8240, Q=6583.2)}
CAP_ASPECT = 1.6
E, NU, RHO_AL = 73.1e9, 0.33, 2850.0
FTY, SF = 393e6, 1.5
RHO_LH2, G = 70.8, 9.81
LOADS = {'Baseline': (0.0, 0.0, 1.0), 'LC1_Maneuver': (0.0, 0.0, 2.5),
         'LC3_Combined': (1.5, 0.0, 2.5), 'LC7_Emergency': (9.0, 3.0, 6.0)}
H_WET, H_DRY, T_LH2 = 1000.0, 1.46, 20.3
K_AL = 116.0


def sigma_cr(R, t):
    s_cl = E * t / (R * math.sqrt(3.0 * (1.0 - NU * NU)))
    phi = math.sqrt(R / t) / 16.0
    gam = 1.0 - 0.901 * (1.0 - math.exp(-phi))
    return gam * s_cl, gam, s_cl


def cap_area(R, h):
    e = math.sqrt(1.0 - (h / R) ** 2)
    return 2.0 * math.pi * R * R * (1.0 + (1.0 - e * e) / e * math.atanh(e))


def rows(pattern):
    out = []
    for f in sorted(g for d in (RP.PRODUCTION, RP.STABILITY, RP.CONVERGENCE)
                    for g in glob.glob(os.path.join(d, pattern))):
        for r in csv.DictReader(open(f)):
            out.append(r)
    return out


def at(rr, t, p):
    keep = {}
    for r in rr:
        try:
            if abs(float(r['t_wall_mm']) - t) < 1e-6 and abs(float(r['p_gauge_bar']) - p) < 1e-4:
                keep[(r['load_case'], round(float(r['fill']), 2))] = r
        except (KeyError, ValueError):
            pass
    return keep


def f(r, k):
    v = r[k]
    return float('inf') if v in ('inf', 'Infinity') else float(v)


def depth_at_fill(R, L, h, ghat, fill, n=90):
    """
    Liquid depth along ghat below the free surface, for a given fill, by
    direct integration over the tank volume (cylinder plus two half
    spheroids). Same relative-equilibrium physics as the MASTER script:
    the free surface is the plane normal to ghat that encloses fill*V.
    Returns (depth at fill, full extent of the tank along ghat).
    """
    gx, gy, gz = ghat
    pts = []
    # sample points in a bounding box, keep those inside the tank
    nx, nr = 3 * n, n
    xs = [(-h) + (L + 2 * h) * (i + 0.5) / nx for i in range(nx)]
    ys = [-R + 2 * R * (j + 0.5) / nr for j in range(nr)]
    for x in xs:
        if x < 0.0:
            rloc = R * math.sqrt(max(0.0, 1.0 - (x / h) ** 2))
        elif x > L:
            rloc = R * math.sqrt(max(0.0, 1.0 - ((x - L) / h) ** 2))
        else:
            rloc = R
        for y in ys:
            for z in ys:
                if y * y + z * z <= rloc * rloc:
                    pts.append(x * gx + y * gy + z * gz)
    pts.sort()
    k = int(round((1.0 - fill) * len(pts)))
    s_surf = pts[min(len(pts) - 1, k)]
    # the extremes are taken analytically (support function of the capsule),
    # because a sampled grid misses the outermost half cell
    ux = abs(gx)
    half = 0.5 * L * ux + math.sqrt(R * R * (1.0 - gx * gx) + h * h * gx * gx)
    ctr = 0.5 * L * gx
    hi, lo = ctr + half, ctr - half
    return hi - s_surf, hi - lo


def unit(n):
    # ghat components in model axes: X forward, Y up, Z lateral; the liquid
    # settles along the inertial force (+X forward, -Y down, lateral).
    nx, ny, nz = n
    m = math.sqrt(nx * nx + ny * ny + nz * nz)
    return (nx / m, -nz / m, -ny / m), m


D = {'inputs': {'E_Pa': E, 'nu': NU, 'rho_Al': RHO_AL, 'Fty_Pa': FTY, 'SF': SF,
                'sigma_allow_MPa': FTY / SF / 1e6, 'rho_LH2': RHO_LH2,
                'h_wet': H_WET, 'h_dry': H_DRY, 'T_LH2': T_LH2,
                'p_design_bar': P_DES, 'p_min_bar': P_MIN}}
walls = dict((r['aircraft'], r) for r in csv.DictReader(open(RP.SIZED_WALLS)))

for ac in AC:
    g = GEO[ac]
    R, L = g['R'], g['L']
    t = float(walls[ac]['t_wall_mm'])
    tm = t / 1000.0
    h = R / CAP_ASPECT
    A = 2 * math.pi * R * L + cap_area(R, h)
    V = math.pi * R * R * L + 4.0 / 3.0 * math.pi * R * R * h
    scr, gam, scl = sigma_cr(R, tm)
    bl = math.sqrt(R * tm)
    d = {'R_m': R, 'L_barrel_m': L, 'h_cap_m': h, 'L_total_m': L + 2 * h,
         'slenderness_L_2R': L / (2 * R), 't_wall_mm': t,
         'governed_by': walls[ac]['governed_by'], 'R_over_t': R / tm,
         'area_m2': A, 'volume_m3': V, 'Q_W': g['Q'], 'q_Wm2': g['Q'] / A,
         'dT_q_over_hdry_K': g['Q'] / A / H_DRY,
         'fin_length_m': math.sqrt(K_AL * tm / H_DRY),
         'shell_mass_kg': RHO_AL * A * tm,
         'LH2_mass_full_kg': RHO_LH2 * V,
         'boundary_layer_mm': bl * 1e3, 'seed_mm': 0.5 * bl * 1e3,
         'sigma_cl_MPa': scl / 1e6, 'knockdown_gamma': gam, 'sigma_cr_MPa': scr / 1e6,
         'hoop_pR_t_design_MPa': P_DES * 1e5 * R / tm / 1e6,
         'axial_pR_2t_design_MPa': P_DES * 1e5 * R / (2 * tm) / 1e6,
         'axial_pR_2t_min_MPa': P_MIN * 1e5 * R / (2 * tm) / 1e6,
         'sizing_record': walls[ac]}

    # free-surface depth and hydrostatic head, per load case and fill
    heads = {}
    for lc in LCS:
        gh, nmag = unit(LOADS[lc])
        tilt = math.degrees(math.acos(min(1.0, abs(gh[1]))))
        for fi in FILLS:
            dep, ext = depth_at_fill(R, L, h, gh, fi, n=60)
            heads['%s|%.2f' % (lc, fi)] = {
                'tilt_deg': tilt, 'n_mag': nmag, 'depth_m': dep, 'extent_m': ext,
                'head_bar': RHO_LH2 * nmag * G * dep / 1e5,
                'head_over_p': RHO_LH2 * nmag * G * dep / (P_DES * 1e5)}
    d['heads'] = heads

    prod = at(rows('%s_production*_results.csv' % ac), t, P_DES)
    stab = at(rows('%s_stability*_results.csv' % ac), t, P_MIN)
    assert len(prod) == 20 and len(stab) == 20, (ac, len(prod), len(stab))
    keys = ['n_el', 'T_min', 'T_max', 'dT', 'hfl_max_Wm2', 'mises_MPa', 'mises_singular_MPa',
            'fos', 'u_max_mm', 'sumRF_N', 'balance_ratio', 'axial_mem_min_MPa',
            'mean_S11_mem_MPa', 'mean_S22_mem_MPa', 'ovality_mm', 'ur_keel_mm',
            'ur_flank_mm', 'RF_aft_kN', 'RF_mid_kN', 'RF_fwd_kN', 'probe_S11_sneg_MPa',
            'probe_S11_spos_MPa', 'probe_S22_sneg_MPa', 'probe_S22_spos_MPa',
            'cap_S11_mem_min_MPa', 'cap_S22_mem_min_MPa', 'buckling_fos']
    d['production'] = dict(('%s|%.2f' % k, dict((kk, f(v, kk)) for kk in keys)) for k, v in prod.items())
    d['stability'] = dict(('%s|%.2f' % k, dict((kk, f(v, kk)) for kk in keys)) for k, v in stab.items())

    # fill sensitivity of the far-field stress, per load case
    sens = {}
    for lc in LCS:
        a, b = prod[(lc, 0.1)], prod[(lc, 0.9)]
        sens[lc] = 100.0 * (f(b, 'mises_MPa') - f(a, 'mises_MPa')) / f(a, 'mises_MPa')
    d['fill_sensitivity_pct'] = sens
    gov = max(prod.items(), key=lambda kv: f(kv[1], 'mises_MPa'))
    d['governing_yield'] = {'case': '%s|%.2f' % gov[0], 'mises_MPa': f(gov[1], 'mises_MPa'),
                            'fos': f(gov[1], 'fos')}
    gb = min(stab.items(), key=lambda kv: f(kv[1], 'buckling_fos'))
    d['governing_buckling'] = {'case': '%s|%.2f' % gb[0], 'buckling_fos': f(gb[1], 'buckling_fos'),
                               'axial_mem_min_MPa': f(gb[1], 'axial_mem_min_MPa')}
    d['min_axial_stab'] = min(f(v, 'axial_mem_min_MPa') for v in stab.values())
    d['worst_balance_pct'] = 100 * max(abs(f(v, 'balance_ratio') - 1) for v in list(prod.values()) + list(stab.values()))

    conv = {}
    for lab in ('pdesign', 'pmin'):
        rr = rows('%s_convergence_%s_results.csv' % (ac, lab))
        rr = [r for r in rr if abs(float(r['t_wall_mm']) - t) < 1e-6]
        rr.sort(key=lambda r: -float(r['seed_mm']))
        conv[lab] = [dict((kk, f(r, kk)) for kk in ['seed_mm'] + keys) for r in rr]
    d['convergence'] = conv
    fin, prd = conv['pdesign'][-1], conv['pdesign'][1]
    d['mesh_error_yield_pct'] = 100 * (prd['mises_MPa'] - fin['mises_MPa']) / fin['mises_MPa']
    d['fos_yield_finest'] = FTY / SF / 1e6 / fin['mises_MPa']
    d['bfos_finest'] = conv['pmin'][-1]['buckling_fos']

    # older sizing files carry a comma inside the stage label, which shifts
    # every later column by one: rejoin the first two fields for those rows
    siz = []
    _p = os.path.join(RP.SIZING, ac, '%s_sizing_results.csv' % ac)
    _rd = csv.reader(open(_p))
    _hd = next(_rd)
    for _r in _rd:
        if len(_r) == len(_hd) + 1:
            _r = [_r[0] + ';' + _r[1]] + _r[2:]
        siz.append(dict(zip(_hd, _r)))
    d['sizing_iterations'] = [dict(stage=r['stage'], t_wall_mm=f(r, 't_wall_mm'),
                                   p_gauge_bar=f(r, 'p_gauge_bar'), mises_MPa=f(r, 'mises_MPa'),
                                   fos=f(r, 'fos'), buckling_fos=f(r, 'buckling_fos'),
                                   axial_mem_min_MPa=f(r, 'axial_mem_min_MPa'),
                                   load_case=r['load_case'], fill=f(r, 'fill')) for r in siz]
    D[ac] = d

with open(os.path.join(OUT, 'final_numbers.json'), 'w') as fh:
    json.dump(D, fh, indent=1, default=str)

# flat tables for the appendix and for checking
with open(os.path.join(OUT, 'production_matrix_sized_walls.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['aircraft', 'load_case', 'fill', 'n_el', 'T_max_K', 'dT_K', 'hfl_max_Wm2',
                'vM_farfield_MPa', 'fos_yield', 'u_max_mm', 'reaction_kN', 'RF_aft_kN',
                'RF_mid_kN', 'RF_fwd_kN', 'ovality_mm', 'balance_ratio'])
    for ac in AC:
        for lc in LCS:
            for fi in FILLS:
                r = D[ac]['production']['%s|%.2f' % (lc, fi)]
                w.writerow([ac, lc, fi, int(r['n_el']), '%.3f' % r['T_max'], '%.3f' % r['dT'],
                            '%.1f' % r['hfl_max_Wm2'], '%.2f' % r['mises_MPa'], '%.3f' % r['fos'],
                            '%.2f' % r['u_max_mm'], '%.1f' % (r['sumRF_N'] / 1e3),
                            '%.1f' % r['RF_aft_kN'], '%.1f' % r['RF_mid_kN'],
                            '%.1f' % r['RF_fwd_kN'], '%.2f' % r['ovality_mm'],
                            '%.4f' % r['balance_ratio']])

for ac in AC:
    d = D[ac]
    print('\n%s  t=%.3f mm (%s)  L/2R=%.2f  R/t=%.0f  q=%.2f W/m2  seed=%.1f mm  sigma_cr=%.2f MPa'
          % (ac, d['t_wall_mm'], d['governed_by'], d['slenderness_L_2R'], d['R_over_t'],
             d['q_Wm2'], d['seed_mm'], d['sigma_cr_MPa']))
    print('  governing yield %s  vM %.2f  fos %.3f | buckling %s fos %.3f'
          % (d['governing_yield']['case'], d['governing_yield']['mises_MPa'],
             d['governing_yield']['fos'], d['governing_buckling']['case'],
             d['governing_buckling']['buckling_fos']))
    print('  fill sensitivity %s' % ', '.join('%s %+.1f%%' % (k, v) for k, v in d['fill_sensitivity_pct'].items()))
    print('  mesh error yield %+.2f%%, fos finest %.3f, bfos finest %s, worst balance %.2f%%'
          % (d['mesh_error_yield_pct'], d['fos_yield_finest'], d['bfos_finest'], d['worst_balance_pct']))
    for k in ('Baseline|0.90', 'LC7_Emergency|0.90', 'LC3_Combined|0.90'):
        hh = d['heads'][k]
        print('  %s tilt %.1f depth %.2f m extent %.2f m head %.3f bar (%.0f%% of p)'
              % (k, hh['tilt_deg'], hh['depth_m'], hh['extent_m'], hh['head_bar'], 100 * hh['head_over_p']))
print('\nwritten', os.path.join(OUT, 'final_numbers.json'))
