# -*- coding: utf-8 -*-
"""
REVIEW CHECKS (October 2026)
Closed-form and volume-sampling checks added in response to the supervisors'
review of the thesis report. No Abaqus input is needed; every number printed
here is reproduced in the revised report or in the response to the comments.

  1. Liquid depth and hydrostatic head for a pure 9 g forward load, to state
     the severity of the combined LC7 case as a stress ratio rather than a
     load-factor ratio.
  2. Split of the change in head ratio p_h/p between load cases into a
     magnitude factor |n| and a depth factor d (orientation), so that the
     report does not attribute to direction alone what is a product of both.
  3. Cryogenic correction of the thermal quantities (NIST 6061-T6 surrogate,
     values read from hand_calcs.json): thermal stress bound, fin length,
     waterline flux and the bowing force.
  4. One reserve-factor convention for both checks: RF = allowable / applied,
     with the project margin of 1.5 inside both allowables.
  5. The 24 h no-vent hold: allowable heat leak and the factor by which every
     thermal result of the linear model would scale.
  6. The hold time actually delivered by the adopted heat leak when the tank
     is sealed at the mission fill of 0.90 instead of the ullage limit of
     0.97 at which the budget was derived (CoolProp, same lever rule as
     boiloff_requirement_derivation.py).
  7. Is the far-field peak pressure-proportional? A straight-line fit of the
     far-field von Mises peak of the twenty production cases against the
     local-pressure membrane value X = (1 + p_h/p) pR/t. A component that
     does not scale with pressure (thermal, for instance) appears as the
     intercept b; the slope a measures the pressure-proportional part.
  8. The hydrostatic design pressure of the fuel-tank rule 14 CFR 25.963(d),
     P = K rho g L, for a tank within the fuselage: K = 9 forward and 6
     downward, L the distance to the farthest tank boundary in the loading
     direction (total length forward, diameter downward), tank full. These are
     ultimate pressures; they are compared with the LC7 head at fill 0.90.
  9. Benchmark without circularity. The flux of the original benchmark run
     (8.0 W/m2) came from Bagarello's computed outer temperature of 282.2 K.
     The 1-D network of bagarello_benchmark_1D_check.py derives the flux from
     their two PRESCRIBED temperatures (283 K, 20 K) instead. The thermal
     model is linear in the flux, so the peak liner temperature of the FE run
     scales to that flux exactly: T = T_min + (T_max - T_min) q_1D / q_run.
 10. Reserve factors at the three mesh levels of the convergence study
     (LC7, fill 0.90), from final_numbers.json.
 11. Derived figures quoted in the report text: fill sensitivity divided by
     the head ratio at fill 0.90, fall of dT with fill at rest, error of the
     bowing estimate against the model, lowest buckling RF per load case,
     yield RF range per load case.

The free surface is placed exactly as in E190_MASTER_thermo_structural.py:
points spread uniformly through the tank volume are sorted by depth along
g_hat, and the surface is the depth that leaves the fill fraction below it.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402

RHO, G, P = 70.8, 9.81, 1.5e5
E, NU, K_AL, H_DRY = 73.1e9, 0.33, 116.0, 1.46
D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
HC = json.load(open(os.path.join(RP.PROCESSED, 'hand_calcs.json')))
AC = ['E190', 'A320', 'A350']
rng = np.random.default_rng(1)


def samples(R, L, h, n=600000):
    """Uniform points in a cylinder (x in 0..L) closed by two half-spheroids."""
    pts = []
    xmin, xmax = -h, L + h
    while sum(len(p) for p in pts) < n:
        m = 400000
        x = rng.uniform(xmin, xmax, m)
        y = rng.uniform(-R, R, m)
        z = rng.uniform(-R, R, m)
        r2 = y * y + z * z
        ok = np.where(x < 0, (x / h) ** 2 + r2 / R ** 2 <= 1,
                      np.where(x > L, ((x - L) / h) ** 2 + r2 / R ** 2 <= 1, r2 <= R * R))
        pts.append(np.stack([x[ok], y[ok], z[ok]], 1))
    return np.concatenate(pts)[:n]


def depth(pts, n, fill):
    nx, ny, nz = n
    nn = math.sqrt(nx * nx + ny * ny + nz * nz)
    g = np.array([nx, -nz, -ny]) / nn          # axes: +X fwd, +Y up, +Z lateral
    s = pts @ g
    s_surf = np.quantile(s, 1.0 - fill)        # liquid occupies the deepest fill
    return float(s.max() - s_surf), nn


out = {}
CASES = {'Baseline': (0, 0, 1), 'LC1': (0, 0, 2.5), 'LC3': (1.5, 0, 2.5),
         'LC7': (9, 3, 6), 'Fwd9only': (9, 0, 0)}
for ac in AC:
    d = D[ac]
    R, L, h = d['R_m'], d['L_barrel_m'], d['h_cap_m']
    t = d['t_wall_mm'] / 1e3
    pts = samples(R, L, h)
    o = {'R_m': R, 'L_barrel_m': L, 'L_total_m': L + 2 * h, 'slenderness': L / (2 * R)}
    for name, n in CASES.items():
        dep, nn = depth(pts, n, 0.90)
        ph = RHO * nn * G * dep
        o[name] = {'n_mag': nn, 'depth_m': dep, 'head_kPa': ph / 1e3, 'head_over_p': ph / P,
                   'hoop_deepest_MPa': (P + ph) * R / t / 1e6}
    # 1. severity of LC7 against the worst single CS-25.561 direction
    o['LC7_vs_Fwd9'] = {'head_ratio': o['LC7']['head_kPa'] / o['Fwd9only']['head_kPa'],
                        'local_pressure_ratio': (P + o['LC7']['head_kPa'] * 1e3) /
                                                (P + o['Fwd9only']['head_kPa'] * 1e3),
                        'load_factor_ratio': o['LC7']['n_mag'] / 9.0}
    # 2. magnitude x depth split of the head ratio
    def split(a, b):
        return {'head_ratio_factor': o[b]['head_over_p'] / o[a]['head_over_p'],
                'magnitude_factor': o[b]['n_mag'] / o[a]['n_mag'],
                'depth_factor': o[b]['depth_m'] / o[a]['depth_m']}
    o['split_Baseline_to_LC1'] = split('Baseline', 'LC1')
    o['split_LC1_to_LC3'] = split('LC1', 'LC3')
    o['split_LC1_to_LC7'] = split('LC1', 'LC7')
    o['split_Baseline_to_LC7'] = split('Baseline', 'LC7')
    # 3. cryogenic correction (surrogate NIST 6061-T6 data in hand_calcs.json)
    cp = HC[ac]['cryogenic_properties']
    th = HC[ac]['thermal']
    a_sec = cp['alpha_secant_6061_per_K']
    dT = th['dT_asymptote_K']
    Lf_c = cp['fin_length_cryo_m']
    o['cryo'] = {'alpha_RT_over_cryo': cp['alpha_ratio_RT_over_cryo'],
                 'kappa_RT_over_cryo': cp['kappa_ratio'],
                 'thermal_stress_bound_RT_MPa': E * 22.3e-6 * dT / (1 - NU) / 1e6,
                 'thermal_stress_bound_cryo_MPa': E * a_sec * dT / (1 - NU) / 1e6,
                 'fin_length_RT_m': th['fin_length_m'], 'fin_length_cryo_m': Lf_c,
                 'waterline_flux_RT_Wm2': th['wall_flux_est_Wm2'],
                 'waterline_flux_cryo_Wm2': th['q_Wm2'] * Lf_c / t,
                 'bowing_force_RT_kN': HC[ac]['thermal_bowing']['F_thermal_kN'],
                 'bowing_force_cryo_kN': cp['F_thermal_cryo_kN'],
                 'mid_ring_est_RT_kN': HC[ac]['thermal_bowing']['F_mid_est_kN'],
                 'mid_ring_est_cryo_kN': cp['F_mid_est_cryo_kN'],
                 'mid_ring_uniform_beam_kN': HC[ac]['thermal_bowing']['F_mid_mechanical_kN']}
    out[ac] = o

# 4. reserve factors, production and finest mesh (convergence table values)
conv = {'E190': (256.34, 256.78, -5.95, -5.90), 'A320': (256.46, 257.76, None, None),
        'A350': (146.02, 147.48, -20.85, -21.65)}
for ac in AC:
    scr = D[ac]['sigma_cr_MPa']
    pv, fv, pa, fa = conv[ac]
    out[ac]['RF'] = {'yield_production': 262.0 / D[ac]['production']['LC7_Emergency|0.90']['mises_MPa'],
                     'yield_finest': 262.0 / fv,
                     'buckling_production': (scr / 1.5) / -pa if pa else None,
                     'buckling_finest': (scr / 1.5) / -fa if fa else None}
out['A350']['RF']['buckling_coarsest'] = (D['A350']['sigma_cr_MPa'] / 1.5) / 21.83

# 5. 24 h no-vent hold (same bulk-energy rise as the 14.68 h derivation)
for ac in AC:
    Q = D[ac]['Q_W']
    out[ac]['hold24'] = {'Q_14p7h_W': Q, 'Q_24h_W': Q * 14.68 / 24.0, 'scale_factor': 14.68 / 24.0,
                         'dT_24h_K': HC[ac]['thermal']['dT_asymptote_K'] * 14.68 / 24.0}

# 6. hold time at the mission fill of 0.90 with the adopted heat leak
try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    '2_requirements_and_sizing'))
    import boiloff_requirement_derivation as BR    # noqa: E402
    for ac in AC:
        d = D[ac]
        V = BR.tank_volume(d['R_m'], d['L_barrel_m'], 1.6)
        res = {}
        for f in (0.97, 0.90):
            m_tot, _, _, rho_b = BR.lockup_state(V, f, BR.P_FILL)
            u1, _ = BR.internal_energy_at(rho_b, BR.P_FILL)
            u2, _ = BR.internal_energy_at(rho_b, BR.P_VENT)
            res[f] = m_tot * (u2 - u1)                      # J to reach venting
        Q = D[ac]['Q_W']
        out[ac]['hold_at_fill'] = {'t_hold_f097_h': res[0.97] / Q / 3600.0,
                                   't_hold_f090_h': res[0.90] / Q / 3600.0,
                                   'Q_for_14p68h_at_f090_W': res[0.90] / (14.68 * 3600.0)}
except ImportError as exc:                                  # CoolProp missing
    print('hold-at-fill check skipped:', exc)

# 7. pressure-proportional and constant parts of the far-field peak
for ac in AC:
    d = D[ac]
    s_th = d['hoop_pR_t_design_MPa']
    X, Y = [], []
    for k, v in d['production'].items():
        X.append((1 + d['heads'][k]['head_over_p']) * s_th)
        Y.append(v['mises_MPa'])
    X, Y = np.array(X), np.array(Y)
    a, b = np.polyfit(X, Y, 1)
    res = Y - (a * X + b)
    out[ac]['peak_fit'] = {'slope_a': float(a), 'intercept_b_MPa': float(b),
                           'max_abs_residual_MPa': float(np.abs(res).max()),
                           'ratio_to_membrane_vm_0866': float(a / 0.866),
                           'X_range_MPa': [float(X.min()), float(X.max())], 'n_cases': int(len(X))}

# 8. 14 CFR 25.963(d) hydrostatic design pressure, tank inside the fuselage
for ac in AC:
    o = out[ac]
    Lt, R = o['L_total_m'], o['R_m']
    fwd = 9.0 * RHO * G * Lt
    down = 6.0 * RHO * G * 2 * R
    o['far25_963d'] = {'P_forward_K9_kPa': fwd / 1e3, 'P_down_K6_kPa': down / 1e3,
                       'LC7_head_f090_kPa': o['LC7']['head_kPa'],
                       'forward_over_LC7_head': fwd / (o['LC7']['head_kPa'] * 1e3),
                       'local_pressure_ratio_forward_over_LC7': (P + fwd) / (P + o['LC7']['head_kPa'] * 1e3),
                       # stress follows the local pressure (check 7): proportional estimate
                       'RF_y_estimate_under_K9': o['RF']['yield_production'] /
                       ((P + fwd) / (P + o['LC7']['head_kPa'] * 1e3))}

# 9. benchmark peak liner temperature at the independently derived flux
import contextlib, io, csv  # noqa: E402
with contextlib.redirect_stdout(io.StringIO()):
    import bagarello_benchmark_1D_check as BB    # noqa: E402  (runs the 1-D network)
q1d = BB.Q / BB.A_in_total
row = next(csv.DictReader(open(os.path.join(RP.RESULTS, 'benchmark_bagarello',
                                            'BAGARELLO_production_Baseline_f75_results.csv'))))
Tmin, Tmax, q_run = float(row['T_min']), float(row['T_max']), 8.0
out['benchmark'] = {'q_run_Wm2': q_run, 'q_1D_Wm2': q1d, 'T_max_run_K': Tmax,
                    'T_max_at_q1D_K': Tmin + (Tmax - Tmin) * q1d / q_run,
                    'T_max_paper_K': 25.3,
                    'diff_at_q1D_pct': 100 * (Tmin + (Tmax - Tmin) * q1d / q_run - 25.3) / 25.3,
                    'BOR_1D_kg_h': BB.BOR, 'BOR_paper_kg_h': 12.03}

# 10. reserve factors at every mesh level of the convergence study
for ac in AC:
    c = D[ac]['convergence']
    scr = D[ac]['sigma_cr_MPa']
    rows = []
    for rd, rm in zip(c['pdesign'], c['pmin']):
        ax_min = rm['axial_mem_min_MPa']
        rows.append({'seed_mm': rd['seed_mm'], 'n_el': rd['n_el'], 'mises_MPa': rd['mises_MPa'],
                     'RF_y': 262.0 / rd['mises_MPa'], 'axial_min_pmin_MPa': ax_min,
                     'RF_b': (scr / 1.5) / -ax_min if ax_min < 0 else None})
    out[ac]['convergence_RF'] = rows

# 11. derived figures quoted in the text
for ac in AC:
    d = D[ac]
    sens = d['fill_sensitivity_pct']
    ratio = {lc: sens[lc] / (100 * d['heads'][lc + '|0.90']['head_over_p']) for lc in sens}
    pb = d['production']
    dT_fall = 100 * (pb['Baseline|0.10']['dT'] - pb['Baseline|0.90']['dT']) / pb['Baseline|0.10']['dT']
    tb = HC[ac]['thermal_bowing']
    bow_err = abs(tb['F_mid_est_kN'] - tb['F_mid_fe_kN']) / abs(tb['F_mid_fe_kN'])
    rfb = {}
    rfy = {}
    for lc in ('Baseline', 'LC1_Maneuver', 'LC3_Combined', 'LC7_Emergency'):
        vals = [(d['sigma_cr_MPa'] / 1.5) / -v['axial_mem_min_MPa'] for k, v in d['stability'].items()
                if k.startswith(lc) and v['axial_mem_min_MPa'] < 0]
        rfb[lc] = min(vals) if vals else None
        ys_ = [262.0 / v['mises_MPa'] for k, v in pb.items() if k.startswith(lc)]
        rfy[lc] = [min(ys_), max(ys_)]
    out[ac]['derived'] = {'sensitivity_over_head_ratio': ratio, 'dT_fall_rest_pct': dT_fall,
                          'bowing_estimate_rel_error': bow_err, 'min_RF_b_by_case': rfb,
                          'RF_y_range_by_case': rfy}

with open(os.path.join(RP.PROCESSED, 'review_checks.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for ac in AC:
    print('\n=====', ac)
    for k, v in out[ac].items():
        print(' ', k, json.dumps(v) if isinstance(v, dict) else v)
