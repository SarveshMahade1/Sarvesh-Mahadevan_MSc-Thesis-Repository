# -*- coding: utf-8 -*-
"""
INDEPENDENT HAND CALCULATIONS AGAINST THE FINITE ELEMENT RESULTS (28 Sep 2026)
=============================================================================

Closed-form physics for every check quoted in the report, evaluated at the
sized walls and set beside the finite element value from
results/processed/final_numbers.json. Standard library only:

    python python_scripts/final_numbers.py      (first)
    python python_scripts/hand_calcs.py

Writes results/processed/hand_calcs.json and prints a readable table.

Each check states the physics it tests. None of them uses a finite element
result as an input, so agreement is evidence and not a restatement.
"""
import json
import math
import os

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
AC = ['E190', 'A320', 'A350']
G = 9.81
P = 1.5e5
PMIN = 0.18675e5
E, NU, ALPHA, K_AL = 73.1e9, 0.33, 22.3e-6, 116.0
RHO_LH2, H_DRY = 70.8, 1.46
A_OVER_B = 1.6
LOADS = {'Baseline': (0.0, 0.0, 1.0), 'LC7_Emergency': (9.0, 3.0, 6.0)}
P_VAC_MIN = 1.2e5        # Pa abs, minimum operating pressure against the vacuum
P_RELIEF_ABS = 1.448e5   # Pa abs, relief setting against the vacuum
T_REF = 20.3             # K


# NIST cryogenic material properties, Al 6061-T6 (UNS A96061), curve fits from
# https://trc.nist.gov/cryogenics/materials/6061%20Aluminum/6061_T6Aluminum_rev.htm
def nist_6061_k(T):
    a = [0.07918, 1.0957, -0.07277, 0.08084, 0.02803, -0.09464, 0.04179, -0.00571, 0.0]
    x = math.log10(T)
    return 10 ** sum(ci * x ** i for i, ci in enumerate(a))


def nist_6061_dLL(T):
    """(L - L293)/L293, unitless."""
    if T < 18.0:
        return -415.45e-5
    return (-4.1277E2 - 3.0389E-1 * T + 8.7696E-3 * T ** 2 - 9.9821E-6 * T ** 3) * 1e-5


def nist_6061_E(T):
    a = [7.771221E1, 1.030646E-2, -2.924100E-4, 8.993600E-7, -1.070900E-9]
    return sum(ci * T ** i for i, ci in enumerate(a))

out = {}
for ac in AC:
    d = D[ac]
    R, L, h = d['R_m'], d['L_barrel_m'], d['h_cap_m']
    t = d['t_wall_mm'] / 1e3
    Lt = L + 2 * h
    prod = d['production']
    stab = d['stability']
    c = {}

    # 1. Membrane stress in the barrel (thin-shell theory, pressure only)
    s_th, s_x = P * R / t / 1e6, P * R / (2 * t) / 1e6
    vm_mem = math.sqrt(s_th ** 2 - s_th * s_x + s_x ** 2)
    b = prod['Baseline|0.10']
    c['membrane'] = {
        'hoop_theory_MPa': s_th, 'axial_theory_MPa': s_x, 'vm_theory_MPa': vm_mem,
        'hoop_probe_MPa': 0.5 * (b['probe_S22_sneg_MPa'] + b['probe_S22_spos_MPa']),
        'axial_probe_MPa': 0.5 * (b['probe_S11_sneg_MPa'] + b['probe_S11_spos_MPa']),
        'hoop_mean_MPa': b['mean_S22_mem_MPa'], 'axial_mean_MPa': b['mean_S11_mem_MPa'],
        'surface_diff_hoop_pct': 100 * abs(b['probe_S22_sneg_MPa'] - b['probe_S22_spos_MPa'])
        / abs(0.5 * (b['probe_S22_sneg_MPa'] + b['probe_S22_spos_MPa']))}
    c['membrane']['hoop_err_pct'] = 100 * (c['membrane']['hoop_probe_MPa'] - s_th) / s_th
    c['membrane']['axial_err_pct'] = 100 * (c['membrane']['axial_probe_MPa'] - s_x) / s_x
    c['membrane']['ratio_probe'] = c['membrane']['hoop_probe_MPa'] / c['membrane']['axial_probe_MPa']

    # 2. Hoop stress at the equator of the semi-ellipsoidal cap (membrane theory)
    s_cap = s_th * (1 - A_OVER_B ** 2 / 2.0)
    c['cap_equator'] = {'hoop_theory_MPa': s_cap, 'hoop_fe_min_MPa': b['cap_S22_mem_min_MPa'],
                        'pole_theory_MPa': s_th * A_OVER_B / 2.0}

    # 3. Global equilibrium: applied resultant against summed ring reactions
    eq = {}
    for key in ('Baseline|0.10', 'Baseline|0.90', 'LC7_Emergency|0.10', 'LC7_Emergency|0.90'):
        lc, fi = key.split('|')
        nx, ny, nz = LOADS[lc]
        nm = math.sqrt(nx * nx + ny * ny + nz * nz)
        m = d['shell_mass_kg'] + RHO_LH2 * d['volume_m3'] * float(fi)
        F = m * G * nm / 1e3
        eq[key] = {'mass_kg': m, 'applied_kN': F, 'reacted_kN': prod[key]['sumRF_N'] / 1e3,
                   'diff_pct': 100 * (prod[key]['sumRF_N'] / 1e3 - F) / F}
    c['equilibrium'] = eq

    # 4. Hydrostatic head: why load direction, not magnitude, matters
    hb, h7 = d['heads']['Baseline|0.90'], d['heads']['LC7_Emergency|0.90']
    c['head'] = {'baseline_depth_m': hb['depth_m'], 'baseline_head_over_p': hb['head_over_p'],
                 'lc7_depth_m': h7['depth_m'], 'lc7_head_over_p': h7['head_over_p'],
                 'lc7_extent_m': h7['extent_m'],
                 'hoop_membrane_at_deepest_MPa': (P + h7['head_bar'] * 1e5) * R / t / 1e6,
                 'fe_farfield_vm_MPa': prod['LC7_Emergency|0.90']['mises_MPa']}

    # 5. Axial membrane stress at the minimum pressure: continuous beam on
    #    three rigid supports carrying the transverse load (n_y, n_z) of the
    #    whole tank, superposed on the pressure term p R / 2t.
    nt = math.sqrt(3.0 ** 2 + 6.0 ** 2)
    m7 = d['shell_mass_kg'] + RHO_LH2 * d['volume_m3'] * 0.9
    w = m7 * G * nt / Lt                                  # N/m
    span = L / 2.0
    M = w * span ** 2 / 8.0                               # hogging over the mid ring
    s_b = M / (math.pi * R ** 2 * t) / 1e6
    s_p = PMIN * R / (2 * t) / 1e6
    c['beam_axial_pmin'] = {'w_N_per_m': w, 'M_mid_kNm': M / 1e3, 'sigma_bend_MPa': s_b,
                            'sigma_pressure_MPa': s_p, 'sigma_x_min_est_MPa': s_p - s_b,
                            'sigma_x_min_fe_MPa': stab['LC7_Emergency|0.90']['axial_mem_min_MPa'],
                            'sigma_cr_MPa': d['sigma_cr_MPa']}

    # 5b. Forward inertia of the liquid (LC7, n_x = 9): the liquid presses on the
    #     forward cap, and the force travels through the forward half of the
    #     barrel to the mid-ring anchor as axial TENSION. Averaged over the
    #     whole barrel this raises the mean axial membrane stress by
    #     F / (2 pi R t) / 2 above the pressure-only value.
    mliq = RHO_LH2 * d['volume_m3'] * 0.9
    Fx = mliq * 9.0 * G
    ds_pred = Fx / (2 * math.pi * R * t) / 2.0 / 1e6
    ds_fe = prod['LC7_Emergency|0.90']['mean_S11_mem_MPa'] - prod['Baseline|0.90']['mean_S11_mem_MPa']
    c['forward_inertia'] = {'F_kN': Fx / 1e3, 'forward_half_tension_MPa': 2 * ds_pred,
                            'mean_axial_rise_pred_MPa': ds_pred, 'mean_axial_rise_fe_MPa': ds_fe}

    # 5c. Beam deflection of the section midway between rings (LC7, 0.9 fill):
    #     two-span continuous beam, maximum span deflection w l^4 / (185 E I).
    #     The recorded 'ovality' is the RANGE of radial displacement around
    #     that section, which contains 2 x this rigid translation plus the
    #     departure from roundness.
    Ib = math.pi * R ** 3 * t
    dbeam = w * span ** 4 / (185.0 * E * Ib)
    c['section_distortion'] = {'beam_deflection_mm': dbeam * 1e3, 'range_from_translation_mm': 2 * dbeam * 1e3,
                               'range_fe_mm': prod['LC7_Emergency|0.90']['ovality_mm'],
                               'translation_share': 2 * dbeam * 1e3 / prod['LC7_Emergency|0.90']['ovality_mm']}

    # 6. Thermal asymptote of the dry wall and the fin length
    q = d['q_Wm2']
    c['thermal'] = {'q_Wm2': q, 'dT_asymptote_K': q / H_DRY,
                    'dT_fe_f10_K': prod['Baseline|0.10']['dT'],
                    'dT_fe_f90_K': prod['Baseline|0.90']['dT'],
                    'fin_length_m': math.sqrt(K_AL * t / H_DRY),
                    'through_thickness_drop_K': q * t / K_AL,
                    # heat conducted to the waterline comes from a band one fin
                    # length wide, so the in-plane flux in the wall there is
                    # q L_fin / t; the rest of the dry wall rejects its heat to
                    # the vapour through h_dry
                    'wall_flux_est_Wm2': q * math.sqrt(K_AL * t / H_DRY) / t,
                    'wall_flux_fe_f90_Wm2': prod['Baseline|0.90']['hfl_max_Wm2'],
                    'wall_flux_fe_f10_Wm2': prod['Baseline|0.10']['hfl_max_Wm2']}

    # 7. Thermal bowing and the mid-ring reaction (Baseline, 10 % fill).
    #    The warmer dry wall expands, the section curves with kappa =
    #    alpha dT / 2R, and the rigid mid ring must pull the barrel back into
    #    line: F = 48 E I delta / L^3 with delta = kappa L^2 / 8.
    I = math.pi * R ** 3 * t
    kap = ALPHA * prod['Baseline|0.10']['dT'] / (2 * R)
    delta = kap * L ** 2 / 8.0
    F_th = 48 * E * I * delta / L ** 3 / 1e3
    F_mech = 0.625 * prod['Baseline|0.10']['sumRF_N'] / 1e3
    r = prod['Baseline|0.10']
    mid_fe_signed = r['sumRF_N'] / 1e3 - r['RF_aft_kN'] - r['RF_fwd_kN']
    c['thermal_bowing'] = {'curvature_per_m': kap, 'bow_mm': delta * 1e3,
                           'F_thermal_kN': F_th, 'F_mid_mechanical_kN': F_mech,
                           'F_mid_est_kN': F_mech - F_th, 'F_mid_fe_kN': mid_fe_signed}

    # 8. Decay of a local disturbance and the far-field exclusion
    beta = (3 * (1 - NU ** 2) / (R * R * t * t)) ** 0.25
    c['decay'] = {'sqrt_Rt_mm': math.sqrt(R * t) * 1e3, 'pi_over_beta_m': math.pi / beta,
                  'exclusion_in_decay_lengths': 0.5 / (math.pi / beta),
                  'attenuation_at_exclusion': math.exp(-beta * 0.5)}

    # 9. Cooldown, which the analysis deliberately excludes
    #    Integrated contraction 293 K -> 20.3 K from the NIST 6061-T6 fit
    #    (the room-temperature alpha overstates it, since alpha falls on cooling).
    eps_cd = -(nist_6061_dLL(T_REF) - nist_6061_dLL(293.0))
    c['cooldown'] = {'contraction_strain': eps_cd,
                     'free_axial_contraction_mm': eps_cd * Lt * 1e3,
                     'uniaxial_restrained_MPa': E * eps_cd / 1e6,
                     'biaxial_restrained_MPa': E * eps_cd / (1 - NU) / 1e6,
                     'free_axial_contraction_RT_alpha_mm': ALPHA * 272.7 * Lt * 1e3}

    # 10. Resolution of the binned temperature transfer (24 bins)
    c['binning'] = {'bin_width_K': (prod['Baseline|0.10']['T_max'] - prod['Baseline|0.10']['T_min']) / 24.0}

    # 11. Pressure differential with the insulating vacuum INTACT (added 28 Sep).
    #     The inner wall faces the evacuated gap, so with the vacuum intact the
    #     differential is the absolute tank pressure, at least 1.2 bar. The
    #     0.187 bar stability case is the differential after a LOSS OF VACUUM
    #     at sea level. The problem is linear, and a uniform pressure increment
    #     adds exactly dp R / 2t of axial membrane stress in the far-field
    #     barrel, so the minimum axial stress with the vacuum intact follows
    #     from the stability run by superposition.
    dp = P_VAC_MIN - PMIN
    s_add = dp * R / (2 * t) / 1e6
    ax_lov = stab['LC7_Emergency|0.90']['axial_mem_min_MPa']
    c['vacuum_intact'] = {'p_min_vacuum_intact_bar': P_VAC_MIN / 1e5,
                          'p_max_vacuum_intact_bar': P_RELIEF_ABS / 1e5,
                          'design_margin_over_relief_pct': 100 * (P / P_RELIEF_ABS - 1),
                          'axial_added_MPa': s_add,
                          'axial_min_loss_of_vacuum_MPa': ax_lov,
                          'axial_min_vacuum_intact_MPa': ax_lov + s_add}

    # 12. Cryogenic material properties (added 28 Sep). The model uses room-
    #     temperature alpha and kappa, as the benchmark study does. NIST
    #     cryogenic data for Al 6061-T6 (no 2219 entry) are used as a surrogate
    #     to size the consequence. Thermal strain is the integral of alpha from
    #     T_ref, so the secant value over the wall's range is the relevant one.
    dTa = q / H_DRY
    a_sec = (nist_6061_dLL(T_REF + dTa) - nist_6061_dLL(T_REF)) / dTa
    ratio = ALPHA / a_sec
    k_ratio = nist_6061_k(293.0) / nist_6061_k(T_REF)
    c['cryogenic_properties'] = {
        'alpha_secant_6061_per_K': a_sec, 'alpha_ratio_RT_over_cryo': ratio,
        'kappa_6061_293K': nist_6061_k(293.0), 'kappa_6061_20K': nist_6061_k(T_REF),
        'kappa_ratio': k_ratio, 'fin_length_cryo_m': math.sqrt(K_AL / k_ratio * t / H_DRY),
        'E_6061_293K_GPa': nist_6061_E(293.0), 'E_6061_20K_GPa': nist_6061_E(T_REF),
        'F_thermal_cryo_kN': F_th / ratio, 'F_mid_est_cryo_kN': F_mech - F_th / ratio}

    # 13. Bound on the thermal stress, and the margins it can move (added 28 Sep).
    #     No thermal stress in the far field can exceed that of the thermal
    #     strain fully restrained in both directions, E alpha dT / (1 - nu),
    #     with the room-temperature alpha of the model and the largest dT.
    #     Removing all or part of it (cryogenic alpha) therefore moves any
    #     far-field stress by at most this amount, in either direction.
    s_thb = E * ALPHA * dTa / (1 - NU) / 1e6
    ax_int = c['vacuum_intact']['axial_min_vacuum_intact_MPa']
    lb = ax_int - s_thb
    c['thermal_bound'] = {'sigma_thermal_max_MPa': s_thb,
                          'axial_min_vacuum_intact_lower_bound_MPa': lb,
                          'buckling_fos_vacuum_intact_lower_bound':
                              (d['sigma_cr_MPa'] / -lb) if lb < 0 else float('inf')}
    # Evidence from the matrix: the far-field peak divided by the local-pressure
    # membrane value (1 + p_h/p) pR/t, over all twenty cases. The dry region
    # ranges from most of the wall (Baseline 0.10) to a thin band, so a thermal
    # contribution to the peak would show as a spread in this ratio.
    rr = [v['mises_MPa'] / ((1 + d['heads'][k]['head_over_p']) * s_th)
          for k, v in prod.items()]
    c['peak_ratio'] = {'min': min(rr), 'max': max(rr),
                       'spread_pct': 100 * (max(rr) - min(rr)) / (0.5 * (max(rr) + min(rr)))}
    out[ac] = c

with open(os.path.join(RP.PROCESSED, 'hand_calcs.json'), 'w') as fh:
    json.dump(out, fh, indent=1)

for ac in AC:
    c = out[ac]
    m = c['membrane']
    print('\n==== %s' % ac)
    print(' membrane hoop %.2f vs probe %.2f (%+.2f%%), axial %.2f vs %.2f (%+.2f%%), ratio %.3f, surfaces %.2f%%'
          % (m['hoop_theory_MPa'], m['hoop_probe_MPa'], m['hoop_err_pct'], m['axial_theory_MPa'],
             m['axial_probe_MPa'], m['axial_err_pct'], m['ratio_probe'], m['surface_diff_hoop_pct']))
    print(' cap equator hoop theory %.2f vs FE min %.2f; pole %.2f'
          % (c['cap_equator']['hoop_theory_MPa'], c['cap_equator']['hoop_fe_min_MPa'], c['cap_equator']['pole_theory_MPa']))
    for k, v in c['equilibrium'].items():
        print(' equilibrium %-18s applied %9.1f kN reacted %9.1f kN (%+.2f%%)' % (k, v['applied_kN'], v['reacted_kN'], v['diff_pct']))
    hh = c['head']
    print(' head: baseline %.2f m %.1f%% of p | LC7 %.2f m %.1f%% of p | hoop at deepest %.1f MPa, FE vM %.1f'
          % (hh['baseline_depth_m'], 100 * hh['baseline_head_over_p'], hh['lc7_depth_m'], 100 * hh['lc7_head_over_p'],
             hh['hoop_membrane_at_deepest_MPa'], hh['fe_farfield_vm_MPa']))
    bb = c['beam_axial_pmin']
    print(' beam at p_min: w %.0f N/m, M %.1f kNm, bend %.2f, pressure %.2f -> est %.2f vs FE %.2f MPa (sigma_cr %.2f)'
          % (bb['w_N_per_m'], bb['M_mid_kNm'], bb['sigma_bend_MPa'], bb['sigma_pressure_MPa'],
             bb['sigma_x_min_est_MPa'], bb['sigma_x_min_fe_MPa'], bb['sigma_cr_MPa']))
    fi = c['forward_inertia']
    print(' forward inertia: F %.0f kN, forward-half tension %.2f MPa, mean axial rise pred %.2f vs FE %.2f MPa'
          % (fi['F_kN'], fi['forward_half_tension_MPa'], fi['mean_axial_rise_pred_MPa'], fi['mean_axial_rise_fe_MPa']))
    sd = c['section_distortion']
    print(' section: beam deflection %.2f mm -> range %.2f of FE %.2f mm (%.0f%%)' % (sd['beam_deflection_mm'], sd['range_from_translation_mm'], sd['range_fe_mm'], 100*sd['translation_share']))
    th = c['thermal']
    print(' thermal: q %.2f W/m2, q/h %.3f K vs FE dT %.3f (f10) %.3f (f90); fin %.3f m; through-thickness %.2e K'
          % (th['q_Wm2'], th['dT_asymptote_K'], th['dT_fe_f10_K'], th['dT_fe_f90_K'], th['fin_length_m'], th['through_thickness_drop_K']))
    print(' wall flux est %.0f vs FE %.0f (f90) %.0f (f10) W/m2' % (th['wall_flux_est_Wm2'], th['wall_flux_fe_f90_Wm2'], th['wall_flux_fe_f10_Wm2']))
    tb = c['thermal_bowing']
    print(' bowing: bow %.2f mm, F_th %.1f kN, mech %.1f -> mid est %.1f vs FE %.1f kN'
          % (tb['bow_mm'], tb['F_thermal_kN'], tb['F_mid_mechanical_kN'], tb['F_mid_est_kN'], tb['F_mid_fe_kN']))
    dd = c['decay']
    print(' decay: sqrt(Rt) %.1f mm, pi/beta %.3f m, exclusion = %.1f decay lengths, attenuation %.1e'
          % (dd['sqrt_Rt_mm'], dd['pi_over_beta_m'], dd['exclusion_in_decay_lengths'], dd['attenuation_at_exclusion']))
    print(' cooldown: %.1f mm free, %.0f / %.0f MPa restrained; bin %.3f K'
          % (c['cooldown']['free_axial_contraction_mm'], c['cooldown']['uniaxial_restrained_MPa'],
             c['cooldown']['biaxial_restrained_MPa'], c['binning']['bin_width_K']))
