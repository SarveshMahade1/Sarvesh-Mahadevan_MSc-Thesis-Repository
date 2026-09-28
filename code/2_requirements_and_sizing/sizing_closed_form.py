# -*- coding: utf-8 -*-
"""
================================================================================
CLOSED-FORM WALL SIZING - CORRECTED SIZING BASIS (27 Sep 2026)
TU Delft MSc thesis - Sarvesh Mahadevan
================================================================================

WHY THIS EXISTS
    The Abaqus sizing runs of 27 Aug (E190) and 26 Sep (A320, A350) applied the
    design pressure to SIDE1 of their shell surfaces, which on those meshes is
    the OUTER face: the 1.5 bar acted as external pressure. The buckling screen
    then read the resulting hoop compression (-pR/t) as axial compression. The
    run logs show it directly: the "compressive axial membrane stress" equals
    pR/t to within 0-17 % at every thickness tried, in every load case. The
    wall thicknesses 4.785 / 6.382 / 9.948 mm are therefore not sized by any
    valid criterion. This script gives the corrected closed-form sizing; the
    finite-element model then verifies the result on campus.

PHYSICS (every term is a standard closed form; sources in brackets)
    Barrel far field, rings at x = 0, L/2, L (aft sliding, mid anchor, fwd
    sliding). Axial membrane stress, tension positive:

      pressure      s_p = p R / (2 t)                          [thin-shell membrane]
      inertia       s_i = n_x g m_half / (2 pi R t)            [all axial inertia enters
                    m_half = m_shell/2 + rho V min(phi, 0.5)     the mid anchor; the heavier
                                                                 half is in compression]
      bending       s_b = M / (pi R^2 t),  M = w (L/2)^2 / 8   [continuous beam, two equal
                    w = n_t g m_tot / L,  n_t = sqrt(n_y^2+n_z^2) spans, moment at centre]
      minimum       s_min = s_p - s_i - s_b

    Axial buckling allowable, NASA SP-8007 (1968), unpressurised knockdown:
      s_cl  = E t / (R sqrt(3 (1 - nu^2)))
      phi   = (1/16) sqrt(R/t)
      gamma = 1 - 0.901 (1 - exp(-phi))
      s_cr  = gamma s_cl
    Pressure stabilisation is NOT credited (conservative).

    Strength (yield) check at the maximum pressure, with the hydrostatic head
    measured along the tilted resultant:
      p_hyd = rho |n| g d,  d = extent of liquid along g_hat
      s_hoop = (p + p_hyd) R / t
      s_ax   = s_p +/- (s_i + s_b)      (both signs checked)
      von Mises (plane stress), compared with F_ty / SF = 262 MPa.

DESIGN BASIS (pressures from Chapter 3)
    STABILITY  at the MINIMUM differential the tank can hold when the load
               occurs: fill / minimum operating pressure 1.2 bar absolute
               (Winnefeld et al. after Verstraete et al.) against ISA sea-level
               ambient 1.01325 bar -> 0.18675 bar gauge. Pressure stabilises
               the shell against axial buckling, so the minimum is the
               critical value. The emergency landing occurs at ground level.
    STRENGTH   at the MAXIMUM design differential, 1.5 bar gauge (Chapter 3).
    Target     FOS 1.5 against s_cr; von Mises <= 262 MPa (FOS 1.5 on F_ty).

This is a SIZING ESTIMATE. The finite-element model verifies it at the chosen
thickness (E190_MASTER_thermo_structural.py, RUN_MODE='SIZING_CHECK').
================================================================================
"""
import math

E, NU, RHO_AL, RHO_LH2, G = 73.1e9, 0.33, 2850.0, 70.8, 9.81
F_TY, SF = 393e6, 1.5
P_MAX = 150000.0                       # Pa gauge, strength
P_MIN = 120000.0 - 101325.0            # Pa gauge, stability
FOS_BUCKLING = 1.5
LOAD_CASES = {'Baseline': (0, 0, 1), 'LC1_Maneuver': (0, 0, 2.5),
              'LC3_Combined': (1.5, 0, 2.5), 'LC7_Emergency': (9, 3, 6)}
FILLS = (0.1, 0.3, 0.5, 0.75, 0.9)
AIRCRAFT = {'A320': (1.7080, 10.0120), 'E190': (1.3015, 13.9675),
            'A350': (2.5771, 33.8240)}      # R_mid, L_barrel  (AIRCRAFT_DATA)
OLD_T = {'A320': 6.382, 'E190': 4.785, 'A350': 9.948}


def geometry(R, L):
    c = R / 1.6
    e = math.sqrt(1 - (c / R) ** 2)
    A = 2 * math.pi * R * L + 2 * math.pi * R * R * (1 + (1 - e * e) / e * math.atanh(e))
    V = math.pi * R * R * L + 4.0 / 3.0 * math.pi * R * R * c
    return A, V, L + 2 * c


def sigma_cr(R, t):
    s_cl = E * t / (R * math.sqrt(3 * (1 - NU * NU)))
    phi = math.sqrt(R / t) / 16.0
    return (1 - 0.901 * (1 - math.exp(-phi))) * s_cl


def axial_terms(R, L, t, n, phi_f, p):
    A, V, Lt = geometry(R, L)
    m_sh = A * t * RHO_AL
    m_tot = m_sh + RHO_LH2 * V * phi_f
    area = 2 * math.pi * R * t
    s_p = p * R / (2 * t)
    m_half = m_sh / 2 + RHO_LH2 * V * min(phi_f, 0.5)
    s_i = abs(n[0]) * G * m_half / area
    n_t = math.hypot(n[1], n[2])
    w = n_t * G * m_tot / L
    s_b = (w * (L / 2) ** 2 / 8) / (math.pi * R * R * t)
    return s_p, s_i, s_b


def buckling_fos(R, L, t):
    worst = (float('inf'), None)
    for lc, n in LOAD_CASES.items():
        for f in FILLS:
            s_p, s_i, s_b = axial_terms(R, L, t, n, f, P_MIN)
            s_min = s_p - s_i - s_b
            fos = sigma_cr(R, t) / -s_min if s_min < 0 else float('inf')
            if fos < worst[0]:
                worst = (fos, (lc, f, s_min))
    return worst


def yield_fos(R, L, t):
    A, V, Lt = geometry(R, L)
    worst = (float('inf'), None)
    for lc, n in LOAD_CASES.items():
        nm = math.sqrt(sum(x * x for x in n))
        ax_frac, tr_frac = abs(n[0]) / nm, math.hypot(n[1], n[2]) / nm
        d_full = Lt * ax_frac + 2 * R * tr_frac      # extent of the tank along g_hat
        for f in FILLS:
            s_p, s_i, s_b = axial_terms(R, L, t, n, f, P_MAX)
            p_h = RHO_LH2 * nm * G * d_full * f       # liquid column ~ fill x extent
            s_h = (P_MAX + p_h) * R / t
            for s_a in (s_p + s_i + s_b, s_p - s_i - s_b):
                vm = math.sqrt(s_h ** 2 - s_h * s_a + s_a ** 2)
                fos = (F_TY / SF) / vm
                if fos < worst[0]:
                    worst = (fos, (lc, f, vm))
    return worst


def size(R, L):
    lo, hi = 0.5e-3, 30e-3
    for _ in range(60):                               # bisection on buckling FOS
        mid = 0.5 * (lo + hi)
        if buckling_fos(R, L, mid)[0] >= FOS_BUCKLING:
            hi = mid
        else:
            lo = mid
    t = hi
    return t, buckling_fos(R, L, t), yield_fos(R, L, t)


if __name__ == '__main__':
    print('Stability at p = %.4f bar gauge, strength at p = %.2f bar gauge\n'
          % (P_MIN / 1e5, P_MAX / 1e5))
    print('%-5s  %9s  %9s  %-28s  %10s  %-26s' % ('', 't_new', 't_old', 'buckling FOS 1.50 at',
                                               'yield FOS', 'yield governing at'))
    for ac, (R, L) in AIRCRAFT.items():
        t, (bf, bat), (yf, yat) = size(R, L)
        print('%-5s  %6.3f mm  %6.3f mm  %-28s  %10.2f  %-26s'
              % (ac, t * 1e3, OLD_T[ac], '%s f=%.2f (%.1f MPa)' % (bat[0], bat[1], bat[2] / 1e6),
                 yf, '%s f=%.2f (%.1f MPa)' % (yat[0], yat[1], yat[2] / 1e6)))
    print('\nAt the OLD thicknesses, with the corrected physics:')
    for ac, (R, L) in AIRCRAFT.items():
        t = OLD_T[ac] / 1e3
        bf, bat = buckling_fos(R, L, t)
        yf, yat = yield_fos(R, L, t)
        print('  %-5s t=%.3f mm   buckling FOS %.2f (%s f=%.2f)   yield FOS %.2f'
              % (ac, t * 1e3, bf, bat[0], bat[1], yf))
