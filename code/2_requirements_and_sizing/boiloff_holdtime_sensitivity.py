"""
================================================================================
SENSITIVITY OF THE INSULATION DESIGN TO THE REQUIRED NO-VENT HOLD TIME
Pure Python + CoolProp + matplotlib - TU Delft MSc thesis

Author: Sarvesh Mahadevan
Date: September 2026

--------------------------------------------------------------------------------
WHY THIS SCRIPT EXISTS
--------------------------------------------------------------------------------
boiloff_requirement_derivation.py establishes that the heat-leak budget follows
from a required no-vent hold time rather than from an assumed boil-off
percentage. That leaves one design choice: how long the tank must be able to
sit sealed.

This script sweeps that choice and shows what it costs. For each aircraft it
propagates the hold time through the same steady-state resistance network used
in the *_vacuum_mli_gradient.py scripts, and reports the required MLI
thickness, the resulting outer-jacket diameter and the insulation mass.

This is the parametric study the thesis previously lacked: it converts a single
assumed number into a demonstrated trade.

GENERATES: latex/figures/17_boiloff_holdtime_sensitivity.png

--------------------------------------------------------------------------------
FORMULAS USED
--------------------------------------------------------------------------------
FORMULA 1  Constant-volume energy balance    Q * t_hold = m * (u2 - u1)
FORMULA 2  Allowable heat leak               Q_allow = m * delta_u / t_hold
FORMULA 3  Cylindrical conduction resistance R = ln(r_out/r_in) / (2*pi*k*L)
FORMULA 4  Convective surface resistance     R = 1 / (h * A)
FORMULA 5  Series resistance network         Q = (T_amb - T_liq) / R_total
FORMULA 6  Required MLI outer radius         r2 = r1 * exp(R_MLI * 2*pi*k_eff*L)
FORMULA 7  Insulation mass                   m = rho * pi * (r2^2 - r1^2) * L
================================================================================
"""

import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
import thesis_figstyle as tfs
tfs.apply_style()
import numpy as np

import CoolProp.CoolProp as CP

from boiloff_requirement_derivation import (
    FLUID, P_FILL, P_VENT, FILL_FRACTION, H_FG, AIRCRAFT,
    tank_volume, tank_length, lockup_state, internal_energy_at,
)


# --- Thermal model constants, identical to the *_vacuum_mli_gradient.py scripts ---
T_AMBIENT = 293.0       # K
T_LIQUID = 20.0         # K
K_INNER = 116.0         # W/m-K, Al 2219-T87
K_OUTER = 116.0         # W/m-K
K_EFF_MLI = 0.0003      # W/m-K, installed (non-ideal) blanket - SOURCE STILL OPEN
H_CONV = 10.0           # W/m2-K, external natural convection
T_WALL_INNER = 0.002    # m
T_WALL_OUTER = 0.001    # m
RHO_MLI = 173.0         # kg/m3, aluminised-mylar/scrim blanket bulk density

# The design point adopted in the thesis.
T_HOLD_DESIGN_H = 14.68


def required_mli_thickness(R_mid, L_tank, Q_target):
    """
    Invert the steady-state resistance network for the MLI thickness that
    yields exactly Q_target.

    The network is, in series from the cold face outwards:
        inner wall conduction -> MLI conduction -> outer wall conduction
        -> external convection

    FORMULA 5 gives the total resistance the network must have:
        R_total = (T_amb - T_liq) / Q_target

    The wall and convection terms are known, so the MLI term is what remains,
    and FORMULA 6 converts that resistance back into a radius.
    """
    r0 = R_mid - T_WALL_INNER / 2.0
    r1 = R_mid + T_WALL_INNER / 2.0

    R_total_required = (T_AMBIENT - T_LIQUID) / Q_target
    R_inner_wall = math.log(r1 / r0) / (2.0 * math.pi * K_INNER * L_tank)

    # The outer-wall and convection resistances depend on the (unknown) outer
    # radius, so solve by fixed-point iteration. Both terms are tiny compared
    # with the MLI term, so this converges in a handful of passes.
    mli_t = 0.015
    for _ in range(200):
        r2 = r1 + mli_t
        r3 = r2 + T_WALL_OUTER
        R_outer_wall = math.log(r3 / r2) / (2.0 * math.pi * K_OUTER * L_tank)
        R_conv = 1.0 / (H_CONV * 2.0 * math.pi * r3 * L_tank)      # FORMULA 4

        R_MLI_required = R_total_required - R_inner_wall - R_outer_wall - R_conv
        if R_MLI_required <= 0:
            return float('nan'), float('nan')

        # FORMULA 6
        r2_new = r1 * math.exp(R_MLI_required * 2.0 * math.pi * K_EFF_MLI * L_tank)
        mli_t_new = r2_new - r1
        if abs(mli_t_new - mli_t) < 1e-10:
            mli_t = mli_t_new
            break
        mli_t = mli_t_new

    r2 = r1 + mli_t
    # FORMULA 7
    mass = RHO_MLI * math.pi * (r2**2 - r1**2) * L_tank
    return mli_t, mass


def main():
    hold_h = np.linspace(2.0, 36.0, 300)

    results = {}
    for name, g in AIRCRAFT.items():
        V = tank_volume(g['R_mid'], g['L_bar'], g['cap_ar'])
        L = tank_length(g['R_mid'], g['L_bar'], g['cap_ar'])
        m_total, _, _, rho_bulk = lockup_state(V, FILL_FRACTION, P_FILL)
        u1, _ = internal_energy_at(rho_bulk, P_FILL)
        u2, _ = internal_energy_at(rho_bulk, P_VENT)
        du = u2 - u1

        Q = m_total * du / (hold_h * 3600.0)                       # FORMULA 2
        pct = 100.0 * (Q / H_FG) * 3600.0 / (V * 70.8)             # FORMULA 3

        tk, ms = [], []
        for q in Q:
            t_, m_ = required_mli_thickness(g['R_mid'], L, q)
            tk.append(t_ * 1000.0)     # mm
            ms.append(m_)

        # The adopted design point
        Q_d = m_total * du / (T_HOLD_DESIGN_H * 3600.0)
        t_d, m_d = required_mli_thickness(g['R_mid'], L, Q_d)
        pct_d = 100.0 * (Q_d / H_FG) * 3600.0 / (V * 70.8)

        results[name] = dict(Q=Q, pct=pct, t=np.array(tk), m=np.array(ms),
                             R=g['R_mid'], L=L, V=V, m_total=m_total,
                             Q_d=Q_d, t_d=t_d * 1000.0, m_d=m_d, pct_d=pct_d)

    # ---------------- console table ----------------
    print("=" * 96)
    print("INSULATION SENSITIVITY TO THE REQUIRED NO-VENT HOLD TIME")
    print("=" * 96)
    print("Design point adopted: t_hold = %.2f h "
          "(equivalent to the 0.1 %%/hr figure previously used)" % T_HOLD_DESIGN_H)
    print("k_eff_MLI = %s W/m-K   rho_MLI = %s kg/m3" % (K_EFF_MLI, RHO_MLI))
    print("")
    for name, r in results.items():
        print("  %-10s  Q = %8.1f W   MLI = %6.2f mm   "
              "blanket mass = %7.1f kg   (%.4f %%/hr)"
              % (name, r['Q_d'], r['t_d'], r['m_d'], r['pct_d']))
    print("")
    hdr = "  %10s | " % 't_hold [h]'
    for n in results:
        hdr += "%16s | " % (n + ' MLI [mm]')
    print(hdr)
    for th in [3, 6, 9, 12, 14.68, 18, 24, 30, 36]:
        idx = int(np.argmin(np.abs(hold_h - th)))
        row = "  %10.2f | " % th
        for n in results:
            row += "%16.2f | " % results[n]['t'][idx]
        print(row)
    print("=" * 96)

    # ---------------- figure ----------------
    NAVY, RED, GREEN, GREY = '#1F3864', '#C00000', '#2E7D32', '#595959'
    colours = {'E190': NAVY, 'A320': GREEN, 'A350-900': RED}
    NAVY = 'black'

    fig, axes = plt.subplots(1, 3, figsize=(455.24 / 72.27, 2.9))

    # (a) allowable heat leak
    ax = axes[0]
    for n, r in results.items():
        ax.plot(hold_h, r['Q'] / 1000.0, lw=1.3, color=colours[n], label=n)
        ax.plot([T_HOLD_DESIGN_H], [r['Q_d'] / 1000.0], 'o', ms=4,
                color=colours[n], zorder=5)
    ax.set_yscale('log')
    ax.set_xlabel('Hold time $t_{hold}$ [h]')
    ax.set_ylabel('Allowable heat leak $Q$ [kW]')
    ax.set_title('(a) Heat-leak budget', loc='left', fontsize=9.5)
    ax.grid(alpha=0.3, which='both')
    pass

    # (b) required MLI thickness
    ax = axes[1]
    for n, r in results.items():
        ax.plot(hold_h, r['t'], lw=1.3, color=colours[n], label=n)
        ax.plot([T_HOLD_DESIGN_H], [r['t_d']], 'o', ms=4, color=colours[n], zorder=5)
    ax.set_xlabel('Hold time $t_{hold}$ [h]')
    ax.set_ylabel('MLI thickness [mm]')
    ax.set_title('(b) MLI thickness', loc='left', fontsize=9.5)
    ax.grid(alpha=0.3)

    # (c) blanket mass
    ax = axes[2]
    for n, r in results.items():
        ax.plot(hold_h, r['m'], lw=1.3, color=colours[n], label=n)
        ax.plot([T_HOLD_DESIGN_H], [r['m_d']], 'o', ms=4, color=colours[n], zorder=5)
    ax.set_yscale('log')
    ax.set_xlabel('Hold time $t_{hold}$ [h]')
    ax.set_ylabel('MLI mass [kg]')
    ax.set_title('(c) MLI mass', loc='left', fontsize=9.5)
    ax.grid(alpha=0.3, which='both')

    for ax in axes:
        ax.axvline(T_HOLD_DESIGN_H, color=GREY, ls='--', lw=1.2, zorder=1)
        pass

# --- moved to the LaTeX caption; the figure carries labelling only
#     fig.suptitle('Sensitivity of the Insulation Design to the Required No-Vent Hold Time',
#                  fontsize=15, fontweight='bold', color=NAVY, y=0.975)
# --- moved to the LaTeX caption; the figure carries labelling only
#     fig.text(0.5, 0.865,
#              'The heat-leak budget follows from a constant-volume energy balance between '
#              'the fill pressure (1.2 bar) and the venting pressure (1.448 bar).\n'
#              'The adopted hold time of 14.68 h reproduces the 0.1 %/hr figure used '
#              'throughout the thesis, and is identical for all three aircraft.',
#              ha='center', fontsize=9.5, color=GREY)

    out = os.path.join(RP.FIGURES, '17_boiloff_holdtime_sensitivity.png')
    d = os.path.dirname(out)
    if not os.path.isdir(d):
        os.makedirs(d)
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], color=colours[n], lw=1.3, marker='o', ms=4) for n in results] + \
         [Line2D([], [], color=GREY, ls='--', lw=1.2)]
    fig.legend(hs, list(results) + ['Adopted hold time, %.1f h (dots: design point)' % T_HOLD_DESIGN_H],
               loc='lower center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0), fontsize=9)
    fig.tight_layout(rect=(0, 0.12, 1, 1), w_pad=0.8)
    tfs.save(fig, out)
    print("\nFigure written to %s" % out)


if __name__ == '__main__':
    main()
