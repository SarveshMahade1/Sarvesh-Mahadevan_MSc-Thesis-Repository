"""
================================================================================
SENSITIVITY OF THE INSULATION DESIGN TO THE REQUIRED NO-VENT HOLD TIME
Pure Python + CoolProp + matplotlib - TU Delft MSc thesis

Author: Sarvesh Mahadevan

RUN      python code/2_requirements_and_sizing/boiloff_holdtime_sensitivity.py
INPUTS   code/design_basis.json (pressures, ambient and liquid temperature,
         effective MLI conductivity, adopted hold time) and
         results/processed/final_numbers.json (heat leak budgets); needs CoolProp
OUTPUTS  thesis_figures/17_boiloff_holdtime_sensitivity.png and a printed table

--------------------------------------------------------------------------------
PURPOSE
--------------------------------------------------------------------------------
boiloff_requirement_derivation.py establishes that the heat-leak budget follows
from a required no-vent hold time rather than from an assumed boil-off
percentage. That leaves one design choice: how long the tank must be able to
sit sealed.

This script sweeps that choice (2 h to 36 h) and shows what it costs. For each
aircraft it propagates the hold time through a steady-state resistance
network (FORMULAS 3 to 6) and reports the required MLI thickness, the
resulting outer-jacket diameter and the insulation mass, which turns the
single hold-time requirement into a demonstrated trade.

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
import design_basis as DB  # noqa: E402  (single source of the design inputs)
import thesis_figstyle as tfs
tfs.apply_style()
import numpy as np

import CoolProp.CoolProp as CP

from boiloff_requirement_derivation import (
    FLUID, P_FILL, P_VENT, FILL_FRACTION, H_FG, AIRCRAFT,
    tank_volume, tank_length, lockup_state, internal_energy_at,
)


# --- Thermal model constants, identical to the *_vacuum_mli_gradient.py scripts ---
T_AMBIENT = DB.TH['T_ambient_K']     # K
T_LIQUID = 20.0         # K, end temperature of the insulation network (thesis Section 3.4)
K_INNER = 116.0         # W/m-K, Al 2219-T87
K_OUTER = 116.0         # W/m-K
K_EFF_MLI = DB.TH['k_mli_effective_W_per_mK']   # W/m-K, installed blanket, adopted engineering value
H_CONV = 10.0           # W/m2-K, external natural convection
T_WALL_INNER = 0.002    # m
T_WALL_OUTER = 0.001    # m
RHO_MLI = 173.0         # kg/m3, aluminised-mylar/scrim blanket bulk density

# The design point adopted in the thesis.
T_HOLD_DESIGN_H = DB.BUDGET['hold_time_h']


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

    plot_per_aircraft(results, hold_h)


def plot_per_aircraft(results, hold_h):
    """Figure 17: one column per aircraft, each with its own
    scales. Top row: allowable heat leak Q against the hold time. Bottom row:
    the MLI thickness that delivers exactly that Q. The points at 6, 14.7, 24
    and 36 h carry their values; the note gives the constant energy Q t_hold
    that the contents can absorb between 1.2 and 1.448 bar, and the blanket
    mass at the adopted and at the 24 h hold time."""
    from matplotlib.ticker import MaxNLocator
    colours = {'E190': tfs.AC_COLOR['E190'], 'A320': tfs.AC_COLOR['A320'], 'A350-900': tfs.AC_COLOR['A350']}
    marks = [6.0, T_HOLD_DESIGN_H, 24.0, 36.0]
    import json
    FN = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
    QBUDGET = {'E190': FN['E190']['Q_W'], 'A320': FN['A320']['Q_W'], 'A350-900': FN['A350']['Q_W']}
    fig, axs = plt.subplots(2, 3, figsize=(455.24 / 72.27, 4.5))
    for k, (n, r) in enumerate(results.items()):
        col = colours[n]
        # The budget of the report (Table 3.4, final_numbers.json) anchors the
        # curve: Q t_hold is constant, so Q(t) = Q_budget * 14.68 / t. The
        # CoolProp state of this script gives a budget 0.03 % different.
        Qd = QBUDGET[n]
        Qk = Qd * T_HOLD_DESIGN_H / hold_h
        Qm = [Qd * T_HOLD_DESIGN_H / m for m in marks]
        tm = [1000.0 * required_mli_thickness(r['R'], r['L'], q)[0] for q in Qm]
        top, bot = axs[0, k], axs[1, k]
        top.plot(hold_h, Qk, color=col, lw=1.4)
        tk = np.array([1000.0 * required_mli_thickness(r['R'], r['L'], q)[0] for q in Qk])
        bot.plot(hold_h, tk, color=col, lw=1.4)
        top.plot(marks, Qm, 'o', color=col, ms=4)
        bot.plot(marks, tm, 'o', color=col, ms=4)
        tfs.label_points(top, marks, Qm, col, '%.0f', side='ur')
        tfs.label_points(bot, marks, tm, col, '%.1f', side='lr')
        for ax in (top, bot):
            ax.axvline(T_HOLD_DESIGN_H, color='0.4', ls='--', lw=0.9, zorder=1)
            ax.axvline(24.0, color='0.4', ls=':', lw=0.9, zorder=1)
            ax.grid(True)
            ax.yaxis.set_major_locator(MaxNLocator(5))
        E_MJ = Qd * T_HOLD_DESIGN_H * 3600.0 / 1e6
        tfs.note(top, 'Heat absorbed\nbefore venting: %.1f MJ' % E_MJ, loc='ur')
        tfs.note(bot, 'MLI mass %.0f kg at 14.7 h,\n%.0f kg at 24 h' % (required_mli_thickness(r['R'], r['L'], Qd)[1], required_mli_thickness(r['R'], r['L'], Qd * T_HOLD_DESIGN_H / 24.0)[1]), loc='ul')
        top.set_ylim(0, 1.6 * Qd * T_HOLD_DESIGN_H / 4.5)
        top.set_xlim(2, 38)
        bot.set_ylim(0, 1.4 * tk.max())
        top.set_title('(%s) %s' % ('abc'[k], n), loc='left', fontsize=9.5)
        bot.set_title('(%s) %s' % ('def'[k], n), loc='left', fontsize=9.5)
        for ax in (top, bot):     # every panel carries both axis labels
            ax.set_xlabel(r'No-vent hold time $t_{\mathrm{hold}}$ [h]')
            ax.set_xticks([6, 14.68, 24, 36])
            ax.set_xticklabels(['6', '14.7', '24', '36'])
            ax.set_xlim(3, 38)
        top.set_ylabel('Allowable heat leak, $Q$ [W]')
        bot.set_ylabel('Required MLI insulation\nthickness [mm]')
    fig.tight_layout(h_pad=1.0, w_pad=1.2)
    out = os.path.join(RP.FIGURES, '17_boiloff_holdtime_sensitivity.png')
    tfs.save(fig, out)


if __name__ == '__main__':
    main()
