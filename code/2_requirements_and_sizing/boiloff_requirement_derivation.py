"""
================================================================================
DERIVATION OF THE BOIL-OFF / HEAT-LEAK REQUIREMENT FROM A NO-VENT HOLD TIME
Pure Python + CoolProp (no Abaqus needed) - TU Delft MSc thesis

Author: Sarvesh Mahadevan
Date: September 2026

--------------------------------------------------------------------------------
WHY THIS SCRIPT EXISTS
--------------------------------------------------------------------------------
The insulation in this thesis was previously sized against an ASSUMED boil-off
target of 0.1 %/hr of tank capacity, adopted from the literature and then
locked. That number was never derived, and the citation attached to it was
incorrect.

This script replaces the assumption with a derivation. It follows the approach
of Winnefeld et al. (2018), who do not assume a %/hr figure at all: they size
the insulation from a physically meaningful operational requirement, namely
that the tank must be able to sit sealed for a required hold time without the
pressure reaching the relief-valve setting.

The %/hr figure then becomes an OUTPUT of the derivation, not an input.

--------------------------------------------------------------------------------
THE PHYSICS
--------------------------------------------------------------------------------
A sealed ("locked-up") tank contains a saturated two-phase mixture of liquid
and vapour hydrogen. Its total mass and total volume are both fixed, so the
overall (bulk) density rho = m_total / V_tank is a constant.

Heat leaking in cannot escape, so it accumulates as internal energy. The state
therefore moves along a line of CONSTANT DENSITY in the two-phase region, and
the saturation pressure rises. When the pressure reaches the relief-valve
setting, the tank must vent, and fuel is lost.

The energy balance over the hold period, at constant volume and constant mass:

    FORMULA 1 (constant-volume energy balance)
        Q * t_hold = m_total * (u2 - u1)

where u1 is the specific internal energy of the mixture at the fill pressure
and u2 at the venting pressure, both evaluated at the SAME bulk density.

Rearranged, the allowable steady heat leak is

    FORMULA 2 (allowable heat leak)
        Q_allow = m_total * (u2 - u1) / t_hold

For comparison against the older convention, this is converted to an equivalent
boil-off percentage:

    FORMULA 3 (equivalent boil-off rate)
        mdot_boiloff = Q_allow / h_fg
        boiloff_%/hr = 100 * mdot_boiloff * 3600 / m_LH2_capacity

--------------------------------------------------------------------------------
DATA SOURCES (every number traced)
--------------------------------------------------------------------------------
1. Hydrogen properties: CoolProp, which implements the reference equation of
   state of Leachman, Jacobsen, Penoncello & Lemmon (2009), J. Phys. Chem. Ref.
   Data 38(3) 721-748 - already cited as ref14 in the thesis. ParaHydrogen is
   used because stored LH2 is catalysed to the para form to avoid the
   ortho-para conversion heat release.

2. Fill / minimum pressure p_fill = 1.2 bar. Taken from Winnefeld et al. (2018),
   who take it from Verstraete et al. (2010) - ref12 and ref10 respectively.
   The tank is held above ambient so that a leak vents outward rather than
   admitting air, which would freeze and form an explosive mixture.

3. Venting pressure p_vent = 1.448 bar. Winnefeld et al. (2018) state this is
   the operating pressure derived from Brewer, Hydrogen Aircraft Technology -
   ref12 and ref6 respectively.

4. Maximum fill fraction 0.97. Winnefeld et al. (2018), after Verstraete - an
   ullage volume is required so that the tank is never liquid-full.

5. Tank geometry: from the three MASTER thermo-structural scripts.

--------------------------------------------------------------------------------
KEY RESULT
--------------------------------------------------------------------------------
The 0.1 %/hr figure used throughout the thesis corresponds to a no-vent hold
time of 14.68 HOURS - an aircraft standing fuelled overnight without venting.
That hold time is IDENTICAL for all three aircraft, because the bulk density
and hence the internal energy change depend only on the fill fraction and the
two operating pressures, not on tank size. It is therefore a legitimate common
design requirement across a regional, a narrow-body and a wide-body aircraft.
================================================================================
"""

import math

try:
    import CoolProp.CoolProp as CP
except ImportError:
    raise SystemExit("CoolProp is required:  pip install CoolProp")

FLUID = 'ParaHydrogen'

# ============================================================================
# INPUTS
# ============================================================================

# --- Operating pressures (see DATA SOURCES 2 and 3) ---
P_FILL = 1.2e5      # Pa, fill / minimum pressure          [Winnefeld, Verstraete]
P_VENT = 1.448e5    # Pa, relief-valve / venting pressure  [Winnefeld, after Brewer]

# --- Maximum fill fraction (see DATA SOURCE 4) ---
FILL_FRACTION = 0.97

# --- Latent heat used by the existing thermal scripts, kept for consistency ---
H_FG = 446000.0     # J/kg

# --- The previously assumed target, retained only for back-comparison ---
LEGACY_BOILOFF_PCT_PER_HR = 0.1

# --- Tank geometry, matching the MASTER thermo-structural scripts ---
#     R_mid  : mid-surface radius of the inner vessel            [m]
#     L_bar  : cylindrical barrel length                         [m]
#     cap_ar : end-cap aspect ratio a/b (semi-ellipsoidal)        [-]
AIRCRAFT = {
    # E190 L_bar is back-solved from the capacity hard-coded in the E190
    # vacuum/MLI sizing script (m_LH2_capacity = 5670.9 kg at rho = 70.8
    # kg/m3, i.e. V = 80.098 m3). The resulting total length, 15.594 m,
    # reproduces L_tank in that script exactly - a useful check.
    'E190':     dict(R_mid=1.3015, L_bar=13.9675, cap_ar=1.6),
    'A320':     dict(R_mid=1.7080, L_bar=10.0120, cap_ar=1.6),
    'A350-900': dict(R_mid=2.5771, L_bar=33.8240, cap_ar=1.6),
}

# Hold times reported in the sensitivity sweep [hours]
HOLD_TIMES_H = [1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 36, 48]


# ============================================================================
# GEOMETRY
# ============================================================================

def tank_volume(R_mid, L_bar, cap_ar):
    """
    FORMULA 4: internal volume of a cylinder closed by two semi-ellipsoidal
    caps of height h_cap = R / (a/b).

        V = pi*R^2*L_barrel  +  (4/3)*pi*R^2*h_cap

    The second term is the volume of the two caps combined, which together
    form one full ellipsoid of semi-axes (R, R, h_cap).
    """
    h_cap = R_mid / cap_ar
    return math.pi * R_mid**2 * L_bar + (4.0 / 3.0) * math.pi * R_mid**2 * h_cap


def tank_length(R_mid, L_bar, cap_ar):
    """Total external length of the inner vessel, barrel plus two caps."""
    return L_bar + 2.0 * (R_mid / cap_ar)


# ============================================================================
# THERMODYNAMICS
# ============================================================================

def lockup_state(V_tank, fill_fraction, p_fill):
    """
    Establish the sealed-tank initial state.

    The tank is filled to `fill_fraction` by VOLUME with saturated liquid at
    p_fill; the remaining ullage holds saturated vapour at the same pressure.
    The bulk density that results is then fixed for the whole hold period,
    because neither the mass nor the volume changes once the tank is sealed.
    """
    rho_liq = CP.PropsSI('D', 'P', p_fill, 'Q', 0, FLUID)
    rho_vap = CP.PropsSI('D', 'P', p_fill, 'Q', 1, FLUID)

    m_liq = rho_liq * V_tank * fill_fraction
    m_vap = rho_vap * V_tank * (1.0 - fill_fraction)
    m_total = m_liq + m_vap

    rho_bulk = m_total / V_tank            # constant for the whole hold
    return m_total, m_liq, m_vap, rho_bulk


def internal_energy_at(rho_bulk, p):
    """
    FORMULA 5: specific internal energy of the two-phase mixture at the given
    saturation pressure, evaluated at the fixed bulk density.

        x = (1/rho_bulk - v_f) / (v_g - v_f)      (vapour quality, by volume)
        u = u_f + x * (u_g - u_f)

    This is the standard two-phase lever rule. It is used rather than calling
    CoolProp with (D, P) directly because inside the vapour dome that pair is
    not independent for a pure fluid - pressure and temperature are linked by
    the saturation line, so the quality must be obtained from the density.
    """
    v_f = 1.0 / CP.PropsSI('D', 'P', p, 'Q', 0, FLUID)
    v_g = 1.0 / CP.PropsSI('D', 'P', p, 'Q', 1, FLUID)
    u_f = CP.PropsSI('U', 'P', p, 'Q', 0, FLUID)
    u_g = CP.PropsSI('U', 'P', p, 'Q', 1, FLUID)

    v_bulk = 1.0 / rho_bulk
    x = (v_bulk - v_f) / (v_g - v_f)
    return u_f + x * (u_g - u_f), x


def allowable_heat_leak(m_total, delta_u, t_hold_s):
    """FORMULA 2: Q_allow = m_total * delta_u / t_hold."""
    return m_total * delta_u / t_hold_s


def equivalent_boiloff_pct_per_hr(Q, m_capacity):
    """FORMULA 3: convert an allowable heat leak to the legacy %/hr figure."""
    mdot = Q / H_FG                                  # kg/s
    return 100.0 * mdot * 3600.0 / m_capacity        # %/hr


# ============================================================================
# MAIN
# ============================================================================

def main():
    T_fill = CP.PropsSI('T', 'P', P_FILL, 'Q', 0, FLUID)
    T_vent = CP.PropsSI('T', 'P', P_VENT, 'Q', 0, FLUID)

    print("=" * 78)
    print("DERIVATION OF THE HEAT-LEAK REQUIREMENT FROM A NO-VENT HOLD TIME")
    print("=" * 78)
    print("Fluid                 : %s (Leachman EOS via CoolProp, ref14)" % FLUID)
    print("Fill pressure         : %.3f bar  -> T_sat = %.3f K" % (P_FILL / 1e5, T_fill))
    print("Venting pressure      : %.3f bar  -> T_sat = %.3f K" % (P_VENT / 1e5, T_vent))
    print("Saturation rise       : %.3f K" % (T_vent - T_fill))
    print("Fill fraction         : %.2f by volume" % FILL_FRACTION)
    print("=" * 78)

    summary = {}

    for name, g in AIRCRAFT.items():
        V = tank_volume(g['R_mid'], g['L_bar'], g['cap_ar'])
        L = tank_length(g['R_mid'], g['L_bar'], g['cap_ar'])

        m_total, m_liq, m_vap, rho_bulk = lockup_state(V, FILL_FRACTION, P_FILL)

        u1, x1 = internal_energy_at(rho_bulk, P_FILL)
        u2, x2 = internal_energy_at(rho_bulk, P_VENT)
        delta_u = u2 - u1

        # Capacity convention used by the existing MLI scripts: rho = 70.8
        # kg/m3 at f = 1.0. Retained so the %/hr numbers are directly
        # comparable with the previously reported values.
        m_capacity_legacy = V * 70.8

        # Hold time that reproduces the legacy 0.1 %/hr target.
        Q_legacy = (LEGACY_BOILOFF_PCT_PER_HR / 100.0) * m_capacity_legacy / 3600.0 * H_FG
        t_equiv_h = m_total * delta_u / Q_legacy / 3600.0

        print("")
        print("-" * 78)
        print(name)
        print("-" * 78)
        print("  R_mid = %.4f m, L_barrel = %.4f m, L_total = %.3f m"
              % (g['R_mid'], g['L_bar'], L))
        print("  V_tank                     = %10.3f m3" % V)
        print("  Sealed mass (f = %.2f)      = %10.1f kg (liquid %.1f + vapour %.1f)"
              % (FILL_FRACTION, m_total, m_liq, m_vap))
        print("  Bulk density (constant)    = %10.3f kg/m3" % rho_bulk)
        print("  Vapour quality  fill/vent  = %.6f / %.6f" % (x1, x2))
        print("  u at fill  / u at vent     = %.3f / %.3f kJ/kg" % (u1 / 1000, u2 / 1000))
        print("  delta_u                    = %10.4f kJ/kg" % (delta_u / 1000))
        print("  Capacity (legacy 70.8 conv)= %10.1f kg" % m_capacity_legacy)
        print("")
        print("  >> The previously assumed %.2f %%/hr (Q = %.1f W) is equivalent to a"
              % (LEGACY_BOILOFF_PCT_PER_HR, Q_legacy))
        print("     no-vent hold time of %.2f HOURS." % t_equiv_h)
        print("")
        print("  %11s | %12s | %12s | %8s" % ('t_hold [h]', 'Q_allow [W]', 'equiv. %/hr', '%/day'))
        print("  %s-+-%s-+-%s-+-%s" % ('-' * 11, '-' * 12, '-' * 12, '-' * 8))
        for th in HOLD_TIMES_H:
            Q = allowable_heat_leak(m_total, delta_u, th * 3600.0)
            pct = equivalent_boiloff_pct_per_hr(Q, m_capacity_legacy)
            print("  %11d | %12.1f | %12.4f | %8.2f" % (th, Q, pct, pct * 24))

        summary[name] = dict(V=V, m_total=m_total, delta_u=delta_u,
                             t_equiv_h=t_equiv_h, Q_legacy=Q_legacy)

    print("")
    print("=" * 78)
    print("SUMMARY: hold time equivalent to the 0.1 %/hr target previously used")
    print("=" * 78)
    print("  %-12s | %9s | %13s | %8s | %10s"
          % ('Aircraft', 'V [m3]', 'm_sealed [kg]', 'Q [W]', 't_hold [h]'))
    print("  %s-+-%s-+-%s-+-%s-+-%s"
          % ('-' * 12, '-' * 9, '-' * 13, '-' * 8, '-' * 10))
    for name, s in summary.items():
        print("  %-12s | %9.2f | %13.1f | %8.1f | %10.2f"
              % (name, s['V'], s['m_total'], s['Q_legacy'], s['t_equiv_h']))
    print("")
    print("  Because delta_u and the capacity convention are identical for all")
    print("  three tanks, the equivalent hold time is a property of the fill and")
    print("  venting pressures alone and is therefore the SAME for every")
    print("  aircraft. This is what makes it a valid common design requirement.")
    print("=" * 78)


if __name__ == '__main__':
    main()
