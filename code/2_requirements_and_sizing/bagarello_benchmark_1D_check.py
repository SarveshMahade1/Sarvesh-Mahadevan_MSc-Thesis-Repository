"""
================================================================================
ANALYTICAL PRE-CHECK OF THE BAGARELLO BENCHMARK
1-D steady-state thermal resistance network, external film coefficient derived
Pure Python - TU Delft MSc thesis

Author: Sarvesh Mahadevan
Date: September 2026

--------------------------------------------------------------------------------
WHY THIS SCRIPT EXISTS
--------------------------------------------------------------------------------
Bagarello, Campagna & Benedetti (2026) report, for their reference tank at a
filling ratio of 0.75, a boil-off rate of 12.03 kg/h and a peak liner
temperature of 25.3 K in the dry (ullage) region. Both are RESULTS of their
analysis, not prescribed values: the only temperatures they fix are the
fuselage compartment at 283 K and the liquid hydrogen at 20 K. Everything
between those two sinks is solved for.

That makes them usable as benchmark quantities. This script reproduces them
with a 1-D resistance network before any finite element work is done, for two
reasons:

  1. If the hand calculation cannot get close to 12.03 kg/h and 25.3 K, the
     FEM will not either, and the problem is in the physics or the property
     data rather than in the mesh. Far cheaper to find out here.

  2. Their external convection coefficient is not quoted in the text. Rather
     than assume one, it is DERIVED here from the Churchill-Chu correlation
     for natural convection about a horizontal cylinder. That keeps the
     benchmark self-contained: the only inputs are the two sink temperatures,
     the geometry and published material properties.

--------------------------------------------------------------------------------
FORMULAS USED
--------------------------------------------------------------------------------
 1  Cylindrical conduction resistance   R = ln(r_o/r_i) / (2*pi*k*L)
 2  Spherical conduction resistance     R = (1/r_i - 1/r_o) / (4*pi*k)
 3  Convective surface resistance       R = 1 / (h*A)
 4  Series network                      Q = (T_hot - T_cold) / R_total
 5  Volumetric expansion coefficient    beta = 1/T_film        (ideal gas)
 6  Grashof number                      Gr = g*beta*dT*D^3 / nu^2
 7  Rayleigh number                     Ra = Gr*Pr
 8  Churchill-Chu, horizontal cylinder
        Nu = { 0.60 + 0.387*Ra^(1/6) / [1 + (0.559/Pr)^(9/16)]^(8/27) }^2
 9  Film coefficient                    h = Nu*k_air / D
10  Boil-off rate                       BOR = Q / lambda_v
11  Dry-wall temperature rise           dT = q / h_int

--------------------------------------------------------------------------------
DATA SOURCES
--------------------------------------------------------------------------------
  Geometry, layer thicknesses, T_fus, T_LH2, h_int, material properties:
      Bagarello, Campagna & Benedetti, "Computational thermo-mechanical
      modelling and design-space exploration of cryogenic hydrogen tanks for
      aviation", Aerospace Science and Technology 168 (2026) 110755.
      Tables 1, 2, 4 and Section 3.1.
  Churchill-Chu correlation:
      Churchill & Chu (1975), as given in Bergman & Lavine, Fundamentals of
      Heat and Mass Transfer. Standard textbook correlation.
  Air properties near 283 K: standard tabulated values.

--------------------------------------------------------------------------------
KEY FINDING
--------------------------------------------------------------------------------
The PUR foam carries about 98.5 percent of the total thermal resistance. Both
parameters that had to be assumed here - the CFRP conductivity and the external
film coefficient - therefore have almost no influence, which means the
benchmark rests on Bagarello's own published foam conductivity and thickness
rather than on anything assumed in this script.
================================================================================
"""

import math

# ============================================================================
# INPUTS - Bagarello reference tank (their Table 4)
# ============================================================================

R_liner_in = 2.12          # m, inner radius of the aluminium liner
L_cyl      = 4.24          # m, cylindrical section length
t_Al       = 0.001         # m, liner thickness
t_PUR      = 0.50          # m, polyurethane foam insulation
t_CFRP     = 0.010         # m, outer composite shell

T_fus      = 283.0         # K, fuselage compartment (their Sec. 3.1)
T_LH2      = 20.0          # K, liquid hydrogen
h_int_dry  = 1.46          # W/m2-K, GH2 natural convection (their Sec. 3.1)
h_int_wet  = 1000.0        # W/m2-K, effectively pins the wall to T_LH2

k_Al       = 121.0         # W/m-K   (their Table 1)
k_PUR      = 0.0163        # W/m-K   (their Table 1)
k_CFRP     = 0.60          # W/m-K   through-thickness, representative
                           #          (not quoted in their Table 2; the foam
                           #           dominates so this barely matters - the
                           #           sensitivity is checked at the end)

lambda_v   = 446000.0      # J/kg, latent heat of vaporisation of LH2
phi_f      = 0.75          # filling ratio at their reported operating point

# Their published results, for comparison
BOR_PAPER   = 12.03        # kg/h
T_MAX_PAPER = 25.3         # K

# --- Air properties near 283 K ---
nu_air = 1.43e-5           # m2/s,  kinematic viscosity
k_air  = 0.0248            # W/m-K, thermal conductivity
Pr_air = 0.71              # -
g      = 9.81              # m/s2


# ============================================================================
# GEOMETRY
# ============================================================================

r0 = R_liner_in                     # liner inner
r1 = r0 + t_Al                      # liner outer / foam inner
r2 = r1 + t_PUR                     # foam outer / CFRP inner
r3 = r2 + t_CFRP                    # CFRP outer, exposed to the compartment

D_out = 2.0 * r3

# Their caps are spherical (aspect ratio 1), so the two caps together form one
# complete sphere of the corresponding radius.
A_in_cyl   = 2.0 * math.pi * r0 * L_cyl
A_in_caps  = 4.0 * math.pi * r0**2
A_in_total = A_in_cyl + A_in_caps

A_out_cyl  = 2.0 * math.pi * r3 * L_cyl
A_out_caps = 4.0 * math.pi * r3**2
A_out_total = A_out_cyl + A_out_caps

V_tank = math.pi * r0**2 * L_cyl + (4.0 / 3.0) * math.pi * r0**3

print('=' * 78)
print('ANALYTICAL PRE-CHECK OF THE BAGARELLO BENCHMARK')
print('=' * 78)
print('Geometry')
print('  liner inner radius r0      = %.4f m' % r0)
print('  foam outer radius  r2      = %.4f m' % r2)
print('  CFRP outer radius  r3      = %.4f m   (D_out = %.3f m)' % (r3, D_out))
print('  cylindrical length         = %.3f m' % L_cyl)
print('  tank volume                = %.2f m3   (their Table 4 states 100 m3)'
      % V_tank)
print('  inner surface area         = %.2f m2' % A_in_total)
print('  outer surface area         = %.2f m2' % A_out_total)
print('')


# ============================================================================
# EXTERNAL FILM COEFFICIENT, DERIVED
# ============================================================================

def churchill_chu(dT, D):
    """
    FORMULAS 5 to 9. Natural convection film coefficient for a horizontal
    cylinder in air.

    beta = 1/T_film for an ideal gas; the film temperature is the mean of the
    surface and the free stream. Returns h in W/m2-K.
    """
    if dT <= 0.0:
        return 1e-6
    T_film = T_fus - dT / 2.0
    beta = 1.0 / T_film                                       # FORMULA 5
    Gr = g * beta * dT * D**3 / nu_air**2                     # FORMULA 6
    Ra = Gr * Pr_air                                          # FORMULA 7
    num = 0.387 * Ra**(1.0 / 6.0)
    den = (1.0 + (0.559 / Pr_air)**(9.0 / 16.0))**(8.0 / 27.0)
    Nu = (0.60 + num / den)**2                                # FORMULA 8
    return Nu * k_air / D                                     # FORMULA 9


def resistances(h_ext, h_int):
    """
    FORMULAS 1 to 3. Series network from the compartment air to the fuel.

    Cylinder and cap are treated as parallel paths, each with its own series
    chain, then combined. The foam dominates by orders of magnitude, which is
    checked in the breakdown printed below.
    """
    # Cylindrical path
    R_conv_out_c = 1.0 / (h_ext * A_out_cyl)
    R_cfrp_c = math.log(r3 / r2) / (2.0 * math.pi * k_CFRP * L_cyl)
    R_pur_c  = math.log(r2 / r1) / (2.0 * math.pi * k_PUR * L_cyl)
    R_al_c   = math.log(r1 / r0) / (2.0 * math.pi * k_Al * L_cyl)
    R_conv_in_c = 1.0 / (h_int * A_in_cyl)
    R_cyl = R_conv_out_c + R_cfrp_c + R_pur_c + R_al_c + R_conv_in_c

    # Spherical cap path
    R_conv_out_s = 1.0 / (h_ext * A_out_caps)
    R_cfrp_s = (1.0 / r2 - 1.0 / r3) / (4.0 * math.pi * k_CFRP)
    R_pur_s  = (1.0 / r1 - 1.0 / r2) / (4.0 * math.pi * k_PUR)
    R_al_s   = (1.0 / r0 - 1.0 / r1) / (4.0 * math.pi * k_Al)
    R_conv_in_s = 1.0 / (h_int * A_in_caps)
    R_sph = R_conv_out_s + R_cfrp_s + R_pur_s + R_al_s + R_conv_in_s

    R_total = 1.0 / (1.0 / R_cyl + 1.0 / R_sph)
    parts = dict(conv_out=R_conv_out_c, cfrp=R_cfrp_c, pur=R_pur_c,
                 al=R_al_c, conv_in=R_conv_in_c, cyl=R_cyl, sph=R_sph)
    return R_total, parts


# The outer wall temperature is unknown, and the film coefficient depends on
# it, so iterate. The foam carries almost all the resistance, so the outer
# wall sits very close to the compartment temperature and this converges in a
# handful of passes.
dT_ext = 1.0
for it in range(200):
    h_ext = churchill_chu(dT_ext, D_out)
    R_tot, parts = resistances(h_ext, h_int_wet)
    Q = (T_fus - T_LH2) / R_tot                               # FORMULA 4
    dT_new = Q * (1.0 / (h_ext * A_out_total))
    if abs(dT_new - dT_ext) < 1e-9:
        dT_ext = dT_new
        break
    dT_ext = 0.5 * dT_ext + 0.5 * dT_new

h_ext = churchill_chu(dT_ext, D_out)
R_tot, parts = resistances(h_ext, h_int_wet)
Q = (T_fus - T_LH2) / R_tot
T_ext_wall = T_fus - dT_ext

print('External film coefficient, DERIVED (Churchill-Chu)')
print('  converged in %d iterations' % (it + 1))
print('  outer wall temperature drop dT_ext = %.4f K' % dT_ext)
print('  outer wall temperature    T_ext    = %.3f K' % T_ext_wall)
print('  film coefficient          h_ext    = %.3f W/m2-K' % h_ext)
print('')

print('Resistance breakdown (cylindrical path, K/W)')
tot_c = parts['cyl']
for nm in ['conv_out', 'cfrp', 'pur', 'al', 'conv_in']:
    print('  %-10s %12.6f   (%5.2f %% of the path)'
          % (nm, parts[nm], 100.0 * parts[nm] / tot_c))
print('  %-10s %12.6f' % ('TOTAL', tot_c))
print('')


# ============================================================================
# HEAT LEAK AND BOIL-OFF
# ============================================================================

BOR = Q / lambda_v * 3600.0                                   # FORMULA 10
q_mean = Q / A_in_total

print('=' * 78)
print('RESULT 1: HEAT LEAK AND BOIL-OFF')
print('=' * 78)
print('  total heat leak Q      = %8.1f W' % Q)
print('  mean inner-wall flux q = %8.3f W/m2' % q_mean)
print('  boil-off rate          = %8.3f kg/h' % BOR)
print('  Bagarello report       = %8.3f kg/h' % BOR_PAPER)
print('  difference             = %+7.1f %%'
      % (100.0 * (BOR - BOR_PAPER) / BOR_PAPER))
print('')


# ============================================================================
# DRY-WALL PEAK TEMPERATURE
# ============================================================================

# In the ullage region the liner receives essentially the same flux from the
# foam but must reject it into hydrogen VAPOUR, whose film coefficient is
# three orders of magnitude lower. The wall therefore runs warmer.
#
# FORMULA 11. This is an upper bound: it ignores conduction ALONG the liner
# from the cold wetted region towards the dry region. Because aluminium
# conducts extremely well (k = 121 W/m-K), that lateral path pulls the real
# peak down. The finite element model captures it; this hand calculation
# does not, so the number below should OVERESTIMATE the peak.

dT_dry = q_mean / h_int_dry                                   # FORMULA 11
T_dry_upper = T_LH2 + dT_dry

print('=' * 78)
print('RESULT 2: PEAK DRY-WALL TEMPERATURE')
print('=' * 78)
print('  flux to reject          = %8.3f W/m2' % q_mean)
print('  h_int (dry, GH2)        = %8.3f W/m2-K' % h_int_dry)
print('  wall rise q/h           = %8.3f K' % dT_dry)
print('  upper-bound T_max       = %8.2f K   (no lateral conduction)' % T_dry_upper)
print('  Bagarello report        = %8.2f K' % T_MAX_PAPER)
print('')
if T_dry_upper >= T_MAX_PAPER:
    print('  CONSISTENT: the hand calculation bounds their FEM result from above,')
    print('  which is exactly what is expected, because lateral conduction along')
    print('  the aluminium liner towards the cold wetted region is neglected here')
    print('  and is included in their model. The FEM should land between %.1f K'
          % T_LH2)
    print('  and %.2f K, and their %.1f K sits in that window.'
          % (T_dry_upper, T_MAX_PAPER))
else:
    print('  INCONSISTENT: the hand calculation falls BELOW their reported peak.')
    print('  Lateral conduction can only reduce the peak, so a bounding estimate')
    print('  that is already too low indicates a problem in the property data,')
    print('  the geometry, or the interpretation of their boundary conditions.')
print('')


# ============================================================================
# SENSITIVITY
# ============================================================================

print('=' * 78)
print('SENSITIVITY')
print('=' * 78)
print('k_CFRP is not quoted in their Table 2 and was assumed. Its influence:')
print('  %-14s %10s %10s %12s' % ('k_CFRP', 'Q (W)', 'BOR', 'T_dry_upper'))
for kc in [0.2, 0.4, 0.6, 1.0, 5.0]:
    k_save = k_CFRP
    globals()['k_CFRP'] = kc
    Rt, _ = resistances(h_ext, h_int_wet)
    Qk = (T_fus - T_LH2) / Rt
    print('  %-14.2f %10.1f %10.3f %12.2f'
          % (kc, Qk, Qk / lambda_v * 3600.0,
             T_LH2 + (Qk / A_in_total) / h_int_dry))
    globals()['k_CFRP'] = k_save
print('')
print('h_ext was derived rather than assumed. Its influence:')
print('  %-14s %10s %10s' % ('h_ext', 'Q (W)', 'BOR'))
for hx in [2.0, h_ext, 10.0, 50.0]:
    Rt, _ = resistances(hx, h_int_wet)
    Qh = (T_fus - T_LH2) / Rt
    print('  %-14.3f %10.1f %10.3f' % (hx, Qh, Qh / lambda_v * 3600.0))
print('')
print('If both sensitivities are weak, the foam controls the problem and the')
print('benchmark rests on their published foam conductivity and thickness')
print('rather than on anything assumed here.')
print('=' * 78)
