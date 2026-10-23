# Analysis guide

This guide describes the analysis of the thesis as a method that can be reused, extended or rebuilt
in another finite element code. For each step it gives the purpose, the physics and equations, the
assumptions, the inputs and outputs, where the step is implemented, and how to add something at that
step. The numbers quoted are those of the thesis (mostly the E190) and come from
`results/processed/`.

Load cases are named by their physical condition. The codes in brackets are used only in file and
job names:

| Name | Code | Load factors (nx forward, ny lateral, nz vertical) [g] | Basis |
|---|---|---|---|
| Baseline condition | `Baseline` | (0, 0, 1) | 1 g, ground or level flight |
| Symmetric manoeuvre | `LC1_Maneuver` | (0, 0, 2.5) | CS 25.337, limit |
| Combined manoeuvre | `LC3_Combined` | (1.5, 0, 2.5) | constructed case, not CS-25 |
| Emergency landing | `LC7_Emergency` | (9, 3, 6) | CS 25.561 factors applied together, ultimate |

## Contents

1. [The chain at a glance](#1-the-chain-at-a-glance)
2. [Tank geometry from the mission](#2-tank-geometry-from-the-mission)
3. [Heat leak budget from a no-vent hold](#3-heat-leak-budget-from-a-no-vent-hold)
4. [Insulation and outer jacket](#4-insulation-and-outer-jacket)
5. [Closed-form starting wall](#5-closed-form-starting-wall)
6. [Finite element model of one design point](#6-finite-element-model-of-one-design-point)
7. [Wall sizing loop](#7-wall-sizing-loop)
8. [Design matrix at two pressures](#8-design-matrix-at-two-pressures)
9. [Mesh convergence](#9-mesh-convergence)
10. [Post-processing and verification](#10-post-processing-and-verification)
11. [Assumptions and limits](#11-assumptions-and-limits)
12. [Recipes for extensions](#12-recipes-for-extensions)
13. [Rebuilding the analysis in another code](#13-rebuilding-the-analysis-in-another-code)

## 1. The chain at a glance

```
design_basis.json
   |
   v
[2] mission -> geometry (R, L, V)          tank_geometry_from_mission.py
[3] hold time -> heat leak budget Q        boiloff_requirement_derivation.py
[4] Q -> insulation, outer jacket          insulation_sizing.py, outer_jacket_closed_form.py
[5] closed-form starting wall              sizing_closed_form.py
   |
   v
[6] FE design point (load case, fill)      E190_MASTER_thermo_structural.py: run_case, extract
      free surface -> heat transfer -> static -> extraction and checks
   |
[7] sizing loop on the governing point     stage_sizing        -> sized_walls.csv
[8] matrix at design and stability p       stage_matrix        -> production, stability results
[9] mesh convergence                       stage_convergence   -> convergence results
   |
   v
[10] merge, closed-form checks, figures    final_numbers.py, hand_calcs.py, review_checks.py, ...
```

Every input of every step is read from `code/design_basis.json` (through `code/design_basis.py`).
Steps 2 to 5 and 10 run in plain Python; `python code/reproduce_all.py` runs them in order. Steps 6
to 9 need Abaqus; `E190_DRY_RUN=1` checks their inputs without a licence.

## 2. Tank geometry from the mission

**Purpose.** Turn the kerosene load of a reference aircraft into the inner-vessel geometry.

**Equations.**

```
m_LH2 = m_JetA * LHV_JetA / LHV_LH2          energy equivalence (43 and 120 MJ/kg)
V     = m_LH2 / (rho_LH2 * phi_mission)      volume at the mission fill (70.8 kg/m3, 0.90)
R     = 0.8648 * D_fuselage / 2              radius from the fuselage diameter
h_cap = R / (a/b)                            semi-ellipsoidal caps, a/b = 1.6
L     = (V - 4/3 pi R^2 h_cap) / (pi R^2)    barrel length
```

E190: m_LH2 = 5103.7 kg, V = 80.1 m3, R = 1.3015 m, L = 13.97 m, total length 15.59 m.

**Assumptions.** Equal stored energy (no credit for different propulsive efficiency); one tank per
aircraft; a fixed ratio of tank to fuselage diameter; caps of fixed aspect ratio.

**Inputs and outputs.** `mission` and `aircraft` sections of the design basis;
`results/processed/tank_geometry.json`.

**Code.** `design_basis.mission_geometry()`, `design_basis.geometry()`,
`2_requirements_and_sizing/tank_geometry_from_mission.py`.

**To extend.** Another packaging rule (two tanks, a different diameter ratio): change
`mission_geometry()`. Another cap shape: `cap_aspect_ratio_a_over_b` (1.0 gives hemispheres; the
model handles that case). A record with `R_mid_m` and `L_barrel_m` set uses that geometry directly.

## 3. Heat leak budget from a no-vent hold

**Purpose.** Derive the allowable heat leak from an operational requirement: the tank must stay
sealed for a required time without reaching the relief pressure.

**Physics.** A sealed tank holds a saturated two-phase mixture at fixed mass and volume, so its bulk
density is constant. Heat that leaks in raises the internal energy along a line of constant density,
and the saturation pressure rises. Between the fill pressure and the venting pressure:

```
Q * t_hold = m_total * (u2 - u1)        constant-volume energy balance
Q_allow    = m_total * (u2 - u1) / t_hold
```

with u1, u2 from the parahydrogen equation of state (CoolProp) at 1.2 bar and 1.448 bar, the tank
sealed at a fill of 0.97. Delta u = 6809 J/kg is the same for every tank, so the hold time for a
given boil-off rate does not depend on size: the 0.1 %/h rate corresponds to a hold of 14.68 h for
all three aircraft. E190: Q = 702.6 W (A320 919.3 W, A350-900 6583.2 W).

**Assumptions.** Thermodynamic equilibrium in the tank (no stratification); steady heat leak.

**Inputs and outputs.** `pressures`, `heat_budget`, aircraft geometry;
`results/processed/heat_budget.json`.

**Code.** `2_requirements_and_sizing/boiloff_requirement_derivation.py`; the sensitivity to the hold
time is in `boiloff_holdtime_sensitivity.py`.

**To extend.** Change `hold_time_h` or the pressures, and set the aircraft's `Q_budget_W` to `null`
so that the rule of the `heat_budget` section is applied.

## 4. Insulation and outer jacket

**Purpose.** Choose an insulation that meets the budget and size the vacuum jacket.

**Physics.** A series resistance network (conduction through each layer, convection outside)
between the compartment at 293 K and the liquid at 20.3 K; three architectures (vacuum with MLI,
foam only, vacuum only by radiation). The outer jacket carries one atmosphere of external pressure:
Bresse collapse for a plain cylinder, Windenburg and Trilling for a ring-stiffened one, with a factor
of 1.5.

**Link to the finite element model.** The insulation is not modelled. It enters only as a uniform
heat flux on the outer face of the vessel, q = Q / A_shell (E190: 5.39 W/m2).

**Code.** `insulation_sizing.py` (writes `results/processed/insulation.json`),
`outer_jacket_closed_form.py`.

## 5. Closed-form starting wall

**Purpose.** A starting wall for the finite element sizing, from membrane and beam theory.

**Equations.** Axial membrane stress in the barrel, tension positive, with rings at 0, L/2 and L:
pressure pR/2t, axial inertia carried by the mid anchor, and beam bending of two equal spans; the
minimum is compared with the NASA SP-8007 allowable

```
sigma_cl = E t / (R sqrt(3 (1 - nu^2)))
gamma    = 1 - 0.901 (1 - exp(-sqrt(R/t) / 16))
sigma_cr = gamma * sigma_cl          (pressure stabilisation not credited)
```

and the von Mises stress at the design pressure with the hydrostatic head with F_ty / 1.5.

**Code.** `2_requirements_and_sizing/sizing_closed_form.py`; the results are the
`t_wall_initial_m` values of the design basis.

## 6. Finite element model of one design point

A design point is one load case at one fill level, at one pressure. It is solved as a steady-state
heat transfer analysis followed by a static analysis that reads its temperatures (sequential
coupling: the deformation does not change the thermal problem at these strains). The whole design
point is the function `run_case`; its results come from `extract`.

### 6.1 Geometry and mesh

- Inner vessel only, as a shell of revolution: aft cap, barrel, forward cap, revolved about the
  global X axis (X forward, Y up, Z lateral). Mid-surface radius R, barrel 0 <= x <= L.
- Thin shell: R/t is 450 to 1100 at the sized walls.
- Two models with the same mesh: diffusive shells DS4/DS3 for heat transfer, S4R/S3 for the static
  analysis. Structured quadrilateral mesh, with a free mesh as fallback.
- Element size from the shell bending length: seed = 0.5 sqrt(R t), two elements per bending
  length, from the convergence study (Section 9). E190: sqrt(R t) = 41.0 mm, seed 20.5 mm.

Code: `build_geometry_and_mesh`. To add a feature (a sump, a nozzle, a stiffening ring), add it to
the sketch or partition the part here; a feature that breaks axial symmetry needs a free mesh.

### 6.2 Liquid free surface

**Physics.** A liquid under steady acceleration reaches relative equilibrium: it behaves as under
gravity, with "down" along the resultant load factor. In model axes the direction is

```
g_hat = (nx, -nz, -ny) / |n|
```

The free surface is the plane normal to g_hat that leaves the fill fraction phi of the tank volume
below it. It has no closed form for a tilted plane in a capped tank, so the volume is sampled with
200 000 uniform points and the plane is the (1 - phi) quantile of their coordinate along g_hat. The
depth of a point is its distance below that plane along g_hat.

E190 at a fill level of 0.90:

| Load case | Tilt from vertical | |n| [g] | Maximum depth [m] | p_h / p at the deepest point |
|---|---|---|---|---|
| Baseline condition | 0 | 1.00 | 2.19 | 0.010 |
| Symmetric manoeuvre | 0 | 2.50 | 2.19 | 0.025 |
| Combined manoeuvre | 31.0 | 2.92 | 7.91 | 0.107 |
| Emergency landing | 57.7 | 11.23 | 11.45 | 0.595 |

A load factor tilted towards the tank axis makes the liquid column much deeper than the tank is
tall, which is why the emergency landing governs.

**The same surface drives both analyses**: it decides which wall is wetted (thermal) and the depth of
liquid (structural). It is computed once per design point.

Code: `effective_down`, `free_surface_offset`, `elem_depth_coord`. A sloshing or transient liquid
would replace these functions.

### 6.3 Thermal problem

**Equation** (steady state, per unit area of a thin shell):

```
div(k t grad T) + q_MLI - h (T - T_fluid) = 0
```

in-plane conduction in the wall, the heat arriving through the insulation, and convection to the
fluid inside. The inner surface is split by the free surface:

| Region | h [W/m2K] | T_fluid [K] | Source |
|---|---|---|---|
| Wetted | 1000 | 20.3 (liquid) | holds the wall at the liquid temperature |
| Dry | 1.46 | 20.3 (vapour) | natural convection of the vapour, Bagarello et al. Sec. 3.1 |

Only the two sink temperatures are prescribed; every wall temperature is solved for. Far from the
waterline the dry wall tends to `dT_inf = q / h_dry` (E190 3.69 K); the transition from the wetted
wall takes place over the conduction length `L_c = sqrt(k t / h_dry)` (E190 0.32 m). Near the
waterline heat flows along the wall into the liquid; the peak wall flux is 1280 to 1800 W/m2.
The resulting wetted-to-dry differences are 3.4 to 7.4 K across the three tanks.

Code: first half of `run_case`; `split_wetted_dry`. To add a heat source, a different film model or
a transient (cooldown, filling), change the step and the loads here; a transient needs an initial
temperature and a time-dependent sink.

### 6.4 Structural problem

**Constitutive law.** Linear thermo-elasticity, `sigma = C : (eps - alpha (T - T_ref))`, aluminium
2219-T87 at room-temperature properties (E = 73.1 GPa, nu = 0.33, alpha = 22.3e-6 1/K).

**Stress-free state.** T_ref = 20.3 K, the cold operating state. The cooldown from 293 K and a crash
load do not act together: the tank cools on the ground, contracting freely on sliding supports, and
is cold and settled when a flight load occurs. The only thermal load is therefore the wetted-to-dry
difference. (With T_ref = 293 K the fully restrained cooldown alone would be 663 MPa.) The cooldown
is a separate ground case.

**Temperatures.** The nodal temperatures of the heat transfer solution are grouped into 24 bins and
applied with a zero through-thickness gradient (the real drop through the wall is of order 1e-4 K).

**Pressure.** On the inner face, in 64 bins of depth:

```
p = p_gauge + rho_LH2 * |n| g * d
```

**Inertia.** The shell mass is loaded by the acceleration (+nx, -nz, -ny) g (CS-25 factors are named
by the direction in which the inertia force acts on the item). The inertia of the liquid is carried
by the scaled hydrostatic pressure, not by a body force.

Code: second half of `run_case`. A new load (a point mass, an external pressure, a pipe load) is
added here.

### 6.5 Supports

Three circumferential rings, at the two cap-to-barrel junctions and at mid-barrel (Bagarello et al.
Sec. 4.4.1, from H2FLY testbed hardware), as an anchor-and-slide scheme:

| Ring | Held | Free |
|---|---|---|
| Aft junction (x = 0) | radial and circumferential (u2, u3) | axial |
| Mid-barrel (x = L/2) | all translations (u1, u2, u3), the axial anchor | none |
| Forward junction (x = L) | radial and circumferential (u2, u3) | axial |

A CS-25 load factor presupposes a load path: the stress arises from the tank being held. The ring
reactions are the loads the support structure must carry (E190, emergency landing at 0.90: 614 kN in
total). A statically determinate 3-2-1 restraint or inertia relief cannot represent this; both are
kept in the script as alternatives for other studies.

Code: `apply_ring_constraints`. Flexible supports (springs or connector elements with the axial and
flexural stiffness of real struts) would replace the displacement conditions here.

### 6.6 Extraction and criteria

**Strength.** The rings are rigid line constraints, so the stress at them is singular and grows with
mesh refinement. The peak von Mises stress is taken outside 0.5 m of each ring (Saint-Venant). The
shell edge disturbance decays over pi/beta, with `beta = [3 (1 - nu^2) / (R t)^2]^(1/4)`: 0.10,
0.13 and 0.30 m for the sized walls, so 0.5 m is 5.0, 3.9 and 1.7 decay lengths. Reserve factor:

```
RF_yield = (F_ty / 1.5) / sigma_vM = 262 MPa / sigma_vM     (requirement RF >= 1)
```

**Stability.** The minimum axial membrane stress (mean of the two section points) in the far-field
barrel, at the stability pressure, against SP-8007:

```
RF_buckling = (sigma_cr / 1.5) / (-sigma_axial,min)          (requirement RF >= 1)
```

The results files record `buckling_fos = sigma_cr / (-sigma_axial,min)`, which must be at least 1.5;
RF_buckling = buckling_fos / 1.5. The margin of 1.5 in both allowables is a project design margin,
not the CS-25 factor of safety. The buckling check is a criterion, not an eigenvalue analysis.

Code: `extract`. A new criterion (fatigue, a fitting stress, a cap buckling check) is computed here,
added to the returned dictionary and its key to `EXTRA_KEYS`.

### 6.7 Verification gates (every design point)

| Check | What it tests | Result column |
|---|---|---|
| Load balance | total reaction against (m_shell + m_liquid) n g; a uniform pressure must cancel on the closed surface, so an imbalance exposes wrong element normals | `balance_ratio` (worst E190 case 1.2 %) |
| Stress orientation | mean hoop membrane stress above mean axial one under pressure, so S11 is axial | `orientation_ok` |
| Membrane probe | dry crown element: both section points at pR/t and pR/2t, no spurious bending | `probe_*` |
| Ovality | out-of-roundness between the rings | `ovality_mm`, `ur_*` |
| Ring reactions | load split between the rings | `RF_aft_kN`, `RF_mid_kN`, `RF_fwd_kN` |
| Cap compression | hoop compression near the equator of an a/b = 1.6 head | `cap_*` |

A sizing run stops if the balance is off by more than 2 % or the orientation check fails.

## 7. Wall sizing loop

The wall is sized on the governing design point (emergency landing at a fill of 0.90) at both
pressures, by whichever requirement needs the most material:

- buckling: RF_buckling >= 1 (buckling_fos >= 1.5) at the stability pressure;
- yield: RF_yield >= 1.0 at the design pressure;
- optionally a minimum gauge (`E190_MIN_GAUGE_MM`).

The buckling reserve scales about as t^2 and the yield reserve about as t, which gives the next wall;
steps are limited to a factor of two, and once a passing and a failing wall are known the loop
bisects. A wall is accepted when the governing requirement is met with at most about 6 % surplus.
The other candidate design points are then checked for buckling at the final wall.

| Aircraft | Sized wall [mm] | Governed by | Lowest RF_yield | Lowest RF_buckling |
|---|---|---|---|---|
| E190 | 1.289 | yield | 1.020 | 1.105 |
| A320 | 1.572 | yield | 1.018 | no axial compression |
| A350-900 | 5.674 | buckling | 1.798 | 1.060 |

Code: `stage_sizing`; record `results/sized_walls.csv`; iterations in
`results/wall_sizing_iterations/`.

## 8. Design matrix at two pressures

Every load case at every fill level (0.10, 0.30, 0.50, 0.75, 0.90), at two pressure differences:

| Pressure | Value (gauge) | Meaning | Used for |
|---|---|---|---|
| Design | 1.5 bar | maximum design difference, above venting (1.448 bar abs) | strength |
| Stability | 0.18675 bar | fill pressure 1.2 bar abs minus sea-level ambient, after a loss of vacuum | buckling |

Pressure stabilises the shell, so the lowest difference is critical for buckling. The thermal field
does not depend on the pressure and is reused between the two matrices.

Code: `stage_matrix`; results in `results/production_1p5bar/` and `results/stability_0p187bar/`.

## 9. Mesh convergence

The governing design point at seeds of 1, 0.5 and 0.25 sqrt(R t), at both pressures. At 0.5 the
far-field stress is within 1 % of the finest mesh for all three tanks (the coarsest overestimates it
by 5 to 11 %), while the peak including the rings keeps rising: the singularity that the far-field
extraction excludes. E190: far field 271.5, 256.3, 256.8 MPa; with the rings 271.5, 256.3,
295.1 MPa; total reaction constant at 614.2 kN.

Code: `stage_convergence`; results in `results/mesh_convergence/`.

## 10. Post-processing and verification

| Script | What it does |
|---|---|
| `final_numbers.py` | Merges every results file at the sized walls and stage pressures into `final_numbers.json`, the single source of the reported numbers |
| `hand_calcs.py` | Closed-form checks of the model: membrane stress, equilibrium, head, beam bending, forward inertia, thermal asymptote q/h, conduction length, bowing |
| `review_checks.py` | Further checks: depth and magnitude split of the head, the emergency landing against a pure 9 g forward load, 14 CFR 25.963(d), cryogenic properties, the 24 h hold, and a fit of the far-field peak against the local membrane stress |
| `bagarello_benchmark_1D_check.py` | One-dimensional check of the benchmark heat flux |

Two results show the physics the model captures:

- **The far-field peak follows the local pressure.** Over the twenty E190 production cases,
  `sigma_vM ~ a (p + p_h) R / t + b` with a = 1.055 and b = 1.5 MPa (largest residual 0.9 MPa): the
  stress is the membrane response to the local pressure, the thermal part is small.
- **Benchmark.** The model of the Bagarello et al. reference tank gives a peak dry-wall temperature
  of 25.1 K at fill 0.75 against their 25.3 K.

## 11. Assumptions and limits

| Assumption | Consequence | Where to change it |
|---|---|---|
| Steady state, quasi-static liquid | no sloshing, no cooldown transient | 6.2, 6.3 |
| Linear elastic, small deflection | no plasticity, no geometric nonlinearity | 6.4 |
| Room-temperature properties | conservative yield; thermal strain overestimated near 20 K (check: `E190_ALPHA`, `E190_KAPPA`) | `material` |
| Rigid rings | local fitting stresses excluded (doubler or wear plate in a detailed design) | 6.5 |
| Insulation as a uniform flux | no local heat paths at supports or pipes | 6.3 |
| SP-8007 criterion | no eigenvalue or post-buckling analysis; pressure stabilisation not credited | 6.6 |
| Load factors applied together in the emergency landing | more severe than CS 25.561, where they act separately | `load_cases` |

## 12. Recipes for extensions

**A new aircraft.** Add a record to `aircraft` in the design basis with `jet_a_mass_kg`,
`fuselage_diameter_m`, a starting wall and `null` geometry and budget; run
`tank_geometry_from_mission.py` and `boiloff_requirement_derivation.py`; check with
`E190_DRY_RUN=1 E190_AIRCRAFT=<key>`; run `E190_RUN_MODE=FINAL`; copy the results into `results/`
and run `python code/reproduce_all.py --no-figures`.

**A new load case.** Add a record to `load_cases` with nx, ny, nz. The free surface, pressure and
inertia follow from the vector. If it should size the wall, set `analysis.governing_case`.

**Another material.** Edit `material` and `criteria`. For temperature-dependent properties, replace
the constant tables in `build_geometry_and_mesh` by temperature tables (Abaqus accepts them
directly).

**Flexible supports.** In `apply_ring_constraints`, replace the displacement conditions by springs
or connector elements to ground with the stiffness of the support structure; keep the reactions in
`extract`.

**A transient analysis.** Change `HeatTransferStep` to `response=TRANSIENT` with an initial
temperature, and keep one temperature field per time step for the static analysis.

**A new result or criterion.** Compute it in `extract`, add it to the returned dictionary and to
`EXTRA_KEYS`; it then appears in every results file. To size against it, add its required thickness
ratio to the requirement `g` in `stage_sizing`.

**Using the results without Abaqus.** All results files are CSV with documented columns (README,
"Results files: columns"), and `results/processed/final_numbers.json` holds every reported number.
A new comparison or figure can start from those files.

## 13. Rebuilding the analysis in another code

The method does not depend on Abaqus. A reimplementation needs these choices, in this order:

1. Shell of revolution of the inner vessel, R and L from the mission, caps a/b = 1.6.
2. Element size 0.5 sqrt(R t); shell elements with conduction for the thermal model.
3. Free surface normal to (nx, -nz, -ny) in X forward, Y up, Z lateral, positioned by volume.
4. Steady heat transfer: wetted h = 1000, dry h = 1.46 W/m2K, both sinks at 20.3 K, uniform flux
   Q / A on the outer face.
5. Static analysis: stress-free at 20.3 K; the temperature field with no through-thickness
   gradient; pressure p_gauge + rho |n| g d on the inner face; shell inertia (+nx, -nz, -ny) g.
6. Three rings: aft and forward junctions held radially and circumferentially, mid-barrel fully held.
7. Far-field extraction outside 0.5 m of the rings; yield RF with 262 MPa; buckling RF with SP-8007
   at the stability pressure; load balance within 2 %.
8. Sizing on the emergency landing at 0.90; matrix of four load cases by five fills at 1.5 and
   0.18675 bar; convergence at 1, 0.5 and 0.25 sqrt(R t).

The results files of this repository are the reference solution for such a reimplementation.
