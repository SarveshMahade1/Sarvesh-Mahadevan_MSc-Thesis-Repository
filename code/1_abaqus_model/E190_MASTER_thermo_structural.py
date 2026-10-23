"""
================================================================================
LH2 INNER VESSEL: SEQUENTIALLY COUPLED THERMO-STRUCTURAL ANALYSIS
E190_MASTER_thermo_structural.py

Abaqus/CAE Python script of the MSc thesis "Digitalized Workflow of a Liquid
Hydrogen Tank Integrating Thermo-Structural Analysis", TU Delft.
Author: Sarvesh Mahadevan

Despite the file name, one script serves every aircraft of the design basis
(in the thesis: E190, A320 and A350-900). Only the geometry, the wall
thickness and the heat leak budget differ between them; every boundary
condition, load case, criterion and extraction step is shared.

--------------------------------------------------------------------------------
WHAT THE SCRIPT DOES
--------------------------------------------------------------------------------
For each design point (one load case at one fill level) the inner vessel is
built as a revolved shell and two analyses are solved in sequence:

  STEP 1  Steady-state heat transfer. The wall below the liquid free surface
          exchanges heat with the liquid (film coefficient h_wetted), the wall
          above it with the ullage gas (h_dry). The heat leak through the
          insulation enters as a uniform flux on the outer face.
          Result: the nodal temperature field T(x).

  STEP 2  Static structural analysis. Loads: the temperature field of step 1,
          the internal pressure plus the hydrostatic head of the liquid, and
          the inertia of the shell. Supports: three circumferential rings.

The liquid free surface is normal to the resultant load factor (relative
equilibrium of a liquid under steady acceleration). The load case therefore
changes both the wetted area, which drives the thermal problem, and the depth
of liquid, which drives the structural problem. This is how the fill level
couples the two analyses.

From each design point the script extracts the far-field von Mises stress,
the axial membrane stress used by the buckling criterion, temperatures,
displacements, support reactions and a set of verification checks, and writes
them to a results file.

The run modes build the workflow of the thesis on this core:

  SIZING       iterate the wall until the yield or the buckling criterion is
               just met in the governing case; records sized_walls.csv
  PRODUCTION   every load case at every fill level, at the design pressure
  STABILITY    the same matrix at the stability pressure
  CONVERGENCE  the governing case over a sweep of mesh densities
  FINAL        SIZING, PRODUCTION, STABILITY and CONVERGENCE in one session

--------------------------------------------------------------------------------
HOW TO RUN
--------------------------------------------------------------------------------
From an empty working folder, which receives the .odb and results files:

    abaqus cae noGUI=<repository>/code/1_abaqus_model/E190_MASTER_thermo_structural.py

The run is controlled by code/design_basis.json and by environment variables,
so nothing in this file needs editing. The most used variables:

    E190_AIRCRAFT     E190 | A320 | A350 | any aircraft added to the design basis
    E190_RUN_MODE     FINAL | SIZING | PRODUCTION | STABILITY | CONVERGENCE
    E190_STAGES       subset of sizing,production,stability,convergence (FINAL)
    E190_DRY_RUN=1    print the resolved inputs and first-order checks, then
                      stop. Runs in plain Python 3, without an Abaqus licence:
                      E190_DRY_RUN=1 python E190_MASTER_thermo_structural.py
    LH2_DESIGN_BASIS  path of an alternative design basis file

The others (E190_ONLY_LC, E190_ONLY_FILL, E190_T_WALL, E190_P_GAUGE,
E190_BENCHMARK, E190_WORKDIR, E190_ALPHA, E190_KAPPA, E190_MIN_GAUGE_MM,
E190_CONV_FRACTIONS) are listed in the README and explained where they are
read below.

--------------------------------------------------------------------------------
INPUTS AND OUTPUTS
--------------------------------------------------------------------------------
Inputs
  code/design_basis.json  aircraft geometry and heat leak budget, material,
                          margins, pressures, thermal boundary conditions,
                          load cases, fill levels, mesh and extraction rules
                          (read through code/design_basis.py)
  sized_walls.csv         sized wall of each aircraft, from the working
                          folder or from results/ of the repository

Outputs (in the working folder; <AC> is the aircraft code)
  <AC>_production_results.csv           design pressure matrix
  <AC>_stability_results.csv            stability pressure matrix
  <AC>_convergence_pdesign_results.csv  mesh study, design pressure
  <AC>_convergence_pmin_results.csv     mesh study, stability pressure
  sizing_fe/<AC>/<AC>_sizing_results.csv, sized_walls.csv
  HT_<tag>.odb and ST_<tag>.odb         heat transfer and static databases
  The columns of the results files are described in the README.

--------------------------------------------------------------------------------
CONVENTIONS
--------------------------------------------------------------------------------
  Axes          Global X forward along the tank axis, Y up, Z lateral
                (right-handed, as the Abaqus triad). The tank runs from the aft
                pole at x = -h_cap to the forward pole at x = L_barrel + h_cap;
                the barrel occupies 0 <= x <= L_barrel.
  Load factors  n = (nx forward, ny lateral, nz vertical) in g, named as in
                CS-25 by the direction in which the inertia force acts on the
                tank. In model axes the acceleration is (+nx, -nz, -ny) g.
  Shell faces   side2 is the inner (wetted) face, side1 the outer face.
  Units         SI throughout (m, kg, s, K, Pa, W). Printed values use mm,
                MPa, kN and bar where stated.
  Load-case codes, used in job and file names (labels from the design basis):
                Baseline        Baseline condition   (nx, ny, nz) = (0, 0, 1)
                LC1_Maneuver    Symmetric manoeuvre  (0, 0, 2.5)
                LC3_Combined    Combined manoeuvre   (1.5, 0, 2.5)
                LC7_Emergency   Emergency landing    (9, 3, 6)

--------------------------------------------------------------------------------
SECTION MAP
--------------------------------------------------------------------------------
  0  Imports, design basis, working folder
  1  USER INPUTS       run mode, aircraft, material, criteria, pressures,
                       thermal boundary conditions, load cases, fill levels,
                       mesh rules (values from the design basis)
  2  DERIVED GEOMETRY  cap height, volume, exact shell area, insulation flux
  3  HELPER FUNCTIONS
       effective_down, free_surface_offset  free surface for any load case
       build_geometry_and_mesh              revolved shell, material, mesh
       apply_ring_constraints               the three ring supports
       split_wetted_dry                     wetted and dry inner surfaces
       job_succeeded                        job completion read from the .log
       Alternative models kept for future studies (not called):
       fill_surface_height, elem_centroid_z, apply_inertia_relief,
       apply_rigid_body_restraint
  4  THE ANALYSIS
       run_case   one design point: heat transfer, then static
       extract    results and verification checks of one design point
  5  DRIVER      stage_sizing, stage_matrix, stage_convergence, the dry run
                 and the selection of the run modes
  6  OUTPUT      results file and summary of a single-mode run

Where to add something
  aircraft, material, load case, fill level   design_basis.json only
  a boundary condition or a load              run_case
  a result column                             extract (returned dict) and EXTRA_KEYS
  a sizing criterion                          extract, then stage_sizing
docs/analysis_guide.md describes every step of the method, its equations and
assumptions, and how to extend it. Development history: DEVELOPMENT_NOTES.md
in this folder.

--------------------------------------------------------------------------------
METHODOLOGY AND SOURCES
--------------------------------------------------------------------------------
The modelling process follows the published approach of

    Bagarello, S., Campagna, D. and Benedetti, I., "Computational
    thermo-mechanical modelling and design-space exploration of cryogenic
    hydrogen tanks for aviation", Aerospace Science and Technology 168 (2026)
    110755.

supplemented on the inner-vessel treatment by

    Biancotto, A., "Design and analysis of LH2 tank structures for aircraft
    retrofit applications", MSc thesis, Delft University of Technology, 2024.

Each step is tagged "PROCESS n" in the code where it is implemented.

  1  Prescribe only the two SINK temperatures, of the compartment outside the
     tank and of the fuel inside it. Every wall temperature is solved for.
                                                        [Bagarello Sec. 3.1]
  2  Steady-state, linear thermo-elastic analysis.        [Bagarello Sec. 3]
  3  The fill level partitions the inner surface into a wetted and a dry
     region.                                            [Bagarello Sec. 3.1]
  4  Two film coefficients that differ by orders of magnitude: the wetted
     wall is held near the fuel temperature, the dry wall rejects heat to the
     vapour by natural convection and runs warmer. This makes the fill level
     a thermal variable and couples the thermal and structural problems.
                                                        [Bagarello Sec. 3.1]
  5  Thermal problem first, then the nodal temperature field carried into the
     structural problem through sigma = C : (eps - alpha*dT).
                                                         [Bagarello Eq. (6)]
  6  Supports: three circumferential rings at the two cap-to-barrel
     junctions and at mid-barrel, taken from H2FLY testbed hardware.
                                                      [Bagarello Sec. 4.4.1]
  7  Mesh convergence against a physical output.          [Bagarello Fig. 4]
  8  Failure criterion matched to the material class: von Mises for an
     isotropic metal, with a project design margin of 1.5 (not the CS-25
     factor of safety; thesis Section 3.6), plus an axial buckling criterion
     (NASA SP-8007).                                    [Bagarello Sec. 4.4]
  9  Realism added in increments whose effects can be attributed.
                                                           [Bagarello Sec. 4]
 10  The inner vessel as the design subject, with dT and dp as the primary
     loads. The outer jacket, insulation and struts are not modelled; the
     insulation enters only as the heat flux on the outer face.
                                                     [Biancotto Sec. 5.8.2]
 11  Supports characterised by axial and flexural stiffness rather than as
     rigid restraints, where they are introduced.     [Biancotto Sec. 5.8.1]

Extensions made in this thesis, which neither source performs:

  A  A CS-25 load matrix coupled to the fill level. Bagarello et al. apply
     thermal load and internal pressure only; here the CS 25.337 manoeuvre
     and CS 25.561 emergency landing factors act on the vessel wall through
     the hydrostatic head and the shell inertia, at every fill level.
  B  Cross-aircraft consistency. One sizing rule and one method applied to a
     regional, a narrow-body and a wide-body aircraft.

The thermal and structural problems are solved sequentially (PROCESS 5)
because the mechanical response does not feed back into the thermal problem
at these strain levels.

--------------------------------------------------------------------------------
FORMULAS (tagged "FORMULA n" at each point of use)
--------------------------------------------------------------------------------
  1  Tank volume           V = pi*R^2*L_barrel + (4/3)*pi*R^2*h_cap
  2  Cap height            h_cap = R / (a/b)
  3  Free surface          plane normal to the load factor that leaves the
                           fraction phi_f of the tank volume below it
                           (horizontal special case: theta - sin(theta) = 2*pi*phi_f)
  4  Convective flux       q = h * (T_wall - T_fluid)
  5  Insulation heat flux  q_MLI = Q_budget / A_shell
  6  Pressure              p = p_gauge + rho_LH2 * |n|*g * depth
  7  Inertial body force   a = n * g on the shell mass
  8  Thermal strain        eps_th = alpha * (T - T_ref)
  9  von Mises stress      sig_eq = sqrt(0.5*[(s1-s2)^2+(s2-s3)^2+(s3-s1)^2])
 10  Allowable stress      sig_allow = F_ty / SF

--------------------------------------------------------------------------------
DATA SOURCES (values in code/design_basis.json)
--------------------------------------------------------------------------------
  Al 2219-T87 properties   thesis Table 3.3 (MIL-HDBK-5J), cross-checked
                           against Bagarello Table 1
  T_LH2 = 20.3 K           Bagarello Sec. 3.1
  h_dry = 1.46 W/m2K       Bagarello Sec. 3.1, their ref. [45]: natural
                           convection of hydrogen vapour at the dry wall
  SF = 1.5                 project design margin (thesis Section 3.6), not
                           the CS 25.303 factor of safety
  Load factors             CS 25.337 and CS 25.561
  Q_budget                 boiloff_requirement_derivation.py, from the
                           14.68 h no-vent hold time (E190: 702.6 W)
  Geometry                 sizing rule of thesis Section 3.1, R = 0.8648 D/2
                           (code/2_requirements_and_sizing/tank_geometry_from_mission.py)

--------------------------------------------------------------------------------
BENCHMARK
--------------------------------------------------------------------------------
Bagarello et al. report for their reference tank at a fill level of 0.75 a
wetted liner at 20 K and a dry region peaking at 25.3 K. E190_BENCHMARK=1
runs their geometry and operating point instead of an aircraft, and the
summary compares the predicted peak wall temperature with 25.3 K.
================================================================================
"""

# ============================================================================
# 0. IMPORTS, DESIGN BASIS, WORKING FOLDER
# ============================================================================

# DRY RUN. With E190_DRY_RUN=1 the script resolves every input, prints them
# with the derived quantities, and stops before building a model. The Abaqus
# imports are then optional, so it runs in any Python 3 as well as in Abaqus
# and a new study can be checked without a licence.
import os as _os_dry
DRY_RUN = _os_dry.environ.get('E190_DRY_RUN', '').strip() not in ('', '0')
try:
    from abaqus import *
    from abaqusConstants import *
except ImportError:
    if not DRY_RUN:
        raise

# THESE STAR IMPORTS ARE REQUIRED FOR HEADLESS RUNS.
# Abaqus attaches most Model methods only when the module that defines them is
# imported. The CAE GUI loads every module at startup; a noGUI run loads only
# the kernel, and a missing module shows up as, for example,
#     AttributeError: 'Model' object has no attribute 'FieldOutputRequest'
# The set below is the one Abaqus writes at the top of every .rpy replay file,
# so the script behaves the same with and without the GUI. It comes before the
# aliased imports underneath, so that abqmesh and abqjob are not shadowed.
try:
    from part import *
    from material import *
    from section import *
    from assembly import *
    from step import *
    from interaction import *
    from load import *
    from mesh import *
    from job import *
    from sketch import *
except ImportError:
    if not DRY_RUN:
        raise

# RESTORE PYTHON'S BUILT-INS. The star imports above bring in the Abaqus
# XY-data operators, which replace sum() (and possibly other built-ins) with
# versions that reject a generator ("found 'generator', expecting a recognized
# type"). The script does not use the XY operators, so the Python built-ins
# are restored for everything below.
import builtins as _builtins
sum = _builtins.sum
any = _builtins.any
all = _builtins.all
max = _builtins.max
min = _builtins.min
abs = _builtins.abs
round = _builtins.round

try:
    import regionToolset
    import mesh as abqmesh
    import job as abqjob
except ImportError:
    if not DRY_RUN:
        raise
import math
import os
import random

# ---------------------------------------------------------------------------
# DESIGN BASIS
# ---------------------------------------------------------------------------
# Every design input (aircraft records, material, criteria, pressures, thermal
# boundary conditions, load cases, fill levels, mesh and extraction rules) is
# read from code/design_basis.json through code/design_basis.py, the same file
# the plain-Python scripts read. To run another study, edit that file or point
# the environment variable LH2_DESIGN_BASIS at a copy of it.
#
# Abaqus does not always define __file__ (not under File > Run Script, and not
# in noGUI mode of Abaqus 2024 on Linux), so the folder code/ is searched in
# this order, taking the first that contains design_basis.py:
#   1. the folder above this script, from __file__
#   2. the folder above this script, from the command line (noGUI=<path>)
#   3. <working folder>/code (the cluster job runs from the repository root)
#   4. the folder of the file named by LH2_DESIGN_BASIS
# The values written further down in this file are the thesis values. They are
# used only when none of these is found, and the script then says so.
DB = None


def _find_code_dir():
    """The folder holding design_basis.py, or None."""
    cands = []
    try:
        cands.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    except NameError:
        pass
    import sys as _sys_cd
    for _arg in _sys_cd.argv:
        _arg = _arg.split('=', 1)[-1]
        if _arg.endswith('.py'):
            cands.append(os.path.dirname(os.path.dirname(os.path.abspath(_arg))))
    cands.append(os.path.join(os.getcwd(), 'code'))
    _env_db = os.environ.get('LH2_DESIGN_BASIS', '').strip()
    if _env_db:
        cands.append(os.path.dirname(os.path.abspath(_env_db)))
    for _c in cands:
        if os.path.exists(os.path.join(_c, 'design_basis.py')):
            return _c
    return None


_CODE_DIR = _find_code_dir()
if _CODE_DIR:
    import sys as _sys_db
    if _CODE_DIR not in _sys_db.path:
        _sys_db.path.insert(0, _CODE_DIR)
    import design_basis as DB          # errors in a user's design basis stop the run
    print('DESIGN BASIS: %s' % DB.PATH)
else:
    print('DESIGN BASIS: not found, using the thesis values written in this file')


def _db(getter, default):
    """Value from the design basis, or the thesis value written here."""
    return default if DB is None else getter(DB)

# ---------------------------------------------------------------------------
# WORKING FOLDER
# ---------------------------------------------------------------------------
# Abaqus writes its .odb, .dat, .msg and .log files and the results .csv files
# into the working folder, so the outputs never mix with the source code.
#   _THESIS_DIR  where the jobs run: E190_WORKDIR if set (the cluster job
#                sets it), otherwise the folder the script was started from
#   _ROOT_DIR    where the results files and sized_walls.csv are kept: the
#                folder the script was started from
# Both are the current folder in the normal use, abaqus cae noGUI=... started
# in an empty working folder.
_THESIS_DIR = os.getcwd()
_SCRIPT_DIR = None
try:
    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    _THESIS_DIR = os.path.dirname(_SCRIPT_DIR)
    # In the repository layout (code/1_abaqus_model/) the analysis runs in the
    # current working folder, not next to the source.
    if os.path.basename(_SCRIPT_DIR) == '1_abaqus_model':
        _THESIS_DIR = os.getcwd()
except NameError:
    pass
_ROOT_DIR = _THESIS_DIR          # results folder, independent of E190_WORKDIR
import csv
_ENV_DIR = os.environ.get('E190_WORKDIR', '').strip()
if _ENV_DIR:
    if not os.path.isdir(_ENV_DIR):
        raise ValueError('E190_WORKDIR is set to a directory that does not '
                         'exist: %s' % _ENV_DIR)
    _THESIS_DIR = _ENV_DIR
# Safeguard: if the results folder does not exist, use the working folder.
if not os.path.isdir(_ROOT_DIR):
    _ROOT_DIR = _THESIS_DIR if os.path.isdir(_THESIS_DIR) else os.getcwd()
    print('Thesis folder not found at the default path; using %s' % _ROOT_DIR)
if os.path.isdir(_THESIS_DIR):
    os.chdir(_THESIS_DIR)
print('Working directory: %s' % os.getcwd())

try:
    from odbAccess import openOdb
except ImportError:
    if not DRY_RUN:
        raise


# ============================================================================
# 1. USER INPUTS
# ============================================================================

# --- Run mode -------------------------------------------------------------
# PRODUCTION, STABILITY, CONVERGENCE, SIZING or FINAL; see the DRIVER section.
RUN_MODE = 'PRODUCTION'

# True runs the reference tank of Bagarello et al. instead of an aircraft.
BENCHMARK_MODE = False

# Both are set from the environment, for example (PowerShell):
#     $env:E190_RUN_MODE = "CONVERGENCE"
#     $env:E190_BENCHMARK = "1"
# Combine with E190_WORKDIR to keep the output of each run in its own folder.
_rm = os.environ.get('E190_RUN_MODE', '').strip().upper()
if _rm:
    if _rm not in ('PRODUCTION', 'CONVERGENCE', 'SIZING', 'STABILITY', 'FINAL'):
        raise ValueError("E190_RUN_MODE must be PRODUCTION, CONVERGENCE, SIZING, "
                         "STABILITY or FINAL, got %r" % _rm)
    RUN_MODE = _rm
    print('RUN MODE OVERRIDE: %s' % RUN_MODE)

_bm = os.environ.get('E190_BENCHMARK', '').strip().lower()
if _bm in ('1', 'true', 'yes', 'on'):
    BENCHMARK_MODE = True
    print('BENCHMARK MODE ACTIVE: Bagarello reference tank, not the aircraft '
          'geometry')
elif _bm in ('0', 'false', 'no', 'off', ''):
    pass
else:
    raise ValueError('E190_BENCHMARK must be 1/0 (or true/false), got %r' % _bm)


# --- WHICH AIRCRAFT ---------------------------------------------------------
# One script serves every aircraft: only the geometry, the wall thickness and
# the heat leak budget differ. The aircraft is chosen with E190_AIRCRAFT:
#     set E190_AIRCRAFT=A320
#     abaqus cae noGUI=code/1_abaqus_model/E190_MASTER_thermo_structural.py
# CASE_TAG follows AIRCRAFT, so each aircraft writes its own results files.
# Several aircraft can run at the same time in separate working folders
# (E190_WORKDIR), because the job names carry the aircraft (JOB NAME PREFIX).
AIRCRAFT = 'E190'            # default; any aircraft of the design basis
AIRCRAFT = os.environ.get('E190_AIRCRAFT', AIRCRAFT).strip() or AIRCRAFT

# Aircraft records:
# R_mid      shell mid-surface radius                                     [m]
# L_barrel   cylindrical barrel length                                    [m]
# cap_aspect a/b of the semi-ellipsoidal end caps                         [-]
# t_wall     starting wall of the sizing loop (closed-form sizing)        [m]
# Q_budget   allowable heat leak, from the 14.68 h no-vent hold           [W]
#            (boiloff_requirement_derivation.py)
# P_gauge    internal gauge pressure for strength                        [Pa]
# The wall used by every analysis after the sizing is the FE-sized wall read
# from sized_walls.csv below. The dictionary holds the thesis values and is
# rebuilt from the design basis whenever that is found.
AIRCRAFT_DATA = {
    'E190': dict(R_mid=1.3015, L_barrel=13.9675, cap_aspect=1.6,
                 t_wall=0.002747, Q_budget=702.6,  P_gauge=150000.0),
    'A320': dict(R_mid=1.7080, L_barrel=10.0120, cap_aspect=1.6,
                 t_wall=0.002694, Q_budget=919.3,  P_gauge=150000.0),
    'A350': dict(R_mid=2.5771, L_barrel=33.8240, cap_aspect=1.6,
                 t_wall=0.008884, Q_budget=6583.2, P_gauge=150000.0),
}
# From the design basis (the same values for the thesis aircraft). A new
# aircraft is added in design_basis.json, not here.
if DB is not None:
    AIRCRAFT_DATA = dict(
        (_k, dict(R_mid=DB.geometry(_k)['R'], L_barrel=DB.geometry(_k)['L'],
                  cap_aspect=DB.geometry(_k)['cap_aspect'],
                  t_wall=DB.record(_k)['t_wall_initial_m'], Q_budget=DB.heat_budget(_k),
                  P_gauge=DB.P_DESIGN))
        for _k in DB.AIRCRAFT)

if not BENCHMARK_MODE:
    if AIRCRAFT not in AIRCRAFT_DATA:
        raise ValueError('AIRCRAFT must be one of %s'
                         % sorted(AIRCRAFT_DATA.keys()))
    _ac        = AIRCRAFT_DATA[AIRCRAFT]
    R_mid      = _ac['R_mid']
    L_barrel   = _ac['L_barrel']
    cap_aspect = _ac['cap_aspect']
    t_wall     = _ac['t_wall']
    Q_budget   = _ac['Q_budget']
    P_gauge    = _ac['P_gauge']
    CASE_TAG   = AIRCRAFT
    # WALL FROM THE FE SIZING. RUN_MODE SIZING or FINAL records the sized wall
    # of each aircraft in sized_walls.csv. Every later run of this script, and
    # the capture, merge and plotting scripts, take the wall from there, so the
    # sized value is never retyped. The file is looked for in the results
    # folder, then in results/ of the repository.
    SIZED_WALLS = os.path.join(_ROOT_DIR, 'sized_walls.csv')
    if not os.path.exists(SIZED_WALLS) and (_SCRIPT_DIR or _CODE_DIR):
        _repo_sw = os.path.join(os.path.dirname(_CODE_DIR) if _CODE_DIR else
                                os.path.dirname(os.path.dirname(_SCRIPT_DIR)),
                                'results', 'sized_walls.csv')
        if os.path.exists(_repo_sw):
            SIZED_WALLS = _repo_sw
    if os.path.exists(SIZED_WALLS):
        for _row in csv.DictReader(open(SIZED_WALLS)):
            if _row.get('aircraft') == AIRCRAFT:
                t_wall = float(_row['t_wall_mm']) / 1000.0
                print('WALL FROM FE SIZING (sized_walls.csv): %.3f mm' % (t_wall * 1e3))
    # Overrides from the environment, for a trial wall or another pressure:
    #     set E190_T_WALL=2.75        (mm)
    #     set E190_P_GAUGE=0.18675    (bar gauge)
    _tw = os.environ.get('E190_T_WALL', '').strip()
    if _tw:
        t_wall = float(_tw) / 1000.0
        print('WALL THICKNESS OVERRIDE: %.3f mm' % (t_wall * 1000.0))
    _pg = os.environ.get('E190_P_GAUGE', '').strip()
    if _pg:
        P_gauge = float(_pg) * 1.0e5
        print('PRESSURE OVERRIDE: %.5f bar gauge' % (P_gauge / 1.0e5))
else:
    # ----------------------------------------------------------------
    # BAGARELLO REFERENCE TANK, their Table 4:
    #
    #     storage volume        V        100 m3
    #     cylinder length       L        4.24 m
    #     radius                R        2.12 m       spherical caps, a/b = 1
    #     inner layer           t_Al     10 mm
    #     middle layer          t_PUR    0.5 m        PUR-96 foam
    #     outer layer           t_CFRP   10 mm
    #     filling ratio         phi_f    0.75
    #
    # They model three solid layers; this model has one shell, the aluminium
    # liner. The foam and the composite enter through the outer heat flux,
    # computed from their stated boundary conditions (below), not fitted to
    # their result.
    # ----------------------------------------------------------------
    R_mid      = 2.12        # m,  their Table 4
    L_barrel   = 4.24        # m,  their Table 4
    cap_aspect = 1.0         # spherical caps
    t_wall     = 0.010       # m,  their t_Al = 10 mm, Table 4
    Q_budget   = None        # not used in benchmark mode; flux set explicitly
    P_gauge    = 170000.0    # Pa, their delta_p = 1.7 bar
    CASE_TAG   = 'BAGARELLO'


# --- JOB NAME PREFIX -------------------------------------------------------
# Job names carry the aircraft, so that the output databases of different
# aircraft never overwrite one another. The E190 keeps bare names, which the
# capture and post-processing scripts of the thesis refer to:
#
#     E190        HT_Baseline_f10.odb
#     A320        HT_A320_Baseline_f10.odb
#     A350        HT_A350_Baseline_f10.odb
#     benchmark   HT_BAGARELLO_Baseline_f10.odb
TAG_PREFIX = '' if CASE_TAG == 'E190' else (CASE_TAG + '_')


# --- Material: Al 2219-T87 ------------------------------------------------
# Yield is taken at ROOM temperature rather than the higher cryogenic value,
# which is conservative. Bagarello et al. use 500 MPa (cryogenic) with the
# same margin, an allowable of 333 MPa against 262 MPa here.
E_mod     = _db(lambda d: d.MAT['E_Pa'], 73.1e9)            # Pa
nu        = _db(lambda d: d.MAT['nu'], 0.33)                # -
rho_shell = _db(lambda d: d.MAT['rho_kg_m3'], 2850.0)       # kg/m3
alpha     = _db(lambda d: d.MAT['alpha_per_K'], 22.3e-6)    # 1/K
k_shell   = _db(lambda d: d.MAT['k_W_per_mK'], 116.0)       # W/m-K
cp_shell  = _db(lambda d: d.MAT['cp_J_per_kgK'], 864.0)     # J/kg-K
Fty_RT    = _db(lambda d: d.MAT['Fty_Pa'], 393e6)           # Pa, room-temperature yield
SF        = _db(lambda d: d.CRIT['yield_margin'], 1.5)      # project design margin, NOT the CS-25 factor of safety (thesis 3.6)
sigma_allow = Fty_RT / SF                                   # FORMULA 10

# Property overrides for the cryogenic-property check. alpha and k_shell
# above are ROOM-TEMPERATURE values, as in the benchmark (Bagarello Table 1).
# Near 20 K the secant expansion coefficient of aluminium is roughly 20 to 30
# times smaller (NIST, 6061-T6). The check brackets the effect without
# changing any default:
#     set E190_ALPHA=0          thermal strain switched off
#     set E190_KAPPA=21.7       conductivity scaled to 20 K (optional)
_al = os.environ.get('E190_ALPHA', '').strip()
if _al:
    alpha = float(_al)
    print('THERMAL EXPANSION OVERRIDE: alpha = %.3e 1/K' % alpha)
_ka = os.environ.get('E190_KAPPA', '').strip()
if _ka:
    k_shell = float(_ka)
    print('CONDUCTIVITY OVERRIDE: k = %.2f W/m-K' % k_shell)

# --- Stability -------------------------------------------------------------
# Axial buckling allowable of the barrel, NASA SP-8007 (1968), unpressurised
# knockdown. Pressure stabilisation is not credited.
#   sigma_cl = E t / (R sqrt(3(1-nu^2))),  phi = sqrt(R/t)/16,
#   gamma = 1 - 0.901 (1 - exp(-phi)),     sigma_cr = gamma sigma_cl
def sigma_cr_sp8007(R, t):
    s_cl = E_mod * t / (R * math.sqrt(3.0 * (1.0 - nu * nu)))
    phi_ = math.sqrt(R / t) / 16.0
    return (1.0 - 0.901 * (1.0 - math.exp(-phi_))) * s_cl

# Design basis (Chapter 3): STABILITY at the minimum differential the tank can
# hold when the load occurs, 1.2 bar abs (Winnefeld after Verstraete) against
# ISA sea-level 1.01325 bar; STRENGTH at the maximum design differential.
P_STABILITY = _db(lambda d: d.P_STABILITY, 120000.0 - 101325.0)   # Pa gauge = 0.18675 bar
FOS_BUCKLING = _db(lambda d: d.CRIT['buckling_margin'], 1.5)
# Optional minimum wall, off by default (the thesis sizes against buckling and
# yield only). To apply the ASME BPVC VIII-1 UG-16(b) floor of 1/16 in.:
#     set E190_MIN_GAUGE_MM=1.59
MIN_GAUGE = float(os.environ.get('E190_MIN_GAUGE_MM', '%r' % (1000.0 * _db(lambda d: d.CRIT['min_gauge_m'], 0.0)))) / 1000.0

# --- Fluid -----------------------------------------------------------------
rho_LH2 = _db(lambda d: d.MISSION['rho_lh2_kg_m3'], 70.8)     # kg/m3
T_LH2   = _db(lambda d: d.TH['T_liquid_K'], 20.3)            # K, Bagarello Sec. 3.1
T_GH2   = _db(lambda d: d.TH['T_ullage_K'], 20.3)            # K, ullage gas bulk temperature
# --- STRESS-FREE REFERENCE TEMPERATURE -------------------------------------
# THIS IS A MODELLING DECISION WITH A LARGE CONSEQUENCE. READ BEFORE CHANGING.
#
# Set to 20.3 K, matching Bagarello et al. Table 1, which lists a reference
# temperature of 20 K for the aluminium. The vessel is therefore defined as
# stress-free in its COLD OPERATING STATE, and the only thermal load the
# structure sees is the wetted-to-dry temperature difference, 3.4 to 7.4 K
# for the sized E190, A320 and A350-900 walls.
#
# WHY NOT 293 K, THE ASSEMBLY TEMPERATURE
# ---------------------------------------
# Using 293 K applies the full 273 K cooldown on top of every load case. For
# the E190 the model then gives 649.5 MPa, 98 percent of the biaxial fully
# restrained value E*alpha*dT/(1-nu) = 663 MPa, against a 262 MPa allowable:
# the rings lock the vessel diameter, so the radial contraction the tank wants
# is resisted entirely by the metal. E190, reference 293 K against 20.3 K
# (largest wetted-to-dry difference 3.68 K):
#
#   excursion from stress-free : 272.7 K        ->  3.68 K
#   radial contraction wanted  :   7.9 mm       ->  0.11 mm
#   axial contraction wanted   :  47.4 mm       ->  0.64 mm   (towards mid-barrel)
#   stress if fully blocked    : 663 MPa        ->  9.0 MPa
#
# A factor of 74 less for the supports to resist.
#
# THE PHYSICAL ARGUMENT
# ---------------------
# A 273 K cooldown and a 9 g crash do not occur together. The tank is filled
# on the ground and cools over hours, contracting freely because real mounts
# use sliding bearings and flexures for exactly that purpose. It then flies.
# A crash happens to a vessel that is already cold and already settled.
# Superimposing the two is not conservative, it combines events separated by
# hours.
#
# WHAT THIS EXCLUDES
# ------------------
# The cooldown itself is not covered here. It belongs in its own analysis:
# a ground condition at 1 g with no crash factors, where the question is
# whether the supports can accommodate 47 mm of contraction. Biancotto
# Sec. 5.8.1 addresses precisely this and finds that supports resisting the
# inner vessel's thermal displacement are overstressed by 1.6 to 3.8 times.
# Recorded as further work rather than silently omitted.
T_ref   = _db(lambda d: d.TH['T_reference_K'], 20.3)         # K, stress-free reference (cold operating state)
T_ASSEMBLY = 293.0       # K, retained for the cooldown case, not used here
g       = 9.81           # m/s2

# --- Thermal boundary conditions ------------------------------------------
# Inner surface, WETTED region: liquid hydrogen against aluminium. The film
# coefficient is high, so the wall is effectively pinned at the fuel
# temperature. Bagarello justify this by "the high thermal conductivity of the
# aluminium alloy and the great transfer at the fluid-structure interface".
h_wetted = _db(lambda d: d.TH['h_wetted_W_per_m2K'], 1000.0)   # W/m2-K, effectively pins the wall to T_LH2

# Inner surface, DRY region: hydrogen vapour by natural convection.
# THIS IS THE KEY NUMBER. It is far smaller than h_wetted, which is why the
# dry wall runs warmer and a gradient forms.
h_dry    = _db(lambda d: d.TH['h_dry_W_per_m2K'], 1.46)       # W/m2-K, Bagarello Sec. 3.1, their ref [45]

# Outer surface: heat arriving through the MLI blanket, applied as a uniform
# inward flux. The inner vessel is the only body modelled, so the insulation
# enters solely through this number.
if not BENCHMARK_MODE:
    MLI_FLUX_MODE = 'FROM_BUDGET'   # q = Q_budget / A_shell      FORMULA 5
else:
    MLI_FLUX_MODE = 'EXPLICIT'
    # ----------------------------------------------------------------
    # BENCHMARK OUTER FLUX, DERIVED FROM THEIR INPUTS, NOT THEIR RESULT.
    #
    # Their Sec. 3.1 states the external boundary condition explicitly:
    #
    #     q_ext = h_c^ext * (T_fus - T_ext)
    #     h_c^ext = 10.00 W/m2K        their Ref. [45]
    #     T_fus   = 283 K              fuselage-side ambient
    #     T_ext   = 282.2 K            outer wall, their reported solution
    #
    #     q = 10.00 * (283 - 282.2) = 8.0 W/m2
    #
    # Cross-check through the insulation, using their Table 1 conductivity
    # and Table 4 thickness, entirely independently of the above:
    #
    #     q = kappa_PUR * dT / t_PUR = 0.0163 * (282.2 - 25) / 0.5
    #       = 8.38 W/m2
    #
    # The two agree to 5 percent, the difference being the CFRP and liner
    # resistances. 8.0 W/m2 is used.
    #
    # WHY THIS IS NOT CIRCULAR. The value comes from their stated convection
    # coefficient and ambient temperature, and is confirmed by their stated
    # foam conductivity and thickness. It is not back-calculated from the
    # 25.3 K result. That number remains the prediction being tested:
    #
    #     predicted dT = q / h_dry = 8.0 / 1.46 = 5.48 K
    #     predicted T_max = 20 + 5.48 = 25.5 K   against their 25.3 K
    # ----------------------------------------------------------------
    q_outer_explicit = 8.0          # W/m2, Bagarello Sec. 3.1


# --- CS-25 load cases ------------------------------------------------------
# Load factors (nx forward, ny lateral, nz vertical) in g; see CONVENTIONS in
# the header for the mapping to the model axes. Codes and names:
#   Baseline       Baseline condition: 1 g, ground or level flight
#   LC1_Maneuver   Symmetric manoeuvre: CS 25.337 limit factor
#   LC3_Combined   Combined manoeuvre: constructed case, not from CS-25
#   LC7_Emergency  Emergency landing: CS 25.561 factors applied together,
#                  which is more severe than the regulation, where they act
#                  separately
load_cases = {
    'Baseline':      {'nx': 0.0, 'ny': 0.0, 'nz': 1.0},
    'LC1_Maneuver':  {'nx': 0.0, 'ny': 0.0, 'nz': 2.5},
    'LC3_Combined':  {'nx': 1.5, 'ny': 0.0, 'nz': 2.5},
    'LC7_Emergency': {'nx': 9.0, 'ny': 3.0, 'nz': 6.0},
}
if DB is not None:
    load_cases = dict((_k, {'nx': DB.LOAD_CASES[_k][0], 'ny': DB.LOAD_CASES[_k][1],
                            'nz': DB.LOAD_CASES[_k][2]}) for _k in DB.LOAD_CASE_ORDER)
# The case that sizes the wall and is used for the mesh study (thesis 3.6, 4.9).
GOVERNING_CASE = _db(lambda d: d.ANALYSIS['governing_case'], 'LC7_Emergency')
GOVERNING_FILL = _db(lambda d: d.ANALYSIS['governing_fill'], 0.9)
ALL_LOAD_CASES = dict(load_cases)       # before any E190_ONLY_LC filter
GOVERNING_FACTORS = ALL_LOAD_CASES[GOVERNING_CASE]

# ----------------------------------------------------------------------------
# OPTIONAL LOAD-CASE FILTER, FOR RUNNING THE MATRIX IN PARALLEL SESSIONS
# ----------------------------------------------------------------------------
# E190_ONLY_LC restricts a session to a subset of the load cases, so the
# matrix can be split across several Abaqus sessions running at the same time.
# On a shared licence pool most of the wall-clock time is spent waiting for a
# solver token, and several sessions hold several places in that queue.
#
# The analysis does not change: every design point is built from scratch,
# meshed at the same seed and solved with the same loads under its own job
# name. The design points are independent, so splitting them across sessions
# cannot alter any number.
#
# Use, one Command Prompt window per load case:
#
#     set E190_ONLY_LC=Baseline
#     abaqus cae noGUI=code/1_abaqus_model/E190_MASTER_thermo_structural.py
#
# then the same with LC1_Maneuver, LC3_Combined, LC7_Emergency. Each session
# writes its own CSV, named after the filter, so nothing is overwritten:
#
#     E190_production_Baseline_results.csv
#     E190_production_LC1_Maneuver_results.csv       ... and so on
#
# In a FINAL run the sessions resume from one another's files; in a
# PRODUCTION run, concatenate the files afterwards. Each session also checks
# out one CAE licence, so keep the number of sessions below the free seats.
_only_lc = os.environ.get('E190_ONLY_LC', '').strip()
LC_FILTER = ''
if _only_lc:
    _wanted = [s.strip() for s in _only_lc.split(',') if s.strip()]
    _unknown = [s for s in _wanted if s not in load_cases]
    if _unknown:
        raise ValueError(
            'E190_ONLY_LC names load cases that do not exist: %s\n'
            'Valid names: %s' % (', '.join(_unknown), ', '.join(load_cases)))
    load_cases = dict((k, v) for k, v in load_cases.items() if k in _wanted)
    LC_FILTER = '_'.join(_wanted)
    print('\nLOAD CASE FILTER ACTIVE: %s' % ', '.join(_wanted))
    print('This session runs %d of the 4 load cases. Run the remainder in'
          % len(load_cases))
    print('parallel sessions and concatenate the CSV files afterwards.')

# --- Fill ratios -----------------------------------------------------------
# Independent static design points, NOT a mission timeline.
fill_ratios = _db(lambda d: list(d.FILLS), [0.1, 0.3, 0.5, 0.75, 0.9])

# OPTIONAL FILL FILTER - the second axis for splitting the matrix.
# Same purpose and same guarantees as E190_ONLY_LC above: it selects which
# points this session computes and changes nothing about how they are
# computed. Together the two filters address any subset of the 4 x 5 matrix,
# so the twenty cases can be spread over up to twenty simultaneous sessions.
#
#     set E190_ONLY_LC=LC7_Emergency
#     set E190_ONLY_FILL=0.75,0.9
#
# EVERY SESSION MUST COVER A DIFFERENT SUBSET. Two sessions given the same
# (load case, fill) pair will write the same job name, the same .odb and the
# same lock file at the same time, and the results of both are then worthless.
# The tag printed at the start of each case is the thing to check.
_only_fill = os.environ.get('E190_ONLY_FILL', '').strip()
FILL_FILTER = ''
if _only_fill:
    _fwanted = []
    for _s in _only_fill.split(','):
        _s = _s.strip()
        if not _s:
            continue
        try:
            _v = float(_s)
        except ValueError:
            raise ValueError('E190_ONLY_FILL is not a number: %r' % _s)
        _match = [f for f in fill_ratios if abs(f - _v) < 1e-9]
        if not _match:
            raise ValueError(
                'E190_ONLY_FILL value %s is not one of the defined fill '
                'ratios %s' % (_s, fill_ratios))
        _fwanted.append(_match[0])
    fill_ratios = sorted(set(_fwanted))
    FILL_FILTER = 'f' + '_'.join('%02d' % int(round(f * 100))
                                 for f in fill_ratios)
    print('FILL RATIO FILTER ACTIVE: %s' % fill_ratios)

# --- Mesh ------------------------------------------------------------------
# The shell bending boundary layer at a ring or a cap-to-barrel junction has a
# characteristic length of sqrt(R*t): 41.0, 51.8 and 120.9 mm for the sized
# E190, A320 and A350-900 walls. The element size must resolve it, so
# seed_size below is a fraction of that length (thesis Equation 4.8).
MESHED_MASS = 0.0        # filled in by build_geometry_and_mesh, kg

# Axial half-width of the zone excluded around each ring when reporting peak
# stress. For the sized walls the decay length pi/beta is 0.10, 0.13 and
# 0.30 m, so 0.5 m spans 5.0, 3.9 and 1.7 decay lengths (thesis Section 4.7).
# See the FAR-FIELD block in extract() for the full justification.
FAR_FIELD_EXCLUSION = _db(lambda d: d.ANALYSIS['far_field_exclusion_m'], 0.5)   # m

boundary_layer = math.sqrt(R_mid * t_wall)
#
# MESH_FRACTION IS SET FROM THE CONVERGENCE STUDY, NOT ASSUMED.
# The governing case was solved with 1, 2 and 4 elements per bending length
# (fractions 1, 0.5 and 0.25). At 0.5 the far-field von Mises stress lies
# within 1 % of the finest mesh for all three tanks; the coarsest mesh
# overestimates it by 5 to 11 % (thesis Section 4.9, Table 4.5).
MESH_FRACTION = _db(lambda d: d.ANALYSIS['mesh_fraction_of_sqrt_Rt'], 0.5)   # elements per boundary layer = 2
seed_size = boundary_layer * MESH_FRACTION

# Densities swept by the convergence study: 1, 2 and 4 elements per boundary
# layer. A finer level (8 per layer) would be 2 to 4 million elements at the
# sized walls, because a thin wall shortens sqrt(R t) and the element count
# grows as 1/t. Other levels from the environment, for example
#     set E190_CONV_FRACTIONS=2,1,0.5,0.25
convergence_fractions = _db(lambda d: list(d.ANALYSIS['convergence_fractions']), [1.0, 0.5, 0.25])
_cf = os.environ.get('E190_CONV_FRACTIONS', '').strip()
if _cf:
    convergence_fractions = [float(x) for x in _cf.split(',') if x.strip()]
    print('CONVERGENCE LEVELS OVERRIDE: %s' % convergence_fractions)

n_bins = _db(lambda d: d.ANALYSIS['hydrostatic_bins'], 64)   # depth bins for the hydrostatic pressure field


# ============================================================================
# 2. DERIVED GEOMETRY
# ============================================================================

h_cap = R_mid / cap_aspect                                   # FORMULA 2
L_total = L_barrel + 2.0 * h_cap
V_tank = (math.pi * R_mid**2 * L_barrel
          + (4.0 / 3.0) * math.pi * R_mid**2 * h_cap)        # FORMULA 1

# Wetted-surface area is computed from the mesh at run time. This analytic
# value is used only for the MLI flux and for reporting.
A_barrel = 2.0 * math.pi * R_mid * L_barrel
# EXACT surface area of the two semi-ellipsoidal caps, which together form one
# complete oblate spheroid of equatorial radius R and polar semi-axis h_cap:
#
#     A = 2*pi*R^2 * [ 1 + ((1-e^2)/e) * artanh(e) ],   e = sqrt(1 - (h/R)^2)
#
# The approximation 2*pi*R*h_cap for the two caps together understates their
# area by 59 percent for the E190, and the shell area by about 8 percent,
# which would appear as an equal imbalance in the load balance check.
#
# SPHERICAL CAPS (cap_aspect = 1, the benchmark tank) are a special case: the
# eccentricity is zero and the formula divides by e = 0. Its limit as e -> 0,
# the sphere A = 4*pi*R^2, is used instead.
_e2 = 1.0 - (h_cap / R_mid) ** 2
if _e2 <= 1e-12:
    # hemispherical caps: the two together form a complete sphere
    A_caps = 4.0 * math.pi * R_mid ** 2
    _ecc = 0.0
else:
    _ecc = math.sqrt(_e2)
    A_caps = 2.0 * math.pi * R_mid ** 2 * (
        1.0 + ((1.0 - _ecc ** 2) / _ecc) * math.atanh(_ecc))
A_shell = A_barrel + A_caps

if MLI_FLUX_MODE == 'FROM_BUDGET':
    q_outer = Q_budget / A_shell                             # FORMULA 5
else:
    q_outer = q_outer_explicit

print('=' * 78)
print('LH2 INNER VESSEL - THERMO-STRUCTURAL ANALYSIS (Bagarello methodology)')
print('=' * 78)
print('Case                 : %s' % CASE_TAG)
print('Run mode             : %s' % RUN_MODE)
print('R_mid                : %.4f m' % R_mid)
print('L_barrel / L_total   : %.4f / %.4f m' % (L_barrel, L_total))
print('Wall thickness       : %.3f mm' % (t_wall * 1000.0))
print('Tank volume          : %.3f m3' % V_tank)
print('Shell area (analytic): %.3f m2' % A_shell)
print('MLI flux on outer    : %.4f W/m2' % q_outer)
print('h_wetted / h_dry     : %.1f / %.2f W/m2-K' % (h_wetted, h_dry))
print('Boundary layer sqrt(R*t) : %.1f mm' % (boundary_layer * 1000.0))
print('Seed size            : %.1f mm  (%.1f elements per boundary layer)'
      % (seed_size * 1000.0, 1.0 / MESH_FRACTION))
print('Allowable stress     : %.1f MPa  (Fty %.0f MPa / SF %.1f)'
      % (sigma_allow / 1e6, Fty_RT / 1e6, SF))
print('=' * 78)


# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

# ----------------------------------------------------------------------------
# TILTED FREE SURFACE
# ----------------------------------------------------------------------------
# A liquid under steady acceleration reaches relative equilibrium: it behaves
# exactly as under gravity, except that "down" is the direction of the
# resultant load factor rather than the vertical. The free surface is the
# plane normal to that direction which leaves the fill fraction of the tank
# volume below it, and the depth is measured along that direction.
#
# For the E190 at a fill level of 0.90 the maximum depth of liquid is
#
#     load case             tilt from vertical   maximum depth
#     Baseline condition         0.0 deg            2.19 m
#     Symmetric manoeuvre        0.0 deg            2.19 m
#     Combined manoeuvre        31.0 deg            7.91 m
#     Emergency landing         57.7 deg           11.45 m
#
# The tank is 15.6 m long and 2.6 m tall, so a load factor tilted toward the
# axis makes the liquid column much deeper than the tank is tall. For a purely
# vertical load factor the treatment reduces to the horizontal free surface.

_VOL_PTS = None          # volume sample points, built once and reused


def _sample_tank_volume(n_pts=200000):
    """
    Uniform points inside the vessel, by rejection sampling in the bounding
    box. Used to locate the free surface for an arbitrary tilt, for which
    there is no closed-form submerged-volume expression.

    Built once and cached, because the point cloud depends only on geometry,
    not on the load case or the fill ratio.
    """
    global _VOL_PTS
    if _VOL_PTS is not None:
        return _VOL_PTS
    random.seed(12345)                       # reproducible
    pts = []
    x_lo, x_hi = -h_cap, L_barrel + h_cap
    RR = R_mid * R_mid
    while len(pts) < n_pts:
        x = random.uniform(x_lo, x_hi)
        y = random.uniform(-R_mid, R_mid)
        z = random.uniform(-R_mid, R_mid)
        rr = y * y + z * z
        if 0.0 <= x <= L_barrel:
            ok = rr <= RR
        elif x < 0.0:
            ok = (x / h_cap) ** 2 + rr / RR <= 1.0
        else:
            ok = ((x - L_barrel) / h_cap) ** 2 + rr / RR <= 1.0
        if ok:
            pts.append((x, y, z))
    _VOL_PTS = pts
    return pts


def effective_down(nx, ny, nz):
    """
    Unit vector in the direction the liquid settles, in global (X, Y, Z) =
    (forward, up, lateral).

    Taken from the body force the model applies (comp1 = +nx*g,
    comp2 = -nz*g, comp3 = -ny*g, see run_case), so the liquid settles in the
    direction in which the structure is loaded, as relative equilibrium
    requires.
    """
    v = (nx, -nz, -ny)
    mag = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    if mag == 0.0:
        return (0.0, -1.0, 0.0)
    return (v[0] / mag, v[1] / mag, v[2] / mag)


def free_surface_offset(g_hat, phi_f):
    """
    The free-surface plane offset s_surf, such that the fraction of the tank
    volume with s = r.g_hat greater than s_surf equals phi_f.

    Because the sample points are uniform in VOLUME, this is simply a quantile
    of the s distribution: sort and index. No bisection is required, and the
    result is exact in the limit of many samples.
    """
    if phi_f <= 0.0:
        return 1e30                          # nothing wetted
    if phi_f >= 1.0:
        return -1e30                         # everything wetted
    pts = _sample_tank_volume()
    s_vals = []
    for (x, y, z) in pts:
        s_vals.append(x * g_hat[0] + y * g_hat[1] + z * g_hat[2])
    s_vals.sort()
    idx = int((1.0 - phi_f) * (len(s_vals) - 1))
    return s_vals[idx]


def elem_depth_coord(el, g_hat):
    """Depth coordinate s = r.g_hat at an element centroid. Larger is deeper."""
    nds = el.getNodes()
    sx = sy = sz = 0.0
    cnt = 0
    for nd in nds:
        c = nd.coordinates
        sx += c[0]
        sy += c[1]
        sz += c[2]
        cnt += 1
    cnt = float(cnt)
    return (sx / cnt) * g_hat[0] + (sy / cnt) * g_hat[1] + (sz / cnt) * g_hat[2]


def fill_surface_height(phi_f):
    """
    ALTERNATIVE MODEL, KEPT FOR VERIFICATION AND FUTURE STUDIES (not called).

    FORMULA 3, horizontal special case. Height of a horizontal liquid free
    surface, measured from the tank axis, for a given fill ratio.

    For a horizontal cylinder the submerged area is a circular segment, and
    the fill ratio relates to the segment half-angle theta through

        theta - sin(theta) = 2*pi*phi_f

    which has no closed-form inverse, so Newton's method is used. The result
    is converted back to a height by y = -R*cos(theta/2).

    It assumes a HORIZONTAL free surface, i.e. a vertical load factor, which
    holds for the Baseline condition and the symmetric manoeuvre but not for
    the combined manoeuvre or the emergency landing. It also treats the tank
    as a cylinder without caps. Using it in the load path would make the
    thermal model ignore the load case.

    Use as a verification: for a purely vertical load factor,
    free_surface_offset() must reproduce -fill_surface_height(phi_f) to within
    the sampling error of the point cloud and the effect of the caps.
    """
    if phi_f <= 0.0:
        return -R_mid
    if phi_f >= 1.0:
        return R_mid
    theta = 2.0 * math.pi * phi_f
    for _ in range(100):
        gval = theta - math.sin(theta) - 2.0 * math.pi * phi_f
        gprime = 1.0 - math.cos(theta)
        if abs(gprime) < 1e-12:
            gprime = 1e-12
        dtheta = gval / gprime
        theta -= dtheta
        if abs(dtheta) < 1e-10:
            break
    return -R_mid * math.cos(theta / 2.0)


def elem_centroid_z(el):
    """
    ALTERNATIVE MODEL, KEPT FOR VERIFICATION AND FUTURE STUDIES (not called).

    Vertical coordinate of an element centroid, the companion of
    fill_surface_height for a horizontal free surface.

    NOTE ON THE ABAQUS API: MeshElement.connectivity returns node INDICES into
    the instance's node array, not node LABELS. Indexing a label-keyed
    dictionary with them raises KeyError, because labels are 1-based and
    indices are 0-based. getNodes() sidesteps the distinction by returning the
    MeshNode objects directly, so it is used here.

    NOTE ON sum(): the Abaqus kernel replaces the builtin sum() with a version
    that rejects generator expressions ("found 'generator', expecting a
    recognized type"). Every accumulation in this script is therefore written
    as an explicit loop rather than sum(... for ... in ...).

    AXIS CONVENTION (fixed throughout this script, see the free body diagram
    figure 11): +X FORWARD, +Y UP, +Z LATERAL. Right-handed, and matching
    the Abaqus viewport triad.

    The meridian is sketched and revolved about the sketch's horizontal
    axis, so the vessel axis lies along GLOBAL X with its circular
    cross-section in the GLOBAL Y-Z plane. The tank spans x = -h_cap at the
    AFT pole to x = L_barrel + h_cap at the FORWARD pole.

    Vertical is therefore Y, which is coordinates[1]. (Using coordinates[2]
    would fill the tank sideways: a von Mises field symmetric top to bottom
    instead of rising from the crown to the keel.)
    """
    nds = el.getNodes()
    tot = 0.0
    for n in nds:
        tot += n.coordinates[1]
    return tot / float(len(nds))


def build_geometry_and_mesh(model_name, seed, analysis_type='STRUCTURAL'):
    """
    Build the axisymmetric tank as a revolved shell and mesh it.

    analysis_type selects the element family: 'THERMAL' gives diffusive shell
    elements for the heat transfer step, 'STRUCTURAL' gives stress elements
    for the static step.

    Shell elements are used rather than solids. The sized walls have R/t of
    about 450 (A350-900) to 1100 (A320), firmly in the thin-shell regime, and
    with k = 116 W/m-K the through-thickness temperature drop is negligible,
    so a shell conduction model is adequate.

    Geometry: the meridian (aft cap, barrel, forward cap) is sketched in the
    X-r plane and revolved 360 degrees about X. The caps are exact quarter
    ellipses with semi-axes h_cap (axial) and R_mid (radial).
    """
    model = mdb.Model(name=model_name)

    # --- Sketch the meridian: half an ellipse, the barrel, half an ellipse --
    s = model.ConstrainedSketch(name='meridian', sheetSize=4.0 * L_total)
    s.ConstructionLine(point1=(0.0, 0.0), point2=(1.0, 0.0))

    # Aft cap: quarter ellipse from the aft pole (x = -h_cap) to the start of
    # the barrel (x = 0). Parametrised, so the profile is exact rather than a
    # spline fit.
    n_cap_pts = 40
    pts = []
    for i in range(n_cap_pts + 1):
        ang = (math.pi / 2.0) * i / float(n_cap_pts)
        x = -h_cap * math.cos(ang)      # from -h_cap up to 0
        r = R_mid * math.sin(ang)
        pts.append((x, r))
    for i in range(len(pts) - 1):
        s.Line(point1=pts[i], point2=pts[i + 1])

    # Barrel
    s.Line(point1=(0.0, R_mid), point2=(L_barrel, R_mid))

    # Forward cap, mirror of the aft one, from x = L_barrel to the forward
    # pole at x = L_barrel + h_cap
    pts2 = []
    for i in range(n_cap_pts + 1):
        ang = (math.pi / 2.0) * i / float(n_cap_pts)
        x = L_barrel + h_cap * math.cos(ang)
        r = R_mid * math.sin(ang)
        pts2.append((x, r))
    pts2.reverse()
    for i in range(len(pts2) - 1):
        s.Line(point1=pts2[i], point2=pts2[i + 1])

    part = model.Part(name='Vessel', dimensionality=THREE_D,
                      type=DEFORMABLE_BODY)
    part.BaseShellRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)

    # --- Material and section ---------------------------------------------
    mat = model.Material(name='Al2219T87')
    mat.Elastic(table=((E_mod, nu), ))
    mat.Density(table=((rho_shell, ), ))
    mat.Expansion(table=((alpha, ), ), zero=T_ref)            # FORMULA 8
    mat.Conductivity(table=((k_shell, ), ))
    mat.SpecificHeat(table=((cp_shell, ), ))

    # NOTE ON THE KEYWORD: the argument that controls how temperature varies
    # through the shell thickness is called `temperature`, NOT
    # `temperatureType`. GRADIENT stores a mid-surface value plus a
    # through-thickness gradient. The applied temperature fields set the
    # gradient to zero explicitly (CONSTANT_THROUGH_THICKNESS in run_case).
    model.HomogeneousShellSection(name='Wall', material='Al2219T87',
                                  thickness=t_wall,
                                  temperature=GRADIENT)
    all_faces = part.Set(name='AllFaces', faces=part.faces[:])
    part.SectionAssignment(region=all_faces, sectionName='Wall')

    # --- Mesh --------------------------------------------------------------
    part.seedPart(size=seed, deviationFactor=0.1, minSizeFactor=0.1)

    # A revolved surface with elliptical caps is not always structured-
    # meshable. Try STRUCTURED for element quality, fall back to FREE so the
    # script does not stop on a meshing technicality.
    # WHICH TECHNIQUE SUCCEEDS MATTERS. A structured mesh on a revolved
    # surface has consistent element normals by construction. The free
    # quad-dominated fallback does not guarantee that, and inconsistent
    # normals would make "side1" point outward on some elements and inward on
    # others. The gauge pressure of 150 kPa acts on a 5.32 m2 cross section of
    # the E190 (798 kN), so a few percent of flipped elements would produce
    # hundreds of kN of spurious resultant on a self-equilibrating load. The
    # technique used is therefore reported, and the load balance check in
    # extract() would expose it.
    mesh_technique = 'STRUCTURED'
    try:
        part.setMeshControls(regions=part.faces[:], elemShape=QUAD,
                             technique=STRUCTURED)
        part.generateMesh()
        if len(part.elements) == 0:
            raise ValueError('structured meshing produced no elements')
    except Exception:
        part.deleteMesh()
        part.setMeshControls(regions=part.faces[:], elemShape=QUAD_DOMINATED,
                             technique=FREE)
        part.generateMesh()
        mesh_technique = 'FREE (fallback) - normals may be inconsistent'
    print('    mesh: %s, %d elements' % (mesh_technique, len(part.elements)))
    # Capture the mass Abaqus actually meshed, for the load balance check.
    # getMassProperties() can return a dict whose 'mass' entry is None rather
    # than raising, so the result is validated before it is stored.
    try:
        global MESHED_MASS
        mp = part.getMassProperties()
        m_ = mp['mass'] if 'mass' in mp else None
        if m_ is not None and m_ > 0.0:
            MESHED_MASS = m_
            print('    meshed mass: %.1f kg (analytic estimate %.1f kg)'
                  % (m_, A_shell * t_wall * rho_shell))
    except Exception:
        pass

    # ELEMENT TYPE MUST MATCH THE ANALYSIS PROCEDURE.
    # A heat transfer step needs diffusive shell elements (DS4/DS3); a static
    # stress step needs structural shell elements (S4R/S3). The same mesh
    # cannot serve both, which is why the thermal and structural models are
    # built separately rather than by swapping the step on one model.
    if analysis_type == 'THERMAL':
        etypes = (abqmesh.ElemType(elemCode=DS4, elemLibrary=STANDARD),
                  abqmesh.ElemType(elemCode=DS3, elemLibrary=STANDARD))
    else:
        etypes = (abqmesh.ElemType(elemCode=S4R, elemLibrary=STANDARD),
                  abqmesh.ElemType(elemCode=S3, elemLibrary=STANDARD))
    part.setElementType(regions=(part.faces[:], ), elemTypes=etypes)

    model.rootAssembly.DatumCsysByDefault(CARTESIAN)
    inst = model.rootAssembly.Instance(name='VesselInst', part=part,
                                       dependent=ON)
    model.rootAssembly.regenerate()
    return model, part, inst


def apply_ring_constraints(model, inst, step_name, seed):
    """
    PROCESS 6.  THREE CIRCUMFERENTIAL RING CONSTRAINTS.

    PROVENANCE
    ----------
    Bagarello et al. Sec. 4.4.1: "Following an approach similar to that used
    by the company H2FLY in testbeds, the component is constrained along three
    circular edges of the outer layer. These correspond to the interfaces
    between the hemispherical and cylindrical regions, as well as the midpoint
    along the cylinder's length."

    So: rings at the two cap-to-barrel junctions and at mid-barrel. This is
    taken from flying demonstrator hardware rather than assumed, which is why
    it is preferred over an arbitrary support position.

    WHY RINGS AND NOT POINTS
    ------------------------
    A point support on a dome apex would concentrate the whole reaction on a
    few elements of a thin shell, giving a spurious local stress of the order
    of a thousand MPa (see apply_rigid_body_restraint). A ring spreads the
    reaction over the full circumference, as real tank hardware does.

    WHY THREE AND NOT TWO
    ---------------------
    The middle ring halves the bending span. Since bending moment scales with
    the square of the span, it removes roughly three quarters of the
    longitudinal bending.

    DEGREES OF FREEDOM - AN ENGINEERING DECISION, STATED OPENLY
    ----------------------------------------------------------
    Bagarello et al. do not specify which degrees of freedom are restrained,
    and it does not greatly matter for them because their stress-free
    reference temperature is 20 K: their model never sees the cooldown.

    This model also uses 20.3 K as the stress-free reference, for the reason
    given where T_ref is defined: cooldown and the CS-25 event are hours
    apart, so the vessel is cold and settled when the load arrives. The only
    thermal load carried here is the wetted-to-dry gradient.

    The anchor-and-slide scheme is nevertheless used rather than fixing all
    three rings. Were the reference ever moved to 293 K, fully
    fixing all three would block roughly 95 mm of axial contraction and
    return the fully restrained thermal stress, which exceeds yield - an
    artefact of over-constraining rather than a physical result. Keeping the
    scheme physical means that choice of reference does not silently change
    the answer.

    Real cryogenic tank supports are designed to accommodate contraction.
    Biancotto's central finding was precisely that supports which resist the
    inner vessel's thermal displacement are overstressed by a factor of
    between 1.6 and 3.8.

    The arrangement used here is the standard anchor-and-slide scheme found in
    pressure vessel and cryogenic piping practice:

        forward ring : radial and circumferential held, FREE to slide axially
        middle  ring : fully held - this is the anchor
        aft     ring : radial and circumferential held, FREE to slide axially

    Every ring therefore carries vertical and lateral load, while the vessel
    remains free to shrink axially toward the middle ring. In global terms,
    with the vessel axis along X, this means u2 and u3 are held at every ring
    and u1 is held only at the middle one.
    """
    x_stations = [
        # AXIS CONVENTION: +X forward, +Y up, +Z lateral (right-handed,
        # matching the Abaqus triad). The tank runs from x = -h_cap at the
        # AFT pole to x = L_barrel + h_cap at the FORWARD pole, so the
        # junction at x = L_barrel is the forward one.
        ('aft', 0.0,             False),   # aft cap-to-barrel junction
        ('mid', L_barrel / 2.0,  True),    # mid-barrel, the axial anchor
        ('fwd', L_barrel,        False),   # forward cap-to-barrel junction
    ]

    # SNAP TO THE NEAREST NODE PLANE, DO NOT USE A FIXED TOLERANCE.
    #
    # The two junction stations sit on a geometric edge, so nodes line up on
    # them exactly. Mid-barrel has no such feature, and a fixed tolerance can
    # catch a partial ring (an asymmetric restraint). Instead, the x coordinate
    # of the node plane nearest the requested station is found, and every node
    # on that plane is taken. This always gives one complete ring, wherever
    # the mesher put the nodes.
    xs = []
    for nd in inst.nodes:
        xs.append(nd.coordinates[0])

    snap_tol = max(seed * 0.05, 1e-4)      # nodes within this are one plane

    for name, x0, is_anchor in x_stations:
        # nearest actual node plane to the requested station
        x_near = min(xs, key=lambda v: abs(v - x0))
        labels = []
        for nd in inst.nodes:
            if abs(nd.coordinates[0] - x_near) <= snap_tol:
                labels.append(nd.label)
        if not labels:
            print('    *** no nodes found for ring %s at x = %.3f' % (name, x0))
            continue
        if abs(x_near - x0) > seed:
            print('    !!! ring %s snapped %.1f mm from the requested station'
                  % (name, (x_near - x0) * 1000.0))
        nset = model.rootAssembly.Set(
            name='Ring_%s' % name,
            nodes=inst.nodes.sequenceFromLabels(tuple(labels)))
        if is_anchor:
            model.DisplacementBC(name='Ring_%s' % name, createStepName='Initial',
                                 region=nset, u1=0.0, u2=0.0, u3=0.0)
        else:
            model.DisplacementBC(name='Ring_%s' % name, createStepName='Initial',
                                 region=nset, u2=0.0, u3=0.0)
        print('    ring %s at x = %7.3f m : %5d nodes, %s'
              % (name, x_near, len(labels),
                 'ANCHOR (u1,u2,u3)' if is_anchor else 'sliding (u2,u3)'))


def apply_inertia_relief(model, step_name):
    """
    ALTERNATIVE MODEL, KEPT FOR FUTURE STUDIES (not called).

    INERTIA RELIEF: the applied loads are balanced by a distributed
    d'Alembert body force proportional to the mass instead of by support
    reactions. Oom Ortiz de Montellano et al. compared a constrained and an
    inertia relief model of the same structure and adopted inertia relief.

    WHY THE THESIS USES RINGS INSTEAD
    ---------------------------------
    The balancing acceleration is distributed over the model's mass, which is
    the shell only: the liquid enters as a pressure field, not as mass. Under
    the emergency landing the whole liquid reaction (562 kN for the E190 at a
    fill level of 0.90) would be smeared over the shell as a fictitious body
    force of tens of g, and it would cancel the shell's own inertia exactly.
    More fundamentally, a CS-25 load factor presupposes a load path: it
    describes the airframe decelerating the tank, and the stress arises from
    being held. The ring supports model that load path, and their reactions
    are the loads the support structure must carry.

    Inertia relief remains useful for a free-body study with a liquid mass
    model, for example with the liquid as distributed point masses. Abaqus
    requires the relieved directions to be unconstrained, so no displacement
    boundary condition may be applied together with it.
    """
    model.InertiaRelief(name='InertiaRelief', createStepName=step_name,
                        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
    print('  inertia relief active in all six directions (no mounts, no BCs)')


def apply_rigid_body_restraint(model, inst):
    """
    ALTERNATIVE MODEL, KEPT FOR FUTURE STUDIES (not called).

    3-2-1 STATICALLY DETERMINATE RESTRAINT. Valid only when the applied loads
    are self-equilibrated, as in a thermal-only or a uniform-pressure-only
    study. Under an inertia load the whole resultant (743 kN for the E190 in
    the emergency landing) would pass through three nodes of the shell and
    give a meaningless local stress of thousands of MPa.

    This is NOT a mount. It removes the six rigid body modes so that the
    stiffness matrix is non-singular, and nothing else. Because it is
    statically determinate it carries zero reaction under a self-equilibrated
    load, and it does not resist thermal contraction: from 293 K the E190
    would shrink freely by alpha*dT*L = 94.9 mm.

    Node A: fix u1, u2, u3      (3 DOF)
    Node B: fix u2, u3          (2 DOF)  - allows axial growth
    Node C: fix u3              (1 DOF)  - allows axial and one lateral
    """
    nodes = inst.nodes

    def nearest(target):
        best, bd = None, 1e30
        for n in nodes:
            c = n.coordinates
            d = ((c[0] - target[0])**2 + (c[1] - target[1])**2
                 + (c[2] - target[2])**2)
            if d < bd:
                bd, best = d, n
        return best

    # A at the forward pole, B on the barrel at the same azimuth, C at 90 deg
    nA = nearest((-h_cap, 0.0, 0.0))
    nB = nearest((L_barrel * 0.5, 0.0, R_mid))
    nC = nearest((L_barrel * 0.5, R_mid, 0.0))

    setA = model.rootAssembly.Set(name='RB_A', nodes=nodes.sequenceFromLabels((nA.label, )))
    setB = model.rootAssembly.Set(name='RB_B', nodes=nodes.sequenceFromLabels((nB.label, )))
    setC = model.rootAssembly.Set(name='RB_C', nodes=nodes.sequenceFromLabels((nC.label, )))

    model.DisplacementBC(name='RB_A', createStepName='Initial', region=setA,
                         u1=0.0, u2=0.0, u3=0.0)
    model.DisplacementBC(name='RB_B', createStepName='Initial', region=setB,
                         u2=0.0, u3=0.0)
    model.DisplacementBC(name='RB_C', createStepName='Initial', region=setC,
                         u3=0.0)
    print('  3-2-1 restraint applied at nodes %d, %d, %d'
          % (nA.label, nB.label, nC.label))


def job_succeeded(jname):
    """
    Did the Abaqus job finish cleanly?

    WHY mdb.Job.status IS NOT USED. Inside the CAE GUI a job monitor keeps
    the status member of the Job object up to date. A noGUI run has no such
    monitor, so status can still read RUNNING after waitForCompletion() has
    returned, and a job that succeeded would be rejected. The job's own .log
    file is the authority: Abaqus writes the line "JOB <name> COMPLETED"
    there itself, identically in GUI and noGUI runs.
    """
    logf = jname + '.log'
    if not os.path.isfile(logf):
        print('    *** no log file for %s' % jname)
        return False
    fh = open(logf)
    txt = fh.read()
    fh.close()
    if ('JOB ' + os.path.basename(jname) + ' COMPLETED') in txt:
        return True
    for marker in ('Abaqus/Analysis exited with error',
                   'exited with an error',
                   'ABORTED'):
        if marker in txt:
            print('    *** %s reports: %s' % (jname, marker))
            return False
    print('    *** %s did not report COMPLETED' % jname)
    return False


def split_wetted_dry(model, inst, g_hat, s_surf, tag):
    """
    Partition the mesh into the region the liquid touches (wetted) and the
    region exposed to ullage gas (dry).

    PROCESS 3.  THIS IS WHERE THE FILL RATIO ENTERS THE THERMAL PROBLEM. As the
    fill ratio falls, the wetted set shrinks, more of the wall is exposed to
    the poorly conducting ullage gas, and the temperature gradient grows.

    THE SPLIT USES THE FREE SURFACE OF THE LOAD CASE, the same g_hat and
    s_surf as the hydrostatic pressure, so the thermal and structural models
    place the liquid in the same place. An element is wetted when its centroid
    lies below the free surface. Under the emergency landing the resultant
    tilts 57.7 degrees and the liquid moves forward: at a fill level of 0.10
    the wetted zone of the E190 is a wedge in the last 1.9 m of the barrel,
    about 7 percent of the wall, not a strip along the bottom. A horizontal
    split (fill_surface_height) would give the same temperature field for
    every load case.

    Returns SURFACES, not sets. Film conditions and surface heat fluxes act on
    a surface, because a shell has two faces and Abaqus has to be told which
    one the load applies to. Passing a plain element Set raises a region error.

    SIDE CONVENTION: side2 is the INNER face on this revolved geometry, so it
    is side2 that carries the film conditions. The pressure loads use the
    same face, which the sign of the load balance confirms (see the note at
    the Pressure() call).
    """
    wet, dry = [], []
    for el in inst.elements:
        s_c = elem_depth_coord(el, g_hat)
        (wet if s_c >= s_surf else dry).append(el.label)

    surfs = {}
    if wet:
        els = inst.elements.sequenceFromLabels(tuple(wet))
        surfs['wet'] = model.rootAssembly.Surface(name='Wet_%s' % tag,
                                                  side2Elements=els)
    if dry:
        els = inst.elements.sequenceFromLabels(tuple(dry))
        surfs['dry'] = model.rootAssembly.Surface(name='Dry_%s' % tag,
                                                  side2Elements=els)
    print('    wetted elements %d, dry elements %d' % (len(wet), len(dry)))
    return surfs


# ============================================================================
# 4. THE ANALYSIS
# ============================================================================

def _usable_ht(path):
    """
    Return `path` if it is a finished heat-transfer job at the CURRENT wall,
    else None. Used to reuse a thermal result instead of solving it again.
    """
    if not path or not os.path.isfile(path + '.odb') or not job_succeeded(path):
        return None
    try:
        o_ = openOdb(path=path + '.odb', readOnly=True)
        tk = None
        for sec in o_.sections.values():
            tk = getattr(sec, 'thickness', None) or tk
        o_.close()
        if tk is not None and abs(float(tk) - t_wall) > 1e-9:
            return None
    except Exception:
        return None
    return path


def run_case(lc_name, factors, phi_f, seed, tag, ht_from=None):
    """
    One design point (load case, fill level). Runs the heat transfer
    analysis, then the static analysis reading its temperature field, and
    returns the dictionary of results from extract(), or None if a job fails.

        lc_name   load-case code, for example 'LC7_Emergency'
        factors   {'nx': ..., 'ny': ..., 'nz': ...}, load factors in g
        phi_f     fill level, 0 to 1
        seed      element size [m]
        tag       job name stem: jobs HT_<tag> and ST_<tag>
        ht_from   optional path (without .odb) of an existing heat transfer
                  job for the SAME wall, mesh, load direction and fill.
                  Pressure does not enter the heat transfer problem, so the
                  temperature field at the stability and at the design
                  pressure is identical, and a usable ht_from is read instead
                  of being solved again. The static model checks that its
                  mesh matches it node for node and stops if it does not.

    Loads and boundary conditions use the module-level inputs of section 1
    (P_gauge, h_wetted, h_dry, q_outer, T_ref, ...). This is the place to add
    a new load or boundary condition.
    """
    reuse_t = _usable_ht(ht_from) if ht_from else None
    nx, ny, nz = factors['nx'], factors['ny'], factors['nz']
    n_mag = math.sqrt(nx * nx + ny * ny + nz * nz)

    # THE FREE SURFACE IS COMPUTED ONCE, HERE, AND SHARED BY BOTH MODELS.
    # The thermal model uses it to decide which wall is wetted; the structural
    # model uses it to decide the hydrostatic head. They must be the same
    # surface, or the analysis is coupled to two different tanks.
    g_hat = effective_down(nx, ny, nz)
    s_surf = free_surface_offset(g_hat, phi_f)

    mname_t = None
    if reuse_t:
        jname_t = reuse_t
        print('    thermal field reused from %s (pressure does not enter the '
              'heat transfer problem)' % os.path.basename(reuse_t))
    else:
        # ------ STEP 1: steady-state heat transfer  (PROCESS 1, 2, 4) ---------
        # Built as its OWN model, with diffusive shell elements. The structural
        # model below is built separately with stress elements; a single mesh
        # cannot carry both element families.
        mname_t = 'MT_%s' % tag
        model_t, part_t, inst_t = build_geometry_and_mesh(mname_t, seed,
                                                          analysis_type='THERMAL')
        surfs = split_wetted_dry(model_t, inst_t, g_hat, s_surf, tag)

        model_t.HeatTransferStep(name='Thermal', previous='Initial',
                                 response=STEADY_STATE, amplitude=RAMP)

        # FORMULA 4 on both inner-surface regions, with very different h.
        # These are the only temperatures prescribed anywhere in the analysis:
        # the two SINKS. Every wall temperature is solved for (PROCESS 1).
        if 'wet' in surfs:
            model_t.FilmCondition(name='Film_wet', createStepName='Thermal',
                                  surface=surfs['wet'], definition=EMBEDDED_COEFF,
                                  filmCoeff=h_wetted, sinkTemperature=T_LH2)
        if 'dry' in surfs:
            model_t.FilmCondition(name='Film_dry', createStepName='Thermal',
                                  surface=surfs['dry'], definition=EMBEDDED_COEFF,
                                  filmCoeff=h_dry, sinkTemperature=T_GH2)

        # FORMULA 5: heat arriving through the insulation, on the OUTER face
        # (side1). The inner face (side2) carries the film conditions above.
        outer_surf = model_t.rootAssembly.Surface(
            name='Outer_%s' % tag, side1Elements=inst_t.elements[:])
        model_t.SurfaceHeatFlux(name='MLI_flux', createStepName='Thermal',
                                region=outer_surf, magnitude=q_outer)

        model_t.FieldOutputRequest(name='F_thermal', createStepName='Thermal',
                                   variables=('NT', 'HFL'))

        jname_t = 'HT_%s' % tag
        jt = mdb.Job(name=jname_t, model=mname_t, type=ANALYSIS)
        jt.submit(consistencyChecking=OFF)
        jt.waitForCompletion()
        if not job_succeeded(jname_t):
            print('    *** heat transfer job failed: %s' % tag)
            return None

    # ------ STEP 2: static structural  (PROCESS 5) -------------------------
    mname = 'MS_%s' % tag
    model, part, inst = build_geometry_and_mesh(mname, seed,
                                                analysis_type='STRUCTURAL')
    all_el = model.rootAssembly.Set(name='AllEl_%s' % tag,
                                    elements=inst.elements[:])

    model.StaticStep(name='Static', previous='Initial', timePeriod=1.0)

    # PROCESS 6. Three circumferential rings after Bagarello Sec. 4.4.1. The
    # applied loads are reacted at the rings, which is what a CS-25 load
    # factor physically describes. (apply_inertia_relief is the free-body
    # alternative, not used.)
    apply_ring_constraints(model, inst, 'Static', seed)

    # INITIAL TEMPERATURE: THE STRESS-FREE STATE. DO NOT REMOVE.
    #
    # Abaqus computes thermal strain from the INITIAL temperature, not only
    # from the Expansion 'zero' value. Without an initial condition the
    # initial temperature defaults to 0, and the model would treat the vessel
    # as warmed from 0 K.
    #
    # T_ref is the COLD operating state, 20.3 K (see the block where it is
    # defined). The structure is therefore stress-free when cold, and the only
    # thermal load is the wetted-to-dry gradient, which is what couples the
    # fill ratio to the structural response.
    all_nodes = model.rootAssembly.Set(name='AllNodes_%s' % tag,
                                       nodes=inst.nodes[:])
    model.Temperature(name='T_initial', createStepName='Initial',
                      region=all_nodes, distributionType=UNIFORM,
                      crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                      magnitudes=(T_ref, ))

    # ------------------------------------------------------------------
    # TEMPERATURE FIELD.  FORMULA 8
    #
    # The nodal temperatures are read from the heat transfer .odb, grouped
    # into 24 bins, and applied as UNIFORM fields with
    # crossSectionDistribution=CONSTANT_THROUGH_THICKNESS. The bin width is
    # (T_max - T_min)/24, 0.14 to 0.31 K for the thesis tanks.
    #
    # WHY NOT distributionType=FROM_FILE. The shell section is declared
    # temperature=GRADIENT, so every temperature supplied to it is a PAIR: a
    # mid-surface value and a through-thickness gradient. A predefined field
    # read FROM_FILE does not reliably deliver a zero gradient, and a spurious
    # through-thickness difference of about 20 K was observed, which roughly
    # doubled the stresses through bending. The real through-thickness drop
    # of a thin aluminium wall carrying a few W/m2 is of order 1e-4 K.
    # CONSTANT_THROUGH_THICKNESS sets the gradient to zero by construction.
    #
    # ACCEPTANCE TEST after any change here: in a dry region far from the
    # rings, both section points (SNEG, SPOS) must give the membrane values,
    # hoop pR/t and axial pR/2t. extract() records this as the probe_*
    # columns of the results file.
    # ------------------------------------------------------------------
    odb_t = openOdb(path=jname_t + '.odb', readOnly=True)
    t_frame = odb_t.steps['Thermal'].frames[-1]
    t_nodal = {}
    for v in t_frame.fieldOutputs['NT11'].values:
        t_nodal[v.nodeLabel] = v.data
    odb_t.close()

    if not t_nodal:
        raise ValueError('NT11 was empty in %s - the thermal job did not '
                         'produce a temperature field.' % (jname_t + '.odb'))

    node_T = []
    n_missing = 0
    for nd in inst.nodes:
        if nd.label in t_nodal:
            node_T.append((nd.label, t_nodal[nd.label]))
        else:
            node_T.append((nd.label, T_ref))
            n_missing += 1
    if n_missing:
        print('    *** %d structural nodes had no thermal counterpart and '
              'were set to T_ref' % n_missing)
    if reuse_t and (n_missing or len(t_nodal) != len(inst.nodes)):
        raise RuntimeError('reused thermal field %s does not match this mesh '
                           '(%d thermal nodes, %d structural, %d missing)'
                           % (reuse_t, len(t_nodal), len(inst.nodes), n_missing))

    T_lo = min(v for (_, v) in node_T)
    T_hi = max(v for (_, v) in node_T)
    n_tbins = 24
    T_width = (T_hi - T_lo) / n_tbins
    if T_width <= 0.0:
        T_width = 1e-9

    tbins = {}
    for (lbl, tv) in node_T:
        idx = int((tv - T_lo) / T_width)
        idx = max(0, min(n_tbins - 1, idx))
        tbins.setdefault(idx, []).append(lbl)

    for idx, labels in tbins.items():
        T_c = T_lo + (idx + 0.5) * T_width
        tset = model.rootAssembly.Set(
            name='TB_%s_%02d' % (tag, idx),
            nodes=inst.nodes.sequenceFromLabels(tuple(labels)))
        model.Temperature(name='T_%s_%02d' % (tag, idx),
                          createStepName='Static', region=tset,
                          distributionType=UNIFORM,
                          crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                          magnitudes=(T_c, ))

    print('    temperature: %.3f to %.3f K applied in %d bins of %.4f K '
          '(constant through thickness)'
          % (T_lo, T_hi, len(tbins), T_width))

    # FORMULA 6: internal pressure plus hydrostatic head, binned by DEPTH
    # ALONG THE RESULTANT LOAD DIRECTION, not by vertical height.
    #
    # g_hat and s_surf are NOT recomputed here. They are taken from the top of
    # run_case, where they were computed once and handed to the thermal model
    # as well, so both models see the same free surface.
    #
    # The elements are grouped into n_bins bins of depth along g_hat, and each
    # bin carries the uniform pressure at its mid-depth.
    s_list = []
    for el in inst.elements:
        s_list.append((el.label, elem_depth_coord(el, g_hat)))
    s_only = [v for (_, v) in s_list]
    s_lo, s_hi = min(s_only), max(s_only)
    width = (s_hi - s_lo) / n_bins
    if width <= 0.0:
        width = 1e-6

    tilt_deg = math.degrees(math.acos(min(1.0, abs(g_hat[1]))))
    print('    free surface: tilt %.1f deg from vertical, max depth %.3f m'
          % (tilt_deg, max(0.0, s_hi - s_surf)))

    bins = {}
    for (lbl, s_v) in s_list:
        idx = int((s_v - s_lo) / width)
        idx = max(0, min(n_bins - 1, idx))
        bins.setdefault(idx, []).append(lbl)

    for idx, labels in bins.items():
        s_c = s_lo + (idx + 0.5) * width
        depth = max(0.0, s_c - s_surf)
        p_tot = P_gauge + rho_LH2 * (n_mag * g) * depth       # FORMULA 6
        bset = model.rootAssembly.Set(
            name='PB_%s_%02d' % (tag, idx),
            elements=inst.elements.sequenceFromLabels(tuple(labels)))
        # side2 IS THE INNER FACE on this revolved geometry. Confirmed by the
        # load balance: on side1 the hydrostatic resultant has the right
        # magnitude but the opposite sign (the liquid pushes inward).
        bsurf = model.rootAssembly.Surface(name='PS_%s_%02d' % (tag, idx),
                                           side2Elements=bset.elements)
        model.Pressure(name='P_%s_%02d' % (tag, idx), createStepName='Static',
                       region=bsurf, magnitude=p_tot,
                       distributionType=UNIFORM)

    # FORMULA 7: inertial body force on the shell's own mass. The inertia of
    # the liquid is carried by the scaled hydrostatic pressure above, not here.
    #
    # AXIS MAPPING. The vessel is revolved about global X, so:
    #     comp1 = global X = FORWARD  -> n_x
    #     comp2 = global Y = UP       -> n_z
    #     comp3 = global Z = LATERAL  -> n_y
    #
    # SIGNS. The CS 25.561 load factors are named by the direction in which
    # the inertia force ACTS ON THE ITEM. "Downward 6.0 g" is a force downward
    # and +Y is up, so comp2 = -nz*g. "Forward 9.0 g" is a force forward and
    # +X is forward, so comp1 = +nx*g. The opposite signs come from the axis
    # naming, not from the physics. The liquid therefore collects against the
    # forward cap in the emergency landing.
    #
    # UNITS. Gravity takes an ACCELERATION and multiplies it by the material
    # density itself. (BodyForce would expect force per unit volume, rho*a.)
    model.Gravity(name='Inertia', createStepName='Static',
                  comp1=nx * g, comp2=-nz * g, comp3=-ny * g,
                  distributionType=UNIFORM)

    model.FieldOutputRequest(name='F_static', createStepName='Static',
                             variables=('S', 'U', 'RF', 'E', 'NT'))

    jname_s = 'ST_%s' % tag
    js = mdb.Job(name=jname_s, model=mname, type=ANALYSIS)
    js.submit(consistencyChecking=OFF)
    js.waitForCompletion()
    if not job_succeeded(jname_s):
        print('    *** static job failed: %s' % tag)
        return None

    res = extract(jname_t, jname_s, lc_name, phi_f, seed, factors)
    _free(mname_t, mname, jname_t, jname_s)
    return res


def _free(*names):
    """
    Remove finished models and jobs from the CAE session. Each design point
    builds two models with a full mesh and hundreds of load sets; twenty fit
    in memory, the roughly 70 of a FINAL session may not. The .odb files on
    disk are untouched.
    """
    for n in names:
        for repo in (mdb.jobs, mdb.models):
            try:
                if n in repo.keys():
                    del repo[n]
            except Exception:
                pass


def extract(jname_t, jname_s, lc_name, phi_f, seed, factors):
    """
    Read the results of one design point from the two output databases and
    return them as a dictionary (one row of the results file; the columns
    are listed in BASE_KEYS and EXTRA_KEYS and described in the README).

    Thermal     T_min, T_max, dT, peak wall heat flux
    Strength    far-field von Mises stress (FORMULA 9, evaluated by the
                solver) and the reserve factor fos = sigma_allow / sigma
    Stability   minimum axial membrane stress of the barrel against the
                SP-8007 allowable: buckling_fos
    Checks      stress orientation, load balance (reaction against applied
                load), membrane probe, ovality, reaction at each ring,
                membrane stress in the caps

    This is the place to add a new result: compute it here, add it to the
    returned dictionary and its key to EXTRA_KEYS.
    """
    odb_t = openOdb(path=jname_t + '.odb')
    fr_t = odb_t.steps['Thermal'].frames[-1]
    temps = [v.data for v in fr_t.fieldOutputs['NT11'].values]
    T_min, T_max = min(temps), max(temps)
    # Peak wall heat flux, so the waterline flux quoted in the thesis is a
    # recorded number, not one read off a colour legend.
    hfl_max = 0.0
    try:
        for v in fr_t.fieldOutputs['HFL'].values:
            d_ = v.data
            m_ = math.sqrt(sum(c * c for c in d_))
            if m_ > hfl_max:
                hfl_max = m_
    except Exception:
        hfl_max = float('nan')
    odb_t.close()

    odb_s = openOdb(path=jname_s + '.odb')
    fr_s = odb_s.steps['Static'].frames[-1]
    inst_s = odb_s.rootAssembly.instances['VESSELINST']

    umag = [v.magnitude for v in fr_s.fieldOutputs['U'].values]
    u_max = max(umag)

    # ------------------------------------------------------------------
    # FAR-FIELD PEAK STRESS.  THIS IS THE NUMBER THE THESIS REPORTS.
    # ------------------------------------------------------------------
    # A rigid displacement boundary condition applied along a line of nodes
    # on a shell is a mathematical idealisation with zero contact area. The
    # stress there is SINGULAR: it grows without limit as the mesh is
    # refined and never converges. In the convergence study of the E190
    # (governing case, seeds 41.0, 20.5 and 10.2 mm) the peak including the
    # rings goes 271.5 -> 256.3 -> 295.1 MPa, rising again at the finest
    # mesh, while the far-field peak goes 271.5 -> 256.3 -> 256.8 MPa and the
    # total reaction is constant at 614.2 kN to 0.03 percent. The divergence
    # is therefore local to the constraint, not a structural response
    # (results/mesh_convergence/).
    #
    # The peak is consequently taken OUTSIDE an exclusion zone around each
    # ring. The justification is Saint-Venant's principle: a local load
    # introduction decays over a distance comparable with the shell decay
    # length pi/beta, 0.10, 0.13 and 0.30 m for the sized E190, A320 and
    # A350-900 walls. The 0.5 m zone adopted spans 5.0, 3.9 and 1.7 decay
    # lengths, where the edge disturbance has fallen to 2e-7, 4e-6 and 5e-3
    # of its value at the ring, so the reported stress is a genuine far-field
    # value governed by global equilibrium rather than by the details of
    # the attachment.
    #
    # This follows Oom Ortiz de Montellano et al., who documented the same
    # artefact for constrained boundary conditions, and is the treatment
    # set out in thesis Section 4.7.
    #
    # WHAT THIS EXCLUDES: the local fitting design at the rings, which in a
    # detailed design would need a doubler or wear plate (a stated
    # limitation of the thesis).
    ring_x = [0.0, L_barrel / 2.0, L_barrel]

    # element label -> centroid x, taken from the odb mesh.
    #
    # PERFORMANCE. The node coordinates are read into a dictionary once.
    # Calling getNodeFromLabel() per node per element would be about a
    # million lookups on the finest mesh, minutes in the Abaqus interpreter.
    ncoord = {}
    nxyz = {}
    for nd in inst_s.nodes:
        ncoord[nd.label] = nd.coordinates[0]
        nxyz[nd.label] = nd.coordinates

    cx = {}
    for el in inst_s.elements:
        xs_ = 0.0
        cnt = 0
        for nl in el.connectivity:
            xs_ += ncoord[nl]
            cnt += 1
        cx[el.label] = xs_ / float(cnt)

    s_field = fr_s.fieldOutputs['S']
    s_all = []
    s_far = []
    for v in s_field.values:
        s_all.append(v.mises)
        x_ = cx.get(v.elementLabel, None)
        if x_ is None:
            continue
        near_ring = False
        for xr in ring_x:
            if abs(x_ - xr) < FAR_FIELD_EXCLUSION:
                near_ring = True
                break
        if not near_ring:
            s_far.append(v.mises)

    s_max_all = max(s_all) if s_all else 0.0
    s_max = max(s_far) if s_far else s_max_all

    # ---------------- AXIAL MEMBRANE STRESS FOR THE BUCKLING CRITERION ----
    # Barrel only (SP-8007 is a cylinder result), far field only (outside
    # FAR_FIELD_EXCLUSION of every ring).
    #
    # S11 IS AXIAL HERE. No material orientation is assigned, so Abaqus takes
    # local direction 1 as the projection of global X onto the shell, and X is
    # the tank axis. (With the tank axis along Z the same default would make
    # S11 the HOOP stress at the crown and keel.) This is checked, not
    # assumed: under internal pressure the mean hoop membrane stress must
    # exceed the mean axial one, so if the mean of S11 exceeds the mean of
    # S22 the orientation is wrong and the run says so (orientation_ok).
    #
    # Membrane = mean of the two section points (SNEG, SPOS); bending, the
    # half-difference, is excluded, as the SP-8007 criterion consumes membrane
    # stress.
    #
    # The caps are recorded separately for a head stability check: under
    # internal pressure an a/b = 1.6 head carries hoop COMPRESSION near its
    # equator, (pR/t)(1 - a^2/2b^2) in membrane theory.
    sp_vals = {}
    cap_vals = {}
    for v in s_field.values:
        x_ = cx.get(v.elementLabel, None)
        if x_ is None:
            continue
        in_barrel = (FAR_FIELD_EXCLUSION <= x_ <= L_barrel / 2.0 - FAR_FIELD_EXCLUSION or
                     L_barrel / 2.0 + FAR_FIELD_EXCLUSION <= x_ <= L_barrel - FAR_FIELD_EXCLUSION)
        if x_ < 0.0 or x_ > L_barrel:
            # end caps: recorded for the head stability check
            spn_ = v.sectionPoint.number if v.sectionPoint is not None else 0
            cap_vals.setdefault(v.elementLabel, []).append((spn_, v.data[0], v.data[1]))
            continue
        if not in_barrel:
            continue
        spn = v.sectionPoint.number if v.sectionPoint is not None else 0
        sp_vals.setdefault(v.elementLabel, []).append((spn, v.data[0], v.data[1]))
    ax_min = float('inf')
    sum11, sum22, n_mem = 0.0, 0.0, 0
    for lbl, pts in sp_vals.items():
        pts.sort()
        lo_, hi_ = pts[0], pts[-1]
        m11 = 0.5 * (lo_[1] + hi_[1])
        m22 = 0.5 * (lo_[2] + hi_[2])
        sum11 += m11
        sum22 += m22
        n_mem += 1
        if m11 < ax_min:
            ax_min = m11
    cap11_min, cap22_min = float('inf'), float('inf')
    for lbl, pts in cap_vals.items():
        pts.sort()
        cap11_min = min(cap11_min, 0.5 * (pts[0][1] + pts[-1][1]))
        cap22_min = min(cap22_min, 0.5 * (pts[0][2] + pts[-1][2]))
    mean11 = sum11 / n_mem if n_mem else 0.0
    mean22 = sum22 / n_mem if n_mem else 0.0
    orient_ok = mean22 > mean11
    s_cr = sigma_cr_sp8007(R_mid, t_wall)
    buck_fos = (s_cr / -ax_min) if ax_min < 0 else float('inf')

    # ---------------- RECORDED CHECKS ------------------------------------
    # (1) MEMBRANE AND THROUGH-THICKNESS PROBE. The crown element nearest
    #     x = L/4 (mid aft bay, at least 1 m from any ring on all three tanks)
    #     and its S11 (axial) and S22 (hoop) at both section points. Under a
    #     vertical load at low fill the crown is dry, so hoop ~ pR/t and
    #     axial ~ pR/2t, and the two surfaces must agree (no fictitious
    #     through-thickness bending).
    ecen = {}
    for el in inst_s.elements:
        xs_, ys_, zs_, c_ = 0.0, 0.0, 0.0, 0
        for nl in el.connectivity:
            q_ = nxyz[nl]
            xs_ += q_[0]; ys_ += q_[1]; zs_ += q_[2]; c_ += 1
        ecen[el.label] = (xs_ / c_, ys_ / c_, zs_ / c_)
    x_probe = L_barrel / 4.0
    probe_lbl = min(ecen, key=lambda k: (ecen[k][0] - x_probe) ** 2 +
                    (ecen[k][1] - R_mid) ** 2 + ecen[k][2] ** 2)
    pr = sorted([(v.sectionPoint.number if v.sectionPoint is not None else 0,
                  v.data[0], v.data[1]) for v in s_field.values
                 if v.elementLabel == probe_lbl])
    if pr:
        probe = dict(probe_S11_sneg_MPa=pr[0][1] / 1e6, probe_S11_spos_MPa=pr[-1][1] / 1e6,
                     probe_S22_sneg_MPa=pr[0][2] / 1e6, probe_S22_spos_MPa=pr[-1][2] / 1e6)
    else:
        probe = dict(probe_S11_sneg_MPa=float('nan'), probe_S11_spos_MPa=float('nan'),
                     probe_S22_sneg_MPa=float('nan'), probe_S22_spos_MPa=float('nan'))

    # (2) OVALITY between the rings. Radial displacement u_r = (u_y y + u_z z)/r
    #     around the section at the middle of each bay (x = L/4 and 3L/4).
    #     Its range is the out-of-roundness; the rings hold it at zero where
    #     they act. Reported for the worse bay, with u_r at the keel (lowest
    #     node) and at the flank (largest |z|).
    u_by_node = {}
    for v in fr_s.fieldOutputs['U'].values:
        u_by_node[v.nodeLabel] = v.data
    tol_x = 1.0 * seed      # a full element either side: always catches a node ring
    oval, ur_keel, ur_flank = 0.0, float('nan'), float('nan')
    for xm in (L_barrel / 4.0, 3.0 * L_barrel / 4.0):
        sec = []
        for lbl, q_ in nxyz.items():
            if abs(q_[0] - xm) < tol_x and lbl in u_by_node:
                r_ = math.sqrt(q_[1] * q_[1] + q_[2] * q_[2])
                if r_ > 0.5 * R_mid:
                    u_ = u_by_node[lbl]
                    sec.append((q_[1], abs(q_[2]), (u_[1] * q_[1] + u_[2] * q_[2]) / r_))
        if len(sec) > 4:
            urs = [c[2] for c in sec]
            rng = max(urs) - min(urs)
            if rng > oval:
                oval = rng
                ur_keel = min(sec, key=lambda c: c[0])[2]
                ur_flank = max(sec, key=lambda c: c[1])[2]

    # (3) REACTION AT EACH RING, so the load split between the three rings is
    #     a recorded result rather than a beam-theory estimate. The total
    #     reaction (sumRF_N) is the load the support structure must carry.
    ring_rf = [[0.0, 0.0, 0.0] for _ in ring_x]
    for v in fr_s.fieldOutputs['RF'].values:
        q_ = nxyz.get(v.nodeLabel)
        if q_ is None:
            continue
        # Only restrained nodes carry a reaction, so every node with a
        # non-zero RF belongs to a ring; assign it to the nearest station.
        # No distance tolerance: the mid ring snaps to the nearest node plane,
        # which can lie up to one element from L/2.
        if max(abs(c) for c in v.data) < 1e-9:
            continue
        k_ = min(range(len(ring_x)), key=lambda i: abs(q_[0] - ring_x[i]))
        for i in range(3):
            ring_rf[k_][i] += v.data[i]
    ring_mag = [math.sqrt(sum(c * c for c in r_)) / 1e3 for r_ in ring_rf]
    n_excl = len(s_all) - len(s_far)
    rf = fr_s.fieldOutputs['RF']
    rf_sum = [0.0, 0.0, 0.0]
    for v in rf.values:
        for i in range(3):
            rf_sum[i] += v.data[i]
    n_el = len(odb_s.rootAssembly.instances['VESSELINST'].elements)
    n_nd = len(odb_s.rootAssembly.instances['VESSELINST'].nodes)
    odb_s.close()

    fos = sigma_allow / s_max if s_max > 0 else float('inf')
    rf_sq = 0.0
    for c in rf_sum:
        rf_sq += c * c
    rf_mag = math.sqrt(rf_sq)

    # ---------------- LOAD BALANCE SELF-CHECK ----------------------------
    # Equilibrium demands that the reaction equal the applied resultant. The
    # applied load has exactly two parts, both computable in closed form:
    #
    #   gravity on the shell  : m_shell * a, along the load-factor direction
    #   hydrostatic resultant : rho_LH2 * |n| * g * V_liquid, along the load
    #                           direction, because depth is measured along it
    #
    # The uniform p_gauge term nets to zero on a closed surface and does not
    # appear. If it DOES appear, the shell normals are inconsistent and some
    # elements are being pressurised on the wrong face, which would corrupt
    # the stress field as well as the balance. A ratio outside 1 +/- 0.02 is
    # printed with its components, and stops a sizing run.
    nx_, ny_, nz_ = factors['nx'], factors['ny'], factors['nz']
    n_mag_ = math.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
    # Use the mass Abaqus actually meshed rather than the analytic estimate,
    # so the check compares the model with itself.
    m_shell = A_shell * t_wall * rho_shell
    if MESHED_MASS is not None and MESHED_MASS > 0.0:
        m_shell = MESHED_MASS
    m_liq = V_tank * rho_LH2 * phi_f

    # Global axes: 1 = X forward, 2 = Y up, 3 = Z lateral.
    #
    # THE LIQUID RESULTANT ACTS ALONG THE LOAD DIRECTION, NOT VERTICALLY.
    # With the liquid settling along g_hat = (-nx, -nz, -ny)/|n| and its
    # resultant being m_liq * |n| * g in that direction, the |n| cancels and
    # the applied load collapses to the obvious result: the total mass times
    # the load factor, component by component.
    #
    #     applied = -(nx, nz, ny) * g * (m_shell + m_liq)
    m_tot_ = m_shell + m_liq
    app = [nx_ * g * m_tot_,
           -nz_ * g * m_tot_,
           -ny_ * g * m_tot_]
    app_sq = 0.0
    for c in app:
        app_sq += c * c
    app_mag = math.sqrt(app_sq)
    ratio = rf_mag / app_mag if app_mag > 0 else float('nan')

    print('    T %.3f - %.3f K (dT %.3f) | vM far-field %.2f MPa | FOS %.2f | '
          'u_max %.2f mm | R_mount %.1f kN | %d el'
          % (T_min, T_max, T_max - T_min, s_max / 1e6, fos,
             u_max * 1000.0, rf_mag / 1000.0, n_el))
    print('      barrel axial membrane: min %.2f MPa | mean S11 %.2f  mean S22 %.2f MPa '
          '| sigma_cr %.2f MPa | buckling FOS %.2f'
          % (ax_min / 1e6, mean11 / 1e6, mean22 / 1e6, s_cr / 1e6, buck_fos))
    if not orient_ok:
        print('      *** ORIENTATION CHECK FAILED: mean S11 exceeds mean S22, so S11')
        print('      *** is not the axial stress. The buckling screen is INVALID.')
    print('      peak at the rings %.2f MPa is SINGULAR and discarded '
          '(%d points excluded within %.2f m)'
          % (s_max_all / 1e6, n_excl, FAR_FIELD_EXCLUSION))
    print('      balance: applied (%8.1f,%8.1f,%8.1f) kN  |R| = %7.1f kN'
          % (app[0] / 1e3, app[1] / 1e3, app[2] / 1e3, app_mag / 1e3))
    print('               reacted (%8.1f,%8.1f,%8.1f) kN  |R| = %7.1f kN'
          % (rf_sum[0] / 1e3, rf_sum[1] / 1e3, rf_sum[2] / 1e3, rf_mag / 1e3))
    if abs(ratio - 1.0) > 0.02:
        print('      *** OUT OF BALANCE by %.1f %%. Component-wise ratios:'
              % (100.0 * (ratio - 1.0)))
        for i, ax in enumerate(['X fwd', 'Y up ', 'Z lat']):
            r = rf_sum[i] / app[i] if abs(app[i]) > 1.0 else float('nan')
            print('          %s applied %9.1f kN   reacted %9.1f kN   '
                  'ratio %6.3f' % (ax, app[i] / 1e3, rf_sum[i] / 1e3, r))
    else:
        print('      balance OK (within 2 %%)')

    return dict(load_case=lc_name, fill=phi_f, seed_mm=seed * 1000.0,
                mises_singular_MPa=s_max_all / 1e6,
                n_el=n_el, n_nd=n_nd, T_min=T_min, T_max=T_max,
                dT=T_max - T_min, mises_MPa=s_max / 1e6, fos=fos,
                u_max_mm=u_max * 1000.0, sumRF_N=rf_mag,
                t_wall_mm=t_wall * 1000.0, p_gauge_bar=P_gauge / 1.0e5,
                axial_mem_min_MPa=ax_min / 1e6, sigma_cr_MPa=s_cr / 1e6,
                buckling_fos=buck_fos, orientation_ok=orient_ok,
                balance_ratio=ratio, hfl_max_Wm2=hfl_max,
                mean_S11_mem_MPa=mean11 / 1e6, mean_S22_mem_MPa=mean22 / 1e6,
                ovality_mm=oval * 1000.0, ur_keel_mm=ur_keel * 1000.0,
                ur_flank_mm=ur_flank * 1000.0, RF_aft_kN=ring_mag[0],
                RF_mid_kN=ring_mag[1], RF_fwd_kN=ring_mag[2],
                cap_S11_mem_min_MPa=cap11_min / 1e6, cap_S22_mem_min_MPa=cap22_min / 1e6,
                **probe)


# ============================================================================
# 5. DRIVER
# ============================================================================
#
# RUN MODES
#   PRODUCTION   every load case at every fill level (thesis: 4 x 5) at the
#                design pressure (strength).
#   STABILITY    the same matrix at the stability pressure (buckling).
#   CONVERGENCE  the governing case (thesis: emergency landing at fill 0.90)
#                over a sweep of mesh densities, at the design pressure AND
#                at the stability pressure.
#   SIZING       finite-element wall sizing, recorded in sized_walls.csv.
#   FINAL        SIZING -> PRODUCTION -> STABILITY -> CONVERGENCE in one
#                session. Results are written after every design point; a
#                stage whose file is already complete at the current wall is
#                skipped, and an interrupted stage resumes where it stopped.
#
# Results files: one row per design point, columns BASE_KEYS + EXTRA_KEYS
# (a sizing file has a leading 'stage' column). A new result is added to the
# dictionary returned by extract() and its key appended to EXTRA_KEYS.

BASE_KEYS = ['load_case', 'fill', 'seed_mm', 'n_el', 'n_nd', 'T_min', 'T_max',
             'dT', 'mises_MPa', 'mises_singular_MPa', 'fos', 'u_max_mm', 'sumRF_N']
EXTRA_KEYS = ['t_wall_mm', 'p_gauge_bar', 'axial_mem_min_MPa', 'sigma_cr_MPa',
              'buckling_fos', 'orientation_ok', 'balance_ratio', 'hfl_max_Wm2',
              'mean_S11_mem_MPa', 'mean_S22_mem_MPa', 'ovality_mm', 'ur_keel_mm',
              'ur_flank_mm', 'RF_aft_kN', 'RF_mid_kN', 'RF_fwd_kN',
              'probe_S11_sneg_MPa', 'probe_S11_spos_MPa',
              'probe_S22_sneg_MPa', 'probe_S22_spos_MPa',
              'cap_S11_mem_min_MPa', 'cap_S22_mem_min_MPa']
SIZED_WALLS = os.path.join(_ROOT_DIR, 'sized_walls.csv')
P_DESIGN = P_gauge                     # strength pressure, as loaded above


def write_rows(path, rows, lead=()):
    """Write result dictionaries as a CSV file with the standard columns."""
    keys = list(lead) + BASE_KEYS + EXTRA_KEYS
    fh = open(path, 'w')
    fh.write(','.join(keys) + '\n')
    for r in rows:
        # commas inside a value (e.g. a stage label) would shift the columns
        fh.write(','.join([str(r.get(k, '')).replace(',', ';') for k in keys]) + '\n')
    fh.close()
    print('\n  written %s (%d rows)' % (path, len(rows)))


def stage_done(path, n_expected):
    """A stage is complete if its CSV holds n_expected rows at the current wall."""
    if not os.path.exists(path):
        return False
    rows = list(csv.DictReader(open(path)))
    good = [r for r in rows if r.get('t_wall_mm') and
            abs(float(r['t_wall_mm']) - t_wall * 1e3) < 1e-6]
    return len(good) >= n_expected


def _typed(row):
    """A CSV row with its numbers and booleans converted back from text."""
    out = {}
    for k, v in row.items():
        if v in ('True', 'False'):
            out[k] = (v == 'True')
        else:
            try:
                out[k] = float(v)
            except (TypeError, ValueError):
                out[k] = v
    return out


def done_rows(path, pressure, pattern=None):
    """
    Rows already computed at the current wall and this pressure, for the load
    cases of this session. `pattern` (a glob) lets a session pick up cases
    finished by ANY window of the same aircraft, not only its own file.
    """
    import glob as _glob
    paths = sorted(_glob.glob(pattern), key=os.path.getmtime) if pattern else [path]
    keep = []
    for pth in paths:
        if os.path.exists(pth):
            for r in csv.DictReader(open(pth)):
                if r.get('load_case') in load_cases or 'load_case' not in r:
                    keep.append(r)
    rows_in = keep
    keep = []
    for r in rows_in:
        try:
            if abs(float(r['t_wall_mm']) - t_wall * 1e3) < 1e-6 and \
                    abs(float(r['p_gauge_bar']) - pressure / 1e5) < 1e-6:
                keep.append(_typed(r))
        except (KeyError, ValueError):
            pass
    return keep


def set_wall(t_new):
    """Change the wall thickness and the mesh seed that follows from it."""
    global t_wall, boundary_layer, seed_size
    t_wall = t_new
    boundary_layer = math.sqrt(R_mid * t_wall)
    seed_size = boundary_layer * MESH_FRACTION


def record_sized_wall(fos_b, fos_y, govern=''):
    """Write or replace this aircraft's row of sized_walls.csv."""
    rows = []
    if os.path.exists(SIZED_WALLS):
        rows = [r for r in csv.DictReader(open(SIZED_WALLS)) if r.get('aircraft') != CASE_TAG]
    rows.append(dict(aircraft=CASE_TAG, t_wall_mm='%.6f' % (t_wall * 1e3),
                     buckling_fos_min='%.4f' % fos_b, yield_fos_lc7_f90='%.4f' % fos_y,
                     governed_by=govern.replace(',', ';'),
                     stability_pressure_bar='%.5f' % (P_STABILITY / 1e5),
                     design_pressure_bar='%.5f' % (P_DESIGN / 1e5)))
    keys = ['aircraft', 't_wall_mm', 'governed_by', 'buckling_fos_min', 'yield_fos_lc7_f90',
            'stability_pressure_bar', 'design_pressure_bar']
    fh = open(SIZED_WALLS, 'w')
    fh.write(','.join(keys) + '\n')
    for r in sorted(rows, key=lambda r: r['aircraft']):
        fh.write(','.join([r.get(k, '') for k in keys]) + '\n')
    fh.close()
    print('  sized wall recorded in %s' % SIZED_WALLS)


def in_dir(sub):
    """Change to a sub-folder of the results folder, creating it if needed."""
    d = os.path.join(_ROOT_DIR, sub)
    if not os.path.isdir(d):
        os.makedirs(d)
    os.chdir(d)
    return d


def stage_sizing():
    """
    PROCESS 8. FINITE-ELEMENT WALL SIZING.

    The wall is sized by WHICHEVER requirement needs the most metal:
      buckling       buckling_fos >= FOS_BUCKLING (1.5) against SP-8007, at
                     the stability pressure
      yield          far-field von Mises <= F_ty/1.5 = 262 MPa (fos >= 1.0),
                     at the design pressure
      minimum gauge  MIN_GAUGE, if set (ASME BPVC VIII-1, UG-16(b))

    Each iteration solves the governing case at both pressures (the second
    reuses the thermal field of the first). The buckling reserve scales about
    as t^2 and the yield reserve about as t, which gives the next wall, with
    the step limited to a factor of two and replaced by bisection once a
    passing and a failing wall are known. Then the other candidate design
    points are checked for buckling at the final wall, and the wall is
    thickened once more if one of them governs. The result is recorded in
    sized_walls.csv; every iteration is written to
    sizing_fe/<AC>/<AC>_sizing_results.csv.

    A new sizing criterion is added here: compute its reserve in extract(),
    and include its required thickness ratio in g below.
    """
    global P_gauge
    rows = []
    in_dir(os.path.join('sizing_fe', CASE_TAG))
    lc7 = GOVERNING_FACTORS              # governing case (thesis: emergency landing, fill 0.90)
    t_try, sized, govern = max(t_wall, MIN_GAUGE), False, ''
    MIN_T = max(MIN_GAUGE, 0.2e-3)       # numerical floor only, never reached in practice
    # ACCEPTANCE. The governing requirement must be met with at most about
    # 6 % surplus (buckling FOS 1.50 to 1.60, yield fos 1.00 to 1.06). The
    # window cannot be much narrower: the mesh changes with the wall, so a
    # small change in t can move the buckling FOS by several percent (for the
    # A350-900, 2.5 % in t moved it from 1.44 to 1.59). In addition, once a
    # failing and a passing wall are within 3 % of each other, the passing
    # one is taken.
    G_LO = 0.94
    t_fail, t_pass = 0.0, None           # largest failing / smallest passing wall
    print('\nFE WALL SIZING %s: buckling at %.5f bar, yield at %.3f bar, '
          'minimum gauge %.2f mm' % (CASE_TAG, P_STABILITY / 1e5, P_DESIGN / 1e5, MIN_GAUGE * 1e3))
    for it in range(10):
        set_wall(t_try)
        ttag = ('%.3f' % (t_wall * 1e3)).replace('.', 'p')
        P_gauge = P_STABILITY
        tag_b = '%sSZ%d_t%s_pmin' % (TAG_PREFIX, it, ttag)
        rb = run_case(GOVERNING_CASE, lc7, GOVERNING_FILL, seed_size, tag_b)
        P_gauge = P_DESIGN
        ry = run_case(GOVERNING_CASE, lc7, GOVERNING_FILL, seed_size,
                      '%sSZ%d_t%s_pdes' % (TAG_PREFIX, it, ttag),
                      ht_from=os.path.join(os.getcwd(), 'HT_' + tag_b))
        if not rb or not ry:
            raise RuntimeError('sizing case failed at t = %.3f mm' % (t_wall * 1e3))
        rb['stage'] = 'iteration %d, stability pressure' % it
        ry['stage'] = 'iteration %d, design pressure' % it
        rows += [rb, ry]
        write_rows('%s_sizing_results.csv' % CASE_TAG, rows, lead=('stage',))
        for r_ in (rb, ry):
            if not r_['orientation_ok']:
                raise RuntimeError('S11 is not the axial stress; sizing stopped.')
            if abs(r_['balance_ratio'] - 1.0) > 0.02:
                raise RuntimeError('load balance off by more than 2 %; sizing stopped.')
        fb, fy = rb['buckling_fos'], ry['fos']
        need_b = 0.0 if fb == float('inf') else FOS_BUCKLING / fb     # <= 1 means met
        need_y = 1.0 / fy
        g = max(need_b, need_y)
        govern = 'buckling' if need_b >= need_y else 'yield'
        print('  iteration %d: t = %.3f mm | buckling FOS %.3f | yield fos %.3f | %s governs (%.3f)'
              % (it, t_wall * 1e3, fb, fy, govern, g))
        if g > 1.0:
            t_fail = max(t_fail, t_wall)
        elif t_pass is None or t_wall < t_pass[0]:
            t_pass = (t_wall, fb, fy, govern)
        if G_LO <= g <= 1.0:
            sized = True
            break
        if t_pass is not None and t_fail > 0 and t_pass[0] / t_fail < 1.03:
            set_wall(t_pass[0])
            fb, fy, govern = t_pass[1], t_pass[2], t_pass[3]
            print('  bracket %.3f (fails) / %.3f mm (passes) within 3 %%: taking %.3f mm'
                  % (t_fail * 1e3, t_pass[0] * 1e3, t_pass[0] * 1e3))
            sized = True
            break
        if g < 0.98 and MIN_GAUGE > 0 and t_wall <= MIN_GAUGE * (1 + 1e-9):
            sized, govern = True, 'minimum gauge'
            break
        t_b = t_wall * math.sqrt(need_b * 1.01) if need_b > 0 else 0.0
        t_y = t_wall * need_y * 1.01
        t_new = max(t_b, t_y, MIN_T)
        t_try = min(max(t_new, 0.5 * t_wall), 2.0 * t_wall)   # at most halve or double
        if t_pass is not None and t_fail > 0:
            # a pass and a fail are known: bisect between them instead
            t_try = math.sqrt(t_fail * t_pass[0])
    if not sized:
        raise RuntimeError('sizing did not converge in 10 iterations (last t = %.3f mm)'
                           % (t_wall * 1e3))
    # Buckling verification of the other candidate design points at the final
    # wall: the governing case at the other fill levels, and the other load
    # cases at the governing fill level.
    checks = [(GOVERNING_CASE, f_) for f_ in fill_ratios if abs(f_ - GOVERNING_FILL) > 1e-9] + \
             [(lc_, GOVERNING_FILL) for lc_ in reversed(list(ALL_LOAD_CASES)) if lc_ != GOVERNING_CASE]
    for lc_name, f_ in checks:
        P_gauge = P_STABILITY
        r = run_case(lc_name, ALL_LOAD_CASES[lc_name], f_, seed_size,
                     '%sSZV_%s_f%02d_pmin' % (TAG_PREFIX, lc_name, int(round(f_ * 100))))
        if r:
            r['stage'] = 'verify, stability pressure'
            rows.append(r)
    write_rows('%s_sizing_results.csv' % CASE_TAG, rows, lead=('stage',))
    stab = [x for x in rows if x['p_gauge_bar'] < 1.0 and abs(x['t_wall_mm'] - t_wall * 1e3) < 1e-6]
    worst = min(stab, key=lambda x: x['buckling_fos'])
    fb_min = worst['buckling_fos']
    if fb_min < FOS_BUCKLING:
        print('  *** %s fill %.2f governs buckling (FOS %.3f); thickening'
              % (worst['load_case'], worst['fill'], fb_min))
        set_wall(t_wall * math.sqrt(FOS_BUCKLING * 1.01 / fb_min))
        P_gauge = P_STABILITY
        r = run_case(worst['load_case'], load_cases[worst['load_case']], worst['fill'], seed_size,
                     '%sSZF_t%s' % (TAG_PREFIX, ('%.3f' % (t_wall * 1e3)).replace('.', 'p')))
        r['stage'] = 'final, governing case'
        rows.append(r)
        fb_min, govern = r['buckling_fos'], 'buckling (%s f%.2f)' % (worst['load_case'], worst['fill'])
        write_rows('%s_sizing_results.csv' % CASE_TAG, rows, lead=('stage',))
    P_gauge = P_DESIGN
    record_sized_wall(fb_min, fy, govern)
    os.chdir(_ROOT_DIR)
    print('\n' + '=' * 78)
    print('FE SIZING RESULT %s: t = %.3f mm | governed by %s | lowest buckling FOS %.3f | '
          'yield fos %.3f' % (CASE_TAG, t_wall * 1e3, govern, fb_min, fy))
    print('=' * 78)
    return rows


def stage_matrix(pressure, label, tag_fmt, out_csv, ht_fmt=None, resume_glob=None):
    """
    Every load case at every fill level at `pressure`. RESUMES: design points
    already in out_csv (or in any file matching resume_glob) at this wall and
    pressure are kept, not re-run. ht_fmt, if given, names the heat transfer
    job to reuse for each design point (the stability matrix reuses the
    production thermal fields, which are identical).
    """
    global P_gauge
    P_gauge = pressure
    rows = done_rows(out_csv, pressure, resume_glob)
    uniq = {}
    for r in rows:
        uniq[(r['load_case'], round(r['fill'], 2))] = r
    rows = list(uniq.values())
    have = set(uniq.keys())
    print('\n%s RUN %s: %.5f bar, seed %.1f mm, load cases %s x %d fills (%d already done)'
          % (label, CASE_TAG, P_gauge / 1e5, seed_size * 1e3, ', '.join(load_cases.keys()),
             len(fill_ratios), len(have)))
    for lc_name, factors in load_cases.items():
        for phi_f in fill_ratios:
            if (lc_name, round(phi_f, 2)) in have:
                continue
            nn = int(round(phi_f * 100))
            tag = tag_fmt % (TAG_PREFIX, lc_name, nn)
            src = (ht_fmt % (TAG_PREFIX, lc_name, nn)) if ht_fmt else None
            print('  %s' % tag)
            r = run_case(lc_name, factors, phi_f, seed_size, tag, ht_from=src)
            if r:
                rows.append(r)
                write_rows(out_csv, rows)          # saved after every case
    P_gauge = P_DESIGN
    return rows


def stage_convergence(pressure, label, out_csv, ht_label=None):
    """
    PROCESS 7. The governing case over the mesh densities
    convergence_fractions (seed = fraction x sqrt(R t)), run at both
    pressures: the design pressure for the far-field von Mises stress and
    the stability pressure for the axial membrane stress of the buckling
    criterion. Resumes per level; the stability-pressure sweep reuses the
    design-pressure thermal field of the same level.
    """
    global P_gauge
    P_gauge = pressure
    rows = done_rows(out_csv, pressure)
    have = set(round(r['seed_mm'], 3) for r in rows)
    print('\nMESH CONVERGENCE %s (%s, %.5f bar): %s fill %.2f, levels %s (%d already done)'
          % (CASE_TAG, label, P_gauge / 1e5, GOVERNING_CASE, GOVERNING_FILL,
             convergence_fractions, len(have)))
    for frac in convergence_fractions:
        sd = boundary_layer * frac
        if round(sd * 1e3, 3) in have:
            continue
        ftag = str(frac).replace('.', 'p')
        tag = '%sconv_%s_%s' % (TAG_PREFIX, label, ftag)
        src = os.path.join(os.getcwd(), 'HT_%sconv_%s_%s' % (TAG_PREFIX, ht_label, ftag)) \
            if ht_label else None
        print('  seed %.1f mm (%.2f elements per boundary layer)' % (sd * 1e3, 1.0 / frac))
        r = run_case(GOVERNING_CASE, GOVERNING_FACTORS, GOVERNING_FILL, sd, tag, ht_from=src)
        if r:
            rows.append(r)
            write_rows(out_csv, rows)
    P_gauge = P_DESIGN
    return rows


# ============================================================================
# DRY RUN: resolved inputs and first-order checks, no Abaqus model
# ============================================================================
# Closed-form values printed for a quick check of a new study: membrane hoop
# stress pR/t, SP-8007 allowable, dry-wall temperature asymptote q/h_dry and
# conduction length sqrt(k t / h_dry), shell edge decay length pi/beta with
# beta = [3(1-nu^2)/(R t)^2]^(1/4), mesh size and number of analyses.
if DRY_RUN:
    _beta = (3.0 * (1.0 - nu * nu) / (R_mid * R_mid * t_wall * t_wall)) ** 0.25
    _gam_phi = math.sqrt(R_mid / t_wall) / 16.0
    _s_cl = E_mod * t_wall / (R_mid * math.sqrt(3.0 * (1.0 - nu * nu)))
    _s_cr = (1.0 - 0.901 * (1.0 - math.exp(-_gam_phi))) * _s_cl
    _n_el = A_shell / (seed_size * seed_size)
    print('\nDRY RUN: inputs resolved, no model built')
    print('-' * 78)
    print('Design basis          : %s' % (DB.PATH if DB is not None else 'thesis values in this file'))
    print('Material              : E %.1f GPa, nu %.2f, alpha %.2e 1/K, k %.0f W/m-K, Fty %.0f MPa'
          % (E_mod / 1e9, nu, alpha, k_shell, Fty_RT / 1e6))
    print('Margins               : yield %.2f -> allowable %.1f MPa; buckling %.2f'
          % (SF, sigma_allow / 1e6, FOS_BUCKLING))
    print('Pressures             : design %.4f bar, stability %.5f bar (gauge)' % (P_gauge / 1e5, P_STABILITY / 1e5))
    print('Thermal               : T_liquid %.1f K, h_wetted %.0f, h_dry %.2f W/m2-K, q_MLI %.3f W/m2'
          % (T_LH2, h_wetted, h_dry, q_outer))
    print('Load cases            : ' + '; '.join('%s (%.1f, %.1f, %.1f)' % (k, v['nx'], v['ny'], v['nz'])
                                                for k, v in load_cases.items()))
    print('Fill levels           : %s' % fill_ratios)
    print('Governing case        : %s at fill %.2f' % (GOVERNING_CASE, GOVERNING_FILL))
    print('-' * 78)
    print('Hoop stress pR/t      : %.1f MPa at the design pressure (von Mises 0.866 pR/t = %.1f MPa)'
          % (P_gauge * R_mid / t_wall / 1e6, 0.866 * P_gauge * R_mid / t_wall / 1e6))
    print('SP-8007 sigma_cr      : %.2f MPa, allowable %.2f MPa (R/t = %.0f)'
          % (_s_cr / 1e6, _s_cr / FOS_BUCKLING / 1e6, R_mid / t_wall))
    print('Dry-wall asymptote    : q/h_dry = %.3f K; conduction length %.3f m'
          % (q_outer / h_dry, math.sqrt(k_shell * t_wall / h_dry)))
    print('Edge decay length     : pi/beta = %.3f m; far-field exclusion %.2f m = %.1f decay lengths'
          % (math.pi / _beta, FAR_FIELD_EXCLUSION, FAR_FIELD_EXCLUSION / (math.pi / _beta)))
    print('Mesh                  : seed %.1f mm, of order %.0f elements (area / seed squared)' % (seed_size * 1e3, _n_el))
    print('FINAL run             : %d production + %d stability + %d convergence analyses'
          % (len(load_cases) * len(fill_ratios), len(load_cases) * len(fill_ratios), 2 * len(convergence_fractions)))
    print('-' * 78)
    raise SystemExit(0)

results = []

# Which FINAL stages to run, for example in a time-limited session:
#     set E190_STAGES=production            production only, then stop
#     set E190_STAGES=stability             stability only
# Default: all four, in order.
STAGES = [x.strip().lower() for x in os.environ.get(
    'E190_STAGES', 'sizing,production,stability,convergence').split(',') if x.strip()]
# Guard: a FINAL session that does not size must use the recorded wall. If
# sized_walls.csv is missing or has no row for this aircraft, stop rather than
# run the matrix at the starting thickness of AIRCRAFT_DATA.
if RUN_MODE == 'FINAL' and 'sizing' not in STAGES:
    _have_row = os.path.exists(SIZED_WALLS) and any(
        r.get('aircraft') == CASE_TAG for r in csv.DictReader(open(SIZED_WALLS)))
    if not _have_row:
        raise RuntimeError('No sized wall for %s in %s: refusing to run the '
                           'final stages at the fallback thickness.'
                           % (CASE_TAG, SIZED_WALLS))
print('FINAL STAGES: %s' % ', '.join(STAGES))

if RUN_MODE in ('SIZING', 'FINAL') and (RUN_MODE == 'SIZING' or 'sizing' in STAGES):
    if os.path.exists(SIZED_WALLS) and RUN_MODE == 'FINAL' and \
            os.path.exists(os.path.join(_ROOT_DIR, 'sizing_fe', CASE_TAG,
                                        '%s_sizing_results.csv' % CASE_TAG)) and \
            any(r.get('aircraft') == CASE_TAG for r in csv.DictReader(open(SIZED_WALLS))):
        print('\nSIZING already done for %s (sized_walls.csv): %.3f mm, skipped'
              % (CASE_TAG, t_wall * 1e3))
    else:
        results += stage_sizing()

if RUN_MODE == 'PRODUCTION' or (RUN_MODE == 'FINAL' and 'production' in STAGES):
    set_wall(t_wall)
    _sfx = ('_' + LC_FILTER) if LC_FILTER else ''
    N_EXPECT = len(load_cases) * len(fill_ratios)
    out = os.path.join(_ROOT_DIR, '%s_production%s_results.csv' % (CASE_TAG, _sfx)) \
        if RUN_MODE == 'FINAL' else None
    if RUN_MODE == 'FINAL' and stage_done(out, N_EXPECT):
        print('\nPRODUCTION already complete for %s at %.3f mm, skipped' % (CASE_TAG, t_wall * 1e3))
    else:
        if RUN_MODE == 'FINAL':
            os.chdir(_ROOT_DIR)
            results += stage_matrix(P_DESIGN, 'PRODUCTION', '%s%s_f%02d', out,
                                    resume_glob=os.path.join(_ROOT_DIR, '%s_production*_results.csv' % CASE_TAG))
        else:
            print('\nPRODUCTION RUN')
            print('Seed %.1f mm, %d load cases x %d fill ratios\n'
                  % (seed_size * 1000.0, len(load_cases), len(fill_ratios)))
            for lc_name, factors in load_cases.items():
                for phi_f in fill_ratios:
                    tag = '%s%s_f%02d' % (TAG_PREFIX, lc_name, int(phi_f * 100))
                    print('  %s' % tag)
                    r = run_case(lc_name, factors, phi_f, seed_size, tag)
                    if r:
                        results.append(r)

if RUN_MODE == 'STABILITY' or (RUN_MODE == 'FINAL' and 'stability' in STAGES):
    set_wall(t_wall)
    _sfx = ('_' + LC_FILTER) if LC_FILTER else ''
    out = os.path.join(_ROOT_DIR, '%s_stability%s_results.csv' % (CASE_TAG, _sfx))
    if stage_done(out, len(load_cases) * len(fill_ratios)):
        print('\nSTABILITY already complete for %s, skipped' % CASE_TAG)
    else:
        in_dir(os.path.join('stability', CASE_TAG))
        stage_matrix(P_STABILITY, 'STABILITY', '%sSTAB_%s_f%02d', out,
                     ht_fmt=os.path.join(_ROOT_DIR, 'HT_%s%s_f%02d'),
                     resume_glob=os.path.join(_ROOT_DIR, '%s_stability*_results.csv' % CASE_TAG))
        os.chdir(_ROOT_DIR)

if RUN_MODE == 'CONVERGENCE' or (RUN_MODE == 'FINAL' and 'convergence' in STAGES
                                  and GOVERNING_CASE in load_cases):
    # (in a FINAL run split by E190_ONLY_LC, only the session that holds
    # the governing case runs it)
    set_wall(t_wall)
    for pressure, label in ((P_DESIGN, 'pdesign'), (P_STABILITY, 'pmin')):
        out = os.path.join(_ROOT_DIR, '%s_convergence_%s_results.csv' % (CASE_TAG, label))
        if stage_done(out, len(convergence_fractions)):
            print('\nCONVERGENCE (%s) already complete for %s, skipped' % (label, CASE_TAG))
            continue
        in_dir(os.path.join('convergence_final', CASE_TAG))
        stage_convergence(pressure, label, out,
                          ht_label='pdesign' if label == 'pmin' else None)
        os.chdir(_ROOT_DIR)

if RUN_MODE == 'FINAL':
    print('\n' + '=' * 78)
    print('FINAL SESSION COMPLETE FOR %s at t = %.3f mm (stages: %s)'
          % (CASE_TAG, t_wall * 1e3, ', '.join(STAGES)))
    print('Next: copy the *_results.csv files into results/ and run '
          'python code/3_results_processing/final_numbers.py')
    print('=' * 78)


# ============================================================================
# 6. OUTPUT
# ============================================================================

out_csv = '%s_%s%s%s_results.csv' % (
    CASE_TAG, RUN_MODE.lower(),
    ('_' + LC_FILTER) if LC_FILTER else '',
    ('_' + FILL_FILTER) if FILL_FILTER else '')
if results and RUN_MODE in ('PRODUCTION', 'SIZING'):
    write_rows(out_csv, results, lead=('stage',) if RUN_MODE == 'SIZING' else ())
if results:
    worst = min(results, key=lambda r: r['fos'])
    print('\n' + '=' * 78)
    print('SUMMARY')
    print('=' * 78)
    print('Governing point : %s at fill %.2f' % (worst['load_case'], worst['fill']))
    print('Peak von Mises  : %.2f MPa' % worst['mises_MPa'])
    print('Allowable       : %.1f MPa' % (sigma_allow / 1e6))
    print('Factor of safety: %.2f  -> %s'
          % (worst['fos'], 'PASS' if worst['fos'] >= 1.0 else 'FAIL'))
    print('Largest wetted-to-dry dT : %.3f K'
          % max([r['dT'] for r in results]))
    rmax = max([r['sumRF_N'] for r in results])
    print('Largest mount reaction : %.1f N  (%.1f kN)' % (rmax, rmax / 1000.0))
    if BENCHMARK_MODE:
        print('\nBENCHMARK: Bagarello et al. report T_max = 25.3 K at fill 0.75.')
        for r in results:
            if abs(r['fill'] - 0.75) < 1e-6:
                print('           This model gives T_max = %.2f K  (error %.1f%%)'
                      % (r['T_max'], 100.0 * (r['T_max'] - 25.3) / 25.3))
    print('=' * 78)
elif RUN_MODE not in ('FINAL', 'STABILITY', 'CONVERGENCE'):
    print('\nNo results produced. Check the job logs above.')
