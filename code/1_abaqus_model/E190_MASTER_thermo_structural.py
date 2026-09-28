"""
================================================================================
E190 LH2 INNER VESSEL - SEQUENTIALLY COUPLED THERMO-STRUCTURAL ANALYSIS
E190_MASTER_thermo_structural.py   (version 3)
Abaqus/CAE Python script - TU Delft MSc thesis

Author: Sarvesh Mahadevan
Date: September 2026

Version 2 (yield + buckling, ring/saddle mounts, uniform 20 K field) is
preserved at _originals/E190_MASTER_thermo_structural_v2_buckling_mounts.py.
See CHANGELOG.md for why it was superseded.

--------------------------------------------------------------------------------
THE ADOPTED METHODOLOGY, AND WHERE EACH STEP COMES FROM
--------------------------------------------------------------------------------
The modelling PROCESS used here is not invented for this thesis. It is the
published, peer-reviewed approach of

    Bagarello, S., Campagna, D. and Benedetti, I., "Computational
    thermo-mechanical modelling and design-space exploration of cryogenic
    hydrogen tanks for aviation", Aerospace Science and Technology 168 (2026)
    110755.

supplemented on the inner-vessel treatment by

    Biancotto, A., "Design and analysis of LH2 tank structures for aircraft
    retrofit applications", MSc thesis, Delft University of Technology, 2024.

Adopting a verified process, rather than assembling one, is what allows the
results of this thesis to be trusted. Each step below is tagged in the code
with "PROCESS n" at the point where it is implemented.

  1  Prescribe only the two SINK temperatures: the compartment outside the
     tank and the fuel inside it. Every wall temperature is unknown and is
     solved for.                                        [Bagarello Sec. 3.1]

  2  Declare the analysis regime explicitly: steady-state, linear
     thermo-elastic.                                      [Bagarello Sec. 3]

  3  Treat the fill ratio as a free parameter that partitions the inner
     surface into a wetted region and a dry region.     [Bagarello Sec. 3.1]

  4  Apply two film coefficients that differ by orders of magnitude. The
     wetted wall is effectively pinned to the fuel temperature; the dry wall
     rejects heat into hydrogen vapour by natural convection and therefore
     runs warmer. THIS IS THE STEP THAT MAKES FILL RATIO A THERMAL VARIABLE,
     and it is the mechanism by which the thermal and structural sides of
     this thesis are coupled.                           [Bagarello Sec. 3.1]

  5  Solve the thermal problem first, then carry the resulting nodal
     temperature field into the structural problem through
     sigma = C : (eps - alpha*dT).                       [Bagarello Eq. (6)]

  6  Stage the mechanical boundary conditions. The baseline removes rigid
     body motion and nothing else; a realistic anchoring layout is a
     SEPARATE, LATER study. Bagarello et al. report that adding realistic
     integration constraints produces a fourfold increase in liner stress,
     so this omission is quantified rather than unstated.
                                                       [Bagarello Sec. 4.4.1]

  7  Converge the mesh against a physical output plotted versus degrees of
     freedom, rather than against element count alone.    [Bagarello Fig. 4]

  8  Match the failure criterion to the material class and apply the
     regulatory safety factor. For an isotropic metallic vessel this is von
     Mises with SF = 1.5 per CS-25.303.                  [Bagarello Sec. 4.4]

  9  Build realism progressively: ideal baseline first, then constraints,
     then geometric features, each increment isolated so its effect can be
     attributed.                                           [Bagarello Sec. 4]

 10  Treat the inner vessel as the design subject in its own right, with
     delta-T and delta-p as the primary load case.       [Biancotto Sec. 5.8.2]

 11  Where supports are eventually introduced, characterise them by axial
     and flexural stiffness rather than as rigid restraints.
                                                        [Biancotto Sec. 5.8.1]

--------------------------------------------------------------------------------
WHAT THIS THESIS DOES THAT THE SOURCE WORK DOES NOT
--------------------------------------------------------------------------------
The process above is adopted. The following extensions are the contribution
of this thesis, and neither source paper performs them:

  A  A FULL CS-25 LOAD MATRIX COUPLED TO FILL RATIO. Bagarello et al. apply
     thermal load and internal pressure only. They do not apply the CS-25.337
     manoeuvre factor or the CS-25.561 emergency landing factors, and they do
     not couple those factors to the fill ratio through the hydrostatic head.
     Biancotto applies CS-25.561 ultimate forces, but to the support
     structure rather than as a stress field on the vessel wall.

  B  CROSS-AIRCRAFT CONSISTENCY. Both sources analyse a single tank. Here one
     sizing rule and one methodology are applied across a regional, a
     narrow-body and a wide-body aircraft, which is what permits a statement
     about how the design scales with mission.

--------------------------------------------------------------------------------
WHAT CHANGED FROM VERSION 2
--------------------------------------------------------------------------------
Four deliberate changes:

  1. INNER VESSEL ONLY (PROCESS 10). The outer jacket, MLI, struts and
     airframe mounts are not modelled. The MLI appears only as a heat flux
     boundary condition on the outer face of the vessel.

  2. A REAL LOAD PATH (PROCESS 6). The vessel is restrained on THREE
     CIRCUMFERENTIAL RINGS, at the two cap-to-barrel junctions and at
     mid-barrel, after Bagarello et al. Sec. 4.4.1, who take the arrangement
     from H2FLY testbed hardware.

     Two earlier attempts were rejected on physical grounds:

       - A 3-2-1 restraint carries zero load only for SELF-EQUILIBRATED
         loads. The CS-25 inertia cases are not: under LC7 the hydrostatic
         pressure field alone has a resultant of 562 kN, which was driven
         through three nodes and returned 4943 MPa.

       - Inertia relief removed that artefact but introduced another. The
         balancing acceleration is distributed in proportion to mass, and
         the liquid has no mass in the model, so the entire liquid reaction
         was smeared over the 1648 kg shell as a fictitious body force of
         up to 34.8 g. It also cancels the shell's own body force exactly,
         so the applied inertia did nothing at all.

     The deeper point is that a CS-25 load factor PRESUPPOSES a load path.
     A freely floating tank does not experience a crash load; it simply
     accelerates. The stress arises from being held. Modelling the restraint
     is therefore not an optional refinement, it is what makes the load case
     meaningful. The reaction at the rings is reported as a result, because
     it is the force the support structure must carry.

     Bagarello et al. report that introducing this anchoring raises liner
     stress roughly fourfold relative to rigid-body-removal alone.

  3. NO BUCKLING (PROCESS 8). The sizing criterion is yield only, assessed
     through von Mises with a safety factor of 1.5 per CS-25.303. Stability assessment
     requires imperfection sensitivity and post-buckling behaviour, which are
     outside the linear scope of this work.

  4. THE THERMAL LOAD IS NOW REAL (PROCESS 3 and 4). The previous model
     applied a uniform 20 K field to a body that was free to contract, which produces exactly zero
     thermal stress. Here the wetted and dry regions of the wall are treated
     separately, exactly as Bagarello et al. do: the wall below the liquid
     surface is pinned at the fuel temperature, while the wall above it
     exchanges heat with the ullage gas by natural convection at a much lower
     film coefficient, and therefore runs warmer. THE FILL RATIO GENERATES THE
     TEMPERATURE GRADIENT, so the thermal and structural sides are coupled
     through the variable this thesis is about.

--------------------------------------------------------------------------------
ANALYSIS CHAIN
--------------------------------------------------------------------------------
  STEP 1  Steady-state heat transfer    -> nodal temperature field T(x)
  STEP 2  Static structural analysis, reading T(x) as a predefined field,
          plus internal pressure, hydrostatic head and CS-25 inertia

This is a SEQUENTIALLY coupled analysis (PROCESS 5). It is valid because the
mechanical response does not feed back into the thermal problem at these strain levels.

--------------------------------------------------------------------------------
FORMULAS USED (tagged "FORMULA n" at each point of use below)
--------------------------------------------------------------------------------
  1  Tank volume            V = pi*R^2*L_barrel + (4/3)*pi*R^2*h_cap
  2  Cap height             h_cap = R / (a/b)
  3  Fill surface height    theta - sin(theta) = 2*pi*phi_f   (circular segment)
  4  Convective flux        q = h * (T_wall - T_fluid)
  5  MLI heat flux          q_MLI = Q_budget / A_shell
  6  Hydrostatic pressure   p(z) = p_gauge + rho_LH2 * (n*g) * depth
  7  Inertial body force    f = n * g   (applied to the shell's own mass)
  8  Thermal strain         eps_th = alpha * (T - T_ref)
  9  von Mises stress       sig_eq = sqrt(0.5*[(s1-s2)^2+(s2-s3)^2+(s3-s1)^2])
 10  Allowable stress       sig_allow = sigma_yield / SF

--------------------------------------------------------------------------------
DATA SOURCES
--------------------------------------------------------------------------------
  Al 2219-T87 properties      : thesis material table (ref7), cross-checked
                                against Bagarello Table 1
  T_LH2 = 20.3 K              : Bagarello Sec. 3.1
  h_int = 1.46 W/m2K          : Bagarello Sec. 3.1, their ref [45], for GH2
                                natural convection against the dry wall
  SF = 1.5                    : CS-25.303 (ref22); Bagarello cite FAR 25.303
  CS-25 load factors          : CS-25.337 and CS-25.561 (ref22)
  Q budget = 702.6 W          : derived in boiloff_requirement_derivation.py
                                from a 14.68 h no-vent hold time
  Geometry                    : thesis sizing rule, R = 0.8648*D_fuselage/2

--------------------------------------------------------------------------------
BENCHMARK TARGET
--------------------------------------------------------------------------------
Bagarello et al. report, for their reference tank at a filling ratio of 0.75,
that the wetted liner sits at 20 K while the dry region peaks at 25.3 K. Set
BENCHMARK_MODE = True to run their geometry and operating point instead of the
E190, and compare the predicted peak wall temperature against that 25.3 K.
================================================================================
"""

from abaqus import *
from abaqusConstants import *

# THESE STAR IMPORTS ARE REQUIRED FOR HEADLESS RUNS. DO NOT REMOVE.
#
# Abaqus attaches most Model methods lazily, when the module that defines them
# is imported. Abaqus/CAE with its GUI loads every module at startup, so a
# script launched through File -> Run Script finds them already present. A
# headless run started with
#
#     abaqus cae noGUI=E190_MASTER_thermo_structural.py
#
# loads only the kernel, and any method whose module was never imported is
# simply absent. The failure looks like a typo rather than a missing import:
#
#     AttributeError: 'Model' object has no attribute 'FieldOutputRequest'
#
# even though the call is correct and runs perfectly inside the GUI. Importing
# the standard set below - the same block Abaqus itself writes at the top of
# every .rpy replay file - makes the script behave identically either way.
#
# These come BEFORE the aliased imports underneath, so that abqmesh and
# abqjob remain bound to the modules and are not shadowed.
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

# RESTORE PYTHON'S BUILT-INS (27 Sep 2026). The star imports above bring in
# the Abaqus XY-data operators, which replace sum() (and may replace other
# built-ins) with versions that reject a generator: "TypeError: arg1; found
# 'generator', expecting a recognized type". This script never uses the XY
# operators, so the Python built-ins are put back for everything below.
import builtins as _builtins
sum = _builtins.sum
any = _builtins.any
all = _builtins.all
max = _builtins.max
min = _builtins.min
abs = _builtins.abs
round = _builtins.round

import regionToolset
import mesh as abqmesh
import job as abqjob
import math
import os
import random

# ---------------------------------------------------------------------------
# WORKING DIRECTORY
# ---------------------------------------------------------------------------
# This script lives in files\python_scripts\ but Abaqus should write its .odb,
# .dat, .msg and the results .csv into the parent files\ folder, so the outputs
# sit alongside the LaTeX and the changelog rather than among the source.
#
# Pinning it here means the script can be launched from anywhere - directly via
# File -> Run Script, or through E190_CHECK_smoketest.py - and the outputs
# always land in the same place.
# Abaqus does NOT define __file__ when a script is launched through
# File -> Run Script, so the automatic path resolution silently fails and the
# outputs land in whatever directory CAE happens to be using (C:\temp by
# default). An explicit fallback is therefore required.
#
# PORTABILITY. The Windows path below is only a last-resort fallback for
# File -> Run Script on Sarvesh's laptop. It must never be the value used on
# another machine, and in particular not on DelftBlue, where it does not
# exist. Resolution order:
#
#   1. E190_WORKDIR environment variable, if set. This is what the SLURM
#      submission script on DelftBlue uses:
#          export E190_WORKDIR=/scratch/<netid>/e190
#   2. The parent of the directory holding this file, when __file__ exists
#      (true for `abaqus cae noGUI=...`, false for File -> Run Script).
#   3. The hard-coded Windows path.
#   4. The current directory, if none of the above is a real directory. On a
#      cluster this is the SLURM submission directory, which is correct.
_THESIS_DIR = r'C:\Users\hrish\Thesis_Code\files'
_SCRIPT_DIR = None
try:
    _SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    _THESIS_DIR = os.path.dirname(_SCRIPT_DIR)
    # Repository layout (code/1_abaqus_model/): the analysis runs in the
    # current working folder, which receives the .odb and results files.
    if os.path.basename(_SCRIPT_DIR) == '1_abaqus_model':
        _THESIS_DIR = os.getcwd()
except NameError:
    pass
_ROOT_DIR = _THESIS_DIR          # the thesis folder, whatever E190_WORKDIR says
import csv
_ENV_DIR = os.environ.get('E190_WORKDIR', '').strip()
if _ENV_DIR:
    if not os.path.isdir(_ENV_DIR):
        raise ValueError('E190_WORKDIR is set to a directory that does not '
                         'exist: %s' % _ENV_DIR)
    _THESIS_DIR = _ENV_DIR
# DelftBlue (27 Sep 2026): Abaqus 2024 on Linux does not define __file__ in
# noGUI mode, so _ROOT_DIR fell back to the Windows path, which does not exist
# there. sized_walls.csv was then not found and the wall would silently have
# reverted to the AIRCRAFT_DATA value. Use the working directory instead.
if not os.path.isdir(_ROOT_DIR):
    _ROOT_DIR = _THESIS_DIR if os.path.isdir(_THESIS_DIR) else os.getcwd()
    print('Thesis folder not found at the default path; using %s' % _ROOT_DIR)
if os.path.isdir(_THESIS_DIR):
    os.chdir(_THESIS_DIR)
print('Working directory: %s' % os.getcwd())

from odbAccess import openOdb


# ============================================================================
# 1. USER INPUTS
# ============================================================================

# --- Run mode -------------------------------------------------------------
# 'PRODUCTION'  : single run at the converged mesh, full load-case x fill matrix
# 'CONVERGENCE' : mesh sensitivity study, one representative load case
RUN_MODE = 'PRODUCTION'

# Set True to reproduce Bagarello's reference tank instead of the E190.
BENCHMARK_MODE = False

# Both overridable from the environment, so a convergence study or a benchmark
# run can be launched without editing this file:
#
#     $env:E190_RUN_MODE = "CONVERGENCE"
#     $env:E190_BENCHMARK = "1"
#
# Combine with E190_WORKDIR to keep each run's output in its own folder.
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
# One physics file serves all three tanks. Only the geometry, wall thickness
# and heat-leak budget differ; every boundary condition, load case, failure
# criterion and post-processing step is shared.
#
# This matters more than it looks. Bringing the E190 model up took eight rounds
# of debugging - a transposed axis, a missing initial temperature, a body force
# in the wrong units, a partial constraint ring, and a pressure applied to the
# wrong face. Keeping three separate copies would mean applying every future
# fix three times and hoping none was missed. There is one file.
#
# Set AIRCRAFT here, or use the thin runner scripts that override it:
#     run_A320.py, run_A350.py
AIRCRAFT = 'E190'            # 'E190' | 'A320' | 'A350'

# Overridable from the environment, so the three aircraft can be launched from
# the command line without editing this file:
#
#     set E190_AIRCRAFT=A320
#     abaqus cae noGUI=python_scripts\E190_MASTER_thermo_structural.py
#
# CASE_TAG follows AIRCRAFT, so each aircraft writes its own results CSV.
#
# WARNING - DO NOT RUN TWO AIRCRAFT AT THE SAME TIME YET. The per-case job
# name is built from the load case and fill ratio only, not the aircraft, so
# an A320 run would overwrite HT_Baseline_f10.odb from the E190 run while it
# was still being read. The CSV files are safe; the .odb files are not. Fix
# the job naming before running aircraft concurrently.
AIRCRAFT = os.environ.get('E190_AIRCRAFT', AIRCRAFT).strip() or AIRCRAFT

# R_mid      shell mid-surface radius                                     [m]
# L_barrel   cylindrical barrel length                                    [m]
# cap_aspect a/b for the semi-ellipsoidal end caps                        [-]
# t_wall     wall thickness                                               [m]
# Q_budget   allowable heat leak, from the ~14.7 h no-vent hold           [W]
#            (boiloff_requirement_derivation.py)
# P_gauge    internal gauge pressure                                      [Pa]
# WALL THICKNESSES - CORRECTED 27 Sep 2026.
# The earlier values (E190 4.785, A320 6.382, A350 9.948 mm) came from Abaqus
# sizing runs that put the design pressure on SIDE1 of their shell surfaces,
# the OUTER face on those meshes, so the tank was squeezed rather than
# pressurised and the buckling screen sized against a hoop compression of
# -pR/t that does not exist in service. See python_scripts/sizing_closed_form.py.
# The values below are the corrected closed-form sizing (buckling at the
# minimum pressure governs, FOS 1.50). RUN_MODE='SIZING' verifies and, if
# needed, adjusts them with THIS model; update them here with its result.
AIRCRAFT_DATA = {
    'E190': dict(R_mid=1.3015, L_barrel=13.9675, cap_aspect=1.6,
                 t_wall=0.002747, Q_budget=702.6,  P_gauge=150000.0),
    'A320': dict(R_mid=1.7080, L_barrel=10.0120, cap_aspect=1.6,
                 t_wall=0.002694, Q_budget=919.3,  P_gauge=150000.0),
    'A350': dict(R_mid=2.5771, L_barrel=33.8240, cap_aspect=1.6,
                 t_wall=0.008884, Q_budget=6583.2, P_gauge=150000.0),
}

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
    # WALL FROM THE FE SIZING (added 27 Sep 2026). RUN_MODE SIZING / FINAL
    # records the sized wall of each aircraft in <thesis>/sized_walls.csv.
    # Every later run of this script, and the capture, merge and plotting
    # scripts, take the wall from there, so the sized value flows through
    # without anybody retyping it.
    SIZED_WALLS = os.path.join(_ROOT_DIR, 'sized_walls.csv')
    # Repository layout: fall back to results/sized_walls.csv of the repository.
    if not os.path.exists(SIZED_WALLS) and _SCRIPT_DIR:
        _repo_sw = os.path.join(os.path.dirname(os.path.dirname(_SCRIPT_DIR)), 'results', 'sized_walls.csv')
        if os.path.exists(_repo_sw):
            SIZED_WALLS = _repo_sw
    if os.path.exists(SIZED_WALLS):
        for _row in csv.DictReader(open(SIZED_WALLS)):
            if _row.get('aircraft') == AIRCRAFT:
                t_wall = float(_row['t_wall_mm']) / 1000.0
                print('WALL FROM FE SIZING (sized_walls.csv): %.3f mm' % (t_wall * 1e3))
    # Overrides from the environment, so a trial wall or a different pressure
    # can be run without editing this file:
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
    # BAGARELLO REFERENCE TANK - their Table 4, read from the paper.
    #
    #     storage volume        V        100 m3
    #     cylinder length       L        4.24 m
    #     radius                R        2.12 m       spherical caps, a/b = 1
    #     inner layer           t_Al     10 mm
    #     middle layer          t_PUR    0.5 m        PUR-96 foam
    #     outer layer           t_CFRP   10 mm
    #     filling ratio         phi_f    0.75
    #
    # CORRECTED 19 September. t_wall was previously set to 0.001 m with the
    # comment "their t_Al = 1 mm". Their Table 4 states 10 mm. The tenfold
    # error put the liner hoop stress at pR/t = 360 MPa instead of 36 MPa and
    # made the benchmark's structural output meaningless.
    #
    # THEY MODEL THREE SOLID LAYERS; THIS MODEL HAS ONE SHELL. Only the
    # aluminium liner is represented here. The foam and composite enter
    # through the outer heat flux, computed from their stated boundary
    # conditions below - not fitted to their answer.
    # ----------------------------------------------------------------
    R_mid      = 2.12        # m,  their Table 4
    L_barrel   = 4.24        # m,  their Table 4
    cap_aspect = 1.0         # spherical caps
    t_wall     = 0.010       # m,  their t_Al = 10 mm, Table 4
    Q_budget   = None        # not used in benchmark mode; flux set explicitly
    P_gauge    = 170000.0    # Pa, their delta_p = 1.7 bar
    CASE_TAG   = 'BAGARELLO'


# --- JOB NAME PREFIX -------------------------------------------------------
# Abaqus job names are built from the load case and fill ratio alone, which
# means an A320 run would write HT_Baseline_f10.odb on top of the E190 file of
# the same name. The CSV files were always safe - CASE_TAG is in their names -
# but the .odb files were not.
#
# E190 KEEPS ITS BARE NAMES. The E190 matrix of 19 September is already on
# disk as HT_Baseline_f10.odb and so on, and renaming it would invalidate
# every reference to those files. Only the other cases take a prefix:
#
#     E190        HT_Baseline_f10.odb
#     A320        HT_A320_Baseline_f10.odb
#     A350        HT_A350_Baseline_f10.odb
#     benchmark   HT_BAGARELLO_Baseline_f10.odb
#
# Combined with a separate output directory per aircraft (E190_WORKDIR), the
# three runs cannot touch one another's results.
TAG_PREFIX = '' if CASE_TAG == 'E190' else (CASE_TAG + '_')


# --- Material: Al 2219-T87 ------------------------------------------------
# Yield is taken at ROOM temperature rather than the higher cryogenic value.
# This is deliberate and conservative: not every point of the wall reaches
# 20 K. Note that Bagarello use 500 MPa (the cryogenic value) with the same
# safety factor, giving 333 MPa allowable against our 262 MPa.
E_mod     = 73.1e9       # Pa
nu        = 0.33         # -
rho_shell = 2850.0       # kg/m3
alpha     = 22.3e-6      # 1/K
k_shell   = 116.0        # W/m-K
cp_shell  = 864.0        # J/kg-K
Fty_RT    = 393e6        # Pa, room-temperature yield
SF        = 1.5          # CS-25.303
sigma_allow = Fty_RT / SF                                   # FORMULA 10

# Property overrides for the cryogenic-property check (added 28 Sep 2026).
# alpha and k_shell above are ROOM-TEMPERATURE values, as in the benchmark
# (Bagarello Table 1). Near 20 K the secant expansion coefficient of
# aluminium is roughly 20 to 30 times smaller (NIST, 6061-T6). The check
# brackets the effect without changing any default:
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

# --- Stability (added 27 Sep 2026) -----------------------------------------
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
P_STABILITY = 120000.0 - 101325.0     # Pa gauge = 0.18675 bar
FOS_BUCKLING = 1.5
# Optional minimum wall. OFF by default: the agreed basis is buckling + yield
# only. To apply the ASME BPVC VIII-1 UG-16(b) floor of 1/16 in. (1.5 mm):
#     set E190_MIN_GAUGE_MM=1.5
MIN_GAUGE = float(os.environ.get('E190_MIN_GAUGE_MM', '0')) / 1000.0

# --- Fluid -----------------------------------------------------------------
rho_LH2 = 70.8           # kg/m3
T_LH2   = 20.3           # K, Bagarello Sec. 3.1
T_GH2   = 20.3           # K, ullage gas bulk temperature
# --- STRESS-FREE REFERENCE TEMPERATURE -------------------------------------
# THIS IS A MODELLING DECISION WITH A LARGE CONSEQUENCE. READ BEFORE CHANGING.
#
# Set to 20.3 K, matching Bagarello et al. Table 1, which lists a reference
# temperature of 20 K for the aluminium. The vessel is therefore defined as
# stress-free in its COLD OPERATING STATE, and the only thermal load the
# structure sees is the wetted-to-dry gradient of roughly 2.6 to 3.9 K.
#
# WHY NOT 293 K, THE ASSEMBLY TEMPERATURE
# ---------------------------------------
# Using 293 K applies the full 273 K cooldown on top of every load case. That
# was tried and produced 649.5 MPa, which is 98 percent of the biaxial fully
# restrained value E*alpha*dT/(1-nu) = 663 MPa, against a 262 MPa allowable.
# The ring constraints lock the vessel diameter, so the 7.9 mm of radial
# contraction the tank wants is resisted entirely by the metal.
#
#   excursion from stress-free : 272.7 K        ->  2.58 K
#   radial contraction wanted  :   7.9 mm       ->  0.075 mm
#   axial contraction wanted   :  47.4 mm       ->  0.449 mm
#   stress if fully blocked    : 663 MPa        ->  6.3 MPa
#
# A factor of 106 less for the supports to resist.
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
# The cooldown itself is no longer covered. It belongs in its own analysis:
# a ground condition at 1 g with no crash factors, where the question is
# whether the supports can accommodate 47 mm of contraction. Biancotto
# Sec. 5.8.1 addresses precisely this and finds that supports resisting the
# inner vessel's thermal displacement are overstressed by 1.6 to 3.8 times.
# Recorded as further work rather than silently omitted.
T_ref   = 20.3           # K, stress-free reference (cold operating state)
T_ASSEMBLY = 293.0       # K, retained for the cooldown case, not used here
g       = 9.81           # m/s2

# --- Thermal boundary conditions ------------------------------------------
# Inner surface, WETTED region: liquid hydrogen against aluminium. The film
# coefficient is high, so the wall is effectively pinned at the fuel
# temperature. Bagarello justify this by "the high thermal conductivity of the
# aluminium alloy and the great transfer at the fluid-structure interface".
h_wetted = 1000.0        # W/m2-K, effectively pins the wall to T_LH2

# Inner surface, DRY region: hydrogen vapour by natural convection.
# THIS IS THE KEY NUMBER. It is far smaller than h_wetted, which is why the
# dry wall runs warmer and a gradient forms.
h_dry    = 1.46          # W/m2-K, Bagarello Sec. 3.1, their ref [45]

# Outer surface: heat arriving through the MLI blanket, applied as a uniform
# inward flux. The inner vessel is the only body modelled, so the insulation
# enters solely through this number.
if not BENCHMARK_MODE:
    MLI_FLUX_MODE = 'FROM_BUDGET'   # q = Q_budget / A_shell      FORMULA 5
else:
    MLI_FLUX_MODE = 'EXPLICIT'
    # ----------------------------------------------------------------
    # BENCHMARK OUTER FLUX - DERIVED FROM THEIR INPUTS, NOT THEIR ANSWER.
    #
    # The previous value of 22.0 W/m2 was a placeholder annotated
    # "representative of 0.5 m PUR foam". It was not Bagarello's number, and
    # it drove the dry wall to 35.4 K against their reported 25.3 K.
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
# Sign convention: +x forward, +y lateral, +z up.
#   Baseline  : 1 g, level flight
#   LC1       : CS-25.337 symmetric manoeuvre
#   LC3       : constructed combined case (see thesis Methodology)
#   LC7       : CS-25.561 emergency landing, three largest factors applied
#               simultaneously. This is more severe than the regulation, which
#               specifies the factors acting separately.
load_cases = {
    'Baseline':      {'nx': 0.0, 'ny': 0.0, 'nz': 1.0},
    'LC1_Maneuver':  {'nx': 0.0, 'ny': 0.0, 'nz': 2.5},
    'LC3_Combined':  {'nx': 1.5, 'ny': 0.0, 'nz': 2.5},
    'LC7_Emergency': {'nx': 9.0, 'ny': 3.0, 'nz': 6.0},
}

# ----------------------------------------------------------------------------
# OPTIONAL LOAD-CASE FILTER - FOR RUNNING THE MATRIX IN PARALLEL SESSIONS
# ----------------------------------------------------------------------------
# WHY THIS EXISTS
# ---------------
# Measured on the 18 September production run: 32 solver jobs consumed 259
# seconds of actual computation in 79 minutes of wall clock. Ninety-five per
# cent of the run was spent queueing for an Abaqus/Standard licence token,
# with the campus pool reporting
#
#     <0 out of 500 licenses remain available>
#
# The script submits one job at a time and waits. With one request
# outstanding, it sits at the back of a shared queue forty times over.
#
# WHAT THIS CHANGES, AND WHAT IT DOES NOT
# ---------------------------------------
# Setting E190_ONLY_LC restricts this session to a subset of the load cases,
# so the matrix can be split across several Abaqus sessions running at the
# same time. Four sessions hold four queue positions instead of one, and are
# served roughly four times as often.
#
# NOTHING ABOUT THE ANALYSIS CHANGES. Each case is still built from scratch,
# meshed at the same seed, solved by the same solver with the same loads, and
# written to its own job name. The cases are independent of one another and
# always were - running them in four sessions rather than one cannot alter a
# single number. Only the order in which they reach the licence server
# differs.
#
# HOW TO USE - four Command Prompt windows, one line each:
#
#     set E190_ONLY_LC=Baseline
#     abaqus cae noGUI=python_scripts\E190_MASTER_thermo_structural.py
#
# then the same with LC1_Maneuver, LC3_Combined, LC7_Emergency. Each session
# writes its own CSV, named after the filter, so nothing is overwritten:
#
#     E190_production_Baseline_results.csv
#     E190_production_LC1_Maneuver_results.csv       ... and so on
#
# Concatenate the four afterwards to obtain the full twenty-row matrix. Each
# session also checks out one cae licence, so keep the number of sessions
# below the number of free cae seats.
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
fill_ratios = [0.1, 0.3, 0.5, 0.75, 0.9]

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
# The shell bending boundary layer at a cap-to-barrel junction has a
# characteristic length of sqrt(R*t). For the E190 that is 78.9 mm. The
# element size must resolve it, which the previous model's 341 mm elements did
# not. seed_size below is chosen as a fraction of that length.
MESHED_MASS = 0.0        # filled in by build_geometry_and_mesh, kg

# Axial half-width of the zone excluded around each ring when reporting peak
# stress. 0.5 m is 6.3 shell decay lengths sqrt(R*t) for the E190. See the
# FAR-FIELD block in extract() for the full justification.
FAR_FIELD_EXCLUSION = 0.5   # m

boundary_layer = math.sqrt(R_mid * t_wall)
#
# MESH_FRACTION IS SET FROM THE CONVERGENCE STUDY, NOT ASSUMED.
# Running E190_CHECK_smoketest.py with MODE = 'CONVERGENCE' gave, for LC7 at
# fill 0.9:
#
#     el/BL   elements   peak vM      change    deviation from finest
#      0.5       5,263   65.85 MPa        -            -3.09 %
#      1.0      18,989   67.73 MPa    +2.86 %          -0.33 %
#      2.0      64,372   68.03 MPa    +0.44 %          +0.11 %
#      4.0     259,944   67.96 MPa    -0.11 %           0.00 %
#
# 0.5 is adopted: two elements across the bending boundary layer, within
# 0.11 % of the finest mesh at a quarter of the element count.
MESH_FRACTION = 0.5                       # elements per boundary layer = 2
seed_size = boundary_layer * MESH_FRACTION

# Densities swept when RUN_MODE = 'CONVERGENCE'.
convergence_fractions = [1.0, 0.5, 0.25]
# 27 Sep 2026: 1, 2 and 4 elements per boundary layer. At the re-sized walls
# the boundary layer sqrt(R t) is short, and the old finest level (8 per layer)
# would be 2 to 4 million elements. The earlier studies at the thicker walls
# covered 0.5 to 8 elements per layer.
# Optional: drop levels from the environment, e.g.  set E190_CONV_FRACTIONS=2,1,0.5,0.25
# (a thinner wall shortens sqrt(R t), so the finest level grows as 1/t).
_cf = os.environ.get('E190_CONV_FRACTIONS', '').strip()
if _cf:
    convergence_fractions = [float(x) for x in _cf.split(',') if x.strip()]
    print('CONVERGENCE LEVELS OVERRIDE: %s' % convergence_fractions)

n_bins = 64   # circumferential/height bins for the hydrostatic pressure field


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
# The obvious approximation 2*pi*R*h_cap understates this by 59 percent, which
# made the whole shell area 7.9 percent low. That error showed up as a
# systematic 1.079 ratio in the X and Z components of the load balance check
# in every single production case - the model was correct, the reference value
# was not.
#
# SPHERICAL CAPS ARE A SPECIAL CASE AND MUST BE HANDLED SEPARATELY.
# When cap_aspect = 1.0 the caps are hemispheres, h_cap = R_mid, and the
# eccentricity is exactly zero. The oblate-spheroid formula then divides by
# e = 0. This is not hypothetical: the Bagarello benchmark tank uses
# spherical caps, so BENCHMARK_MODE = True hits it on the first run.
# The limit of the formula as e -> 0 is the sphere, A = 4*pi*R^2, which is
# what is used below.
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
print('E190 INNER VESSEL - THERMO-STRUCTURAL ANALYSIS (Bagarello methodology)')
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
# exactly as under gravity, except that "down" is the direction of the applied
# inertial body force rather than the vertical.
#
# The earlier version measured depth VERTICALLY while scaling the pressure by
# the magnitude of the resultant load factor. That is exact for Baseline and
# LC1, whose load factors are purely vertical, but wrong for LC3 and LC7. The
# standalone check in check_tilted_free_surface.py quantified the error:
#
#     load case   tilt      depth now    depth correct    p_hyd ratio
#     Baseline     0.0 deg    2.188 m       2.188 m          x1.00
#     LC1          0.0 deg    2.188 m       2.188 m          x1.00
#     LC3         31.0 deg    2.188 m       7.903 m          x3.6
#     LC7         57.7 deg    2.188 m      11.441 m          x5.2
#
# The tank is 15.6 m long and 2.6 m tall, so tilting the effective gravity
# toward the axis makes the liquid column far deeper than the tank is tall.
# At f = 0.90 the LC7 hoop stress at the deepest point rises 43 percent.
#
# Baseline and LC1 return a ratio of exactly 1.00, so this is a strict
# GENERALISATION of the previous treatment, not a replacement: those two load
# cases are unaffected.

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

    Taken directly from the body force the model applies - comp1 = -nx*g,
    comp2 = -nz*g, comp3 = -ny*g - so the liquid settles in the same direction
    the structure is pushed, which is what relative equilibrium requires.
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
    FORMULA 3. Height of the horizontal liquid free surface, measured from the
    tank axis, for a given fill ratio.

    For a horizontal cylinder the submerged area is a circular segment, and
    the fill ratio relates to the segment half-angle theta through

        theta - sin(theta) = 2*pi*phi_f

    which has no closed-form inverse, so Newton's method is used. The result
    is converted back to a height by y = -R*cos(theta/2).

    NO LONGER CALLED. RETAINED FOR THE ANALYTICAL CROSS-CHECK ONLY.
    ---------------------------------------------------------------
    This assumes the free surface is HORIZONTAL, i.e. that the acceleration
    acts vertically. That is true for Baseline and LC1 and false for LC3 and
    LC7. It was used by the thermal wetted/dry split until the split was moved
    onto the tilted surface; using it anywhere in the load path again would
    reintroduce a thermal model that ignores the load case.

    It remains useful as a verification: for a purely vertical load factor,
    free_surface_offset() must reproduce -fill_surface_height(phi_f) to within
    the sampling error of the point cloud.
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
    Vertical coordinate of an element centroid.

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

    Vertical is therefore Y, which is coordinates[1].

    An earlier version used coordinates[2] (Z, the horizontal direction across
    the tank). The symptom was a von Mises field that was symmetric top to
    bottom with a longitudinal band down the middle, instead of increasing
    steadily from the crown to the keel. In effect the tank was being filled
    sideways.
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

    Shell elements are used rather than solids. The wall is 4.785 mm on a
    1.3015 m radius (R/t = 272), which is firmly in the thin-shell regime, and
    with k = 116 W/m-K the through-thickness temperature drop is negligible,
    so a shell conduction model is adequate.
    """
    model = mdb.Model(name=model_name)

    # --- Sketch the meridian: half an ellipse, the barrel, half an ellipse --
    s = model.ConstrainedSketch(name='meridian', sheetSize=4.0 * L_total)
    s.ConstructionLine(point1=(0.0, 0.0), point2=(1.0, 0.0))

    # Forward cap: quarter ellipse from the pole to the barrel start.
    # Parametrised so the profile is exact rather than a spline fit.
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

    # Aft cap, mirror of the forward one
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
    # through-thickness gradient, which is what a predefined temperature field
    # read from a heat transfer .odb supplies.
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
    # others. Since the gauge pressure alone is 150 kPa over a 5.32 m2 cross
    # section, i.e. 798 kN, even a few percent of flipped elements produces
    # hundreds of kN of spurious resultant on what should be a
    # self-equilibrating load. The technique used is therefore reported.
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
    PROCESS 6 (revised).  THREE CIRCUMFERENTIAL RING CONSTRAINTS.

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
    few elements of a 4.785 mm shell, giving a spurious local stress of the
    order of a thousand MPa. That is the same punch-through seen earlier with
    the 3-2-1 restraint. A ring spreads the reaction over the full
    circumference, which is what real tank hardware does.

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

    The anchor-and-slide scheme is nevertheless retained rather than fixing
    all three rings. Were the reference ever moved back to 293 K, fully
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
    # them exactly. Mid-barrel has no such feature, and a fixed tolerance
    # caught a ragged partial set: an earlier run found 124 nodes at each
    # junction but only 31 at mid-barrel, against roughly 52 for a complete
    # circumferential ring. The axial anchor was therefore gripping only part
    # of the circumference, which is an asymmetric restraint.
    #
    # Instead, find the x coordinate of the node plane nearest the requested
    # station, then take every node on that plane. This always yields one
    # complete ring regardless of where the mesher happened to put nodes.
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
    PROCESS 6.  INERTIA RELIEF - the equilibrium device used instead of mounts.

    WHY A 3-2-1 RESTRAINT IS NOT SUFFICIENT HERE
    --------------------------------------------
    A statically determinate 3-2-1 restraint carries zero reaction only when
    the applied loads are SELF-EQUILIBRATED. Thermal expansion and uniform
    internal pressure satisfy that. Inertia does not.

    Under LC7 the hydrostatic pressure field has a net resultant of

        m_LH2 * n * g  =  5104 kg * 11.225 * 9.81  =  562 kN

    and the shell body force adds a further 181 kN. With no mount, nothing
    balances them, so the entire 743 kN passes through three nodes of a
    4.785 mm shell. An earlier run of this model did exactly that and returned
    4943 MPa and a 6.9 m displacement: not a physical result, simply a point
    load punching through a thin plate.

    THE UNDERLYING POINT
    --------------------
    A CS-25 load factor presupposes a restraint. It describes the airframe
    decelerating the tank. A genuinely unrestrained body does not experience a
    load factor at all; it simply accelerates. "No mounts" and "CS-25 inertia"
    are therefore in tension, and the tension has to be resolved explicitly.

    Inertia relief is the standard resolution: the applied loads are balanced
    by a distributed d'Alembert body force proportional to the mass, rather
    than by a concentrated reaction. Equilibrium is satisfied, no artificial
    restraint is introduced, and no stress concentration is created. This is
    the approach taken by Oom Ortiz de Montellano et al., who tested both a
    constrained model and an inertia relief model of the same structure and
    adopted inertia relief for all subsequent analyses.

    ACCEPTED APPROXIMATION
    ----------------------
    The balancing acceleration is distributed over the model's mass, which
    here is the shell only: the liquid appears as a pressure field, not as
    mass. The reaction is therefore spread over the wall rather than placed
    where the liquid actually sits. Because the pressure also acts over the
    wall, the mismatch is in the distribution rather than the magnitude, and
    it is bounded. It is recorded as a limitation rather than hidden.

    Note that Abaqus requires the free directions under inertia relief to be
    unconstrained, so no displacement boundary condition is applied.
    """
    model.InertiaRelief(name='InertiaRelief', createStepName=step_name,
                        u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
    print('  inertia relief active in all six directions (no mounts, no BCs)')


def apply_rigid_body_restraint(model, inst):
    """
    3-2-1 STATICALLY DETERMINATE RESTRAINT.

    RETAINED FOR REFERENCE ONLY - not called in the current analysis chain.
    See apply_inertia_relief above for why. This is valid only for load cases
    whose applied loads are self-equilibrated, such as a thermal-only or a
    uniform-pressure-only study.

    This is NOT a mount. It removes the six rigid body modes so that the
    stiffness matrix is non-singular, and nothing else. Because it is
    statically determinate it carries zero reaction under a self-equilibrated
    load, and it does not resist the thermal contraction: the vessel shrinks
    freely by alpha*dT*L, which for the E190 is 94.9 mm.

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

    WHY THIS EXISTS, AND WHY mdb.Job.status IS NOT USED.
    ----------------------------------------------------
    The obvious test is

        jt.waitForCompletion()
        if jt.status != COMPLETED: ...

    and it works inside Abaqus/CAE, because the GUI runs a job monitor that
    keeps the Job object's status member up to date. A headless run started
    with

        abaqus cae noGUI=E190_MASTER_thermo_structural.py

    has no such monitor, so status can still read RUNNING - or nothing useful
    at all - after waitForCompletion() has returned. The test then rejects a
    job that in fact succeeded.

    That is precisely what happened on the first full headless attempt: all
    twenty heat transfer jobs solved correctly, each reporting

        Abaqus JOB HT_<tag> COMPLETED
        THE ANALYSIS HAS COMPLETED SUCCESSFULLY
        0 ERROR MESSAGES

    in twelve seconds apiece, and all twenty were discarded by the status
    check before the structural half could run. Twenty good thermal .odb files
    on disk, no results, and a misleading "job failed" message for each.

    The job's own .log file is the authority: Abaqus writes the COMPLETED line
    there itself, and it is identical in GUI and headless runs. So that is
    what is read here.
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

    THE SPLIT USES THE TILTED FREE SURFACE, NOT A HORIZONTAL ONE.
    DO NOT REVERT TO elem_centroid_z AND fill_surface_height.
    --------------------------------------------------------------------
    An earlier version partitioned by vertical centroid height against
    fill_surface_height(phi_f), which takes no load case at all. The thermal
    model therefore placed the liquid in a horizontal band along the bottom of
    the tank for EVERY load case, while the structural model placed it against
    the tilted resultant. The two halves of a coupled analysis disagreed about
    where the fuel was.

    The signature in the results was plain once looked for: the wetted-to-dry
    dT came out bit-identical across Baseline, LC1, LC3 and LC7 at every fill
    ratio (3.62878036499023 K at phi_f = 0.10 for all four). That cannot be
    physical. Under LC7 the resultant tilts 57.7 degrees and the fuel surges
    forward: at phi_f = 0.10 the real wetted zone is a wedge in the last 1.9 m
    of the barrel, roughly 7 per cent of the wall, not a strip running its
    whole length. Different wetted footprint, different temperature field.

    Returns SURFACES, not sets. Film conditions and surface heat fluxes act on
    a surface, because a shell has two faces and Abaqus has to be told which
    one the load applies to. Passing a plain element Set raises a region error.

    SIDE CONVENTION: side2 is the INNER face on this revolved geometry, so it
    is side2 that carries the film conditions. This is the same convention the
    pressure loads use, and it was verified there by the sign of the load
    balance - see the note at the Pressure() call.
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
    One (load case, fill ratio) point. Runs the heat transfer analysis, then
    the static analysis reading its temperature field.

    ht_from (added 27 Sep 2026): path of an existing heat-transfer job for
    the SAME wall, mesh, load direction and fill. Pressure does not enter
    the heat-transfer problem, so the temperature field at 0.187 bar and at
    1.5 bar is identical; when ht_from is usable it is read instead of being
    solved again. The structural half checks that its mesh matches it node
    for node and stops if it does not.
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

        # FORMULA 5: heat arriving through the MLI, on the OUTER face. The inner
        # face carries the film conditions above, so this goes on side2.
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

    # PROCESS 6. Three circumferential rings after Bagarello Sec. 4.4.1.
    # Inertia relief is no longer used: with a real load path the applied
    # loads are reacted at the rings, which is what a CS-25 load factor
    # physically describes. apply_inertia_relief() is retained below for the
    # free-body reference case.
    apply_ring_constraints(model, inst, 'Static', seed)

    # INITIAL TEMPERATURE - THE STRESS-FREE STATE. DO NOT REMOVE.
    #
    # Abaqus computes thermal strain as alpha * (T_current - T_INITIAL), NOT
    # relative to the Expansion 'zero' value. Without an initial condition the
    # initial temperature defaults to 0 K, and the model then believes the
    # vessel is being warmed from 0 K, which is wrong by any reckoning. An
    # earlier run showed exactly that: u_max of 3.22 mm where the thermal load
    # was effectively absent.
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
    # THE FIELD IS READ AND RE-APPLIED EXPLICITLY. DO NOT REVERT TO
    # distributionType=FROM_FILE. THE REASON IS MEASURED, NOT ASSUMED.
    #
    # The shell section is declared `temperature=GRADIENT`, which means every
    # temperature supplied to it is a PAIR: a mid-surface value and a
    # through-thickness gradient. A predefined field read FROM_FILE off the
    # heat transfer .odb does not reliably deliver a zero second term, so the
    # wall ends up carrying a through-thickness temperature difference that
    # does not exist in the heat transfer solution at all.
    #
    # WHAT IT COST. Measured on a mid-barrel element far from every ring,
    # ST_Baseline_f10.odb:
    #
    #     SNEG  hoop 85.99 MPa   axial 64.56 MPa
    #     SPOS  hoop 36.75 MPa
    #
    #     membrane hoop = (85.99 + 36.75)/2 = 61.37 MPa   correct: 40.80
    #     bending  hoop = (85.99 - 36.75)/2 = 24.62 MPa   correct: ZERO
    #
    # A cylinder under uniform internal pressure has NO through-thickness
    # bending at mid-span. The 24.62 MPa implies a through-thickness dT of
    # 24.62 / [E*alpha/(2(1-nu))] = 20.2 K. The wall is 4.785 mm of aluminium
    # at k = 116 W/m-K carrying about 9 W/m2; the real through-thickness drop
    # is of order 1e-4 K. The gradient is entirely an artefact of the
    # coupling, and it inflated every reported stress by roughly a factor of
    # two. It is invisible to the load-balance check because bending through
    # the thickness is self-equilibrated.
    #
    # WHAT IS DONE INSTEAD. The nodal temperatures are read straight out of
    # the heat transfer .odb, binned, and applied as UNIFORM fields with
    # crossSectionDistribution=CONSTANT_THROUGH_THICKNESS. That keyword forces
    # the gradient term to zero by construction, so the artefact cannot recur.
    # This mirrors the pressure binning used below, which is already proven.
    #
    # ACCEPTANCE TEST after any change here. Probe a mid-barrel element at
    # least 1 m from any ring. BOTH section points must agree, at hoop
    # 40.8 MPa and axial 20.4 MPa, von Mises 35.3 MPa. If SNEG and SPOS
    # disagree at mid-span, the coupling is wrong and the results are not
    # usable.
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
    # as well. Recomputing them in two places is how the thermal and
    # structural halves drifted apart in the first place.
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
        # side2 IS THE INNER FACE on this revolved geometry. Verified by the
        # load balance: with side1 the hydrostatic reaction came out at
        # -561.4 kN against an applied +562.0 kN - correct magnitude, opposite
        # sign - meaning the liquid was pushing inward instead of outward.
        bsurf = model.rootAssembly.Surface(name='PS_%s_%02d' % (tag, idx),
                                           side2Elements=bset.elements)
        model.Pressure(name='P_%s_%02d' % (tag, idx), createStepName='Static',
                       region=bsurf, magnitude=p_tot,
                       distributionType=UNIFORM)

    # FORMULA 7: inertial body force on the shell's own mass. The liquid's
    # inertia is carried by the scaled hydrostatic pressure above, not here.
    # AXIS MAPPING AND SIGNS. The vessel is revolved about global X, so:
    #     comp1 = global X = FORWARD  -> n_x
    #     comp2 = global Y = UP       -> n_z
    #     comp3 = global Z = LATERAL  -> n_y
    #
    # The vertical and lateral components were transposed in an earlier
    # version, which applied 6 g sideways and 3 g downwards under LC7.
    #
    # SIGNS. The CS-25.561 load factors are named by the direction in which the
    # inertia force ACTS ON THE ITEM. "Downward 6.0 g" is a force downward, and
    # +Y is up, so comp2 = -nz*g. "Forward 9.0 g" is a force FORWARD, and +X is
    # forward, so comp1 = +nx*g. The two therefore carry opposite signs, which
    # looks odd but is correct: the asymmetry comes from the axis naming, not
    # from the physics.
    #
    # An earlier version used -nx*g, which threw the contents AFT. The stress
    # magnitudes were identical because the problem is mirror-symmetric, but
    # the liquid piled against the wrong cap and every contour showed the peak
    # at the wrong end of the tank.
    # UNITS. BodyForce (*DLOAD type BX/BY/BZ) expects force per unit VOLUME,
    # i.e. rho*a, not an acceleration. Passing the acceleration directly made
    # the shell inertia 2850 times too small: 51 N instead of 181.5 kN, which
    # is why the mount reaction came back as the liquid resultant alone.
    #
    # Gravity takes an ACCELERATION and multiplies by the material density
    # itself, so it cannot be got wrong this way. It is used in preference.
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
    Remove finished models and jobs from the CAE session (added 27 Sep 2026).
    Each case builds two models with a full mesh and hundreds of load sets.
    Twenty cases fit in memory; the ~70 of a chained FINAL session may not.
    The ODBs on disk are untouched.
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
    Pull the quantities the thesis reports.

    FORMULA 9 (von Mises) is evaluated by the solver; this reads the maximum.
    Reaction forces are summed as a check that the 3-2-1 restraint really is
    carrying nothing.
    """
    odb_t = openOdb(path=jname_t + '.odb')
    fr_t = odb_t.steps['Thermal'].frames[-1]
    temps = [v.data for v in fr_t.fieldOutputs['NT11'].values]
    T_min, T_max = min(temps), max(temps)
    # Peak wall heat flux (added 27 Sep 2026), so the waterline flux quoted in
    # the report is a recorded number, not one read off a colour legend.
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
    # refined and never converges. The convergence study showed exactly
    # this behaviour, with the peak rising 108.0 -> 108.0 -> 122.1 -> 140.6
    # MPa on successive refinements, the increments GROWING rather than
    # shrinking, and a fitted exponent of sigma ~ h^-0.19.
    #
    # Meanwhile the global response converged cleanly over the same range:
    # the mount reaction to 0.2 percent and the peak displacement to 0.9
    # percent. The divergence is therefore local to the constraint and is
    # not a structural response.
    #
    # The peak is consequently taken OUTSIDE an exclusion zone around each
    # ring. The justification is Saint-Venant's principle: a local load
    # introduction decays over a distance comparable with the shell decay
    # length sqrt(R*t), which is 78.9 mm here. The 0.5 m zone adopted is
    # 6.3 decay lengths, so the reported stress is a genuine far-field
    # value governed by global equilibrium rather than by the details of
    # the attachment.
    #
    # This follows Oom Ortiz de Montellano et al., who documented the same
    # artefact for constrained boundary conditions, and is the treatment
    # already set out in section_boundary_conditions.tex.
    #
    # WHAT THIS EXCLUDES: the local fitting design at the rings, which in a
    # detailed design would need a doubler or wear plate. Stated as a
    # limitation, not concealed.
    ring_x = [0.0, L_barrel / 2.0, L_barrel]

    # element label -> centroid x, taken from the odb mesh.
    #
    # PERFORMANCE. getNodeFromLabel() is a lookup, and calling it once per
    # node per element is roughly a million calls on the finest mesh, which
    # takes minutes in the Abaqus interpreter. The node coordinates are
    # therefore hoisted into a dictionary once, turning the inner loop into
    # a constant-time dictionary access.
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

    # ---------------- AXIAL MEMBRANE STRESS FOR THE BUCKLING SCREEN -------
    # Added 27 Sep 2026. Barrel only (SP-8007 is a cylinder result), far
    # field only (outside FAR_FIELD_EXCLUSION of every ring).
    #
    # S11 IS AXIAL HERE. No material orientation is assigned, so Abaqus takes
    # local direction 1 as the projection of global X onto the shell, and X is
    # the tank axis. (The superseded sizing scripts were Z-axial, where the
    # same default made S11 the HOOP stress at the crown and keel.) This is
    # checked, not assumed: under internal pressure the mean hoop membrane
    # stress must exceed the mean axial one, so if the mean of S11 exceeds
    # the mean of S22 the orientation is wrong and the run says so.
    #
    # Membrane = mean of the two section points (SNEG, SPOS); bending, the
    # half-difference, is excluded, as the SP-8007 screen consumes membrane
    # stress.
    sp_vals = {}
    cap_vals = {}
    for v in s_field.values:
        x_ = cx.get(v.elementLabel, None)
        if x_ is None:
            continue
        in_barrel = (FAR_FIELD_EXCLUSION <= x_ <= L_barrel / 2.0 - FAR_FIELD_EXCLUSION or
                     L_barrel / 2.0 + FAR_FIELD_EXCLUSION <= x_ <= L_barrel - FAR_FIELD_EXCLUSION)
        if x_ < 0.0 or x_ > L_barrel:
            # End caps (heads), recorded for a head-stability check: under
            # internal pressure an a/b = 1.6 head carries hoop COMPRESSION
            # near its equator, (pR/t)(1 - a^2/2b^2) in membrane theory.
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

    # ---------------- RECORDED CHECKS (added 27 Sep 2026) -----------------
    # (1) MEMBRANE AND THROUGH-THICKNESS PROBE. The crown element nearest
    #     x = L/4 (mid aft bay, at least 1 m from any ring on all three tanks)
    #     and its S11 (axial) and S22 (hoop) at both section points. Under a
    #     vertical load at low fill the crown is dry, so hoop ~ pR/t and
    #     axial ~ pR/2t, and the two surfaces must agree (no fictitious
    #     through-thickness bending). Previously done by hand; now recorded.
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
    #     a recorded result rather than a beam-theory estimate.
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
    #   hydrostatic resultant : rho_LH2 * |n| * g * V_liquid, straight DOWN,
    #                           because the model measures depth vertically
    #
    # The uniform p_gauge term nets to zero on a closed surface and does not
    # appear. If it DOES appear, the shell normals are inconsistent and some
    # elements are being pressurised on the wrong face - which would corrupt
    # the stress field as well as the balance.
    nx_, ny_, nz_ = factors['nx'], factors['ny'], factors['nz']
    n_mag_ = math.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_)
    # Use the mass Abaqus actually meshed, not the analytic estimate. The
    # analytic A_shell approximates the cap area as 2*pi*R*h_cap and
    # underestimates by about 7.7 percent, which showed up as a spurious
    # 7.7 percent imbalance in X and Z.
    m_shell = A_shell * t_wall * rho_shell
    if MESHED_MASS is not None and MESHED_MASS > 0.0:
        m_shell = MESHED_MASS
    m_liq = V_tank * rho_LH2 * phi_f

    # Global axes: 1 = X forward, 2 = Y up, 3 = Z lateral.
    #
    # THE LIQUID RESULTANT ACTS ALONG THE TILTED DIRECTION, NOT VERTICALLY.
    # An earlier version of this check placed the whole liquid resultant on
    # the Y axis. That was correct while the model measured depth vertically,
    # but wrong once the tilted free surface was introduced, and it produced
    # component ratios of 3.8 and 0.6 on cases that were in fact balanced to
    # better than one percent.
    #
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
#   PRODUCTION   4 load cases x 5 fills at the design pressure (strength).
#   STABILITY    the same 20 cases at the stability pressure (buckling).
#   CONVERGENCE  LC7 fill 0.9 over a sweep of mesh densities, at the design
#                pressure AND at the stability pressure.
#   SIZING       finite-element wall sizing, recorded in sized_walls.csv.
#   FINAL        SIZING -> PRODUCTION -> STABILITY -> CONVERGENCE in one
#                session (added 27 Sep 2026). Each stage writes its CSV the
#                moment it finishes; a stage whose CSV already exists at the
#                current wall is skipped, so a crash costs one stage, not all.

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
    global t_wall, boundary_layer, seed_size
    t_wall = t_new
    boundary_layer = math.sqrt(R_mid * t_wall)
    seed_size = boundary_layer * MESH_FRACTION


def record_sized_wall(fos_b, fos_y, govern=''):
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
    d = os.path.join(_ROOT_DIR, sub)
    if not os.path.isdir(d):
        os.makedirs(d)
    os.chdir(d)
    return d


def stage_sizing():
    # PROCESS 11 (revised 27 Sep 2026, after the first campus iteration).
    # The wall is sized by WHICHEVER of three requirements needs most metal:
    #   buckling  FOS >= 1.5 against SP-8007, at the stability pressure
    #   yield     far-field von Mises <= F_ty/1.5 = 262 MPa (fos >= 1.0),
    #             at the design pressure
    #   minimum gauge MIN_GAUGE (ASME BPVC VIII-1, UG-16(b): 1/16 in.)
    # Each iteration runs LC7 fill 0.9 at both pressures. Buckling FOS scales
    # about as t^2 and the yield fos about as t, which gives the next wall.
    # Converged when the governing requirement is met within 2 %, or when the
    # wall sits at the minimum gauge with both others satisfied.
    global P_gauge
    rows = []
    in_dir(os.path.join('sizing_fe', CASE_TAG))
    lc7 = load_cases['LC7_Emergency']
    t_try, sized, govern = max(t_wall, MIN_GAUGE), False, ''
    MIN_T = max(MIN_GAUGE, 0.2e-3)       # numerical floor only, never reached in practice
    # ACCEPTANCE (revised on campus, 27 Sep). The governing requirement must be
    # met with at most ~6 % surplus (buckling FOS 1.50-1.60, yield fos
    # 1.00-1.06). A narrower window made the A350 oscillate: the mesh changes
    # with the wall, so a 2.5 % change in t moved its buckling FOS from 1.44
    # to 1.59, straddling a 1.50-1.53 window forever. In addition, once a
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
        rb = run_case('LC7_Emergency', lc7, 0.9, seed_size, tag_b)
        P_gauge = P_DESIGN
        ry = run_case('LC7_Emergency', lc7, 0.9, seed_size,
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
    # Buckling verification on the other candidate cases at the final wall.
    checks = [('LC7_Emergency', f_) for f_ in (0.1, 0.3, 0.5, 0.75)] + \
             [('LC3_Combined', 0.9), ('LC1_Maneuver', 0.9), ('Baseline', 0.9)]
    for lc_name, f_ in checks:
        P_gauge = P_STABILITY
        r = run_case(lc_name, load_cases[lc_name], f_, seed_size,
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
    4 load cases x 5 fills at `pressure`. RESUMES: cases already in out_csv at
    this wall and pressure are kept, not re-run. ht_fmt, if given, names the
    heat-transfer job to reuse for each case (the stability matrix reuses
    the production thermal fields, which are identical).
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
    # PROCESS 7. LC7 fill 0.9 over a sweep of mesh densities, run at both
    # pressures: the design pressure for the far-field von Mises stress and
    # the stability pressure for the axial membrane stress the buckling
    # screen consumes. Resumes per level; the stability-pressure sweep reuses
    # the design-pressure thermal field of the same level.
    global P_gauge
    P_gauge = pressure
    rows = done_rows(out_csv, pressure)
    have = set(round(r['seed_mm'], 3) for r in rows)
    print('\nMESH CONVERGENCE %s (%s, %.5f bar): LC7 fill 0.9, levels %s (%d already done)'
          % (CASE_TAG, label, P_gauge / 1e5, convergence_fractions, len(have)))
    for frac in convergence_fractions:
        sd = boundary_layer * frac
        if round(sd * 1e3, 3) in have:
            continue
        ftag = str(frac).replace('.', 'p')
        tag = '%sconv_%s_%s' % (TAG_PREFIX, label, ftag)
        src = os.path.join(os.getcwd(), 'HT_%sconv_%s_%s' % (TAG_PREFIX, ht_label, ftag)) \
            if ht_label else None
        print('  seed %.1f mm (%.2f elements per boundary layer)' % (sd * 1e3, 1.0 / frac))
        r = run_case('LC7_Emergency', load_cases['LC7_Emergency'], 0.9, sd, tag, ht_from=src)
        if r:
            rows.append(r)
            write_rows(out_csv, rows)
    P_gauge = P_DESIGN
    return rows


results = []

# Which FINAL stages to run (added 27 Sep 2026, for a time-limited session):
#     set E190_STAGES=production            production only, then stop
#     set E190_STAGES=stability             stability only (needs production first)
# Default: all four, in order.
STAGES = [x.strip().lower() for x in os.environ.get(
    'E190_STAGES', 'sizing,production,stability,convergence').split(',') if x.strip()]
# Guard (27 Sep 2026): a FINAL session that does not size must use the recorded
# wall. If sized_walls.csv is missing or has no row for this aircraft, stop
# rather than run the matrix at the AIRCRAFT_DATA fallback thickness.
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
                                  and 'LC7_Emergency' in load_cases):
    # (in a split FINAL session, only the window that owns LC7 runs it)
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
    print('Now run:  abaqus python python_scripts\\check_campus_outputs.py')
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
