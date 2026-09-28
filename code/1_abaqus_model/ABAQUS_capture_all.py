# -*- coding: utf-8 -*-
"""
================================================================================
AUTOMATED CONTOUR CAPTURE FOR EVERY PRODUCTION ODB
TU Delft MSc thesis - Sarvesh Mahadevan

Opens every production ODB in turn and writes contour images at three views,
with locked contour limits, white background and no title block, so that every
frame in a set is directly comparable with every other.

RUN IT WITH (from the thesis root, the folder containing E190\\ A320\\ A350\\):

    abaqus cae noGUI=python_scripts\\ABAQUS_capture_all.py

Do NOT use "abaqus python" - that has no graphics context and cannot print
images. It must be "abaqus cae noGUI=..." (or "abaqus viewer noGUI=...").

--------------------------------------------------------------------------------
WHY THIS EXISTS RATHER THAN CLICKING
--------------------------------------------------------------------------------
Every frame in a comparison set must share a contour scale and a camera. Done by
hand across 120 ODBs and several sittings, the limits drift and the comparison
silently stops meaning anything: a colour in one frame no longer denotes the same
stress as that colour in the next. Here the limits and the viewpoint are the same
lines of code for every frame, so consistency is guaranteed rather than hoped for.

It also costs one licence checkout for the whole run instead of one per session.

--------------------------------------------------------------------------------
WHAT IT PRODUCES
--------------------------------------------------------------------------------
    abaqus_captures\\<AIRCRAFT>\\<AIRCRAFT>_<LOADCASE>_f<NN>_<field>_<view>.png

e.g.  abaqus_captures\\E190\\E190_LC7_Emergency_f90_mises_iso.png

fields : nt11, hfl   (from the HT odb)     mises  (from the ST odb)
views  : front, iso, sect

Full run = 3 aircraft x 4 load cases x 5 fills x (2 thermal + 1 structural)
           fields x 3 views = 540 images, roughly 1.5-2.5 h unattended.

--------------------------------------------------------------------------------
BEFORE THE FULL RUN
--------------------------------------------------------------------------------
Set QUICK_TEST = True and run it once. That does a single ODB (about 9 images)
so you can check the camera, the limits and the crop before committing hours.
Then set it back to False.
================================================================================
"""

from abaqus import *
from abaqusConstants import *
import visualization
# Restore Python's built-ins: the star imports above replace sum() and others
# with the Abaqus XY-data operators, which reject generators.
import builtins as _builtins
sum, any, all, max, min, abs, round = (_builtins.sum, _builtins.any, _builtins.all,
                                       _builtins.max, _builtins.min, _builtins.abs,
                                       _builtins.round)
import os
import sys
import re
import csv
import glob
import math

# ==============================================================================
#  REVISION 27 SEP 2026 - read this before running
# ==============================================================================
# 1. ODB SOURCE CORRECTED. The previous version read A320\ and A350\ sub-
#    folders, which hold the 19 Sep ODBs at the OLD wall thicknesses (7.440 and
#    15.250 mm). The corrected runs (6.382 and 9.948 mm) are in the thesis root.
#    find_odb() now looks there, and every ODB is checked against the wall
#    thickness in AIRCRAFT_DATA of the master script before it is used. An ODB
#    at the wrong thickness is SKIPPED and logged, never plotted.
# 2. CONTOUR WINDOWS FROM THE RESULTS. The thinner walls raise the stresses
#    (A350 LC7 now reaches 87.5 MPa against an old window top of 60), so the old
#    hand-typed windows would have saturated. The top of each window is now
#    taken from the production CSVs at run time. See contour_windows().
# 3. LEGEND NEVER COVERS THE TANK. For the views in SIDE_LEGEND_VIEWS each frame
#    is rendered twice, legend on and legend off, into <OUT_ROOT>\_raw. Then
#        python python_scripts\compose_side_legend.py
#    cuts the legend out (it is exactly where the two renders differ), cuts the
#    tank out of the legend-off render, and places them side by side. Overlap
#    is impossible by construction, whatever the tank shape.
# 4. Output goes to a NEW folder, abaqus_captures_final, so nothing rendered from
#    the old ODBs can be picked up by mistake.
# ==============================================================================

# ==============================================================================
#  CONFIGURATION - everything you might want to change is in this block
# ==============================================================================

QUICK_TEST = False        # True = one ODB only, to check camera and legend

# 'THESIS' = only the frames the report uses (a few minutes).
# 'ALL'    = every aircraft x load case x fill x field x view (1.5-2.5 h).
FIGURE_SET = 'THESIS'

OUT_ROOT = 'abaqus_captures_resized'   # images of the re-sized walls (27 Sep)
                            # canvas size is set per view, see VIEWS below

# Views whose legend is moved beside the tank by compose_side_legend.py.
# The front view is the one where the legend lands on the tank; add 'iso' or
# 'sect' here if you ever want the same treatment for those.
SIDE_LEGEND_VIEWS = ['front', 'bottom']

# Views rendered for the contour frames. Only the front view carries a legend;
# the section view is read against the front view of the same case, so it is
# rendered without one. The iso view added little and is no longer rendered.
CONTOUR_VIEWS = ['front', 'sect', 'bottom']
# 'bottom' (added 27 Sep 2026) looks up at the keel. At low fill the
# out-of-round deformation between the rings sits on the underside, where the
# front view cannot show it.

AIRCRAFT = ['E190', 'A320', 'A350']
LOADCASES = ['Baseline', 'LC1_Maneuver', 'LC3_Combined', 'LC7_Emergency']
FILLS = [10, 30, 50, 75, 90]

# Thermal fields are identical for Baseline and LC1 because both resultants are
# purely vertical, so the liquid sits in the same place. Chapter 5 establishes
# this. Set True to skip the duplicate thermal frames and save 90 images.
SKIP_REDUNDANT_THERMAL = False

# ------------------------------------------------------------------ views -----
# If you have saved a user view in CAE with one of these names, the script uses
# it and ignores the vector. Otherwise it falls back to the viewVector.
# 'size'   canvas in px. The front view of a 6:1 tank wants a wide, short frame;
#          the iso views run diagonally and need the height.
# 'legend' top-left corner of the legend box as (x%, y%) of the viewport. The
#          default (2, 98) puts it straight on top of the tank, which is what
#          the first test run did. Low-left sits in the dead space instead.
# 'zoom'   applied after fitView, so the model does not touch the frame edge.
VIEWS = {
    'front': {'saved': 'FRONT',
              'vector': (0.0, 0.0, 1.0), 'up': (0.0, 1.0, 0.0), 'cut': False,
              'size': (2400, 760), 'legend': (1, 98), 'zoom': 0.86},
    'iso':   {'saved': 'ISO45',
              'vector': (0.707, 1.0, 0.707), 'up': (0.0, 1.0, 0.0), 'cut': False,
              'size': (2400, 1100), 'legend': (1, 98), 'zoom': 0.86},
    'bottom': {'saved': None,
              'vector': (0.0, -1.0, 0.0), 'up': (0.0, 0.0, 1.0), 'cut': False,
              'size': (2400, 760), 'legend': (1, 98), 'zoom': 0.86},
    'sect':  {'saved': 'ISO45',
              'vector': (0.707, 1.0, 0.707), 'up': (0.0, 1.0, 0.0), 'cut': True,
              'size': (2400, 1100), 'legend': (1, 98), 'zoom': 0.86},
}

# Render style of the contour frames (27 Sep, 21:45). FILLED is flat and
# unlit, which reads as bright and flashy; SHADED lights the surface, so the
# colours sit slightly darker and the curvature of the shell is visible.
CONTOUR_RENDER = SHADED

# ------------------------------------------------------ zoomed mesh frames ---
# At the re-sized walls an element is 2 to 4 px in a full view, so a mesh-on
# full view is black. These frames zoom into the middle of the barrel (the mid
# ring, where the stress field has structure) far enough that the elements are
# about 15 to 25 px, and are shot twice, with and without edges, so the mesh
# and the contour it carries can both be shown.
ZOOM_SHOT = True
ZOOM_FACTOR = 5.0
ZOOM_CASES = [('LC7_Emergency', 90, 'mises'), ('LC7_Emergency', 90, 'nt11')]
ZOOM_VIEWS = ['front', 'iso']

# ALL, EXTERIOR, FEATURE, FREE or NONE.
# The first test run came back speckled: at 2400 px an E190 element is about
# 10 px, so 64372 element edges read as grey noise laid over the contours
# rather than as a mesh. NONE gives the cleanest contour figure; FEATURE keeps
# the silhouette and the cap-to-barrel junctions. Use FEATURE if NONE looks
# too floaty, and ALL only for the deliberate mesh figure in Chapter 4.
EDGE_STYLE = NONE

# ------------------------------------------------------------- mesh figure ----
# Chapter 4 needs one figure that DOES show the mesh, to document the element
# size derived from the bending boundary layer. That is the only place the
# edges belong: everywhere else they swamp the contours.
#
# Shot as a plain undeformed mesh with no contour and no legend, because the
# subject is the discretisation, not a result. The mesh is identical across
# load cases and fills, so any ODB will do.
MESH_SHOT = True         # plain mesh of every aircraft, for the 'three models' figure
MESH_SOURCE = {'loadcase': 'Baseline', 'fill': 50}
MESH_VIEWS = ['front', 'iso']   # front for the to-scale comparison, iso as before

# ------------------------------------------------------- mesh-on comparison ---
# One case shot twice, identical in every respect except the element edges, so
# you can see for yourself what the mesh costs a contour figure before
# committing to mesh-off for the other 540. Produces
#   <AC>_COMPARE_f90_mises_front_mesh.png
#   <AC>_COMPARE_f90_mises_front_nomesh.png
COMPARE_SHOT = False
COMPARE_SOURCE = {'loadcase': 'LC7_Emergency', 'fill': 90, 'field': 'mises'}
COMPARE_VIEWS = ['front']

# ------------------------------------------------- MESH-ON VARIANT RUN -------
# Set MESH_VARIANT = True to shoot a SECOND copy of the report figures with the
# element edges left ON, so the mesh-on look can be judged in the report rather
# than argued about.
#
# It writes to a DIFFERENT folder (abaqus_captures_mesh), so the finished
# mesh-off set is never touched, and it restricts itself to the cases that
# actually appear in the thesis, so it takes a couple of minutes rather than
# rerunning all 549.
#
# Drop the two sets side by side in Chapter 5, decide, and keep one.
MESH_VARIANT = False      # 27 Sep: contour figures use outline edges only.
                          # The full mesh is shown once, in the MESH figure
                          # (MESH_SHOT above), where it is the subject.
EDGE_STYLE = NONE         # 27 Sep, 21:40: edges OFF. At the re-sized walls the
                          # mesh is 3 to 4 times finer (up to 238,169 elements)
                          # and the edges turned the contours black. The mesh
                          # itself is still shot once, in MESH_SHOT.
# (previous setting, 27 Sep afternoon:)
# EDGE_STYLE = ALL        # 27 Sep decision: element edges ON for the report
                          # contours (the colours stay readable). NB: FEATURE
                          # gives the same result on this model - the free mesh
                          # has inconsistent normals, so every edge counts as a
                          # feature edge. Use NONE for edge-free contours.

if MESH_VARIANT:
    EDGE_STYLE = ALL

# (aircraft, load case, fill) triples used by figures in the thesis. Taken from
# the \includegraphics{abaqus/...} lines of Chapters 4 and 5, plus the A320 and
# A350 LC7 f10 frames so the fill effect can be shown for all three aircraft.
THESIS_CASES = []
# The SAME cases for every aircraft, so each observation in Chapter 5 can be
# shown for all three tanks side by side:
#   Baseline f10 / f50 / f90   thermal field and its dependence on fill
#   LC1 f10                    out-of-round deformation between the rings
#   LC7 f10 / f50 / f90        emergency landing, the governing case
for _ac in ('E190', 'A320', 'A350'):
    for _lc, _f in (('Baseline', 10), ('Baseline', 50), ('Baseline', 90),
                    ('LC1_Maneuver', 10),
                    ('LC7_Emergency', 10), ('LC7_Emergency', 50),
                    ('LC7_Emergency', 90)):
        THESIS_CASES.append((_ac, _lc, _f))
CASE_FILTER = THESIS_CASES if FIGURE_SET == 'THESIS' else []


# --------------------------------------------------------------- limits -------
# Per aircraft, because the three tanks have different wall thicknesses and
# therefore different stress ranges. A single window would render two of the
# three nearly flat. Values below span the far-field range from the production
# CSV with a little headroom, and deliberately exclude the ring singularity
# (E190 68.2, A320 55.0, A350 64.8 MPa) which does not converge and is discarded
# by the far-field convention.
#
NT11_LIMITS = {
    'E190': (20.3, 24.0),
    'A320': (20.3, 25.0),
    'A350': (20.3, 27.7),
}

# Von Mises needs a window per LOAD CASE, not just per aircraft. The first test
# run proved why: a 35-65 MPa window on the Baseline case rendered almost the
# whole tank below the floor, in the out-of-range colour, because Baseline sits
# at 38.7-39.0 MPa while LC7 reaches 61.6. One window cannot serve both.
#
# Within a load case the window is shared across all five fills - that is the
# comparison the fill sweep makes, and it must not move between frames.
# Across load cases the windows differ, but those are separate figures.
#
# Values in MPa. Since 27 Sep these are the BASE windows: the floor is used as
# given, and the top is raised at run time to the next 5 MPa above the largest
# far-field stress of that (aircraft, load case) in the production CSVs, so a
# window can never saturate after a re-run. NT11 is treated the same way. The
# windows actually used are printed to the capture log.
MISES_WINDOW_BASE = {
    ('E190', 'Baseline'):      (20, 40),
    ('E190', 'LC1_Maneuver'):  (20, 40),
    ('E190', 'LC3_Combined'):  (20, 45),
    ('E190', 'LC7_Emergency'): (20, 65),
    ('A320', 'Baseline'):      (15, 35),
    ('A320', 'LC1_Maneuver'):  (15, 35),
    ('A320', 'LC3_Combined'):  (15, 35),
    ('A320', 'LC7_Emergency'): (15, 50),
    ('A350', 'Baseline'):      (10, 30),
    ('A350', 'LC1_Maneuver'):  (10, 30),
    ('A350', 'LC3_Combined'):  (10, 30),
    ('A350', 'LC7_Emergency'): (15, 60),
}

# HFL peak varies with the heat-leak budget (E190 702.6 W, A320 919.3 W,
# A350 6583.2 W) and with wetted area, so it is not known in advance. The script
# auto-ranges HFL on the FIRST odb it opens for each aircraft - the f10 Baseline,
# which has the largest dry area and therefore the highest flux concentration at
# the waterline - then locks that value for every other frame of that aircraft.
# The value used is printed in the log. To override, put it here:
HFL_FIXED = {}          # e.g.  {'E190': (0.0, 750.0)}

NUM_INTERVALS = 12

# Filled in by main() from the BASE tables and the production CSVs.
MISES_WINDOW = dict(MISES_WINDOW_BASE)
NT11_USED = dict(NT11_LIMITS)

# ==============================================================================
#  END OF CONFIGURATION
# ==============================================================================

ROOT = os.getcwd()
LOG = []


def log(msg):
    """
    Print and append, then flush the whole log to disk.

    Written out on every call rather than once at the end, because a two-hour
    run that dies at ODB 90 must still leave a record of how far it got and
    what it was doing. Rewriting a few kilobytes 600 times costs nothing
    against the run time.
    """
    print(msg)
    sys.__stdout__.flush()
    LOG.append(msg)
    try:
        d = os.path.join(ROOT, OUT_ROOT)
        if not os.path.isdir(d):
            os.makedirs(d)
        fh = open(os.path.join(d, 'capture_log.txt'), 'w')
        fh.write('\n'.join(LOG))
        fh.close()
    except Exception:
        pass        # logging must never be what breaks the run


def open_odb(path):
    """
    Open an ODB, tolerating the signature differences between Abaqus releases.

    session.openOdb's first parameter is `name`, not `path`. Calling it as
    session.openOdb(path=..., readOnly=True) leaves `name` unfilled and fails
    with "not all required arguments specified; expected 1, got 0" on Abaqus
    2025. The positional form works, but the others are kept as fallbacks so
    this runs on older releases too.

    If the ODB is already open in the session it is reused rather than
    reopened, which would otherwise raise.
    """
    for key in session.odbs.keys():
        if os.path.normcase(key) == os.path.normcase(path):
            return session.odbs[key]

    attempts = (
        ('positional', lambda: session.openOdb(path, readOnly=True)),
        ('positional-only', lambda: session.openOdb(path)),
        ('visualization', lambda: visualization.openOdb(path=path,
                                                        readOnly=True)),
        ('name+path', lambda: session.openOdb(name=path, path=path,
                                              readOnly=True)),
    )
    last = None
    for label, fn in attempts:
        try:
            return fn()
        except Exception as e:
            last = '%s -> %s' % (label, e)
    raise Exception(last)


# Where each aircraft's CURRENT ODBs live, in order of preference.
#   E190 : E190\ (19 Sep) or the root (26 Sep re-run; bit-identical results).
#   A320, A350 : the ROOT ONLY. Their sub-folders hold the superseded
#                old-thickness ODBs and must never be read.
ODB_DIRS = {'E190': ['.'], 'A320': ['.'], 'A350': ['.']}   # 27 Sep: all current ODBs are in the root

MASTER_SCRIPT = os.path.join('python_scripts', 'E190_MASTER_thermo_structural.py')
try:    # repository layout: the master script sits next to this one
    MASTER_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'E190_MASTER_thermo_structural.py')
except NameError:
    pass
EXPECTED_T_FALLBACK = {'E190': 0.002747, 'A320': 0.002694, 'A350': 0.008884}
EXPECTED_T_NOW = dict(EXPECTED_T_FALLBACK)


def expected_thickness():
    """
    Wall thickness per aircraft, read from AIRCRAFT_DATA in the master script
    so the check follows the model rather than a number copied into this file.
    Falls back to the values above only if the master script cannot be parsed.
    """
    out = dict(EXPECTED_T_FALLBACK)
    try:
        txt = open(os.path.join(ROOT, MASTER_SCRIPT)).read()
        for ac, t in re.findall(
                r"'(E190|A320|A350)'\s*:\s*dict\([^)]*?t_wall\s*=\s*([0-9.eE+-]+)",
                txt, re.S):
            out[ac] = float(t)
        src = MASTER_SCRIPT
    except Exception as e:
        src = 'fallback values (%s)' % e
    # The FE sizing (RUN_MODE SIZING / FINAL) records the walls actually run
    # in sized_walls.csv; those take precedence over AIRCRAFT_DATA.
    sw = os.path.join(ROOT, 'sized_walls.csv')
    if os.path.exists(sw):
        for r in csv.DictReader(open(sw)):
            if r.get('aircraft') in out:
                out[r['aircraft']] = float(r['t_wall_mm']) / 1000.0
        src = 'sized_walls.csv'
    log('  wall thickness check against %s: %s' % (
        src, ', '.join('%s %.3f mm' % (k, out[k] * 1e3) for k in sorted(out))))
    return out


def odb_name(aircraft, kind, loadcase, fill):
    """
    Path of the ODB for one case. E190 files carry no aircraft tag; A320 and
    A350 do:  HT_Baseline_f10.odb,  HT_A320_Baseline_f10.odb
    Returns the first existing candidate in ODB_DIRS, or the preferred path if
    none exists (the caller reports it as missing).
    """
    tag = '' if aircraft == 'E190' else (aircraft + '_')
    base = '%s_%s%s_f%02d.odb' % (kind, tag, loadcase, fill)
    cands = [os.path.normpath(os.path.join(ROOT, d, base))
             for d in ODB_DIRS[aircraft]]
    for c in cands:
        if os.path.exists(c):
            return c
    return cands[0]


def odb_thickness(odb):
    """Shell thickness stored in the ODB's section definition, in metres."""
    try:
        for s in odb.sections.values():
            t = getattr(s, 'thickness', None)
            if t:
                return float(t)
    except Exception:
        pass
    return None


def thickness_ok(odb, aircraft, path, expected):
    """
    Refuse any ODB whose wall thickness is not the current one. This is the
    guard against plotting a superseded model: the image would look entirely
    plausible and be wrong.
    """
    t = odb_thickness(odb)
    if t is None:
        log('  WARNING %s: thickness not readable from the ODB; accepted on '
            'location only' % os.path.basename(path))
        return True
    if abs(t - expected[aircraft]) > 1e-7:
        log('  REJECTED %s: wall %.3f mm, current model is %.3f mm'
            % (os.path.basename(path), t * 1e3, expected[aircraft] * 1e3))
        return False
    return True


def read_results():
    """
    Far-field results per (aircraft, load case, fill) from the CURRENT CSVs.

    E190   : E190_production*_results.csv in the root (re-run at the corrected wall).
    A320/50: <AC>_production*_results.csv in the ROOT only - this picks up the
             full run and the single-case re-runs. The A320/A350 rows inside
             ALL_AIRCRAFT and the files in the sub-folders are the superseded
             old-thickness results and are deliberately not read.
    Later files override earlier ones for the same case.
    """
    rows = {}

    def take(ac, path):
        try:
            for r in csv.DictReader(open(path)):
                if r.get('aircraft', ac) != ac:
                    continue
                key = (ac, r['load_case'], int(round(float(r['fill']) * 100)))
                tw = r.get('t_wall_mm')
                if not tw or abs(float(tw) - EXPECTED_T_NOW[ac] * 1e3) > 1e-4:
                    continue        # superseded wall, or a file from before 27 Sep
                rows[key] = r
        except Exception as e:
            log('  NOTE: could not read %s (%s)' % (os.path.basename(path), e))

    for ac in ('E190', 'A320', 'A350'):
        files = glob.glob(os.path.join(ROOT, '%s_production*_results.csv' % ac))
        for p in sorted(files, key=os.path.getmtime):
            take(ac, p)
    return rows


def contour_windows(results):
    """
    Final contour limits. Floor from the BASE tables; top raised to the next
    5 MPa (von Mises) or 0.1 K (NT11) above the largest far-field value of the
    case group, so the plotted range always contains the tabulated result.
    """
    mises = {}
    for (ac, lc), (lo, hi) in MISES_WINDOW_BASE.items():
        vals = [float(r['mises_MPa']) for (a, l, f), r in results.items()
                if a == ac and l == lc]
        top = 5.0 * math.ceil(max(vals) / 5.0) if vals else hi
        mises[(ac, lc)] = (lo, max(hi, top))
        log('  window %s %-14s von Mises %5.1f - %5.1f MPa  (far-field max %s)'
            % (ac, lc, lo, mises[(ac, lc)][1],
               ('%.2f' % max(vals)) if vals else 'n/a'))
    nt11 = {}
    for ac, (lo, hi) in NT11_LIMITS.items():
        vals = [float(r['T_max']) for (a, l, f), r in results.items() if a == ac]
        top = math.ceil(max(vals) * 10.0) / 10.0 if vals else hi
        nt11[ac] = (lo, max(hi, top))
        log('  window %s NT11 %.1f - %.1f K' % (ac, lo, nt11[ac][1]))
    return mises, nt11


# ------------------------------------------------------------------------------
#  Viewport, once for the whole run
# ------------------------------------------------------------------------------
def make_viewport():
    vp = session.Viewport(name='Capture', origin=(0.0, 0.0),
                          width=320, height=130)
    vp.makeCurrent()
    # NOT maximize(). A maximised viewport ignores later setValues(width=,
    # height=), so the per-view aspect ratio set in set_canvas() was being
    # silently discarded: the viewport stayed near-square, fitView fitted the
    # tank to a square, and the render into a wide PNG then letterboxed the
    # front views and clipped the iso ones. Size is set per view instead.

    session.graphicsOptions.setValues(backgroundStyle=SOLID,
                                      backgroundColor='#FFFFFF')

    # Legend stays: it carries the scale, and the figure is unreadable without
    # it. Everything else is caption material and only wastes frame.
    vp.viewportAnnotationOptions.setValues(
        triad=OFF, title=OFF, state=OFF, compass=OFF, legend=ON,
        legendBackgroundStyle=MATCH, legendBox=ON)
    # An OPAQUE white legend background, so that where the legend sits over
    # the model in the legend-on render it hides the model rather than letting
    # it show through. compose_side_legend.py relies on that to cut a clean
    # legend. MATCH is already white and opaque on Abaqus 2025; this makes it
    # explicit where the release supports it.
    try:
        vp.viewportAnnotationOptions.setValues(
            legendBackgroundStyle=OTHER, legendBackgroundColor='#FFFFFF')
    except Exception:
        pass

    session.printOptions.setValues(vpDecorations=OFF, vpBackground=OFF,
                                   compass=OFF, reduceColors=False)
    return vp


def style_contours(vp, lo, hi):
    """Locked limits and render style. Edges are applied later, see note."""
    vp.odbDisplay.contourOptions.setValues(
        minAutoCompute=OFF, minValue=lo,
        maxAutoCompute=OFF, maxValue=hi,
        numIntervals=NUM_INTERVALS)
    vp.odbDisplay.commonOptions.setValues(renderStyle=CONTOUR_RENDER)


def apply_edges(vp):
    """
    Set the edge style immediately before printing.

    The first test run came back with the mesh visible despite visibleEdges
    being set in style_contours(), so it is reapplied here, after the field,
    the view and the cut are all in place. Cheap insurance: whatever reset it,
    this runs last.
    """
    try:
        vp.odbDisplay.commonOptions.setValues(visibleEdges=EDGE_STYLE)
    except Exception as e:
        log('    NOTE: could not set edge style (%s)' % e)


def set_section_results(vp, mode):
    """
    Choose which shell section point the contour reports.

    This is VIEWPORT state and it PERSISTS between fields, which caused a real
    failure: USE_ENVELOPE was set once for von Mises and then leaked into the
    next case's HFL frames. NT11 survived it, being a nodal field with no
    section points, but HFL is an integration-point vector and came back drawn
    in flat white - no values plotted at all.

    So it is set explicitly for every field rather than assumed.
    """
    try:
        vp.odbDisplay.basicOptions.setValues(sectionResults=mode)
        return True
    except Exception as e:
        log('    NOTE: could not set section results (%s)' % e)
        return False


def set_field(vp, field):
    """Select the primary variable. Returns True on success."""
    try:
        if field == 'nt11':
            # Nodal, no section points involved.
            set_section_results(vp, USE_BOTTOM)
            vp.odbDisplay.setPrimaryVariable(
                variableLabel='NT11', outputPosition=NODAL)
        elif field == 'hfl':
            # Must NOT be the envelope, see set_section_results above.
            set_section_results(vp, USE_BOTTOM)
            vp.odbDisplay.setPrimaryVariable(
                variableLabel='HFL', outputPosition=INTEGRATION_POINT,
                refinement=(INVARIANT, 'Magnitude'))
        elif field == 'mises':
            # The master script takes its max over every value in the S field,
            # both section points included, so the tabulated stress is an
            # envelope. Plot the envelope too, or the figures would report a
            # different quantity from the tables.
            set_section_results(vp, USE_ENVELOPE)
            vp.odbDisplay.setPrimaryVariable(
                variableLabel='S', outputPosition=INTEGRATION_POINT,
                refinement=(INVARIANT, 'Mises'))
        return True
    except Exception as e:
        log('    SKIP field %s: %s' % (field, e))
        return False


def hfl_peak(odb):
    """
    Largest heat-flux magnitude in the last frame, read straight from the ODB.

    The viewport's autoMinValue/autoMaxValue were used for this at first and
    returned 3.4e38 - FLT_MAX, the sentinel Abaqus reports when the auto range
    has not been recomputed since the limits were last set explicitly. It
    happened to work on the very first run and broke as soon as another
    capture ran ahead of it and left specified limits in place.

    Reading the field data is deterministic and depends on no viewport state.
    """
    peak = 0.0
    for step in odb.steps.values():
        try:
            fo = step.frames[-1].fieldOutputs['HFL']
        except Exception:
            continue
        for v in fo.values:
            try:
                m = v.magnitude
            except Exception:
                d = v.data
                m = (d[0] * d[0] + d[1] * d[1]) ** 0.5
            if m > peak:
                peak = m
    return peak


def set_view(vp, spec):
    """Saved user view if it exists, else the viewpoint vector."""
    name = spec.get('saved')
    if name and name in session.views.keys():
        vp.view.setValues(session.views[name])
    else:
        vp.view.setViewpoint(viewVector=spec['vector'],
                             cameraUpVector=spec['up'])
    vp.view.setProjection(projection=PARALLEL)     # perspective distorts a
    vp.view.fitView()                              # 15 m tank badly
    try:
        vp.view.zoom(spec.get('zoom', 1.0))        # keep off the frame edge
    except Exception:
        pass


def set_canvas(vp, spec):
    """
    Image size, viewport aspect and legend placement, all per view.

    MUST be called BEFORE set_view(), because fitView fits to the VIEWPORT's
    aspect ratio, not the PNG's. The test run had a viewport of 320x130
    (2.46:1) rendered into a 2400x760 PNG (3.16:1), so Abaqus letterboxed the
    result and the tank came out small with a band of white above and below.
    Matching the viewport shape to the image shape makes fitView fill the frame.
    """
    w, h = spec.get('size', (2400, 900))
    try:
        vp.restore()            # a maximised viewport ignores width/height
    except Exception:
        pass
    try:
        vp.setValues(width=320.0, height=320.0 * float(h) / float(w))
    except Exception as e:
        log('    NOTE: could not resize viewport (%s)' % e)

    # Derive the image size from what the viewport ACTUALLY became rather than
    # from what was asked for. If Abaqus clamps or adjusts the viewport, an
    # imageSize computed independently would no longer match its aspect ratio,
    # and the render would letterbox or crop again. Reading it back makes the
    # two agree by construction.
    try:
        aspect = float(vp.height) / float(vp.width)
        session.pngOptions.setValues(imageSize=(w, max(2, int(round(w * aspect)))))
    except Exception:
        session.pngOptions.setValues(imageSize=(w, h))
    try:
        vp.viewportAnnotationOptions.setValues(
            legendPosition=spec.get('legend', (1, 98)))
    except Exception as e:
        log('    NOTE: could not move legend (%s)' % e)


def set_cut(vp, on):
    """Z-plane section cut, keeping the half below the plane."""
    try:
        if on:
            vp.odbDisplay.setValues(viewCutNames=('Z-Plane', ), viewCut=ON)
            vp.odbDisplay.viewCuts['Z-Plane'].setValues(
                showModelBelowCut=True, showModelOnCut=True,
                showModelAboveCut=False)
        else:
            vp.odbDisplay.setValues(viewCut=OFF)
        return True
    except Exception as e:
        log('    NOTE: view cut failed (%s)' % e)
        return False


# ------------------------------------------------------------------------------
#  Main
# ------------------------------------------------------------------------------
def capture_mesh(vp, aircraft_list):
    """
    The Chapter 4 mesh figure, one set per aircraft.

    Runs BEFORE the contour sweep so that aborting the long run still leaves
    these on disk. Element counts are logged because the caption needs them
    and reading them off the model beats trusting a number typed months ago.
    """
    n = 0
    for ac in aircraft_list:
        out_dir = os.path.join(ROOT, OUT_ROOT, ac)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        path = odb_name(ac, 'ST', MESH_SOURCE['loadcase'], MESH_SOURCE['fill'])
        if not os.path.exists(path):
            log('  MESH: no odb at %s' % os.path.basename(path))
            continue
        try:
            odb = open_odb(path)
        except Exception as e:
            log('  MESH: could not open %s (%s)' % (os.path.basename(path), e))
            continue

        vp.setValues(displayedObject=odb)
        # Undeformed shape, not a contour: the subject is the discretisation.
        vp.odbDisplay.display.setValues(plotState=(UNDEFORMED, ))
        vp.viewportAnnotationOptions.setValues(legend=OFF)

        try:
            inst = odb.rootAssembly.instances
            n_el = 0
            n_nd = 0
            for k in inst.keys():
                n_el = n_el + len(inst[k].elements)
                n_nd = n_nd + len(inst[k].nodes)
            log('  MESH %s: %d elements, %d nodes' % (ac, n_el, n_nd))
        except Exception as e:
            log('  MESH %s: could not read element count (%s)' % (ac, e))

        for vname in MESH_VIEWS:
            spec = VIEWS[vname]
            set_cut(vp, spec['cut'])
            set_canvas(vp, spec)
            set_view(vp, spec)
            # Every element edge, which is the whole point of this one figure.
            vp.odbDisplay.commonOptions.setValues(visibleEdges=ALL,
                                                  renderStyle=FILLED)
            fname = '%s_MESH_%s' % (ac, vname)
            try:
                session.printToFile(fileName=os.path.join(out_dir, fname),
                                    format=PNG, canvasObjects=(vp, ))
                n += 1
            except Exception as e:
                log('  MESH FAILED %s: %s' % (fname, e))

        set_cut(vp, False)
        vp.viewportAnnotationOptions.setValues(legend=ON)
        odb.close()
    return n


def capture_compare(vp, aircraft_list):
    """
    The same contour frame with and without element edges.

    The governing case is used, because that is the frame with the most
    structure in it and therefore the one where the mesh does the most damage
    to legibility. Everything else about the two images is identical, so the
    comparison isolates the edges alone.
    """
    lc = COMPARE_SOURCE['loadcase']
    fill = COMPARE_SOURCE['fill']
    field = COMPARE_SOURCE['field']
    kind = 'ST' if field == 'mises' else 'HT'
    n = 0

    for ac in aircraft_list:
        out_dir = os.path.join(ROOT, OUT_ROOT, ac)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        path = odb_name(ac, kind, lc, fill)
        if not os.path.exists(path):
            log('  COMPARE: no odb at %s' % os.path.basename(path))
            continue
        try:
            odb = open_odb(path)
        except Exception as e:
            log('  COMPARE: could not open %s (%s)'
                % (os.path.basename(path), e))
            continue

        vp.setValues(displayedObject=odb)
        vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF, ))
        if not set_field(vp, field):
            odb.close()
            continue

        if field == 'mises':
            lo, hi = MISES_WINDOW[(ac, lc)]
            lo, hi = lo * 1.0e6, hi * 1.0e6
        elif field == 'nt11':
            lo, hi = NT11_USED[ac]
        else:
            peak = hfl_peak(odb)
            lo, hi = 0.0, 50.0 * (int(peak / 50.0) + 1)
        style_contours(vp, lo, hi)

        for vname in COMPARE_VIEWS:
            spec = VIEWS[vname]
            set_cut(vp, spec['cut'])
            set_canvas(vp, spec)
            set_view(vp, spec)

            for label, edges in (('mesh', ALL), ('nomesh', NONE)):
                try:
                    vp.odbDisplay.commonOptions.setValues(visibleEdges=edges)
                    fname = '%s_COMPARE_f%02d_%s_%s_%s' % (
                        ac, fill, field, vname, label)
                    session.printToFile(
                        fileName=os.path.join(out_dir, fname),
                        format=PNG, canvasObjects=(vp, ))
                    n += 1
                except Exception as e:
                    log('  COMPARE FAILED %s/%s: %s' % (vname, label, e))

        set_cut(vp, False)
        odb.close()
    return n


def capture_zoom(vp, aircraft_list):
    """
    Zoomed frames of the governing case, mesh on and mesh off (ZOOM_SHOT).
    Front views go through print_frame, so they get the side legend like
    every other front view; iso views carry no legend.
    """
    n = 0
    for ac in aircraft_list:
        out_dir = os.path.join(ROOT, OUT_ROOT, ac)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        for lc, fill, field in ZOOM_CASES:
            kind = 'ST' if field == 'mises' else 'HT'
            path = odb_name(ac, kind, lc, fill)
            if not os.path.exists(path):
                log('  ZOOM: no odb at %s' % os.path.basename(path))
                continue
            try:
                odb = open_odb(path)
            except Exception as e:
                log('  ZOOM: could not open %s (%s)' % (os.path.basename(path), e))
                continue
            vp.setValues(displayedObject=odb)
            vp.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF, ))
            if not set_field(vp, field):
                odb.close()
                continue
            if field == 'mises':
                lo, hi = MISES_WINDOW[(ac, lc)]
                lo, hi = lo * 1.0e6, hi * 1.0e6
            else:
                lo, hi = NT11_USED[ac]
            style_contours(vp, lo, hi)
            for vname in ZOOM_VIEWS:
                spec = dict(VIEWS[vname])
                spec['zoom'] = spec.get('zoom', 1.0) * ZOOM_FACTOR
                set_cut(vp, False)
                set_canvas(vp, spec)
                set_view(vp, spec)
                for label, edges in (('mesh', ALL), ('nomesh', NONE)):
                    try:
                        vp.odbDisplay.commonOptions.setValues(visibleEdges=edges)
                        fname = '%s_ZOOM_%s_f%02d_%s_%s_%s' % (ac, lc, fill, field, vname, label)
                        n += print_frame(vp, out_dir, fname, vname)
                    except Exception as e:
                        log('  ZOOM FAILED %s %s %s: %s' % (ac, vname, label, e))
            odb.close()
            log('  zoom  %s' % os.path.basename(path))
    return n


def print_frame(vp, out_dir, fname, vname):
    """
    Write one frame. For a view in SIDE_LEGEND_VIEWS write two raw renders,
    identical except for the legend, for compose_side_legend.py to assemble.
    Returns the number of files written.
    """
    if vname not in SIDE_LEGEND_VIEWS:
        # No legend on this view: it is read against the front view.
        vp.viewportAnnotationOptions.setValues(legend=OFF)
        try:
            session.printToFile(fileName=os.path.join(out_dir, fname),
                                format=PNG, canvasObjects=(vp, ))
        finally:
            vp.viewportAnnotationOptions.setValues(legend=ON)
        return 1
    raw_dir = os.path.join(ROOT, OUT_ROOT, '_raw', os.path.basename(out_dir))
    if not os.path.isdir(raw_dir):
        os.makedirs(raw_dir)
    try:
        vp.viewportAnnotationOptions.setValues(legend=ON)
        try:
            # Opaque white legend box: in the zoomed frames the tank fills the
            # whole image, so a transparent legend would be lifted out of
            # compose_side_legend.py with tank colours behind its text.
            vp.viewportAnnotationOptions.setValues(
                legendBox=ON, legendBackgroundStyle=OTHER,
                legendBackgroundColor='#FFFFFF')
        except Exception:
            pass
        session.printToFile(fileName=os.path.join(raw_dir, fname + '__L'),
                            format=PNG, canvasObjects=(vp, ))
        vp.viewportAnnotationOptions.setValues(legend=OFF)
        session.printToFile(fileName=os.path.join(raw_dir, fname + '__N'),
                            format=PNG, canvasObjects=(vp, ))
    finally:
        vp.viewportAnnotationOptions.setValues(legend=ON)
    return 2


def main():
    global MISES_WINDOW, NT11_USED
    vp = make_viewport()
    hfl_limits = dict(HFL_FIXED)
    n_written, n_missing, n_failed, n_rejected = 0, 0, 0, 0

    log('')
    log('--- checks and contour windows ---')
    expected = expected_thickness()
    global EXPECTED_T_NOW
    EXPECTED_T_NOW = expected
    results = read_results()
    for ac in AIRCRAFT:
        n = len([k for k in results if k[0] == ac])
        log('  results available for %s: %d of 20 cases%s'
            % (ac, n, '' if n == 20 else '   <-- INCOMPLETE'))
    MISES_WINDOW, NT11_USED = contour_windows(results)

    aircraft_list = AIRCRAFT[:1] if QUICK_TEST else AIRCRAFT
    lc_list = LOADCASES[:1] if QUICK_TEST else LOADCASES
    fill_list = FILLS[2:3] if QUICK_TEST else FILLS

    if MESH_SHOT:
        log('')
        log('--- mesh figure (Chapter 4) ---')
        n_written += capture_mesh(vp, aircraft_list)

    if ZOOM_SHOT:
        log('')
        log('--- zoomed frames, mesh on and off ---')
        n_written += capture_zoom(vp, aircraft_list)

    if COMPARE_SHOT:
        log('')
        log('--- mesh on/off comparison ---')
        n_written += capture_compare(vp, aircraft_list)

    for ac in aircraft_list:
        out_dir = os.path.join(ROOT, OUT_ROOT, ac)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        log('')
        log('=' * 70)
        log('  %s   ->  %s' % (ac, out_dir))
        log('=' * 70)

        for lc in lc_list:
            for fill in fill_list:

                for kind, fields in (('HT', ['nt11', 'hfl']),
                                     ('ST', ['mises'])):

                    if (kind == 'HT' and SKIP_REDUNDANT_THERMAL
                            and lc == 'LC1_Maneuver'):
                        continue

                    if CASE_FILTER and (ac, lc, fill) not in CASE_FILTER:
                        continue

                    path = odb_name(ac, kind, lc, fill)
                    if not os.path.exists(path):
                        log('  MISSING  %s' % os.path.basename(path))
                        n_missing += 1
                        continue

                    try:
                        odb = open_odb(path)
                    except Exception as e:
                        log('  FAILED to open %s: %s'
                            % (os.path.basename(path), e))
                        n_failed += 1
                        continue

                    if not thickness_ok(odb, ac, path, expected):
                        odb.close()
                        n_rejected += 1
                        continue

                    vp.setValues(displayedObject=odb)
                    vp.odbDisplay.display.setValues(
                        plotState=(CONTOURS_ON_UNDEF, ))

                    for field in fields:
                        if not set_field(vp, field):
                            continue

                        # ---- limits for this field -----------------------
                        if field == 'mises':
                            lo, hi = MISES_WINDOW[(ac, lc)]
                            lo, hi = lo * 1.0e6, hi * 1.0e6   # MPa -> Pa
                        elif field == 'nt11':
                            lo, hi = NT11_USED[ac]
                        else:
                            if ac not in hfl_limits:
                                amax = hfl_peak(odb)
                                hi = 50.0 * (int(amax / 50.0) + 1)
                                hfl_limits[ac] = (0.0, hi)
                                log('  HFL range for %s locked at 0 - %.0f '
                                    'W/m2 (field peak %.1f on this odb)'
                                    % (ac, hi, amax))
                            lo, hi = hfl_limits[ac]

                        style_contours(vp, lo, hi)

                        for vname in CONTOUR_VIEWS:
                            spec = VIEWS[vname]
                            set_cut(vp, spec['cut'])
                            set_canvas(vp, spec)
                            set_view(vp, spec)
                            apply_edges(vp)

                            fname = '%s_%s_f%02d_%s_%s' % (
                                ac, lc, fill, field, vname)
                            try:
                                n_written += print_frame(vp, out_dir,
                                                         fname, vname)
                            except Exception as e:
                                log('  FAILED %s: %s' % (fname, e))
                                n_failed += 1

                        set_cut(vp, False)

                    odb.close()
                    log('  done  %s' % os.path.basename(path))

    log('')
    log('=' * 70)
    log('  written %d   missing odb %d   failed %d   rejected (wrong '
        'thickness) %d' % (n_written, n_missing, n_failed, n_rejected))
    if SIDE_LEGEND_VIEWS:
        log('')
        log('  NEXT: the %s frames are still raw. Assemble them with' %
            '/'.join(SIDE_LEGEND_VIEWS))
        log('      python python_scripts\\compose_side_legend.py')
    log('  output root: %s' % os.path.join(ROOT, OUT_ROOT))
    if QUICK_TEST:
        log('')
        log('  QUICK_TEST was True - this was a single ODB.')
        log('  Check the images, then set QUICK_TEST = False and rerun.')
    log('=' * 70)

    try:
        fh = open(os.path.join(ROOT, OUT_ROOT, 'capture_log.txt'), 'w')
        fh.write('\n'.join(LOG))
        fh.close()
    except Exception:
        pass


main()
