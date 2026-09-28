# -*- coding: utf-8 -*-
"""
Folder layout of this repository, in one place. Every plain-Python script
imports this module, so the folders can be renamed here without editing
the scripts.
"""
import os
import sys

CODE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(CODE)

RESULTS = os.path.join(REPO, 'results')
SIZED_WALLS = os.path.join(RESULTS, 'sized_walls.csv')
PRODUCTION = os.path.join(RESULTS, 'production_1p5bar')        # <AC>_production_*_results.csv
STABILITY = os.path.join(RESULTS, 'stability_0p187bar')        # <AC>_stability_results.csv
CONVERGENCE = os.path.join(RESULTS, 'mesh_convergence')        # <AC>_convergence_*_results.csv
SIZING = os.path.join(RESULTS, 'wall_sizing_iterations')       # <AC>/<AC>_sizing_results.csv
BENCHMARK = os.path.join(RESULTS, 'benchmark_bagarello')
DELFTBLUE = os.path.join(RESULTS, 'solver_check_delftblue')
PROCESSED = os.path.join(RESULTS, 'processed')                 # final_numbers.json, hand_calcs.json, insulation.json

CAPTURES = os.path.join(REPO, 'abaqus_contour_captures')
FIGURES = os.path.join(REPO, 'thesis_figures')
APPENDIX = os.path.join(REPO, 'thesis_appendix')
DRAWINGS = os.path.join(REPO, 'design_drawings')

# the scripts import one another across the code sub-folders
for _sub in ('2_requirements_and_sizing', '3_results_processing', '4_figures_and_drawings'):
    _p = os.path.join(CODE, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)
