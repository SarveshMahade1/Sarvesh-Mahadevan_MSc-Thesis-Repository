# -*- coding: utf-8 -*-
"""
REPRODUCE EVERYTHING THAT DOES NOT NEED ABAQUS, IN ORDER

Runs the plain-Python chain of the workflow from the design basis and the
stored finite element results to every processed number, document, appendix
table and thesis figure:

    1  tank geometry from the mission          (design basis -> R, L, V)
    2  heat leak budget from the hold time     (needs CoolProp)
    3  insulation and outer jacket
    4  closed-form starting walls
    5  merge of the finite element results     -> results/processed/final_numbers.json
    6  closed-form checks and review checks
    7  documents and appendix tables
    8  thesis figures                          (contours need Tesseract OCR)

Run from anywhere:
    python code/reproduce_all.py               everything
    python code/reproduce_all.py --no-figures  numbers and documents only

On the files of this repository the outputs are identical to the committed
ones (git status shows no change, apart from the last floating-point digit of
one material property on some platforms). With a modified design basis or new
Abaqus results, the same command regenerates the whole chain.
"""
import os
import subprocess
import sys
import time

CODE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(CODE)

STEPS = [
    ('2_requirements_and_sizing/tank_geometry_from_mission.py', []),
    ('2_requirements_and_sizing/boiloff_requirement_derivation.py', []),
    ('2_requirements_and_sizing/insulation_sizing.py', []),
    ('3_results_processing/final_numbers.py', []),
    ('3_results_processing/hand_calcs.py', []),
    ('2_requirements_and_sizing/outer_jacket_closed_form.py', []),
    ('2_requirements_and_sizing/sizing_closed_form.py', []),
    ('2_requirements_and_sizing/bagarello_benchmark_1D_check.py', []),
    ('3_results_processing/review_checks.py', []),
    ('3_results_processing/make_appendix_traceability.py', []),
    ('3_results_processing/make_review_docs.py', []),
]
FIGURES = [
    ('2_requirements_and_sizing/boiloff_holdtime_sensitivity.py', []),
    ('4_figures_and_drawings/plot_final_figures.py', []),
    ('4_figures_and_drawings/draw_workflow.py', ['thesis_figures']),
    ('4_figures_and_drawings/draw_report_diagrams.py', []),
    ('4_figures_and_drawings/make_report_contours.py', ['abaqus_contour_captures_report', 'thesis_figures']),
    ('3_results_processing/verify_figure_annotations.py', []),
]


def main():
    steps = STEPS + ([] if '--no-figures' in sys.argv else FIGURES)
    t0 = time.time()
    for i, (script, args) in enumerate(steps, 1):
        print('[%2d/%d] %s %s' % (i, len(steps), script, ' '.join(args)))
        r = subprocess.run([sys.executable, os.path.join(CODE, script)] + args, cwd=REPO,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        if r.returncode != 0:
            print(r.stdout[-3000:])
            sys.exit('stopped: %s failed' % script)
    print('done in %.0f s; compare with the committed files using  git status  or  git diff' % (time.time() - t0))


if __name__ == '__main__':
    main()
