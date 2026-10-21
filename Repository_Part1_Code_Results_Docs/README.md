# Digitalized Workflow of a Liquid Hydrogen Tank Integrating Thermo-Structural Analysis

Model, scripts and complete results of the MSc thesis
*Digitalized Workflow of a Liquid Hydrogen Tank Integrating Thermo-Structural Analysis*,
Sarveshwaran Mahadevan, MSc Sustainable Energy Technology, Delft University of Technology, 2026.
Supervisors: Dr. Francesco Orefice and Dr. Othon Moultos.

One parametric Abaqus script builds, meshes, sizes, solves and post-processes a liquid hydrogen
tank for three aircraft classes: the Embraer E190, the Airbus A320neo and the Airbus A350-900.
Each tank is sized to carry the hydrogen that delivers the same energy as the reference kerosene
load; the aircraft define the range, not an installation. Its aluminium 2219-T87 inner wall is
sized so that its reserve factors are at least 1 against yield at 1.5 bar and against a NASA SP-8007
buckling screen at 0.187 bar, with a project design margin of 1.5 inside both allowables
(RF = allowable / applied). It is then analysed under four static load cases (one gravity, the
CS 25.337 manoeuvre, a non-regulatory combined case and the CS 25.561 emergency landing factors
applied together) at five fill levels, at both pressures:

- 120 production and stability cases, each a steady-state heat transfer analysis followed by a
  linear static analysis;
- 18 mesh-convergence cases.

Every number in the thesis can be traced to a file in this repository and regenerated from it.

## Load-case names

The thesis names the four load cases in words. The scripts and result files use short codes,
which appear only inside file names and column values:

| Code in file names | Name in the thesis | Basis | n_x, n_y, n_z |
|---|---|---|---|
| `Baseline` | Baseline condition | 1 g, ground or level flight | 0, 0, 1.0 |
| `LC1_Maneuver` | Symmetric manoeuvre | CS 25.337, limit | 0, 0, 2.5 |
| `LC3_Combined` | Combined manoeuvre | assumed, not CS-25 | 1.5, 0, 2.5 |
| `LC7_Emergency` | Emergency landing | CS 25.561 factors combined, ultimate | 9.0, 3.0, 6.0 |

The figure files `fig5_case_lc3` and `fig5_case_lc7` carry the same codes in their names.

## Main results

| | E190 | A320 | A350-900 |
|---|---|---|---|
| Barrel slenderness L/2R | 5.37 | 2.93 | 6.56 |
| Sized wall [mm] | 1.289 | 1.572 | 5.674 |
| Wall governed by | yield | yield | buckling |
| Governing far-field von Mises, emergency landing, fill 0.90 [MPa] | 256.9 | 257.3 | 145.7 |
| Lowest yield RF, production mesh (finest mesh) | 1.020 (1.020) | 1.018 (1.016) | 1.798 (1.777) |
| Lowest buckling RF, production mesh (finest mesh) | 1.105 (1.113) | no compression | 1.060 (1.021) |
| Fill sensitivity under the emergency landing [%] | +46.4 | +36.0 | +100.9 |
| Wetted-to-dry temperature difference, fill 0.10 [K] | 3.683 | 4.644 | 7.364 |

The material properties are room-temperature values. Quantities that depend on the expansion
coefficient or the conductivity are corrected to cryogenic values by scaling in
`docs/detailed_calculations.md`; the reversal of the mid-ring reaction seen in the model does not
survive that correction. Two points are open: the location of the far-field peak of the E190 and
A320 (about 20 % above the membrane stress, pressure-proportional) and the A350-900 buckling margin
with cryogenic properties.

## What is where

```
code/
  1_abaqus_model/              the finite element model (Abaqus/CAE Python)
  2_requirements_and_sizing/   heat leak budget, insulation, outer jacket, closed-form walls
  3_results_processing/        merges the results, closed-form checks, appendix tables
  4_figures_and_drawings/      every graph, contour panel and drawing of the thesis
  repo_paths.py                the folder layout, in one place

results/
  sized_walls.csv              wall of each aircraft from the FE sizing
  production_1p5bar/           production matrix at 1.5 bar, 20 cases per aircraft
  stability_0p187bar/          the same matrix at 0.187 bar
  mesh_convergence/            governing case (emergency landing, fill 0.90), three meshes, both pressures
  wall_sizing_iterations/      every iteration of the FE wall sizing, per aircraft
  benchmark_bagarello/         the reference tank of Bagarello et al. run through the method
  solver_check_delftblue/      one case solved on the DelftBlue cluster (Abaqus 2024)
  processed/                   final_numbers.json, hand_calcs.json, insulation.json,
                               production_matrix_sized_walls.csv

abaqus_contour_captures_report/  the Abaqus captures (legend and tank in one image) from which
                               make_report_contours.py builds every contour and mesh figure of the thesis
abaqus_contour_captures/       earlier draft only: a partial set of 15 side-legend captures, kept for the record
docs/                          full_results_matrices.md (all 120 cases with reserve factors) and
                               detailed_calculations.md (the checks behind the final report)
thesis_figures/                the 29 figures of the final report (PDF as used in the thesis, PNG preview)
thesis_appendix/               appendix_calculation_traceability.tex, the condensed appendix as printed in
                               the thesis, and appendix_calculation_traceability_full.tex, the complete
                               tables written by make_appendix_traceability.py
```

Two folders are created only when the corresponding scripts are run: `design_drawings/` (supporting
drawings that are not figures of the thesis) and `earlier_draft_figures/` (output of the earlier-draft
contour scripts).

In every results file name, `<AC>` is `E190`, `A320` or `A350`. The Abaqus output databases
(.odb) are not included, because they amount to tens of gigabytes; running the model regenerates
them.

## The code

### 1_abaqus_model (needs an Abaqus licence)

| Script | Purpose |
|---|---|
| `E190_MASTER_thermo_structural.py` | The digitalized workflow. Despite the name it serves all three aircraft. Geometry, mesh (seed 0.5 sqrt(Rt)), free surface normal to the load factor, steady-state heat transfer, banded temperature transfer, static solution with pressure, hydrostatic head and inertia, three ring supports, far-field extraction, verification gates, FE wall sizing. |
| `ABAQUS_capture_all.py` | Captures the contour images from the output databases. |
| `run_final_delftblue.sh` | SLURM job script used on the DelftBlue cluster. |

Run the model from an empty working folder; the .odb and results files are written there. Copy
`results/sized_walls.csv` into that folder first, so that every run uses the sized walls. The model
is controlled by environment variables, so nothing in the file needs editing:

| Variable | Values |
|---|---|
| `E190_AIRCRAFT` | `E190`, `A320`, `A350` |
| `E190_RUN_MODE` | `FINAL` (sizing, production, stability, convergence), `SIZING`, `PRODUCTION`, `STABILITY`, `CONVERGENCE` |
| `E190_STAGES` | subset of `sizing,production,stability,convergence` for a `FINAL` run |
| `E190_ONLY_LC`, `E190_ONLY_FILL` | run part of the matrix, for example `LC7_Emergency` and `0.9` |
| `E190_T_WALL`, `E190_P_GAUGE` | override the wall [mm] or the pressure [bar] |
| `E190_BENCHMARK=1` | run the Bagarello et al. reference tank instead |

Example, one aircraft, production and stability at the sized wall (Windows PowerShell):

```
cd <empty working folder>
copy <repository>\results\sized_walls.csv .
$env:E190_AIRCRAFT = "A320"
$env:E190_RUN_MODE = "FINAL"
$env:E190_STAGES   = "production,stability"
abaqus cae noGUI=<repository>\code\1_abaqus_model\E190_MASTER_thermo_structural.py
```

### 2_requirements_and_sizing (plain Python)

| Script | Purpose |
|---|---|
| `boiloff_requirement_derivation.py` | Heat leak budget from the no-vent hold (constant-volume energy balance). |
| `boiloff_holdtime_sensitivity.py` | Hold time swept from 2 h to 36 h, and its effect on the insulation (uses CoolProp). |
| `insulation_sizing.py` | Foam, vacuum and vacuum-with-MLI options against the budget. |
| `outer_jacket_closed_form.py` | External-pressure collapse of the outer jacket, plain and ring-stiffened. |
| `sizing_closed_form.py` | Closed-form starting walls for the FE sizing. |
| `bagarello_benchmark_1D_check.py` | One-dimensional check of the benchmark heat flux. |

### 3_results_processing (plain Python)

| Script | Purpose |
|---|---|
| `final_numbers.py` | Merges every results file at the sized walls into `results/processed/final_numbers.json`. |
| `hand_calcs.py` | Closed-form checks against the model: membrane stress, equilibrium, head, beam bending, forward inertia, thermal asymptote, fin length, bowing. |
| `make_appendix_traceability.py` | Writes the appendix tables of the thesis from the processed results. |
| `review_checks.py` | Checks added after the supervisors' review: magnitude and depth split of the hydrostatic pressure, the emergency landing against a pure 9 g forward load, the 14 CFR 25.963(d) design pressure, cryogenic rescaling of the thermal quantities, reserve factors at every mesh level, the 24 h hold and the hold at fill 0.90, the pressure-proportionality fit of the far-field peak, and the benchmark without circularity. Writes `results/processed/review_checks.json`. |
| `make_review_docs.py` | Writes `docs/full_results_matrices.md` and `docs/detailed_calculations.md` from the processed results. |
| `compare_delftblue_test.py` | Compares the laptop and cluster solutions of one case. |

### 4_figures_and_drawings (plain Python, plus one Node.js script)

| Script | Purpose |
|---|---|
| `plot_final_figures.py` | Graphs of Chapters 3 and 5. |
| `make_report_contours.py` | Contour and mesh figures of the revised report: panel labels only, vector colour bars read from the Abaqus legends (levels by OCR, colours sampled), values in MPa, K and W/m2. Needs `tesseract`. |
| `draw_workflow.py` | The workflow diagram of Chapter 4. |
| `make_contour_panels.py`, `compose_side_legend.py` | Earlier draft of the contour figures, superseded by `make_report_contours.py` and kept for the record. They write to `earlier_draft_figures/` and never to `thesis_figures/`. |
| `draw_report_diagrams.py` | Resistance network, wall stack, free surface. |
| `draw_pressure_vessel_cuts.py` | Supporting drawing of the two cuts that expose the membrane stresses, written to `design_drawings/`. |
| `draw_tank_design_details.py` | Supporting tank general arrangement and details, written to `design_drawings/`. |
| `make_digital_thread_detailed_pptx.js` | The process flowchart as an editable PowerPoint: `node make_digital_thread_detailed_pptx.js <output.pptx>` (needs `npm install pptxgenjs`). |
| `thesis_figstyle.py` | Shared figure style. |

## Reproducing the numbers and figures

Install the plain-Python dependencies once, then run from anywhere:

```
pip install -r requirements.txt
python code/3_results_processing/final_numbers.py
python code/3_results_processing/hand_calcs.py
python code/2_requirements_and_sizing/insulation_sizing.py
python code/3_results_processing/make_appendix_traceability.py
python code/3_results_processing/review_checks.py
python code/3_results_processing/make_review_docs.py
python code/4_figures_and_drawings/plot_final_figures.py
python code/4_figures_and_drawings/make_report_contours.py abaqus_contour_captures_report thesis_figures
python code/4_figures_and_drawings/draw_workflow.py thesis_figures
python code/4_figures_and_drawings/draw_report_diagrams.py
python code/4_figures_and_drawings/draw_pressure_vessel_cuts.py
python code/2_requirements_and_sizing/boiloff_holdtime_sensitivity.py
```

Run from the repository root on the files in this repository (October 2026), these reproduce
`results/processed/`, `docs/` and `thesis_appendix/appendix_calculation_traceability_full.tex`
(identical apart from the last floating-point digit of one property), and regenerate all 29 figures
of `thesis_figures/` identically. They also write a few supporting graphs that are not in the thesis
(for example `fig5_convergence`, `fig5_heatflux_three`, `fig3_resistance_network`).
`make_report_contours.py` calls `tesseract` to read the Abaqus legend levels, so Tesseract OCR must be
installed and on the path. The model scripts run in the Python interpreter bundled with Abaqus and do
not use `requirements.txt`.

## Figure index

| Thesis figure | File in `thesis_figures/` | Script |
|---|---|---|
| 3.1 | `17_boiloff_holdtime_sensitivity` | `2_requirements_and_sizing/boiloff_holdtime_sensitivity.py` |
| 3.2 | `fig3_wall_stack` | `4_figures_and_drawings/draw_report_diagrams.py` |
| 3.3 | `fig3_sizing_history` | `4_figures_and_drawings/plot_final_figures.py` |
| 4.1 | `fig4_workflow` | `4_figures_and_drawings/draw_workflow.py` |
| 4.2 | `fig4_models` | `4_figures_and_drawings/make_report_contours.py` |
| 4.3 | `fig4_free_surface` | `4_figures_and_drawings/draw_report_diagrams.py` |
| 5.1 | `fig5_e190_fields` | `4_figures_and_drawings/make_report_contours.py` |
| 5.2 | `fig5_e190_matrix` | `4_figures_and_drawings/plot_final_figures.py` |
| 5.3 to 5.5 | `fig5_case_baseline`, `fig5_case_lc3`, `fig5_case_lc7` | `4_figures_and_drawings/plot_final_figures.py` |
| 5.6 to 5.9 | `fig5_svp_baseline`, `fig5_svp_symmetric`, `fig5_svp_combined`, `fig5_svp_emergency` | `4_figures_and_drawings/plot_final_figures.py` |
| 5.10 | `fig5_sensitivity` | `4_figures_and_drawings/plot_final_figures.py` |
| 5.11 | `fig5_emergency_three` | `4_figures_and_drawings/make_report_contours.py` |
| 5.12 | `fig5_thermal` | `4_figures_and_drawings/plot_final_figures.py` |
| 5.13 | `fig5_thermal_three` | `4_figures_and_drawings/make_report_contours.py` |
| 5.14 | `fig5_margins` | `4_figures_and_drawings/plot_final_figures.py` |
| B.1 to B.9 | `appx_contours_<AC>_temperature`, `_heatflux`, `_emergency` | `4_figures_and_drawings/make_report_contours.py` |

## Known differences from the printed thesis

`outer_jacket_closed_form.py` was corrected after submission. It now reads the sized inner walls from
`results/sized_walls.csv` and the MLI thickness of each aircraft from `results/processed/insulation.json`;
the version used for the thesis fell back to provisional walls when `sized_walls.csv` was not in the
working folder and applied the E190 MLI thickness to all three aircraft. The required jacket thickness
changes by at most 0.6 % and no conclusion is affected:

| Outer jacket [mm] | E190 | A320 | A350-900 |
|---|---|---|---|
| Unstiffened, thesis Section 3.4 | 25.7 | 33.7 | 50.6 |
| Unstiffened, this repository | 25.7 | 33.5 | 50.4 |
| Ring spacing 0.25 m, thesis Section 3.4 | 3.45 | 3.98 | 4.88 |
| Ring spacing 0.25 m, this repository | 3.45 | 3.97 | 4.87 |

The load-case labels in `docs/` and in `appendix_calculation_traceability_full.tex` were renamed from
the short codes to the names used in the thesis; the numbers are unchanged.

## Verification

Each check tests a different piece of physics (thesis Section 4.9).

| Check | Agreement |
|---|---|
| Equilibrium, applied load against ring reactions, all 120 cases | within 1.2 % |
| Barrel hoop stress against pR/t | within 0.12 % |
| Mean barrel axial stress against pR/2t | within 0.01 % |
| Dry-wall temperature against q/h_dry | within 0.2 % |
| Axial inertia load path | within 2.3 % |
| Mesh convergence of the far-field stress, production against finest mesh | within 1 % |
| Mesh scatter of the A350-900 buckling RF | about 4 % (1.012, 1.060, 1.021) |
| Abaqus 2025 against Abaqus 2024, one case: field quantities, far-field peak | to printed precision, 0.31 % |
| Peak liner temperature of the Bagarello et al. benchmark, flux from their prescribed temperatures | within 3.7 % |

## Computing environment

Sizing and production matrix: Abaqus/Standard 2025 on a Windows workstation. Stability matrix and
mesh convergence: Abaqus/Standard 2024 on the DelftBlue cluster of TU Delft. Each stage ran on one
version only.

## Citation and licence

Cite the thesis (see `CITATION.cff`, or the "Cite this repository" button on GitHub).
The code is released under the MIT licence (`LICENSE`).
