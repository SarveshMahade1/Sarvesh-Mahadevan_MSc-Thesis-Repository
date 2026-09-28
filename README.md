# Digital thread for the thermo-structural analysis of liquid hydrogen tanks

Model, scripts and complete results of the MSc thesis
*Development of a Digital Thread of a Liquid Hydrogen Tank Integrating Thermo-Structural Analysis*,
Sarveshwaran Mahadevan, MSc Sustainable Energy Technology, Delft University of Technology, 2026.

One parametric Abaqus script builds, meshes, sizes, solves and post-processes a liquid hydrogen
tank for three aircraft classes: the Embraer E190, the Airbus A320neo and the Airbus A350-900.
Each tank is sized to carry the hydrogen that delivers the same energy as the reference kerosene
load. Its aluminium 2219-T87 inner wall is sized so that it does not yield at 1.5 bar, and does
not buckle at 0.187 bar (factor of safety 1.5 on the NASA SP-8007 allowable). It is then analysed
under four CS-25 load conditions at five fill ratios, at both pressures:

- 120 production and stability cases, each a steady-state heat transfer analysis followed by a
  linear static analysis;
- 18 mesh-convergence cases.

Every number in the thesis can be traced to a file in this repository and regenerated from it.

## Main results

| | E190 | A320 | A350-900 |
|---|---|---|---|
| Barrel slenderness L/2R | 5.37 | 2.93 | 6.56 |
| Sized wall [mm] | 1.289 | 1.572 | 5.674 |
| Wall governed by | yield | yield | buckling |
| Governing far-field von Mises, LC7, fill 0.90 [MPa] | 256.9 | 257.3 | 145.7 |
| Lowest yield factor of safety (target 1.0 on 262 MPa) | 1.020 | 1.018 | 1.798 |
| Lowest buckling factor of safety (target 1.5) | 1.657 | no compression | 1.590 |
| Fill sensitivity under LC7 [%] | +46.4 | +36.0 | +100.9 |
| Wetted-to-dry temperature difference, fill 0.10 [K] | 3.683 | 4.644 | 7.364 |

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
  mesh_convergence/            governing case (LC7, fill 0.90), three meshes, both pressures
  wall_sizing_iterations/      every iteration of the FE wall sizing, per aircraft
  benchmark_bagarello/         the reference tank of Bagarello et al. run through the method
  solver_check_delftblue/      one case solved on the DelftBlue cluster (Abaqus 2024)
  processed/                   final_numbers.json, hand_calcs.json, insulation.json,
                               production_matrix_sized_walls.csv

abaqus_contour_captures/       the 60 Abaqus images the thesis contour and mesh figures are built from
thesis_figures/                every figure of the thesis, as generated
thesis_appendix/               the calculation-traceability appendix, as generated
design_drawings/               tank general arrangement and details, process flowchart (editable .pptx)
```

In every results file name, `<AC>` is `E190`, `A320` or `A350`. The Abaqus output databases
(.odb) are not included, because they amount to tens of gigabytes; running the model regenerates
them.

## The code

### 1_abaqus_model (needs an Abaqus licence)

| Script | Purpose |
|---|---|
| `E190_MASTER_thermo_structural.py` | The digital thread. Despite the name it serves all three aircraft. Geometry, mesh (seed 0.5 sqrt(Rt)), free surface normal to the load factor, steady-state heat transfer, banded temperature transfer, static solution with pressure, hydrostatic head and inertia, three ring supports, far-field extraction, verification gates, FE wall sizing. |
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
| `compare_delftblue_test.py` | Compares the laptop and cluster solutions of one case. |

### 4_figures_and_drawings (plain Python, plus one Node.js script)

| Script | Purpose |
|---|---|
| `plot_final_figures.py` | Graphs of Chapters 3 and 5. |
| `make_contour_panels.py`, `compose_side_legend.py` | Contour and mesh figures from the captures. |
| `draw_report_diagrams.py` | Resistance network, wall stack, free surface. |
| `draw_pressure_vessel_cuts.py` | The two cuts that expose the membrane stresses. |
| `draw_tank_design_details.py` | Tank general arrangement and details. |
| `make_digital_thread_detailed_pptx.js` | The process flowchart as an editable PowerPoint (needs `npm install pptxgenjs`). |
| `thesis_figstyle.py` | Shared figure style. |

## Reproducing the numbers and figures

Install the plain-Python dependencies once, then run from anywhere:

```
pip install -r requirements.txt
python code/3_results_processing/final_numbers.py
python code/3_results_processing/hand_calcs.py
python code/2_requirements_and_sizing/insulation_sizing.py
python code/3_results_processing/make_appendix_traceability.py
python code/4_figures_and_drawings/plot_final_figures.py
python code/4_figures_and_drawings/make_contour_panels.py
python code/4_figures_and_drawings/draw_report_diagrams.py
python code/4_figures_and_drawings/draw_pressure_vessel_cuts.py
python code/2_requirements_and_sizing/boiloff_holdtime_sensitivity.py
```

Run on the files in this repository, these reproduce `results/processed/` exactly and regenerate
`thesis_figures/`. The model scripts run in the Python interpreter bundled with Abaqus and do not
use `requirements.txt`.

## Verification

Each check tests a different piece of physics (thesis Section 5.12).

| Check | Agreement |
|---|---|
| Equilibrium, applied load against ring reactions, all 120 cases | within 1.2 % |
| Barrel hoop stress against pR/t | within 0.12 % |
| Mean barrel axial stress against pR/2t | within 0.01 % |
| Dry-wall temperature against q/h_dry | within 0.2 % |
| Axial inertia load path | within 2.3 % |
| Mesh convergence, production against finest mesh | within 1 % |
| Abaqus 2025 against Abaqus 2024, one case | within 0.31 % |
| Peak liner temperature of the Bagarello et al. benchmark | within 0.9 % |

## Computing environment

Sizing and production matrix: Abaqus/Standard 2025 on a Windows workstation. Stability matrix and
mesh convergence: Abaqus/Standard 2024 on the DelftBlue cluster of TU Delft. Each stage ran on one
version only.

## Citation and licence

Cite the thesis (see `CITATION.cff`, or the "Cite this repository" button on GitHub).
The code is released under the MIT licence (`LICENSE`).
