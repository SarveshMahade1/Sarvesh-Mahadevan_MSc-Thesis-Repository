# Development notes: E190_MASTER_thermo_structural.py

How the finite element model reached its final form, and which modelling errors were found and
corrected on the way. The script itself describes only the final method; these notes keep the
record of why some of its choices are made, so that a later user does not reintroduce an error that
has already been found. The numbers are those reported by the model at the time (mostly E190, at
earlier wall thicknesses than the sized walls of the thesis).

Load-case codes as in the file names: `Baseline` Baseline condition, `LC1_Maneuver` symmetric
manoeuvre, `LC3_Combined` combined manoeuvre, `LC7_Emergency` emergency landing.

## 1. Model versions

Version 2 of the model sized the wall against yield and buckling, used ring and saddle mounts and
applied a uniform 20 K temperature field. Version 3, the model of the thesis, made four deliberate
changes:

1. **Inner vessel only** (PROCESS 10). The outer jacket, insulation, struts and airframe mounts are
   not modelled. The insulation appears only as a heat flux on the outer face of the vessel.
2. **A real load path** (PROCESS 6). The vessel is held on three circumferential rings, at the two
   cap-to-barrel junctions and at mid-barrel (Bagarello et al. Sec. 4.4.1, from H2FLY testbed
   hardware). Two earlier support models were rejected; see Section 2.
3. **Two sizing criteria** (PROCESS 8). Yield (far-field von Mises at most F_ty/1.5 at the design
   pressure) and buckling (minimum axial membrane stress against the NASA SP-8007 allowable divided
   by 1.5 at the stability pressure). The buckling check is a criterion, not an eigenvalue or
   post-buckling analysis.
4. **A thermal load that depends on the fill level** (PROCESS 3 and 4). A uniform 20 K field on a
   body free to contract produces exactly zero thermal stress. The wetted and dry regions are now
   treated separately, so the fill level generates the temperature gradient.

## 2. Supports and load path

- **3-2-1 restraint.** A statically determinate restraint carries zero reaction only for
  self-equilibrated loads. The CS-25 inertia cases are not: in the emergency landing the hydrostatic
  pressure field alone has a resultant of 562 kN (5104 kg x 11.225 x 9.81 m/s2) and the shell body
  force adds 181 kN. The whole 743 kN passed through three nodes of the 4.785 mm shell and gave
  4943 MPa and a 6.9 m displacement. Kept as `apply_rigid_body_restraint`, for self-equilibrated
  studies only.
- **Inertia relief.** It removed that artefact but introduced another. The balancing acceleration
  is distributed in proportion to mass, and the liquid has no mass in the model, so the whole liquid
  reaction was smeared over the 1648 kg shell as a fictitious body force of up to 34.8 g. It also
  cancelled the shell's own body force exactly, so the applied inertia had no effect. Kept as
  `apply_inertia_relief`.
- **Conclusion.** A CS-25 load factor presupposes a load path: a freely floating tank does not
  experience a crash load, it accelerates. The stress arises from being held, so modelling the
  restraint is what makes the load case meaningful. Bagarello et al. report that realistic anchoring
  raises liner stress about fourfold relative to rigid-body removal alone.
- **Partial mid-barrel ring.** A fixed node-selection tolerance found 124 nodes at each junction
  but only 31 at mid-barrel, against about 52 for a complete ring, so the axial anchor held only part
  of the circumference. The ring now snaps to the nearest node plane.

## 3. Thermal coupling

- **Through-thickness gradient from FROM_FILE.** With the temperature field read FROM_FILE onto a
  `temperature=GRADIENT` shell section, a mid-barrel element of `ST_Baseline_f10.odb` showed hoop
  stresses of 85.99 MPa (SNEG) and 36.75 MPa (SPOS): a membrane hoop stress of 61.37 MPa instead of
  40.80 MPa and a bending stress of 24.62 MPa where a pressurised cylinder has none at mid-span. The
  bending implied a through-thickness difference of 20.2 K; the real one is of order 1e-4 K. It
  roughly doubled every stress and was invisible to the load balance, because through-thickness
  bending is self-equilibrated. The field is now applied in bins with
  `CONSTANT_THROUGH_THICKNESS`, and the membrane probe in the results file checks it.
- **Missing initial temperature.** Without an initial temperature condition Abaqus takes 0 as the
  initial temperature; one run gave u_max = 3.22 mm with the thermal load effectively absent.
- **Stress-free reference at 293 K.** Applying the full cooldown on top of every load case gave
  649.5 MPa, 98 percent of the fully restrained value of 663 MPa. The reference is the cold
  operating state, 20.3 K; the cooldown is a separate ground case (further work).
- **Horizontal wetted/dry split.** The thermal split once used a horizontal free surface for every
  load case while the structural model used the tilted one. The signature was a wetted-to-dry dT
  identical to all digits across the four load cases at each fill (3.62878036499023 K at fill 0.10).
  Both models now use the same free surface, computed once per design point.

## 4. Loads

- **Vertical depth.** Depth was once measured vertically while the pressure was scaled by the
  magnitude of the load factor. Exact for the Baseline condition and the symmetric manoeuvre, but at
  fill 0.90 the correct maximum depth is 7.903 m (combined manoeuvre, 31.0 degree tilt) and
  11.441 m (emergency landing, 57.7 degrees) instead of 2.188 m. The emergency landing hoop stress
  at the deepest point rose 43 percent with the correction.
- **Transposed axes.** The vertical and lateral load factors were once swapped, applying 6 g
  sideways and 3 g downwards in the emergency landing.
- **Forward sign.** comp1 = -nx g threw the liquid aft. Magnitudes were the same by symmetry, but
  the peak appeared at the wrong cap.
- **Body force units.** BodyForce expects force per unit volume; passing an acceleration made the
  shell inertia 2850 times too small (51 N instead of 181.5 kN). Gravity is used instead.
- **Pressure face.** With the pressure on side1 the hydrostatic reaction came out at -561.4 kN
  against an applied +562.0 kN: the liquid pushed inward. side2 is the inner face.
- **Earlier wall thicknesses.** Walls of 4.785 mm (E190), 6.382 mm (A320) and 9.948 mm (A350-900)
  came from sizing runs that put the pressure on the outer face, so the buckling check sized against
  a hoop compression that does not exist in service. They were replaced by the corrected
  closed-form starting walls and the FE sizing (`sizing_closed_form.py`).

## 5. Load balance check

- **Cap area.** The analytic shell area once used 2 pi R h_cap for the caps, 59 percent low for the
  caps and 7.9 percent for the shell, which showed as a systematic 1.079 ratio in the X and Z
  components of every case. The exact spheroid area and the meshed mass are now used.
- **Liquid resultant direction.** The check once placed the liquid resultant on the vertical axis,
  correct while depth was vertical but not after the tilted free surface was introduced. It gave
  component ratios of 3.8 and 0.6 for cases balanced to better than one percent.

## 6. Benchmark

- The benchmark liner was once 1 mm thick; Bagarello Table 4 states 10 mm. The tenfold error put
  the hoop stress at 360 MPa instead of 36 MPa.
- A placeholder outer flux of 22.0 W/m2 drove the dry wall to 35.4 K against their 25.3 K. The flux
  is now derived from their stated boundary conditions (8.0 W/m2, cross-checked at 8.38 W/m2
  through the foam), not from their result.

## 7. Running Abaqus

- **Star imports.** In noGUI mode, methods such as `FieldOutputRequest` were missing until the
  standard module set was imported.
- **Built-ins.** The star imports replace `sum()` with a version that rejects generators; the Python
  built-ins are restored.
- **Job status.** In the first full noGUI run all twenty heat transfer jobs completed in about
  twelve seconds each and were all rejected, because `Job.status` is not updated without the GUI.
  Completion is read from the job's .log file.
- **Licence queue.** On the production run of 18 September 2026, 32 jobs used 259 s of computation
  in 79 min of wall-clock time; most of it was spent waiting for a licence token. The load-case and
  fill filters split the matrix over parallel sessions.
- **Job names.** The E190 matrix of 19 September 2026 was written with bare job names, so the E190
  keeps them and the other aircraft take a prefix.
- **Memory.** A FINAL session builds about 70 design points; finished models are deleted from the
  session (`_free`).
- **DelftBlue.** Abaqus 2024 on Linux does not define `__file__` in noGUI mode, so the default
  folders could not be derived from the script location. The working folder, the design basis and
  `sized_walls.csv` are now located without it (see the DESIGN BASIS block of the script).

## 8. Sizing loop

- Revised on 27 September 2026 after the first iteration on the campus machines: three requirements
  (buckling, yield, optional minimum gauge), FE sizing of the governing case at both pressures.
- A narrow acceptance window (buckling FOS 1.50 to 1.53) made the A350-900 oscillate: the mesh
  changes with the wall, so a 2.5 percent change in t moved its buckling FOS from 1.44 to 1.59.
  The window is now about 6 percent, with bisection between a passing and a failing wall.
- Superseded sizing scripts had the tank axis along Z, where the default local direction made S11
  the hoop stress at crown and keel. The model is X-axial and checks the orientation in every run.

## 9. Mesh convergence

- The earlier studies, at the thicker walls, covered 0.5 to 8 elements per boundary layer. The peak
  stress at the rings rose 108.0, 108.0, 122.1, 140.6 MPa on successive refinements (fitted
  sigma ~ h^-0.19) while the reaction converged to 0.2 percent and the displacement to 0.9 percent:
  the singularity that the far-field extraction excludes.
- At the sized walls the study uses 1, 2 and 4 elements per boundary layer; a level of 8 would be 2
  to 4 million elements.

## 10. Timeline

| Date (2026) | Change |
|---|---|
| 18 to 19 Sep | Production matrix of the E190 (version 3), parallel sessions |
| 19 Sep | Benchmark liner thickness corrected |
| 27 Sep | Walls corrected, FE sizing with buckling and yield, FINAL mode, recorded checks, cluster runs |
| 28 Sep | Cryogenic property overrides (E190_ALPHA, E190_KAPPA) |
| Oct | Design basis file, dry run, comments rewritten for readers; design basis and sized walls located without `__file__` |
