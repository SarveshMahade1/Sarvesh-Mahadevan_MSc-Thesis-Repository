#!/bin/bash
# =============================================================================
# REMAINING FINAL STAGES ON DELFTBLUE (stability and mesh convergence)
# TU Delft MSc thesis, Sarvesh Mahadevan, 27 Sep 2026
#
# One SLURM job = one aircraft x one stage. Submit six of them:
#
#   cd /scratch/$USER/thesis
#   for AC in E190 A320 A350; do
#     sbatch --job-name=${AC}_stab --export=ALL,AC=$AC,STAGE=stability   code/1_abaqus_model/run_final_delftblue.sh
#     sbatch --job-name=${AC}_conv --export=ALL,AC=$AC,STAGE=convergence code/1_abaqus_model/run_final_delftblue.sh
#   done
#
#   squeue --me                 status
#   tail -f slurm-<jobid>.out   live output
#
# Quick comparison test against the laptop (one case, about 10 minutes):
#   see the header of code/3_results_processing/compare_delftblue_test.py
#
# Folder layout expected in /scratch/$USER/thesis:
#   sized_walls.csv
#   code/1_abaqus_model/run_final_delftblue.sh
#
# The wall thickness is read from sized_walls.csv, exactly as on the laptop.
# Sizing and production are NOT re-run (E190_STAGES holds only the one stage).
# The laptop thermal results are not available here, so each stability case
# solves its own heat transfer job first (same model, same seed); the MASTER
# does this automatically when no reusable thermal result is found.
#
# Solver version: Abaqus 2024 (the newest release on DelftBlue). All cases of
# a stage run on the same version, so each stage is self-consistent. State
# the version in the report.
# =============================================================================

#SBATCH --partition=compute
#SBATCH --account=education-eemcs-msc-set
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3900M
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -e
WORKDIR=${WORKDIR:-/scratch/$USER/thesis}
: "${AC:?set AC=E190|A320|A350 with --export}"
: "${STAGE:?set STAGE=production|stability|convergence with --export}"

echo "=============================================================="
echo " $AC  $STAGE   job $SLURM_JOB_ID on $(hostname)   $(date)"
echo "=============================================================="

cd "$WORKDIR"
# repository layout: take the sized walls from results/ if none is in the working folder
[ -f sized_walls.csv ] || cp results/sized_walls.csv . 2>/dev/null || true
[ -f sized_walls.csv ] || { echo "ERROR: sized_walls.csv missing in $WORKDIR"; exit 1; }

export E190_WORKDIR="$WORKDIR"
export E190_AIRCRAFT="$AC"
export E190_RUN_MODE=FINAL
export E190_STAGES="$STAGE"

module load abaqus/2024
# CAE links against libraries the compute nodes do not have (checked with ldd
# on a compute node, 27 Sep). /scratch/$USER/libs holds:
#   libjpeg.so.62  copied from /usr/lib64 on the login node
#   libGLU.so.1    link to the 2025 stack module mesa-glu/9.0.2
#   libOSMesa.so.6 link to libOSMesa.so.8 of the 2025 stack module mesa/23.3.6
#   libstdc++.so.6, libgcc_s.so.1  links to the GCC 13.3.0 runtime
#                  (/apps/generic/compiler/2025/...), needed by Mesa's LLVM
export LD_LIBRARY_PATH=/scratch/$USER/libs:${LD_LIBRARY_PATH:-}
LAUNCHER=$(command -v abq2024 || command -v abaqus)
echo "launcher: $LAUNCHER"

$LAUNCHER cae noGUI=code/1_abaqus_model/E190_MASTER_thermo_structural.py

echo "finished $(date)"
ls -la "$WORKDIR"/${AC}_${STAGE}*_results.csv 2>/dev/null || true
