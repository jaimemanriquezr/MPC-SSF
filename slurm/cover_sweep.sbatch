#!/bin/bash -l
# Covered/uncovered driver sweep -- one array task per registry cell.
#
#   sbatch --export=STAGE=A      --array=1-5  slurm/cover_sweep.sbatch
#   sbatch --export=STAGE=B      --array=1-16 slurm/cover_sweep.sbatch
#   sbatch --export=STAGE=Bprime --array=1-4  slurm/cover_sweep.sbatch
#   sbatch --export=STAGE=C      --array=1-6  slurm/cover_sweep.sbatch
#
# Submit FROM THE REPO ROOT -- every addpath below is relative to the submit
# directory. `#!/bin/bash -l` is required on cosmos or `module` is not found.
#
# The array index is an index into coverCases(STAGE), NOT a hand-typed parallel
# bash array. season90.sbatch and theta_sweep.sbatch carry two or three parallel
# arrays that must be kept in sync by hand; that is the maintenance trap this
# avoids. To change the sweep, edit analysis/probes/coverCases.m and adjust
# --array to match numel(coverCases(STAGE)).
#
# Cell counts (2026-08-22): A=5  B=16  Bprime=4  C=6  all=29
# Wall estimates at 100 cells / MaxDt 3e-6, ~0.8 min per simulated day:
#   30 d ~ 25 min, 104 d ~ 85 min, 10 d ~ 8 min (the n300 Bprime arm is slower).
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 04:00:00
#SBATCH -J coverSw
#SBATCH -o slurm/logs/cover_%A_%a.out
#SBATCH -e slurm/logs/cover_%A_%a.err

set -euo pipefail

: "${STAGE:?set STAGE via --export=STAGE=A}"
mkdir -p slurm/logs

echo "=== cover_sweep stage=$STAGE task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
               probeCover($SLURM_ARRAY_TASK_ID, Stage=\"$STAGE\")"

echo "=== done $(date -Is) ==="
