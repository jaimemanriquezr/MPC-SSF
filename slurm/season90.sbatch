#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 08:00:00
#SBATCH -J season90
#SBATCH -o slurm/logs/season90_%a.out
#SBATCH -e slurm/logs/season90_%a.err
#SBATCH --array=0-5
# 90-day seasonal pairs at manuscript resolution (100 cells, MaxDt 3e-6),
# three variants x two seasons. Isolates the three confounds between the
# probe chain (10-30 d, field influent, corrected kinetics -> winter > summer)
# and fig:seasons-results (90 d, Table B.1 influent, published kinetics ->
# summer > winter). Reports BOTH phi_b and mass metrics, raw and
# porosity-weighted. Reference: E1/E3 at 100 cells took ~75 min locally.
module load matlab/2025b
SEASONS=(summer winter summer winter summer winter)
VARIANTS=(manuscript manuscript corrected corrected field field)
S=${SEASONS[$SLURM_ARRAY_TASK_ID]}
V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); probeSeason90(\"$S\", \"$V\")"
