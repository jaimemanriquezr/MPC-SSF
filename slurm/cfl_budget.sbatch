#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cflbud
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/cflbud-%A_%a.out
#SBATCH -e slurm/logs/cflbud-%A_%a.err
#
# Phase 0 of the Bailo investigation (.claude/plans/2026-08-23-bailo-scheme.md).
# Two meshes at kappa=1e-6, matured 3 d so the supernatant carries biofilm --
# a clean start makes 0c a null test (phib_sup ~ 0 => cohesion inert).
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100 ;;
  2) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=$N, Kappa=1e-6, MatureDays=3)"
