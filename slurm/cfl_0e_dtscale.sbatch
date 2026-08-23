#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cfldt
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/cfldt-%A_%a.out
#SBATCH -e slurm/logs/cfldt-%A_%a.err
#
# 0e follow-up: is the Solver A / diagnostic-sum mismatch a first-order operator
# splitting offset (shrinks with dt) or a structural difference between the two
# equations (dt-independent)? dt is 96% cap-bound, so halving MaxDt halves the
# actual dt for almost every step -- a clean lever.
#   structural  => ratio ~ 1   => component elimination is unsafe
#   splitting   => ratio ~ 4   => benign, converges
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) DT=3e-6 ;;
  2) DT=7.5e-7 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=100, Kappa=1e-6, MatureDays=3, MaxDt=$DT)"
