#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J schcmp
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/schcmp-%A_%a.out
#SBATCH -e slurm/logs/schcmp-%A_%a.err
#
# Three-way: (Shin + explicit Solver B), (Shin + implicit dispersion),
# (Bailo + implicit dispersion). Same model, mesh, inflow and 3 d duration, so
# cost, accuracy, what binds, and the free-running Solver A bounds are directly
# comparable. N=100 for a quick read, N=500 for the production answer.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100 ;;
  2) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeSchemeComparison(NCells=$N, Kappa=1e-6, Days=3)"
