#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cflmrg
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-3
#SBATCH -o slurm/logs/cflmrg-%A_%a.out
#SBATCH -e slurm/logs/cflmrg-%A_%a.err
#
# Re-run 0a-0c with the MARGIN recorded: how close each region comes to being
# the binding one. Reporting only the argmax hid whether the matrix region --
# the one v_b actually enters -- lost narrowly or by orders of magnitude.
# Three meshes, because matrix/binding scales ~dz^-2 and the trend is the point.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100 ;;
  2) N=200 ;;
  3) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=$N, Kappa=1e-6, MatureDays=3)"
