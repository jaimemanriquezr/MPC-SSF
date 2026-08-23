#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cflpar
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/cflpar-%A_%a.out
#SBATCH -e slurm/logs/cflpar-%A_%a.err
#
# 0f: let Solver A's phi_b evolve in PARALLEL -- seeded once, never re-seeded
# from Solver B -- and measure the accumulated drift against the component sum.
# 0e re-seeds every step and so can only ever show O(dt); the component-
# elimination design makes Solver A authoritative for the whole run, which is
# what this tests.
#
# A short N=40 check already blew up (phi_b < 0 at t = 4.2e-7 d, then exponential
# growth to 1e276 and NaN by t = 0.17). That was the INACTIVE regime where
# phi_b ~ 0 and the mobility is degenerate, so it needs confirming with biofilm
# actually present.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100 ;;
  2) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=$N, Kappa=1e-6, MatureDays=3)"
