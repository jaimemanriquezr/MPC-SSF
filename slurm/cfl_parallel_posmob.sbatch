#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cflpm
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/cflpm-%A_%a.out
#SBATCH -e slurm/logs/cflpm-%A_%a.err
#
# 0f with the mobility clamped to its positive part, max(zeta_0*u(1-u), 0), to
# pair against the unclamped run (3531603) in the ACTIVE regime.
#
# Short N=40 result: clamping removes the blow-up (unclamped reaches 1e276 and
# NaN by t=0.17; clamped stays in [-0.081, 0.087] indefinitely) but does NOT fix
# the design -- phi_b still goes negative at the same instant, t = 4.19e-7 d, and
# the drift saturates at 5x the phi_b scale. The two trajectories are identical
# to t=0.09, so the clamp only caps a runaway that is already exponential with
# positive mobility.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100 ;;
  2) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=$N, Kappa=1e-6, MatureDays=3, PositiveMobility=true)"
