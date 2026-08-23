#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J bailoc
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -a 1-4
#SBATCH -o slurm/logs/bailoc-%A_%a.out
#SBATCH -e slurm/logs/bailoc-%A_%a.err
#
# Bailo cost and free-running bounds, in the ACTIVE regime (3 d), paired against
# shin at the same mesh. The only Newton count I have so far -- 0.7 it/step -- is
# from a 0.3 d run, i.e. pure startup, where most steps are inactive and Newton
# exits at zero iterations. The standalone prototype gave 1.8-4.2, which is the
# more credible figure; this settles it in the real solver.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=100; SC=shin  ;;
  2) N=100; SC=bailo ;;
  3) N=500; SC=shin  ;;
  4) N=500; SC=bailo ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=$N, Kappa=1e-6, MatureDays=3, CohesionScheme=\"$SC\")"
