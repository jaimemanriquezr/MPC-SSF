#!/bin/bash
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J msRest
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --array=0-5
#SBATCH -o slurm/logs/msRest-%A_%a.out
#SBATCH -e slurm/logs/msRest-%A_%a.err
# Stage 2: remaining experiments in parallel, each reading the cached
# mature30_summer built by manuscript_summer.sbatch.
# Submit: sbatch --dependency=afterok:<jobid-of-stage-1> slurm/manuscript_rest.sbatch
module load matlab/2025b
IDS=(E1 E2 E4 E5 E6 E10)
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  manuscriptExperiments('${IDS[$SLURM_ARRAY_TASK_ID]}', NCells=500, MaxDt=1e-5)"
