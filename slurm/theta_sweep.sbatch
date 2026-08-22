#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 06:00:00
#SBATCH -J thetaSw
#SBATCH -o slurm/theta_%a.out
#SBATCH --array=0-3
module load matlab/2025b
SEASONS=(summer winter summer winter)
THETAS=(1.09 1.09 1.12 1.12)
S=${SEASONS[$SLURM_ARRAY_TASK_ID]}
T=${THETAS[$SLURM_ARRAY_TASK_ID]}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); probeTheta(\"$S\", $T)"
