#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -t 06:00:00
#SBATCH -J massChk
#SBATCH -o slurm/masscheck_%a.out
#SBATCH --array=0-1
module load matlab/2025b
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
  matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); probeMassCheck(\"winter\", 0)"
else
  matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); probeTheta(\"summer\", 2.0)"
fi
