#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J bisect
#SBATCH -t 08:00:00
#SBATCH -c 1
#SBATCH --array=0-3
#SBATCH -o slurm/logs/bisect_%a.out
#SBATCH -e slurm/logs/bisect_%a.err
#
# Find which recent default made the biofilm evacuate the sand bed.
# Score at 20 d against published Manriquez2026 Fig 9(a): bed ~0.55, supernatant ~0.15.
module load matlab/2025b
ARMS=(baseline kappa0 noupwind rawlight)
A=${ARMS[$SLURM_ARRAY_TASK_ID]}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeRegressionBisect(\"$A\", NCells=100, Days=20)"
