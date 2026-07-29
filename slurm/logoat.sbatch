#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 24:00:00
#SBATCH -J logoat
#SBATCH --array=1-3
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/logoat-%A_%a.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/logoat-%A_%a.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
SCEN=(pulse startup flowstep)
S=${SCEN[$((SLURM_ARRAY_TASK_ID-1))]}
cd /home/jrman/MPC-SSF/julia
echo "scenario=$S host=$(hostname) julia=$(julia --version)"
julia --project=. -e "include(\"analysis/log_oat_sensitivity.jl\"); main(disturbance=:$S)"
