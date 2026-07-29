#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J lightsweep
#SBATCH --array=1-5
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/lightsweep-%A_%a.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/lightsweep-%A_%a.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
# Amplitudes bracket the photoinhibition optimum logarithmically. I_eff = 1 at the
# supernatant top needs A ~ 0.018; the default is 0.8, i.e. ~44x into inhibition.
AMPS=(0.018 0.05 0.2 0.8 3.2)
A=${AMPS[$((SLURM_ARRAY_TASK_ID-1))]}
cd /home/jrman/MPC-SSF/julia
echo "amplitude=$A tsim=7 ncells=500 delta=20mm dt=3e-7 host=$(hostname)"
julia --project=. analysis/light_profile_sweep.jl $A 7.0 500 3e-7
