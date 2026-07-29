#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 24:00:00
#SBATCH -J morris
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/morris-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/morris-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64

# Re-screen post-cohesion-fix. The committed results/morris/mu_sigma.csv is from
# 2026-07-21, a week before the fix, so it carries the double-well potential.
#
# 24 trajectories, not the 6 behind the stale file: 6 is far below the
# conventional 10-50 and left the magnitudes indicative at best. 24 x (22+1) =
# 552 runs. The driver is single-threaded, so -c 2 is honest about what it uses;
# there is nothing here for more cores to do.
R_TRAJ=${MORRIS_R:-24}
TSIM=${MORRIS_TSIM:-0.3}
NCELLS=${MORRIS_NCELLS:-30}

cd /home/jrman/MPC-SSF/julia
echo "host=$(hostname) r_traj=$R_TRAJ tsim=$TSIM ncells=$NCELLS runs=$((R_TRAJ*23))"
julia --project=. analysis/morris_screening.jl $R_TRAJ $TSIM $NCELLS
