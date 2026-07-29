#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J mpcssf-pkgtest
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/pkgtest-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/pkgtest-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
cd /home/jrman/MPC-SSF/julia
echo "julia: $(julia --version)  host: $(hostname)"
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
