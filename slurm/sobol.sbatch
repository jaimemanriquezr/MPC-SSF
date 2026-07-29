#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 48
#SBATCH --exclusive
#SBATCH -t 24:00:00
#SBATCH -J sobol
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/sobol-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/sobol-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
# -c 48 is load-bearing: with --exclusive alone Slurm still sets CPUs/Task=1 and
# pins the task to a single core, so Julia's threads all share one CPU.
NT=${SLURM_CPUS_PER_TASK:-48}
export JULIA_NUM_THREADS=$NT
cd /home/jrman/MPC-SSF/julia
N=${SOBOL_N:-256}
echo "host=$(hostname) cpus_per_task=$NT N=$N runs=$((N*10))"
julia --project=. -t $NT -e 'println("julia sees ", Threads.nthreads(), " threads")'
julia --project=. -t $NT analysis/sobol_sensitivity.jl $N 25 1.0
