#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 48:00:00
#SBATCH -J sobolr
#SBATCH -o /home/jrman/MPC-SSF-sobol/slurm/logs/sobolr-%j.out
#SBATCH -e /home/jrman/MPC-SSF-sobol/slurm/logs/sobolr-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64

# Runs from MPC-SSF-sobol, a copy of the tree, so it cannot collide with job
# 3427791 writing into MPC-SSF/julia/analysis/results/sobol/. Both jobs may run
# at once; whichever finishes first is the result.

# 16 cores, NOT 48 exclusive. Two reasons.
#
# Throughput: job 3427791 ran 48 threads and sustained CPULoad 9-14, so the
# surplus threads contend in GC rather than doing work. Fewer threads is the
# fix, not more.
#
# Queueing: lu48 had 147 nodes allocated and 37 partly used when this was
# written, and an --exclusive 48-core request was scheduled 21 h out. A 16-core
# request fits a partly-used node and starts now. Sharing a node risks the
# memory-bandwidth contention seen in the lightsweep array, but a job that
# starts tomorrow afternoon is worse than one that starts today and contends.
NT=${SOBOL_THREADS:-16}
export JULIA_NUM_THREADS=$NT

# 48 h, not 24 h: lu48 allows 7 days, so the old wall was a choice. With the
# resume path below, a job that still runs out continues in the next one
# instead of starting from nothing.
cd /home/jrman/MPC-SSF-sobol/julia
N=${SOBOL_N:-256}
echo "host=$(hostname) threads=$NT of ${SLURM_CPUS_PER_TASK:-?} N=$N runs=$((N*10))"
julia --project=. -t $NT -e 'println("julia sees ", Threads.nthreads(), " threads")'
julia --project=. -t $NT analysis/sobol_sensitivity.jl $N 25 1.0 resume
