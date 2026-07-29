#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J sobolcost
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/sobolcost-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/sobolcost-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
cd /home/jrman/MPC-SSF/julia
# Single-threaded cost of ONE Sobol design evaluation at the running job's
# settings, so elapsed CPU time can be converted into runs completed.
julia --project=. -t 1 -e '
include("analysis/sobol_sensitivity.jl")
m0 = pathogen_model()
t0 = time(); ms = build_mature(m0, 25, 3.0); tb = time()-t0
println("BUILD_MATURE_S=", round(tb, digits=1))
using Random
rng = MersenneTwister(1)
for i in 1:3
    u = rand(rng, length(SPARAMS))
    t1 = time(); y, fl = eval_qoi(u, ms, m0, 25, 1.0, 30); dt1 = time()-t1
    println("EVAL_S=", round(dt1, digits=2), " flag=", fl, " y=", round(y, digits=4))
end' 2>&1
