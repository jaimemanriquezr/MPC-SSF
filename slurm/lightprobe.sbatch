#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -J lightprobe
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/lightprobe-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/lightprobe-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
cd /home/jrman/MPC-SSF/julia
# Short run at the target mesh: is fixed dt=3e-6 stable at dz=2mm, what does it
# cost, and does any biofilm appear above z=0 now that the roughness layer has cells?
julia --project=. -e '
include("analysis/light_profile_sweep.jl")
for nc in (100, 500)
    try
        run_one(0.8; tsim=0.2, ncells=nc, nframes=10)
    catch e
        println("ncells=", nc, " ERROR: ", sprint(showerror, e))
    end
end'
