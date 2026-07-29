#!/bin/bash -l
#SBATCH -A lu2025-7-124
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J lightdiag
#SBATCH -o /home/jrman/MPC-SSF/slurm/logs/lightdiag-%j.out
#SBATCH -e /home/jrman/MPC-SSF/slurm/logs/lightdiag-%j.err
set -euo pipefail
module purge
module load Julia/1.10.4-linux-x86_64
cd /home/jrman/MPC-SSF/julia
julia --project=. -e '
include("analysis/light_profile_sweep.jl")
m = light_profile_model()
# Is the failure dt-sensitive (numerical) or dt-invariant (physics)?
for dt in (3e-6, 1e-6, 3e-7)
    f = light_filter(amplitude=0.8, ncells=500, delta=20e-3)
    r = simulate(State(f, m); inflow_concentrations=INFLUENT_LIGHT, simulation_time=0.30,
                 time_step=dt, n_frames=30, quiet=true, implicit_osmosis=true)
    err = get(r.simulation_data, :error, nothing)
    println("dt=", dt, "  flag=", r.flag, "  t_final=", round(r.time_final, digits=5),
            "  error=", err === nothing ? "none" : err)
end
# Control: does the FULL Lund model survive the same mesh and dt?
println()
f = SandFilter(temperature=15.0, sand_roughness=20e-3); f = addgridpoints(f, 500)
INF_LUND = [1e-3,1e-3,0.0,1e-5,1e-2,1e-2,1e-5,0.0,1e-4]
r2 = simulate(State(f, modelLund()); inflow_concentrations=INF_LUND, simulation_time=0.30,
              time_step=3e-6, n_frames=30, quiet=true, implicit_osmosis=true)
println("CONTROL modelLund dt=3e-6 ncells=500: flag=", r2.flag, " t_final=", round(r2.time_final,digits=5),
        " error=", get(r2.simulation_data, :error, "none"))'
