using MPCSSF, Printf
include("/Users/jaime/Research/SSF/MPC-SSF/julia/analysis/proxy_model.jl")
include("/Users/jaime/Research/SSF/MPC-SSF/julia/analysis/pathogen_repro.jl")
const INFL = Float64[2.68e-3,1.00e-2,0.0,5.36e-3,9.10e-3,6.23e-3,2.00e-5,0.0,1.75e-4]
gz(f) = try; MPCSSF.gridzero(f); catch; findfirst(z->abs(z)<gridsize(f)/2, f.grid.centers); end
f0 = SandFilter(temperature=15.0)
@printf("sand_roughness δ=%.4g m ; sand_porosity ε0=%.3f ; domain height=%.2f depth=%.2f\n\n",
        f0.sand_roughness, f0.sand_porosity, f0.height, f0.depth)
for nc in (20, 40, 80)
    m = pathogen_model()
    f = addgridpoints(SandFilter(temperature=15.0), nc)
    r = run_proxy(State(f, m); simulation_time=1.5, inflow_concentrations=INFL, n_frames=3, quiet=true)
    z = f.grid.centers; dz = gridsize(f); poro = computeporosity(f, z)
    phib = get_volume_fractions(r).biofilm[:, end]; ip = argmax(phib); n0 = gz(f)
    loc = ip == n0 ? "IN 0-cell" : (ip < n0 ? "above/water side" : "below (sand, +$(ip-n0) cells)")
    @printf("ncells=%d  dz=%.4f m  (dz/δ=%.1f)  0-cell idx=%d z0=%+.4f  eps0=%.3f\n",
            nc, dz, dz / f0.sand_roughness, n0, z[n0], poro[n0])
    @printf("   PEAK phib=%.4f at idx=%d z=%+.4f  → %s\n", maximum(phib), ip, z[ip], loc)
    for i in max(1, n0 - 2):min(length(z), n0 + 3)
        @printf("     z=%+.4f  eps=%.3f  sandfrac(1-eps)=%.3f  phib=%.4f  flowPAT=%.2e\n",
                z[i], poro[i], 1 - poro[i], phib[i], concentration(r, "PAT", :flowing)[i, end])
    end
    println()
end
