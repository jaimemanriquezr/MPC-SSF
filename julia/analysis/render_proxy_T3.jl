# Run the faithful proxy (modelLund + implicit osmosis) to T=3, save all frames to
# CSV, and render depth-vs-time profiles/heatmaps. Run with a project that has both
# MPCSSF (dev) and CairoMakie:
#   julia --project=<plotenv> julia/analysis/render_proxy_T3.jl <outdir>
using MPCSSF, CairoMakie, DelimitedFiles, Printf

include(joinpath(@__DIR__, "proxy_model.jl"))

const INFL = Float64[2.68e-3,1.00e-2,0.0,5.36e-3,9.10e-3,6.23e-3,2.00e-5,0.0,1.75e-4]
const SPECIES = ("HET","PHO","POM","PAT","O2","IC","NH4","HPO4","DOM")

outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results", "proxy_T3")
isdir(outdir) || mkpath(outdir)

f = addgridpoints(SandFilter(temperature=15.0), 30)
@info "running proxy to T=3"
r = run_proxy(State(f, modelLund()); simulation_time=3.0, inflow_concentrations=INFL,
              n_frames=30, quiet=true)
@info "done" flag=string(r.flag) t_end=r.time_final steps=length(r.simulation_data[:step_times])-1

# ---- save frames to CSV ------------------------------------------------------
ts = times(r); z = depths(r); vf = get_volume_fractions(r)
writedlm(joinpath(outdir, "times.csv"), ts, ',')
writedlm(joinpath(outdir, "depths.csv"), z, ',')
writedlm(joinpath(outdir, "phi_biofilm.csv"), vf.biofilm, ',')
writedlm(joinpath(outdir, "phi_matrix.csv"), vf.matrix, ',')
writedlm(joinpath(outdir, "phi_enclosed.csv"), vf.enclosed, ',')
for nm in SPECIES, reg in (:matrix, :enclosed, :flowing)
    writedlm(joinpath(outdir, "conc_$(nm)_$(reg).csv"), concentration(r, nm, reg), ',')
end
@info "CSV frames saved" outdir nfiles=length(readdir(outdir))

# ---- render profiles/heatmaps ------------------------------------------------
save(joinpath(outdir, "phib_heatmap.png"),
     let fig=Figure(size=(720,520)); ax=Axis(fig[1,1]; xlabel="time (day)", ylabel="depth z (m)",
        yreversed=true, title="biofilm fraction φ_b (proxy, T=3)");
        hm=heatmap!(ax, ts, z, permutedims(vf.biofilm); colormap=:viridis);
        Colorbar(fig[1,2], hm; label="φ_b"); fig end)
save(joinpath(outdir, "HET_matrix_heatmap.png"), plot_concentration_heatmap(r, "HET", :matrix))
save(joinpath(outdir, "HET_matrix_profiles.png"), plot_concentration(r, "HET", :matrix))
save(joinpath(outdir, "PHO_matrix_profiles.png"), plot_concentration(r, "PHO", :matrix))
save(joinpath(outdir, "volume_fractions_final.png"), plot_volume_fractions(r))
save(joinpath(outdir, "PAT_flowing_heatmap.png"), plot_concentration_heatmap(r, "PAT", :flowing))
@info "render complete" pngs=length(filter(x->endswith(x,".png"), readdir(outdir)))
