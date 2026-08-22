# Cluster-scale reproduction of the Manriquez2026 BIOFILM figures (Fig 7-12) with
# the (corrected) Julia port. These runs are stiff — the biofilm model needs
# dt ~ 1e-7 d (Diehl2025 §3.7), so a 90-day/500-cell run is ~1e9 steps: run on a
# cluster (one experiment per job). Reproduces the paper setup (srun_biofilm.m +
# batch_*.m + Table B.1) EXCEPT it uses the fixed temperature/hydrolysis code, so
# results are the corrected versions (compare against the paper to see the bug
# impact — most visible in the winter/seasons runs).
#
# Usage (one experiment per invocation):
#   julia --project=julia julia/analysis/biofilm_cluster.jl <experiment> <outdir> [ncells] [tsim_days]
#   experiments: summer winter covered uncovered growth scrape0 scrape5 scrape15 scrape25 scrape50
# Writes frames to <outdir>/<experiment>_{times,depths,conc_<name>_<region>,phib}.csv
#
# For the temperature-bug delta, run each with `BUGGY=1` in the env to use the
# buggy θ^(293−T_K) form (via the identity buggy(T)=correct(40−T); see below).

using SSF
using DelimitedFiles, Printf

# ---- reconstructed biofilm model (modelLund + srun_biofilm.m overrides) -------
# modelLund already matches Manriquez2026 App. B (θ=1.047/1.066/1.080, bsand=5.47e2,
# bP=5.47, h20=9e-2, half-sats, ρ=1117). Overrides: detachment 0.14·√(|v|/18),
# ζ0=1e2. Osmosis τ=1e-7 kept (relaxing it does NOT help — φb drift spikes the
# dispersion CFL weight).
function biofilm_model()
    m = modelLund()
    coh = CahnHilliardModel(kappa=m.cohesion_submodel.kappa, zeta_0=1e2,
                            zeta_1=m.cohesion_submodel.zeta_1)
    Model(m.components, m.reactions; cohesion_submodel=coh, water_density=m.water_density,
          biofilm_porosity=m.biofilm_porosity, osmosis_rate=m.osmosis_rate,
          detachment=(v -> 0.14 .* sqrt.(abs.(v) ./ 18.0)))
end

# Table B.1 influent [HET PHO POM PAT | O2 IC NH4 HPO4 DOM] (kg/m³)
const INFLUENT = Float64[2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 1.0e-5, 1.75e-4]

# ---- light forcings (batch_seasons.m / batch_covering.m) ----------------------
light_uncovered(t) = max(0.5*(sin(2π*(t-0.3)) + 1) - 0.2, 0.0)     # summer/base
light_covered(t)   = 0.01 * light_uncovered(t)                     # 1% (covered)
light_summer(t)    = light_uncovered(t)
light_winter(t)    = max(0.5*(sin(2π*(t-0.2)) + 1) - 0.4, 0.0)

# ---- scrape (port of the fixed scrape.m; refill_sand=true variant) ------------
# Zeros biofilm + enclosed water in the scraped sand cells (0 ≤ z ≤ scraped depth)
# and cleans the flowing suspension. Grid is not shifted (refill-with-fresh-sand
# convention); the shift variant would rebuild the grid.
function scrape!(s::State, z_s::Real)
    f = s.filter; dz = gridsize(f); z = f.grid.centers
    scraped_m = ceil(z_s/dz + 1/2) * sign(z_s) * dz
    cells = (z .>= 0) .& (z .<= scraped_m) .& (z_s > 0)
    gc = s.global_concentration
    gc.matrix[cells, :] .= 0; gc.enclosed_particles[cells, :] .= 0
    gc.enclosed_liquids[cells, :] .= 0; s.enclosed_water_volume[cells] .= 0
    gc.flowing_particles .= 0; gc.flowing_liquids .= 0
    return s
end

# ---- run + save --------------------------------------------------------------
function save_run(r, outdir, tag)
    isdir(outdir) || mkpath(outdir)
    writedlm(joinpath(outdir, "$(tag)_times.csv"), times(r), ',')
    writedlm(joinpath(outdir, "$(tag)_depths.csv"), depths(r), ',')
    writedlm(joinpath(outdir, "$(tag)_phib.csv"), get_volume_fractions(r).biofilm, ',')
    for nm in ("HET","PHO","POM","PAT","O2","IC","NH4","HPO4","DOM"), reg in (:matrix,:flowing)
        writedlm(joinpath(outdir, "$(tag)_conc_$(nm)_$(reg).csv"), concentration(r, nm, reg), ',')
    end
    open(joinpath(outdir, "$(tag)_meta.txt"), "w") do io
        println(io, "flag=", r.flag); println(io, "t_final=", r.time_final)
        println(io, "steps=", length(r.simulation_data[:step_times])-1)
    end
end

function main(exp, outdir, ncells, tsim)
    buggy = get(ENV, "BUGGY", "0") == "1"   # buggy θ^(293−T): use temp 40−T
    Tcorr(T) = buggy ? 40 - T : T           # summer 19→21, winter 3→37 under BUGGY
    f = addgridpoints(SandFilter(temperature=20.0), ncells)   # default 20°C, overridden below
    m = biofilm_model()

    if exp == "summer"
        f.light_irradiation = light_summer; f.temperature = Tcorr(19.0)
        s = State(f, m)
    elseif exp == "winter"
        f.light_irradiation = light_winter;  f.temperature = Tcorr(3.0)
        s = State(f, m)
    elseif exp in ("covered", "uncovered", "growth")
        f.light_irradiation = exp == "covered" ? light_covered : light_uncovered
        f.temperature = Tcorr(20.0); s = State(f, m)
    elseif startswith(exp, "scrape")
        depth_cm = parse(Int, exp[7:end])
        # scraping resumes from a mature filter: expects a prior `growth` run's
        # final state saved as JLD, OR regenerate here (long). Placeholder: start
        # from clean + scrape (user should wire the mature IC from a growth run).
        f.light_irradiation = light_uncovered; f.temperature = Tcorr(20.0)
        s = State(f, m); scrape!(s, depth_cm * 1e-2)
    else
        error("unknown experiment $exp")
    end

    @info "biofilm run" experiment=exp ncells tsim temperature=f.temperature buggy
    r = simulate(s; inflow_concentrations=INFLUENT, simulation_time=tsim,
                 time_step=:adaptive, cfl_factor=0.99, adaptive_initial_dt=1e-8,
                 adaptive_max_dt=1e-5, n_frames=200, quiet=false)
    save_run(r, outdir, exp * (buggy ? "_buggy" : ""))
    @info "done" flag=r.flag steps=length(r.simulation_data[:step_times])-1
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: biofilm_cluster.jl <experiment> <outdir> [ncells=500] [tsim_days]")
    exp = ARGS[1]; outdir = ARGS[2]
    ncells = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 500
    tsim = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) :
           (exp in ("summer","winter") ? 90.0 : 30.0)
    main(exp, outdir, ncells, tsim)
end
