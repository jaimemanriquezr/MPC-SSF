# Julia mirror of analysis/manuscriptExperiments.m (matlab-claude): the
# Manriquez2026 manuscript experiments on the corrected model, for
# cross-implementation comparison against the MATLAB runs. Same protocols
# (HEAD-legacy), same CSV layout, results under julia/analysis/results/manuscript/.
#
# Run:  julia --project=julia julia/analysis/manuscript_experiments.jl E3 [E4 ...] \
#              [--ncells=100] [--maxdt=3e-6] [--smoke]
# IDs: E1 seasons · E2 covered · E3 long-term(+X1, builds the mature cache)
#      E4 scraping · E5 pulses · E6 controls · E10 scrape×marker · X1/E7/E8/E9 aliases
#
# Deviations and provenance: analysis/PARAMETERS.md (worktree) and the
# 2026-08-19 manuscript-driver plan. Model: pathogen_model(phototroph_respiration
# = 0.55, pg_excess = true), κ = 1e-7 (article Table 2), detachment
# 0.14·√(v/18) for all families (2026-08-20 decision), T nominal 15 °C
# (seasons 19/3 °C).

using SSF
using DelimitedFiles, Serialization, Printf, Statistics

include(joinpath(@__DIR__, "pathogen_repro.jl"))   # pathogen_model()

const OUTROOT = joinpath(@__DIR__, "results", "manuscript")

light_summer(t) = max(0.5 * (sin(2π * (t - 0.3)) + 1) - 0.2, 0.0)
light_winter(t) = max(0.5 * (sin(2π * (t - 0.2)) + 1) - 0.4, 0.0)
base_influent(pat) = Float64[2.68e-3, 1.00e-2, 0.0, pat, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]

durations(smoke) = smoke ?
    (long=0.2, matureAt=0.1, covered=0.1, postScrape=0.05,
     patStart=0.1, pulseOn=0.12, pulseOff=0.14, patLen=0.08) :
    (long=90.0, matureAt=30.0, covered=30.0, postScrape=30.0,
     patStart=30.0, pulseOn=30.0, pulseOff=32.0, patLen=10.0)

function the_model()
    m = pathogen_model(phototroph_respiration=0.55, pg_excess=true, normalized_light=true)
    coh = CahnHilliardModel(kappa=1e-7, zeta_0=m.cohesion_submodel.zeta_0,
                            zeta_1=m.cohesion_submodel.zeta_1)
    Model(m.components, m.reactions; cohesion_submodel=coh,
          water_density=m.water_density, biofilm_porosity=m.biofilm_porosity,
          osmosis_rate=m.osmosis_rate,
          detachment=(v -> 0.14 .* sqrt.(abs.(v) ./ 18)))
end

function the_filter(ncells, tempC, light)
    f = addgridpoints(SandFilter(temperature=tempC), ncells)
    f.light_irradiation = light
    return f
end

function run_sim(s, tsim, infl, nframes, maxdt)
    r = simulate(s; inflow_concentrations=infl, simulation_time=tsim,
                 time_step=:adaptive, adaptive_initial_dt=1e-8, adaptive_max_dt=maxdt,
                 n_frames=nframes, implicit_osmosis=true, quiet=true)
    r.flag == "OK" || error("run failed: flag=$(r.flag)")
    return r
end

biofilm_fraction(r) = begin
    phiB = concentration(r, "Water", :enclosed) ./ mean(l.density for l in liquids(r.model))
    dP = mean(p.density for p in particles(r.model))
    for p in particles(r.model)
        phiB = phiB .+ (concentration(r, p.name, :matrix) .+ concentration(r, p.name, :enclosed)) ./ dP
    end
    for l in liquids(r.model)
        phiB = phiB .+ concentration(r, l.name, :enclosed) ./ mean(x.density for x in liquids(r.model))
    end
    phiB
end

function total_biomass(r, dz)
    mass = zeros(1, length(times(r)))
    for p in particles(r.model)
        mass .+= sum(concentration(r, p.name, :matrix) .+ concentration(r, p.name, :enclosed), dims=1) .* dz
    end
    for l in liquids(r.model)
        mass .+= sum(concentration(r, l.name, :enclosed), dims=1) .* dz
    end
    vec(mass)
end

nearest_idx(ts, d) = argmin(abs.(ts .- d))
ensure_dir(d) = (isdir(d) || mkpath(d); d)

# ---- E3/E1-summer/X1 + mature cache ------------------------------------------
function run_summer(ncells, maxdt, D, smoke)
    outdir = ensure_dir(joinpath(OUTROOT, "E3_longterm"))
    sentinel = joinpath(outdir, "done.txt")
    isfile(sentinel) && return
    nframes = max(round(Int, 24 * D.long), 12)
    f = the_filter(ncells, 19.0, light_summer)
    r = run_sim(State(f, the_model()), D.long, base_influent(0.0), nframes, maxdt)
    centers = f.grid.centers
    dz = gridsize(f)
    ts = times(r)

    liq = ["O2", "IC", "NH4", "HPO4", "DOM"]
    eff = hcat(ts, (concentration(r, nm, :flowing)[end, :] for nm in liq)...)
    writedlm(joinpath(outdir, "effluent_liquids.csv"), eff, ',')

    keep = 1:max(1, round(Int, 24)):length(ts)
    for nm in liq
        field = concentration(r, nm, :flowing)[:, keep]
        writedlm(joinpath(outdir, "field_$nm.csv"),
                 vcat(hcat(0.0, ts[keep]'), hcat(centers, field)), ',')
    end

    phiB = biofilm_fraction(r)
    snaps = unique(min.([1, 2, 5, 10, 20, (D.matureAt .+ [0, 1, 2, 5, 10, 20])...], D.long))
    ks = [nearest_idx(ts, d) for d in snaps]
    writedlm(joinpath(outdir, "phib_snapshots_summer.csv"), hcat(centers, phiB[:, ks]), ',')
    writedlm(joinpath(outdir, "phib_snapshot_days.csv"), snaps, ',')
    writedlm(joinpath(outdir, "total_biomass_summer.csv"), hcat(ts, total_biomass(r, dz)), ',')

    xdir = ensure_dir(joinpath(OUTROOT, "X1_o2_diel"))
    o2 = concentration(r, "O2", :flowing)
    lv = light_summer.(ts)
    lastday = ts .>= (D.long - 1)
    o2l = o2[:, lastday]; lvl = lv[lastday]
    writedlm(joinpath(xdir, "o2_profiles_noon_midnight.csv"),
             hcat(centers, o2l[:, argmax(lvl)], o2l[:, argmin(lvl)]), ',')
    writedlm(joinpath(xdir, "effluent_o2_diel.csv"), hcat(ts, o2[end, :], lv), ',')

    st = final_state(r; t=D.matureAt)
    ensure_dir(joinpath(OUTROOT, "cache"))
    serialize(joinpath(OUTROOT, "cache", "mature30_summer_N$(ncells)$(smoke ? "_smoke" : "").jls"), st)
    write(sentinel, string(time()))
    return
end

function get_mature(ncells, maxdt, D, smoke)
    ensure_dir(joinpath(OUTROOT, "cache"))
    p = joinpath(OUTROOT, "cache", "mature30_summer_N$(ncells)$(smoke ? "_smoke" : "").jls")
    isfile(p) || run_summer(ncells, maxdt, D, smoke)
    deserialize(p)
end

# Re-home a cached State onto a fresh filter/model at time t0.
rehome(f, m, st, t0) = State(f, m, Float64(t0),
    GlobalConcentration(copy(st.global_concentration.matrix),
                        copy(st.global_concentration.enclosed_particles),
                        copy(st.global_concentration.flowing_particles),
                        copy(st.global_concentration.enclosed_liquids),
                        copy(st.global_concentration.flowing_liquids)),
    copy(st.enclosed_water_volume), Velocity(copy(st.velocity.biofilm), zeros(length(st.enclosed_water_volume) + 1)))

function scrape_state(st, zs, dz, centers)
    gc = GlobalConcentration(copy(st.global_concentration.matrix),
                             copy(st.global_concentration.enclosed_particles),
                             zero(st.global_concentration.flowing_particles),
                             copy(st.global_concentration.enclosed_liquids),
                             zero(st.global_concentration.flowing_liquids))
    scraped = ceil(zs / dz + 1 / 2) * sign(zs) * dz
    cells = zs > 0 ? ((centers .>= 0) .& (centers .<= scraped)) : falses(length(centers))
    sup = centers .< 0
    mask = cells .| sup
    gc.matrix[mask, :] .= 0; gc.enclosed_particles[mask, :] .= 0
    gc.enclosed_liquids[mask, :] .= 0
    phiW = copy(st.enclosed_water_volume); phiW[mask] .= 0
    (gc=gc, phiW=phiW)
end

# ---- E1 winter ----------------------------------------------------------------
function e1_winter(ncells, maxdt, D)
    outdir = ensure_dir(joinpath(OUTROOT, "E1_seasons"))
    isfile(joinpath(outdir, "phib_final_winter.csv")) && return
    f = the_filter(ncells, 3.0, light_winter)
    r = run_sim(State(f, the_model()), D.long, base_influent(0.0),
                max(round(Int, 2 * D.long), 6), maxdt)
    writedlm(joinpath(outdir, "phib_final_winter.csv"),
             hcat(f.grid.centers, biofilm_fraction(r)[:, end]), ',')
    writedlm(joinpath(outdir, "total_biomass_winter.csv"),
             hcat(times(r), total_biomass(r, gridsize(f))), ',')
end

# ---- E2 covered ----------------------------------------------------------------
function e2_covered(ncells, maxdt, D)
    outdir = ensure_dir(joinpath(OUTROOT, "E2_covered"))
    isfile(joinpath(outdir, "phib_final_covered.csv")) && return
    for (tag, scale) in (("uncovered", 1.0), ("covered", 0.01))
        f = the_filter(ncells, 19.0, t -> scale * light_summer(t))
        r = run_sim(State(f, the_model()), D.covered, base_influent(0.0),
                    max(round(Int, 2 * D.covered), 6), maxdt)
        writedlm(joinpath(outdir, "phib_final_$tag.csv"),
                 hcat(f.grid.centers, biofilm_fraction(r)[:, end]), ',')
    end
end

# ---- E4 scraping ----------------------------------------------------------------
function e4_scraping(ncells, maxdt, D, smoke)
    outdir = ensure_dir(joinpath(OUTROOT, "E4_scraping"))
    isfile(joinpath(outdir, "total_biomass.csv")) && return
    st = get_mature(ncells, maxdt, D, smoke)
    f = the_filter(ncells, 19.0, light_summer)
    centers = f.grid.centers; dz = gridsize(f)
    depths = smoke ? [0.0, 0.05] : [0.0, 0.05, 0.15, 0.25, 0.50]
    massT = nothing
    for zs in depths
        sc = scrape_state(st, zs, dz, centers)
        s = State(f, the_model(), 0.0, sc.gc, sc.phiW,
                  Velocity(zeros(length(centers) - 1), zeros(length(centers) + 1)))
        r = run_sim(s, D.postScrape, base_influent(0.0), max(round(Int, 2 * D.postScrape), 6), maxdt)
        ts = times(r); phiB = biofilm_fraction(r)
        ks = [nearest_idx(ts, d) for d in unique(min.([0, 1, 2, 5, 10, 20], D.postScrape))]
        tag = "GP$(round(Int, zs * 100))"
        writedlm(joinpath(outdir, "phib_snapshots_$tag.csv"), hcat(centers, phiB[:, ks]), ',')
        block = hcat(ts, total_biomass(r, dz))
        massT = massT === nothing ? block : hcat(massT, block)
    end
    writedlm(joinpath(outdir, "total_biomass.csv"), massT, ',')
end

# ---- E5/E6 marker feeds ----------------------------------------------------------
function run_marker(ncells, maxdt, D, smoke, dirname, variants)
    outdir = ensure_dir(joinpath(OUTROOT, dirname))
    isfile(joinpath(outdir, "done.txt")) && return
    st = get_mature(ncells, maxdt, D, smoke)
    f = the_filter(ncells, 19.0, light_summer)
    patNom = 5.3616e-3
    for v in variants
        base = base_influent(0.0)
        pulsed = base_influent(v == "BigPulse" ? 100patNom : patNom)
        infl, cref =
            v == "Pulse"        ? ((t -> (D.pulseOn <= t < D.pulseOff) ? pulsed : base), patNom) :
            v == "BigPulse"     ? ((t -> (D.pulseOn <= t < D.pulseOff) ? pulsed : base), 100patNom) :
            v == "Clean"        ? (base, NaN) :
            v == "ConstantFeed" ? (base_influent(patNom), patNom) : error("variant $v")
        s = rehome(f, the_model(), st, D.patStart)
        r = run_sim(s, D.patLen, infl, max(round(Int, 144 * D.patLen), 12), maxdt)
        ts = times(r)
        patOut = concentration(r, "PAT", :flowing)[end, :]
        F = log10.(max(cref, eps()) ./ max.(patOut, 1e-30))
        hetF = log10.(2.68e-3 ./ max.(concentration(r, "HET", :flowing)[end, :], 1e-30))
        phoF = log10.(1.00e-2 ./ max.(concentration(r, "PHO", :flowing)[end, :], 1e-30))
        writedlm(joinpath(outdir, "effluent_$v.csv"), hcat(ts, patOut, F, hetF, phoF), ',')
        keep = 1:max(1, round(Int, length(ts) / 48)):length(ts)
        for vol in (:flowing, :matrix, :enclosed)
            field = concentration(r, "PAT", vol)[:, keep]
            writedlm(joinpath(outdir, "patfield_$(v)_$(vol).csv"),
                     vcat(hcat(0.0, ts[keep]'), hcat(f.grid.centers, field)), ',')
        end
    end
    write(joinpath(outdir, "done.txt"), string(time()))
end

# ---- E10 scrape × marker ----------------------------------------------------------
function e10_scrape_pat(ncells, maxdt, D, smoke)
    outdir = ensure_dir(joinpath(OUTROOT, "E10_scrape_pat"))
    isfile(joinpath(outdir, "done.txt")) && return
    st = get_mature(ncells, maxdt, D, smoke)
    f = the_filter(ncells, 19.0, light_summer)
    centers = f.grid.centers; dz = gridsize(f)
    patNom = 5.3616e-3
    depths = smoke ? [0.0, 0.15] : [0.0, 0.15, 0.25, 0.50]
    for zs in depths
        sc = scrape_state(st, zs, dz, centers)
        s = State(f, the_model(), 0.0, sc.gc, sc.phiW,
                  Velocity(zeros(length(centers) - 1), zeros(length(centers) + 1)))
        r = run_sim(s, D.patLen, base_influent(patNom), max(round(Int, 24 * D.patLen), 12), maxdt)
        ts = times(r)
        patOut = concentration(r, "PAT", :flowing)[end, :]
        writedlm(joinpath(outdir, "logremoval_GP$(round(Int, zs*100)).csv"),
                 hcat(ts, patOut, log10.(patNom ./ max.(patOut, 1e-30))), ',')
    end
    write(joinpath(outdir, "done.txt"), string(time()))
end

# ---- dispatcher --------------------------------------------------------------
function main(ids; ncells=100, maxdt=3e-6, smoke=false)
    D = durations(smoke)
    order = ["E3", "E1", "E2", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "X1"]
    ids = ids == ["all"] ? order : sort(ids, by=x -> findfirst(==(x), order))
    for id in ids
        @info "stage" id
        t0 = time()
        id == "E3" && run_summer(ncells, maxdt, D, smoke)
        id == "X1" && run_summer(ncells, maxdt, D, smoke)
        id == "E1" && (run_summer(ncells, maxdt, D, smoke); e1_winter(ncells, maxdt, D))
        id == "E2" && e2_covered(ncells, maxdt, D)
        id == "E4" && e4_scraping(ncells, maxdt, D, smoke)
        id in ("E5", "E7", "E8", "E9") && run_marker(ncells, maxdt, D, smoke, "E5_pulses", ["Pulse", "BigPulse"])
        id in ("E6", "E9") && run_marker(ncells, maxdt, D, smoke, "E6_controls", ["Clean", "ConstantFeed"])
        id == "E10" && e10_scrape_pat(ncells, maxdt, D, smoke)
        @printf("stage %s done in %.1f min\n", id, (time() - t0) / 60)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    ids = String[a for a in ARGS if !startswith(a, "--")]
    isempty(ids) && (ids = ["all"])
    getopt(name, default) = begin
        hit = findfirst(a -> startswith(a, "--$name="), ARGS)
        hit === nothing ? default : parse(Float64, split(ARGS[hit], "=")[2])
    end
    main(ids; ncells=Int(getopt("ncells", 100)), maxdt=getopt("maxdt", 3e-6),
         smoke=any(==("--smoke"), ARGS))
end
