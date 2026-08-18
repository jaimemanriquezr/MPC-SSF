# Reproduction of the Manriquez2026 PATHOGEN figures (Fig 13-17) with the
# (corrected) Julia port. IMPORTANT: the paper's pathogen model is the PHYSICAL
# model (SDparameters_pathogen.mat = modelLund density 1117 + the marker reactions
# r7/r8/r9), NOT the density-1.2 `thesis_model` (that's a non-physical test model
# used only for the golden master). The physical model is as stiff as the biofilm
# model (τ=1e-7, tiny half-sats) ⇒ dt~1e-7 ⇒ these figures are ALSO CLUSTER-SCALE.
# Run this on a cluster like biofilm_cluster.jl. Reproduces srun_pathogen.m +
# batch_inflow/scraping.m. The marker reaction rates below are placeholders (from
# thesis_model); fill in the physical values from Appendix B Table B.4/B.5.
#
#   Fig 13  in/outflow of the marker (PAT) over time
#   Fig 14  PAT concentration vs depth/time, reference parameters
#   Fig 15  PAT concentration with lower inactivation + bacterivory (predation) rate
#   Fig 16  log-removal of HET/PHO in a mature filter during/after a 2-d ×100 feed pulse
#   Fig 17  log-removal of the marker under constant feed at scrape depths 0/5/15/25/50 cm
#
# Run:  julia --project=julia julia/analysis/pathogen_repro.jl <outdir> [ncells] [mature_days]

using MPCSSF
using DelimitedFiles, Printf

# ---- pathogen model + srun_pathogen.m overrides ------------------------------
# modelPathogen() = thesis_model.mat. Overrides: detachment 1.4e-5·√(|v|/18),
# ζ0=1e2, PAT transfer(transport)_rate ×40. Fig 15 lowers inactivation(6) +
# bacterivory(7) nominal rates.
_remake(p::Particle; kw...) = Particle(; name=p.name, (; density=p.density,
    dispersivity=p.dispersivity, transport_rate=p.transport_rate, attachment_matrix=p.attachment_matrix,
    attachment_sand=p.attachment_sand, attenuation=p.attenuation,
    sand_attachment_factor=p.sand_attachment_factor, kw...)...)
_remake(r::Reaction; kw...) = Reaction(; name=r.name, (; nominal_rate=r.nominal_rate,
    temperature_correction_factor=r.temperature_correction_factor, order=r.order,
    half_saturation_constants=r.half_saturation_constants,
    stoichiometric_coefficients=r.stoichiometric_coefficients, efficiency_biofilm=r.efficiency_biofilm,
    efficiency_flowing=r.efficiency_flowing, is_light_dependent=r.is_light_dependent,
    minimum_light_factor=r.minimum_light_factor, optimal_light_factor=r.optimal_light_factor, kw...)...)

# Physical pathogen model = modelLund (ρ=1117) + the THREE marker reactions from
# Manriquez2026 Table B.4: aerobic growth r7 (μ̂_PAT=0.2 d⁻¹), inactivation r8
# (d̂_PAT=0.02 d⁻¹), bacterivory r9 (p̂_PAT=8.0 d⁻¹). Lund reactions made inert in
# the flowing phase. PAT gets transfer ×40 (srun_pathogen) and sand_pathogen.
# CAVEATS: (1) Table B.4 lists NO K_pred (bacterivory Monod half-sat) — kept the
# thesis_model value 2e-3, VERIFY against Section 4 text. (2) The r7 aerobic-growth
# stoichiometry (which liquids it consumes) is not in the tables — modelled here
# like heterotroph growth (O2/NH4/DOM Monod, O2 consumed); confirm vs Section 4's
# stoichiometric matrix. Rates θ: growth 1.047, inactivation/bacterivory 1.08.
function pathogen_model(; growth_mult=1.0, inact_mult=1.0, bact_mult=1.0,
                        growth_rate=0.2, inact_rate=0.02, bact_rate=8.0, kpred=2e-3,
                        water_factor=1e-3, sand_pathogen=0.0, phototroph_respiration=0.0)
    m = modelLund(; phototroph_respiration)
    comps = Component[c isa Particle && c.name == "PAT" ?
                      _remake(c; transport_rate=40*c.transport_rate, sand_attachment_factor=sand_pathogen) : c
                      for c in m.components]
    rxs = Reaction[_remake(r; efficiency_flowing=0.0) for r in m.reactions]   # Lund rxns inert in flowing
    push!(rxs,
        # r7 — marker aerobic growth (APPROX stoichiometry; see caveat)
        Reaction(name="MarkerGrowth", nominal_rate=growth_mult*growth_rate, temperature_correction_factor=1.047,
                 efficiency_flowing=0.0, order=Dict("PAT"=>1.0),
                 half_saturation_constants=Dict("O2"=>3.0e-3, "NH4"=>4.0e-3, "DOM"=>2.0e-4),
                 stoichiometric_coefficients=Dict("PAT"=>1.0, "O2"=>-1.2317, "DOM"=>-1.5873)),
        # r8 — inactivation
        Reaction(name="Inactivation", nominal_rate=inact_mult*inact_rate, temperature_correction_factor=1.08,
                 efficiency_flowing=0.0, order=Dict("PAT"=>1.0), stoichiometric_coefficients=Dict("PAT"=>-1.0)),
        # r9 — bacterivory (HET-promoted predation)
        Reaction(name="Bacterivory", nominal_rate=bact_mult*bact_rate, temperature_correction_factor=1.08,
                 efficiency_flowing=water_factor, order=Dict("PAT"=>1.0),
                 half_saturation_constants=Dict("HET"=>kpred),
                 stoichiometric_coefficients=Dict("PAT"=>-1.0, "HET"=>1.0)))
    coh = CahnHilliardModel(kappa=m.cohesion_submodel.kappa, zeta_0=1e2, zeta_1=m.cohesion_submodel.zeta_1)
    Model(comps, rxs; cohesion_submodel=coh, water_density=m.water_density,
          biofilm_porosity=m.biofilm_porosity, osmosis_rate=m.osmosis_rate,
          detachment=(v -> 1.4e-5 .* sqrt.(abs.(v) ./ 18.0)))
end

# Table B.1 influent (kg/m³)  [HET PHO POM PAT | O2 IC NH4 HPO4 DOM]
# HET 2.68e-3, PHO 1.00e-2, POM 0, PAT/marker 5.36e-3, O2 9.10e-3, IC 6.23e-3,
# NH4 2.00e-5, HPO4 0 (Table B.1: phosphate influent is exactly 0), DOM 1.75e-4.
const INFLUENT = Float64[2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]
light_default(t) = max(0.5*(sin(2π*(t-0.3)) + 1) - 0.2, 0.0)

# ---- pulse inflow (port of pulse_error.m): pulsed components ×factor over a window
function pulse_inflow(; base=INFLUENT, pvec, start, len, factor)
    b = (1 .- pvec) .* base; p = pvec .* base
    return t -> (start <= t < start + len) ? b .+ factor .* p : b .+ p
end

# ---- scrape (fixed scrape.m, refill-sand variant) ----------------------------
function scrape!(s::State, z_s::Real)
    f = s.filter; dz = gridsize(f); z = f.grid.centers
    scraped_m = ceil(z_s/dz + 1/2) * sign(z_s) * dz
    cells = (z .>= 0) .& (z .<= scraped_m) .& (z_s > 0)
    gc = s.global_concentration
    gc.matrix[cells,:] .= 0; gc.enclosed_particles[cells,:] .= 0
    gc.enclosed_liquids[cells,:] .= 0; s.enclosed_water_volume[cells] .= 0
    gc.flowing_particles .= 0; gc.flowing_liquids .= 0; s
end

# ---- log-removal: log10(C_in / C_out) per component, per frame ----------------
function log_removal(r, name; Cin)
    Cout = concentration(r, name, :flowing)[end, :]          # effluent (bottom cell), per frame
    return log10.(max.(Cin, eps()) ./ max.(Cout, eps()))
end

# ---- generate a mature filter (run to `age` days under constant influent) -----
function generate_mature(; ncells, age, model=pathogen_model())
    f = addgridpoints(SandFilter(temperature=20.0), ncells); f.light_irradiation = light_default
    r = simulate(State(f, model); inflow_concentrations=INFLUENT, simulation_time=age,
                 time_step=:adaptive, adaptive_initial_dt=1e-6, adaptive_max_dt=1e-4,
                 n_frames=100, quiet=true)
    return r
end

save(a, outdir, tag) = (isdir(outdir)||mkpath(outdir); writedlm(joinpath(outdir,tag), a, ','))

function main(outdir, ncells, mature_days)
    isdir(outdir) || mkpath(outdir)
    @info "generating mature filter" ncells mature_days
    rm0 = generate_mature(ncells=ncells, age=mature_days)
    @info "mature filter" flag=rm0.flag maxphib=maximum(get_volume_fractions(rm0).biofilm[:,end])
    save(depths(rm0), outdir, "depths.csv")

    # Fig 14 — PAT vs depth/time under constant marker feed, reference params
    s = final_state(rm0); f = s.filter; f.light_irradiation = light_default
    r14 = simulate(s; inflow_concentrations=INFLUENT, simulation_time=10.0, time_step=:adaptive,
                   adaptive_initial_dt=1e-6, adaptive_max_dt=1e-4, n_frames=100, quiet=true)
    save(times(r14), outdir, "fig14_times.csv"); save(concentration(r14,"PAT",:flowing), outdir, "fig14_PAT_flowing.csv")
    save(concentration(r14,"PAT",:matrix), outdir, "fig14_PAT_matrix.csv")
    @info "Fig14 done" flag=r14.flag

    # Fig 13 — in/outflow of the marker over time (from the Fig 14 run)
    save(concentration(r14,"PAT",:flowing)[1,:], outdir, "fig13_PAT_inflowcell.csv")
    save(concentration(r14,"PAT",:flowing)[end,:], outdir, "fig13_PAT_outflow.csv")

    # Fig 15 — lower inactivation + bacterivory
    m15 = pathogen_model(inact_mult=0.1, bact_mult=0.1)
    rm15 = generate_mature(ncells=ncells, age=mature_days, model=m15)
    s15 = final_state(rm15); s15.filter.light_irradiation = light_default
    r15 = simulate(s15; inflow_concentrations=INFLUENT, simulation_time=10.0, time_step=:adaptive,
                   adaptive_initial_dt=1e-6, adaptive_max_dt=1e-4, n_frames=100, quiet=true)
    save(concentration(r15,"PAT",:flowing), outdir, "fig15_PAT_flowing.csv")
    @info "Fig15 done" flag=r15.flag

    # Fig 16 — log-removal of HET/PHO during/after a 2-d ×100 feed pulse (pulse HET,PHO,DOM)
    s16 = final_state(rm0); s16.filter.light_irradiation = light_default
    pulse = pulse_inflow(pvec=Float64[1,1,0,0,0,0,0,0,1], start=0.0, len=2.0, factor=100.0)
    r16 = simulate(s16; inflow_concentrations=pulse, simulation_time=10.0, time_step=:adaptive,
                   adaptive_initial_dt=1e-6, adaptive_max_dt=1e-4, n_frames=200, quiet=true)
    save(times(r16), outdir, "fig16_times.csv")
    save(log_removal(r16,"HET"; Cin=100*INFLUENT[1]), outdir, "fig16_logrem_HET.csv")
    save(log_removal(r16,"PHO"; Cin=100*INFLUENT[2]), outdir, "fig16_logrem_PHO.csv")
    @info "Fig16 done" flag=r16.flag

    # Fig 17 — log-removal of marker under constant feed, scrape depths 0/5/15/25/50 cm
    open(joinpath(outdir,"fig17_logrem_marker.csv"),"w") do io
        for cm in (0,5,15,25,50)
            s17 = final_state(rm0); s17.filter.light_irradiation = light_default; scrape!(s17, cm*1e-2)
            r17 = simulate(s17; inflow_concentrations=INFLUENT, simulation_time=10.0, time_step=:adaptive,
                           adaptive_initial_dt=1e-6, adaptive_max_dt=1e-4, n_frames=50, quiet=true)
            lr = log_removal(r17,"PAT"; Cin=INFLUENT[4])
            println(io, cm, ",", join(lr, ","))
            @info "Fig17" depth_cm=cm flag=r17.flag final_logremoval=lr[end]
        end
    end
    @info "all pathogen figures written" outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    outdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "pathogen_figs")
    ncells = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 60
    mature_days = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 30.0
    main(outdir, ncells, mature_days)
end
