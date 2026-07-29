# Global sensitivity SCREENING (Morris elementary effects) for the physical
# slow-sand pathogen model, following the UQ plan in
# Sensitivity_Analysis_UQ_SSF_model.pdf and drawing parameter ranges from
# Campos2006 (biological rates, temperature factors, half-saturations) and
# Schijven2013 (attachment / inactivation coefficients and their ranges).
#
# Stage 1 of the plan: a broad Morris screen over ~20-30 uncertain quantities to
# rank importance (µ*) and flag nonlinearity/interaction (σ) BEFORE the expensive
# variance-based (Sobol) stage on the retained 8-12 parameters.
#
# Feasibility: each evaluation uses the fast faithful proxy (implicit_osmosis, see
# proxy_model.jl) so r·(k+1) runs is workstation-scale. For k=24, r=20 that is
# 500 runs; at ~15-30 s/run a full screen is a few hours (run on a cluster or
# overnight for the production screen; the CLI exposes small smoke sizes).
#
# Run:
#   julia --project=julia julia/analysis/morris_screening.jl [r_trajectories] [tsim] [ncells]
# Writes julia/analysis/results/morris/{elementary_effects,mu_sigma}.csv and prints
# a ranked table per quantity of interest.

using MPCSSF
using Printf, DelimitedFiles, Random, Statistics

include(joinpath(@__DIR__, "proxy_model.jl"))
include(joinpath(@__DIR__, "pathogen_repro.jl"))   # physical pathogen_model()

# ============================================================================
# I.  Base case + forcing
# ============================================================================
# Manriquez2026 Table B.1 influent [HET PHO POM PAT | O2 IC NH4 HPO4 DOM].
const INFLUENT_BASE = Float64[2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]
const PAT_IN = INFLUENT_BASE[4]

# ============================================================================
# II. Uncertain parameters
# ============================================================================
# Each parameter maps a unit coordinate u∈[0,1] to (a) a log-uniform value over
# [lo,hi] and (b) a builder that applies it. `scale=:log` samples geometrically,
# `:lin` arithmetically. `source` records provenance (paper/table) for the range.
# Ranges tagged TODO are provisional (factor ~3 around the modelLund nominal) and
# are refined from the Campos2006 / Schijven2013 extracted CSVs.
#
# apply signature: (filter, model, value) -> (filter', model').  Filters are
# mutable (mutate + return); models are immutable (rebuilt via _remake helpers).

struct Param
    name::String
    block::String
    lo::Float64
    hi::Float64
    scale::Symbol            # :log | :lin
    source::String
    apply::Function          # (f, m, v) -> (f, m)
end

_unit_to_val(p::Param, u) = p.scale === :log ?
    p.lo * (p.hi / p.lo)^u : p.lo + (p.hi - p.lo) * u

# ---- model "remake" helpers (immutable structs) ----------------------------
_rx(r::Reaction; kw...) = Reaction(; name=r.name,
    (; nominal_rate=r.nominal_rate, temperature_correction_factor=r.temperature_correction_factor,
       order=r.order, half_saturation_constants=r.half_saturation_constants,
       stoichiometric_coefficients=r.stoichiometric_coefficients,
       efficiency_biofilm=r.efficiency_biofilm, efficiency_flowing=r.efficiency_flowing,
       is_light_dependent=r.is_light_dependent, minimum_light_factor=r.minimum_light_factor,
       optimal_light_factor=r.optimal_light_factor, kw...)...)
_pt(p::Particle; kw...) = Particle(; name=p.name,
    (; density=p.density, dispersivity=p.dispersivity, transport_rate=p.transport_rate,
       attachment_matrix=p.attachment_matrix, attachment_sand=p.attachment_sand,
       attenuation=p.attenuation, sand_attachment_factor=p.sand_attachment_factor, kw...)...)
_mdl(m::Model; kw...) = Model(; (; components=m.components, reactions=m.reactions,
    cohesion_submodel=m.cohesion_submodel, water_density=m.water_density,
    biofilm_porosity=m.biofilm_porosity, osmosis_rate=m.osmosis_rate, detachment=m.detachment, kw...)...)

# find reaction index by name (physical pathogen model reaction set)
_rxi(m, nm) = findfirst(r -> r.name == nm, m.reactions)

# Builders take the sampled value `v` as the ABSOLUTE parameter value (log-uniform
# over [lo,hi], with the nominal inside the range) — cleaner for the later Sobol
# stage than multiplicative steps. -------------------------------------------
b_rate(names) = (f, m, v) -> (f, _mdl(m; reactions=Reaction[
    (r.name in names) ? _rx(r; nominal_rate=v) : r for r in m.reactions]))
b_theta(names) = (f, m, v) -> (f, _mdl(m; reactions=Reaction[
    (r.name in names) ? _rx(r; temperature_correction_factor=v) : r for r in m.reactions]))
b_halfsat(rxname, key) = (f, m, v) -> (f, _mdl(m; reactions=Reaction[
    (r.name == rxname) ? _rx(r; half_saturation_constants=merge(r.half_saturation_constants, Dict(key=>v))) : r
    for r in m.reactions]))
b_attach_sand() = (f, m, v) -> (f, _mdl(m; components=Component[
    c isa Particle ? _pt(c; attachment_sand=v) : c for c in m.components]))
b_sand_pat() = (f, m, v) -> (f, _mdl(m; components=Component[
    (c isa Particle && c.name=="PAT") ? _pt(c; sand_attachment_factor=v) : c for c in m.components]))
b_dispersivity() = (f, m, v) -> (f, _mdl(m; components=Component[
    c isa Particle ? _pt(c; dispersivity=v) : c for c in m.components]))
b_transport(kind) = (f, m, v) -> (f, _mdl(m; components=Component[
    (c isa Particle && kind==:P) ? _pt(c; transport_rate=v) :
    (c isa Liquid   && kind==:L) ? Liquid(; name=c.name, density=c.density, dispersivity=c.dispersivity,
                                          transport_rate=v) : c
    for c in m.components]))
b_beta() = (f, m, v) -> (f, _mdl(m; biofilm_porosity=v))
b_zeta0() = (f, m, v) -> (f, _mdl(m; cohesion_submodel=CahnHilliardModel(
    kappa=m.cohesion_submodel.kappa, zeta_0=v, zeta_1=m.cohesion_submodel.zeta_1)))
b_kappa() = (f, m, v) -> (f, _mdl(m; cohesion_submodel=CahnHilliardModel(
    kappa=v, zeta_0=m.cohesion_submodel.zeta_0, zeta_1=m.cohesion_submodel.zeta_1)))
b_detach_scale(nom) = (f, m, v) -> (f, _mdl(m; detachment = w -> (v/nom) .* m.detachment(w)))
b_velocity() = (f, m, v) -> (f.inflow_velocity = v; (f, m))
b_temp() = (f, m, v) -> (f.temperature = v; (f, m))
b_light_sand() = (f, m, v) -> (f.light_attenuation_sand = v; (f, m))
b_influent_pat() = (f, m, v) -> (f, m)   # PAT influent handled in run (INFLUENT vector)

# The parameter table is built from the running model's nominals so [lo,hi]
# brackets the nominal; ranges/sources tagged from the two reference papers.
function build_params(m0, f0)
    nom(nm) = (i=_rxi(m0, nm); i===nothing ? NaN : m0.reactions[i].nominal_rate)
    P = Particle[c for c in m0.components if c isa Particle]
    part = P[1]
    v0 = f0.inflow_velocity
    ps = Param[]
    add(args...) = push!(ps, Param(args...))

    # Ranges are log-uniform (positive params) or linear, bracketing the model
    # nominal; `source` cites the paper table + value the range is drawn from.
    # Campos2006 rates are h⁻¹ → ×24 for day⁻¹ (the model's unit).

    # --- Block: Operating & forcing --------------------------------------
    add("velocity",      "forcing", 2.16, 21.6, :log, "Schijven2013 T1 filtration-rate spread; nom q=7.2 m/d", b_velocity())
    add("temperature",   "forcing", 3.0,  25.0, :lin, "Campos2006/Manriquez seasonal 3-25 °C", b_temp())
    add("influent_PAT",  "forcing", 1.6e-3, 1.6e-2, :log, "Manriquez2026 Table B.1 PAT=5.36e-3 ±3×", b_influent_pat())

    # --- Block: Transport & retention ------------------------------------
    add("dispersivity",  "transport", 3.6e-3, 3.6e-2, :log, "Schijven2013 T1 dispersivity; nom 0.012 ±3×", b_dispersivity())
    add("transport_P",   "transport", 1.64, 16.4, :log, "Lund phase-transfer; nom 5.47 ±3×", b_transport(:P))
    add("attach_sand",   "transport", 180.0, 1640.0, :log, "Diehl2025/Lund b_sand=547 ±3×", b_attach_sand())
    add("sand_pathogen", "transport", 5e-3, 0.71, :log, "Schijven2013 T4 sticking efficiency α (MS2 5.6e-3–ECWR1 0.71)", b_sand_pat())

    # --- Block: Biological kinetics (Campos2006 Table 3 rates + Table 2 θ) -
    add("mu_HET",        "kinetics", 6e-3, 5.4e-2, :log, "Campos2006 kgmaxa 0.70–1.0/d (rel. spread); nom 0.018", b_rate(["Heterotroph growth"]))
    add("mu_PHO",        "kinetics", 1.65, 16.5, :log, "Campos2006 kgmaxb algae; nom 5.5 ±3×", b_rate(["Phototroph growth"]))
    add("d_HET",         "kinetics", 0.05, 4.08, :log, "Campos2006 T3 kdb 0.0021–0.17/h ×24 = 0.05–4.08/d", b_rate(["Heterotroph death"]))
    add("hydrolysis",    "kinetics", 0.03, 0.15, :log, "Campos2006 T3 kh 0.0028–0.0034/h ×24≈0.067–0.082; nom 0.09", b_rate(["Hydrolysis"]))
    add("theta_growth",  "kinetics", 1.02, 1.09, :lin, "Campos2006 T2 θ_kgb/kgp=1.047", b_theta(["Heterotroph growth","Phototroph growth"]))
    add("theta_death",   "kinetics", 1.03, 1.12, :lin, "Campos2006 T2 θ_kga=1.066, θ_kr*=1.08", b_theta(["Heterotroph death","Phototroph death"]))
    add("K_O2_HET",      "kinetics", 1e-3, 1e-2, :log, "O2 half-sat; nom 3e-3 (Reichert)", b_halfsat("Heterotroph growth","O2"))
    add("K_DOM_HET",     "kinetics", 5e-5, 1e-3, :log, "Campos2006 T3 ksCd 0.80–0.98 mgC/L; nom 2e-4", b_halfsat("Heterotroph growth","DOM"))
    add("K_HPO4_HET",    "kinetics", 1e-8, 5e-5, :log, "Campos2006 T3 ksp 0.001–0.05 mgP/L=1e-6–5e-5; nom 1.4e-8 (P-limiting, influent HPO4=0)", b_halfsat("Heterotroph growth","HPO4"))

    # --- Block: Pathogen/marker kinetics (Schijven2013 + Manriquez B.4) ---
    add("marker_growth", "pathogen", 0.05, 0.6,  :log, "Manriquez B.4 μ̂_PAT=0.2 ±3×", b_rate(["MarkerGrowth"]))
    add("inactivation",  "pathogen", 5e-3, 1.84, :log, "Schijven2013 T2/T3 µl,µs 0–1.84/d; Manriquez B.4 0.02", b_rate(["Inactivation"]))
    add("bacterivory",   "pathogen", 2.0,  20.0, :log, "Manriquez B.4 p̂_PAT=8.0 ±2.5×", b_rate(["Bacterivory"]))

    # --- Block: Biofilm mechanics ----------------------------------------
    add("beta_porosity", "biofilm", 0.95, 0.995, :lin, "Diehl2025/Lund biofilm porosity (0.99)", b_beta())
    add("zeta_0",        "biofilm", 1e1, 1e3, :log, "Cahn-Hilliard cohesion; nom 1e2", b_zeta0())
    add("detach_scale",  "biofilm", 0.3, 3.0, :log, "detachment-law scale (relative to nominal)", b_detach_scale(1.0))

    return ps
end

# ============================================================================
# III. Quantities of interest (UQ plan §1)
# ============================================================================
# L(t)=log10(cPAT_in/cPAT_out); scalar summaries L̄, Lmin; biofilm mass Mb and
# centroid zb; O2 penetration depth; max HET/PHO; total attached PAT.
function qois(r, f)
    z = depths(r); dz = gridsize(f); poro = computeporosity(f, z)
    ts = times(r)
    si = findfirst(>=(0.0), z)                          # first sand cell
    flowPAT = concentration(r, "PAT", :flowing)         # N × nframes
    cout = max.(flowPAT[end, :], 1e-30)                 # effluent (bottom cell)
    Lt = log10.(PAT_IN ./ cout)
    Lmean = mean(Lt); Lmin = minimum(Lt)
    vf = get_volume_fractions(r)
    phib = vf.biofilm[:, end]
    Mb = sum(poro .* phib .* dz)                         # biofilm mass proxy
    zb = sum(z .* phib) / (sum(phib) + 1e-30)           # centroid
    matHET = concentration(r, "HET", :matrix)[:, end]
    o2 = concentration(r, "O2", :flowing)[:, end]
    surf = si === nothing ? 0.0 : matHET[si]
    peni = surf > 0 ? findlast(>(0.01 * surf), matHET) : nothing
    # o2pen is CENSORED whenever the filter stays aerobic to the bottom, which is
    # the normal case here: measured at ncells=30 the profile only falls to 60% of
    # its maximum, so the 1% threshold is crossed in 0 of 62 cells and `findlast`
    # returns the last index for every parameter set. Reported as a depth that
    # would silently read as "oxygen reaches the bottom" data when it is really
    # "the metric never resolved". NaN so it cannot be averaged into a ranking.
    o2max = maximum(o2); o2thr = 0.01 * o2max + 1e-30
    o2pen_i = findlast(>(o2thr), o2)
    o2_censored = o2pen_i === nothing || o2pen_i == length(o2)
    attPAT = sum(poro .* concentration(r, "PAT", :matrix)[:, end] .* dz)
    return (; Lmean, Lmin, Mb, zb,
            pen = peni === nothing ? 0.0 : z[peni],
            o2pen = o2_censored ? NaN : z[o2pen_i],
            # Uncensored companion: the fraction of oxygen consumed across the
            # domain. Varies continuously whether or not the filter goes anoxic,
            # so it carries the signal o2pen was meant to and cannot saturate.
            o2dep = o2max > 0 ? 1 - minimum(o2) / o2max : 0.0,
            maxHET = maximum(matHET), maxPHO = maximum(concentration(r, "PHO", :matrix)[:, end]),
            attPAT, flag = string(r.flag))
end
const QOI_NAMES = ["Lmean","Lmin","Mb","zb","pen","o2pen","o2dep","maxHET","maxPHO","attPAT"]

# ============================================================================
# IV. Model evaluation at a unit-cube point
# ============================================================================
function evaluate(u::AbstractVector, params::Vector{Param}, f0build, m0; tsim, ncells)
    f = f0build()
    m = m0
    infl = copy(INFLUENT_BASE)
    for (p, ui) in zip(params, u)
        v = _unit_to_val(p, ui)
        if p.name == "influent_PAT"
            infl[4] = v
        else
            f, m = p.apply(f, m, v)
        end
    end
    s = State(f, m)
    r = run_proxy(s; simulation_time=tsim, inflow_concentrations=infl,
                  n_frames=12, quiet=true, adaptive_max_dt=3e-6)
    return qois(r, f)
end

# ============================================================================
# V. Morris trajectory design (elementary effects)
# ============================================================================
# Standard Morris (1991) trajectories on a p-level grid, Δ = p/(2(p-1)).
function morris_trajectory(k, p, rng)
    Δ = p / (2 * (p - 1))
    xstar = rand(rng, 0:(p-1), k) ./ (p - 1)            # random base point on grid
    xstar = clamp.(xstar, 0.0, 1.0 - Δ)                 # room to step +Δ
    perm = randperm(rng, k)                             # order of coordinate changes
    dirs = rand(rng, (-1, 1), k)
    B = zeros(k + 1, k)
    B[1, :] = xstar
    x = copy(xstar)
    for (row, j) in enumerate(perm)
        x[j] = clamp(x[j] + dirs[j] * Δ, 0.0, 1.0)
        B[row + 1, :] = x
    end
    return B, perm, Δ, dirs
end

# ============================================================================
# VI. Driver
# ============================================================================
"µ* (mean |EE|) and σ per parameter per QoI, from the effects accumulated so far."
function _mu_sigma_rows(EE, params)
    rows = Vector{Any}[["qoi","param","block","source","mu_star","sigma","n"]]
    for qn in QOI_NAMES, (j, p) in enumerate(params)
        e = EE[qn][j]
        mustar = isempty(e) ? NaN : mean(abs.(e))
        sig = length(e) < 2 ? NaN : std(e)
        push!(rows, [qn, p.name, p.block, p.source, mustar, sig, length(e)])
    end
    return rows
end

function main(; r_traj=8, tsim=0.3, ncells=30, p_levels=6, seed=20260721)
    rng = MersenneTwister(seed)
    m0 = pathogen_model()
    f0build() = addgridpoints(SandFilter(temperature=15.0), ncells)
    params = build_params(m0, f0build())
    k = length(params)
    outdir = joinpath(@__DIR__, "results", "morris"); isdir(outdir) || mkpath(outdir)

    @info "Morris screening" k_params=k r_traj tsim ncells runs=r_traj*(k+1)
    println("Parameters (", k, "), blocks: ", join(unique(getfield.(params, :block)), ", "))

    # EE[qoi][param] accumulates elementary effects across trajectories
    EE = Dict(q => [Float64[] for _ in 1:k] for q in QOI_NAMES)
    nrun = 0; nfail = 0
    for t in 1:r_traj
        B, perm, Δ, dirs = morris_trajectory(k, p_levels, rng)
        ys = Vector{Any}(undef, k + 1)
        for row in 1:(k + 1)
            q = evaluate(B[row, :], params, f0build, m0; tsim, ncells)
            ys[row] = q; nrun += 1
            q.flag == "OK" || (nfail += 1)
        end
        for (row, j) in enumerate(perm)
            y0 = ys[row]; y1 = ys[row + 1]
            step = dirs[j] * Δ
            for qn in QOI_NAMES
                ee = (getfield(y1, Symbol(qn)) - getfield(y0, Symbol(qn))) / step
                isfinite(ee) && push!(EE[qn][j], ee)
            end
        end
        @printf("  trajectory %d/%d done (%d runs, %d non-OK)\n", t, r_traj, nrun, nfail)
        flush(stdout)   # Julia buffers under Slurm redirection; without this the log
                        # stays empty for hours and progress is unreadable.
        # Trajectories are independent, so effects from the completed ones are a
        # valid — if noisier — screen on their own. Written to a SEPARATE file:
        # mu_sigma.csv must keep the last COMPLETE run rather than be replaced by
        # a truncated one. Same reasoning as the Sobol partial writes.
        writedlm(joinpath(outdir, "mu_sigma_partial.csv"), _mu_sigma_rows(EE, params), ',')
    end

    rows = _mu_sigma_rows(EE, params)
    writedlm(joinpath(outdir, "mu_sigma.csv"), rows, ',')
    println("\nwrote ", joinpath(outdir, "mu_sigma.csv"), "  (", nrun, " runs, ", nfail, " non-OK)")

    # ranked table for the headline QoI (mean log-removal)
    for qn in ("Lmean", "Mb")
        println("\nTop parameters by µ* for QoI = $qn:")
        ranked = sort([(p.name, p.block, mean(abs.(EE[qn][j]))) for (j, p) in enumerate(params) if !isempty(EE[qn][j])];
                      by=x->x[3], rev=true)
        for (nm, bl, ms) in ranked[1:min(10, length(ranked))]
            @printf("  µ*=%.3e  %-16s (%s)\n", ms, nm, bl)
        end
    end
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    r_traj = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 8
    tsim   = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.3
    ncells = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 30
    main(; r_traj, tsim, ncells)
end
