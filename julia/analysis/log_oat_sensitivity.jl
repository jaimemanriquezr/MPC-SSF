# Logarithmic one-at-a-time (OAT) sensitivity analysis for the physical slow-sand
# PATHOGEN model, per the updated instructions in
# Logarithmic_OAT_Sensitivity_SSF_model.pdf.
#
# Protocol:
#   * Build ONE mature-filter state (nominal model, run to maturity). Every
#     perturbed run starts from that same state, so we isolate the sensitivity of
#     the RESPONSE to a disturbance (not of the mature state itself).
#   * Disturbance: an influent pathogen PULSE (a challenge/seeding event), like the
#     Schijven2013 experiments. Removal L(t) = log10(C_ref / c_PAT_out(t)).
#   * Perturb each parameter θ_i by ×2 and ×1/2 (log-symmetric). Bounded params
#     (θ≈1.05, β=0.99) are perturbed on their DEVIATION (θ−1, 1−β).
#
# Sensitivity measures (all over the post-disturbance window [t_d, T]):
#   s_i(t)  = (L+_i(t) − L−_i(t)) / (2 ln 2)  ≈ ∂L/∂ln θ_i
#   I_i     = sqrt(mean s_i(t)^2)                        (RMS — primary ranking)
#   I_i^max = max |s_i(t)|                               (max instantaneous)
#   D_i^min = max(L+min,i − Lmin,0 , L−min,i − Lmin,0)   (effect on worst removal)
#   I_i^min = (L+min,i − L−min,i)/(2 ln 2)
#   A_i     = ||L+ + L− − 2 L0|| / ||L+ − L−||           (asymmetry/nonlinearity)
#
# Run:  julia --project=julia julia/analysis/log_oat_sensitivity.jl [tmature] [tpost] [ncells]
# Writes results/log_oat/{measures.csv, curves_<param>.csv, L0.csv} and a ranked table.

using MPCSSF
using Printf, DelimitedFiles, Statistics

include(joinpath(@__DIR__, "proxy_model.jl"))       # run_proxy (implicit osmosis)
include(joinpath(@__DIR__, "pathogen_repro.jl"))    # pathogen_model()

# ---- forcing / disturbance config ------------------------------------------
const INFLUENT_BASE = Float64[2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]
const PAT_IN   = INFLUENT_BASE[4]
# Challenge = 10× ambient PAT. NB: a 100× step pulse creates a front too sharp for
# the capped proxy dt to resolve, spiking φ_b in the surface cell to a spurious
# clog (verified: ×100 CLOGGED at t=0.34, ×10 runs clean to T). 10× is a
# legitimate contamination challenge and keeps every run physical.
const PULSE_FACTOR = 10.0
const PULSE_T0 = 0.1              # pulse start (day) = disturbance time t_d
const PULSE_T1 = 0.3              # pulse end (day)
const TD = PULSE_T0               # analysis window start
const LN2 = log(2)

# influent as a function of time: PAT pulse over [T0,T1], ×pat_mult ambient else.
function challenge_inflow(pat_mult)
    Cref = PULSE_FACTOR * PAT_IN * pat_mult
    function infl(t)
        v = copy(INFLUENT_BASE)
        v[4] = (PULSE_T0 <= t < PULSE_T1) ? Cref : PAT_IN * pat_mult
        return v
    end
    return infl, Cref
end

# ---- immutable-struct remake helpers ---------------------------------------
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

# ---- parameter definition ---------------------------------------------------
# `apply(f, m, v)` sets the parameter to absolute value `v` and returns (f, m).
# `nominal` is the value perturbed by ×2 and ×1/2. `pat_mult` params are special-
# cased by name in `evaluate` (they scale the influent, not (f,m)).
struct P
    name::String; block::String; nominal::Float64; source::String; apply::Function
end

set_rate(names) = (f, m, v) -> (f, _mdl(m; reactions=Reaction[
    (r.name in names) ? _rx(r; nominal_rate=v) : r for r in m.reactions]))
set_theta(names) = (f, m, dev) -> (f, _mdl(m; reactions=Reaction[    # dev = θ−1
    (r.name in names) ? _rx(r; temperature_correction_factor=1 + dev) : r for r in m.reactions]))
set_halfsat(rxname, key) = (f, m, v) -> (f, _mdl(m; reactions=Reaction[
    (r.name == rxname) ? _rx(r; half_saturation_constants=merge(r.half_saturation_constants, Dict(key=>v))) : r
    for r in m.reactions]))
set_attach_sand() = (f, m, v) -> (f, _mdl(m; components=Component[
    c isa Particle ? _pt(c; attachment_sand=v) : c for c in m.components]))
set_sand_pat() = (f, m, v) -> (f, _mdl(m; components=Component[
    (c isa Particle && c.name=="PAT") ? _pt(c; sand_attachment_factor=v) : c for c in m.components]))
set_disp() = (f, m, v) -> (f, _mdl(m; components=Component[
    c isa Particle ? _pt(c; dispersivity=v) : c for c in m.components]))
set_transport_P() = (f, m, v) -> (f, _mdl(m; components=Component[
    c isa Particle ? _pt(c; transport_rate=v) : c for c in m.components]))
set_beta_gap() = (f, m, gap) -> (f, _mdl(m; biofilm_porosity=1 - gap))    # gap = 1−β
set_zeta0() = (f, m, v) -> (f, _mdl(m; cohesion_submodel=CahnHilliardModel(
    kappa=m.cohesion_submodel.kappa, zeta_0=v, zeta_1=m.cohesion_submodel.zeta_1)))
set_detach(nom) = (f, m, v) -> (f, _mdl(m; detachment = w -> (v/nom) .* m.detachment(w)))
set_velocity() = (f, m, v) -> (f.inflow_velocity = v; (f, m))
set_temp() = (f, m, v) -> (f.temperature = v; (f, m))
noop() = (f, m, v) -> (f, m)   # for pat-mult params (handled in evaluate)

# The parameter set. Ranges/provenance cited from Campos2006 & Schijven2013 (the
# OAT only needs the nominal; ranges live in morris_screening.jl / the papers).
const PARAMS = P[
    P("velocity",      "forcing",   7.2,     "Schijven2013 T1 filtration rate",           set_velocity()),
    P("temperature",   "forcing",   15.0,    "Campos2006/Manriquez seasonal",             set_temp()),
    P("influent_PAT",  "forcing",   1.0,     "Manriquez B.1 PAT (challenge multiplier)",  noop()),
    P("dispersivity",  "transport", 0.012,   "Schijven2013 T1 dispersivity",              set_disp()),
    P("transport_P",   "transport", 5.47,    "Lund phase transfer",                       set_transport_P()),
    P("attach_sand",   "transport", 547.0,   "Diehl2025/Lund b_sand",                     set_attach_sand()),
    P("sand_pathogen", "transport", 0.06,    "Schijven2013 T4 sticking α (geo-mean)",     set_sand_pat()),
    P("mu_HET",        "kinetics",  0.0181,  "Campos2006 kgmaxa",                          set_rate(["Heterotroph growth"])),
    P("mu_PHO",        "kinetics",  5.5,     "Campos2006 kgmaxb",                          set_rate(["Phototroph growth"])),
    P("d_HET",         "kinetics",  2.0,     "Campos2006 T3 kdb",                          set_rate(["Heterotroph death"])),
    P("hydrolysis",    "kinetics",  0.09,    "Campos2006 T3 kh",                           set_rate(["Hydrolysis"])),
    P("theta_growth",  "kinetics",  0.047,   "Campos2006 T2 θ_growth=1.047 (dev)",         set_theta(["Heterotroph growth","Phototroph growth"])),
    P("theta_death",   "kinetics",  0.066,   "Campos2006 T2 θ_death=1.066 (dev)",          set_theta(["Heterotroph death","Phototroph death"])),
    P("K_O2_HET",      "kinetics",  3.0e-3,  "O2 half-sat (Reichert)",                     set_halfsat("Heterotroph growth","O2")),
    P("K_DOM_HET",     "kinetics",  2.0e-4,  "Campos2006 T3 ksCd",                         set_halfsat("Heterotroph growth","DOM")),
    P("K_HPO4_HET",    "kinetics",  1.4e-8,  "Campos2006 T3 ksp (P-limiting)",             set_halfsat("Heterotroph growth","HPO4")),
    P("marker_growth", "pathogen",  0.2,     "Manriquez B.4 μ̂_PAT",                        set_rate(["MarkerGrowth"])),
    P("inactivation",  "pathogen",  0.02,    "Schijven2013 µl,µs; Manriquez B.4",          set_rate(["Inactivation"])),
    P("bacterivory",   "pathogen",  8.0,     "Manriquez B.4 p̂_PAT",                        set_rate(["Bacterivory"])),
    P("beta_porosity", "biofilm",   0.01,    "Lund β=0.99 (gap 1−β)",                      set_beta_gap()),
    P("zeta_0",        "biofilm",   1.0e2,   "Cahn-Hilliard cohesion",                     set_zeta0()),
    P("detach_scale",  "biofilm",   1.0,     "detachment-law scale",                       set_detach(1.0)),
]

# ---- build the shared mature state (nominal model) --------------------------
function build_mature(m0, ncells, tmature)
    f = addgridpoints(SandFilter(temperature=15.0), ncells)
    r = run_proxy(State(f, m0); simulation_time=tmature, inflow_concentrations=INFLUENT_BASE,
                  n_frames=5, quiet=true)
    @info "mature state" flag=string(r.flag) t=r.time_final phib=maximum(get_volume_fractions(r).biofilm)
    return final_state(r)
end

# ---- one challenge run: returns (times, L(t), flag) -------------------------
function run_challenge(ms, m0, p::Union{P,Nothing}, value, ncells, tpost, nframes)
    f2 = addgridpoints(SandFilter(temperature=15.0), ncells)
    m2 = m0
    pat_mult = 1.0
    if p !== nothing
        if p.name == "influent_PAT"
            pat_mult = value
        else
            f2, m2 = p.apply(f2, m2, value)
        end
    end
    # re-home the mature concentrations onto the (possibly perturbed) f2/m2 at t=0
    s = State(f2, m2, 0.0, deepcopy(ms.global_concentration),
              copy(ms.enclosed_water_volume), deepcopy(ms.velocity))
    infl, Cref = challenge_inflow(pat_mult)
    r = run_proxy(s; simulation_time=tpost, inflow_concentrations=infl,
                  n_frames=nframes, quiet=true)
    ts = times(r)
    cout = max.(concentration(r, "PAT", :flowing)[end, :], 1e-30)
    return ts, log10.(Cref ./ cout), string(r.flag)
end

# ---- driver -----------------------------------------------------------------
function main(; tmature=3.0, tpost=1.5, ncells=30, nframes=60)
    m0 = pathogen_model()
    outdir = joinpath(@__DIR__, "results", "log_oat"); isdir(outdir) || mkpath(outdir)
    @info "Log-OAT sensitivity" nparams=length(PARAMS) tmature tpost ncells runs=2*length(PARAMS)+1

    ms = build_mature(m0, ncells, tmature)

    # baseline L0(t)
    t0, L0, flag0 = run_challenge(ms, m0, nothing, 0.0, ncells, tpost, nframes)
    win = t0 .>= TD                       # post-disturbance window mask
    Lmin0 = minimum(L0[win])
    writedlm(joinpath(outdir, "L0.csv"), hcat(t0, L0), ',')
    @printf("baseline: flag=%s  window pts=%d  Lmean=%.3f  Lmin=%.3f\n\n",
            flag0, count(win), mean(L0[win]), Lmin0)

    rows = Vector{Any}[["param","block","I_rms","I_max","D_min","I_min","asymmetry","flag+","flag-","source"]]
    ranking = Tuple{String,String,Float64,Float64,Float64}[]
    for p in PARAMS
        tp, Lp, fp = run_challenge(ms, m0, p, 2 * p.nominal,   ncells, tpost, nframes)  # ×2
        tm, Lm, fm = run_challenge(ms, m0, p, 0.5 * p.nominal, ncells, tpost, nframes)  # ×1/2
        # all runs share the frame grid (same t0, tpost, nframes) → align on window
        s  = (Lp .- Lm) ./ (2LN2)
        sw = s[win]
        Irms = sqrt(mean(sw .^ 2))
        Imax = maximum(abs.(sw))
        Lminp = minimum(Lp[win]); Lminm = minimum(Lm[win])
        Dmin = max(Lminp - Lmin0, Lminm - Lmin0)
        Imin = (Lminp - Lminm) / (2LN2)
        dp = Lp[win] .- L0[win]; dm = Lm[win] .- L0[win]
        asym = sqrt(sum((dp .+ dm) .^ 2)) / (sqrt(sum((dp .- dm) .^ 2)) + 1e-30)
        push!(rows, [p.name, p.block, Irms, Imax, Dmin, Imin, asym, fp, fm, p.source])
        writedlm(joinpath(outdir, "curves_$(p.name).csv"),
                 vcat(["t" "dL_plus" "dL_minus" "s"], hcat(tp, Lp .- L0, Lm .- L0, s)), ',')
        push!(ranking, (p.name, p.block, Irms, Imax, Dmin))
        @printf("  %-14s [%-9s] I_rms=%.3e  I_max=%.3e  D_min=%+.3e  asym=%.2f  (%s/%s)\n",
                p.name, p.block, Irms, Imax, Dmin, asym, fp, fm)
    end
    writedlm(joinpath(outdir, "measures.csv"), rows, ',')

    println("\n── Ranking by RMS log-sensitivity I_i (primary) ──")
    for (nm, bl, ir, im, dm) in sort(ranking; by=x->x[3], rev=true)
        @printf("  I_rms=%.3e  %-14s (%s)\n", ir, nm, bl)
    end
    println("\n── Ranking by D_min (effect on worst removal) ──")
    for (nm, bl, ir, im, dm) in sort(ranking; by=x->abs(x[5]), rev=true)[1:min(8,length(ranking))]
        @printf("  |D_min|=%.3e  %-14s (%s)\n", abs(dm), nm, bl)
    end
    println("\nwrote ", outdir)
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    tmat = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 3.0
    tpost = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.5
    ncells = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 30
    main(; tmature=tmat, tpost=tpost, ncells=ncells)
end
