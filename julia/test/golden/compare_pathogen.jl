# Golden-master comparison for the PATHOGEN model: run the pathogen simulation
# in Julia and diff it against the MATLAB reference exported by
# export_pathogen_reference.m (which drives the authoritative
# slow-sand-filtration @SDfilter/run_pathogen.m + thesis_model.mat).
#
# Standalone (detailed diff table):
#     julia --project=julia julia/test/golden/compare_pathogen.jl <reference_dir>
# From the test suite: runtests.jl includes this file and calls
#     golden_compare_pathogen(refdir; verbose=true)
#
# The run spec MUST match export_pathogen_reference.m exactly. Confounders are
# removed so the two codebases are numerically comparable: 20 C (Julia) / 293 K
# (MATLAB) => mu = nominal; dark_respiration = 0 => identical light factor;
# detachment sqrt(|v|/7.2) on both sides; grid addgridpoints(20) == add_cells(20).

using SSF
using DelimitedFiles
using Printf

# --- shared run spec (keep in sync with export_pathogen_reference.m) --------
const PGOLD_N       = 20
const PGOLD_TSIM    = 1e-5
const PGOLD_DT      = 1e-6
const PGOLD_NFRAMES = 5
#                            HET   PHO   POM  PAT   O2    IC    NH4   HPO4 DOM
const PGOLD_INFLOW = Float64[1e-3, 1e-3, 0.0, 1e-5, 1e-2, 1e-2, 1e-5, 0.0, 1e-4] ./ 10
# uniform mature-biofilm seed (matches ic.biofilm in export_pathogen_reference.m)
const PGOLD_SEED_MATRIX = Float64[0.05, 0.05, 0.02, 0.01]   # HET PHO POM PAT
const PGOLD_SEED_ENCL_L = 1e-3                              # enclosed liquids
const PGOLD_SEED_WATER  = 0.1                               # enclosed water volume

_patreadmat(refdir, name) = readdlm(joinpath(refdir, name), ',', Float64)

function _patparse_meta(refdir)
    meta = Dict{String,String}()
    for line in eachline(joinpath(refdir, "meta.txt"))
        isempty(strip(line)) && continue
        k, v = split(line, '=', limit=2)
        meta[k] = v
    end
    return meta
end

"""
    golden_run_pathogen_julia(; dark_respiration=0.0, light_const=nothing,
                                seed_matrix=PGOLD_SEED_MATRIX,
                                seed_enclosed_liquids=PGOLD_SEED_ENCL_L) -> Results

Run the pathogen golden-master simulation in Julia. `light_const` (a number) sets
a constant light forcing; `nothing` keeps the default diel forcing. Keep the
keyword defaults in sync with `export_pathogen_reference.m`'s option defaults.
"""
function golden_run_pathogen_julia(; dark_respiration=0.0, light_const=nothing,
                                   seed_matrix=PGOLD_SEED_MATRIX,
                                   seed_enclosed_liquids=PGOLD_SEED_ENCL_L)
    filter = light_const === nothing ? SandFilter(temperature=20.0) :
             SandFilter(temperature=20.0, light_irradiation=(_t -> float(light_const)))
    filter = addgridpoints(filter, PGOLD_N)      # 20 C -> theta^0 = 1 (mu = nominal)
    model  = modelPathogen(water_factor=1e-3, sand_pathogen=0.1,
                           dark_respiration=dark_respiration)
    state  = State(filter, model)
    # Seed the same uniform mature biofilm as the MATLAB reference: matrix
    # particles + enclosed liquids + enclosed water (enclosed particles stay 0).
    state.global_concentration.matrix .= reshape(seed_matrix, 1, :)
    state.global_concentration.enclosed_liquids .= seed_enclosed_liquids
    state.enclosed_water_volume .= PGOLD_SEED_WATER
    return simulate(state; inflow_concentrations=PGOLD_INFLOW,
                    simulation_time=PGOLD_TSIM, time_step=PGOLD_DT,
                    n_frames=PGOLD_NFRAMES, clogging_fraction=0.99, quiet=true)
end

function _patdiffs(a, b; floor=1e-12)
    size(a) == size(b) || error("shape mismatch $(size(a)) vs $(size(b))")
    absdiff = abs.(a .- b)
    maxabs = maximum(absdiff; init=0.0)
    sig = abs.(b) .> floor
    maxrel = any(sig) ? maximum(absdiff[sig] ./ abs.(b[sig]); init=0.0) : 0.0
    return maxabs, maxrel
end

const PGOLD_DEFAULT_REF = joinpath(@__DIR__, "reference_pathogen")
const PGOLD_LIGHT_REF   = joinpath(@__DIR__, "reference_pathogen_light")

# Run kwargs for the light+dark variant (see run_export_pathogen.m). Constant
# light + dark_respiration>0 with a phototroph-heavy seed exercises the
# dark-respiration light floor max(fdark, I·e^{1-I}).
const PGOLD_LIGHT_RUN = (dark_respiration=0.1, light_const=0.8,
                         seed_matrix=Float64[0.05, 0.1, 0.02, 0.01],
                         seed_enclosed_liquids=0.05)

"""
    golden_compare_pathogen(refdir=PGOLD_DEFAULT_REF; run_kwargs=(;),
                            rtol=1e-6, atol=1e-7, verbose=false)
        -> (match, worst_abs, worst_rel, flags_agree)

Run the Julia pathogen simulation (with `run_kwargs` forwarded to
[`golden_run_pathogen_julia`](@ref)) and compare every exported field against the
MATLAB reference in `refdir`. A field passes if its max abs OR max rel difference
is within tolerance.
"""
function golden_compare_pathogen(refdir=PGOLD_DEFAULT_REF; run_kwargs=(;),
                                 rtol=1e-6, atol=1e-7, verbose=false)
    meta = _patparse_meta(refdir)
    names = split(meta["names"], ',')
    res = golden_run_pathogen_julia(; run_kwargs...)
    flags_agree = String(strip(get(meta, "flag", "?"))) == res.flag

    verbose && @printf("MATLAB flag: %s | Julia flag: %s\n\n",
                       get(meta, "flag", "?"), res.flag)

    worst_abs = 0.0
    worst_rel = 0.0
    fail = false
    function report(label, a, b)
        ma, mr = _patdiffs(a, b)
        worst_abs = max(worst_abs, ma)
        worst_rel = max(worst_rel, mr)
        ok = ma <= atol || mr <= rtol
        fail |= !ok
        verbose && @printf("  %-28s maxabs=%.3e  maxrel=%.3e  %s\n",
                           label, ma, mr, ok ? "ok" : "FAIL")
    end

    verbose && println("field                         differences (Julia vs MATLAB)")
    report("depths", depths(res), vec(_patreadmat(refdir, "depths.csv")))
    report("times",  times(res),  vec(_patreadmat(refdir, "times.csv")))

    regions = [(:matrix, "matrix"), (:enclosed, "enclosed"), (:flowing, "flowing")]
    for nm in names
        for (rsym, rstr) in regions
            ref = _patreadmat(refdir, "conc_$(nm)_$(rstr).csv")
            jl  = concentration(res, String(nm), rsym)
            report("$(nm)/$(rstr)", jl, ref)
        end
    end
    report("Water/enclosed", concentration(res, "Water", :enclosed),
           _patreadmat(refdir, "conc_Water_enclosed.csv"))

    if verbose
        @printf("\nworst over all fields: maxabs=%.3e  maxrel=%.3e\n", worst_abs, worst_rel)
        println(fail ? "RESULT: MISMATCH beyond tolerance." :
                       "RESULT: MATCH — Julia reproduces the MATLAB pathogen reference.")
    end
    return (match=!fail, worst_abs=worst_abs, worst_rel=worst_rel, flags_agree=flags_agree)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: compare_pathogen.jl <reference_dir>")
    r = golden_compare_pathogen(ARGS[1]; verbose=true)
    exit((r.match && r.flags_agree) ? 0 : 1)
end
