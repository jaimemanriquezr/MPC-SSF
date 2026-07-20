# Golden-master comparison: run the shared SimpleModel simulation in Julia and
# diff it against the MATLAB reference exported by export_reference.m.
#
# Two ways to use it:
#   * Standalone (prints a detailed diff table):
#       julia --project=julia julia/test/golden/compare.jl <reference_dir>
#   * From the test suite: runtests.jl includes this file and calls
#     `golden_compare(refdir; verbose=true)` inside a testset when the env var
#     MPCSSF_GOLDEN_REF points at a reference directory.
#
# The run spec MUST match export_reference.m exactly.

using MPCSSF
using DelimitedFiles
using Printf

# --- shared run spec (keep in sync with export_reference.m) -----------------
const GOLDEN_N       = 20
const GOLDEN_TSIM    = 1e-3
const GOLDEN_DT      = 1e-5
const GOLDEN_NFRAMES = 5
#                              Microorganism  POM   Nutrient
const GOLDEN_INFLOW  = Float64[1e-2,          0.0,  1.0]

_readmat(refdir, name) = readdlm(joinpath(refdir, name), ',', Float64)

function _parse_meta(refdir)
    meta = Dict{String,String}()
    for line in eachline(joinpath(refdir, "meta.txt"))
        isempty(strip(line)) && continue
        k, v = split(line, '=', limit=2)
        meta[k] = v
    end
    return meta
end

"Run the shared golden-master simulation in Julia and return its Results."
function golden_run_julia()
    filter = SandFilter()
    filter = addgridpoints(filter, GOLDEN_N)
    model  = simpleModel()
    state  = State(filter, model)
    return simulate(state; inflow_concentrations=GOLDEN_INFLOW,
                    simulation_time=GOLDEN_TSIM, time_step=GOLDEN_DT,
                    n_frames=GOLDEN_NFRAMES, clogging_fraction=0.99, quiet=true)
end

# max absolute and max relative difference between two arrays. The relative diff
# only counts entries whose reference magnitude exceeds `floor`; entries at/below
# it are judged by the absolute diff alone (a near-zero reference makes relative
# error meaningless — e.g. 1e-9 / 1e-49).
function _diffs(a, b; floor=1e-12)
    size(a) == size(b) || error("shape mismatch $(size(a)) vs $(size(b))")
    absdiff = abs.(a .- b)
    maxabs = maximum(absdiff; init=0.0)
    sig = abs.(b) .> floor
    maxrel = any(sig) ? maximum(absdiff[sig] ./ abs.(b[sig]); init=0.0) : 0.0
    return maxabs, maxrel
end

# Committed reference (produced by export_reference.m); used when no dir is given.
const GOLDEN_DEFAULT_REF = joinpath(@__DIR__, "reference")

"""
    golden_compare(refdir=GOLDEN_DEFAULT_REF; rtol=1e-6, atol=1e-7, verbose=false)
        -> (match::Bool, worst_abs, worst_rel, flags_agree::Bool)

Run the Julia simulation and compare every exported field against the MATLAB
reference in `refdir`. A field passes if its max abs OR max rel difference is
within tolerance.

Tolerances are `atol=1e-7`, `rtol=1e-6`. Every field except the flowing phase
matches to ~1e-15 or better; the flowing-phase concentrations agree to ~9e-9
absolute (~6.5e-8 relative), which is accumulated floating-point difference over
100 nonlinear steps — MATLAB and Julia use different BLAS/sparse-solver
libraries, so bit-identical agreement below ~1e-8 is not achievable. The
micro/nutrient flowing errors are equal-and-opposite (micro+nutrient conserved),
confirming the residual is the Growth reaction's last-bit split, not a port bug.
See README.md.
"""
function golden_compare(refdir=GOLDEN_DEFAULT_REF; rtol=1e-6, atol=1e-7, verbose=false)
    meta = _parse_meta(refdir)
    names = split(meta["names"], ',')
    res = golden_run_julia()
    flags_agree = String(strip(get(meta, "flag", "?"))) == res.flag

    verbose && @printf("MATLAB flag: %s | Julia flag: %s\n\n",
                       get(meta, "flag", "?"), res.flag)

    worst_abs = 0.0
    worst_rel = 0.0
    fail = false
    function report(label, a, b)
        ma, mr = _diffs(a, b)
        worst_abs = max(worst_abs, ma)
        worst_rel = max(worst_rel, mr)
        ok = ma <= atol || mr <= rtol
        fail |= !ok
        verbose && @printf("  %-28s maxabs=%.3e  maxrel=%.3e  %s\n",
                           label, ma, mr, ok ? "ok" : "FAIL")
    end

    verbose && println("field                         differences (Julia vs MATLAB)")
    report("depths", depths(res), vec(_readmat(refdir, "depths.csv")))
    report("times",  times(res),  vec(_readmat(refdir, "times.csv")))

    regions = [(:matrix, "matrix"), (:enclosed, "enclosed"), (:flowing, "flowing")]
    for nm in names
        for (rsym, rstr) in regions
            ref = _readmat(refdir, "conc_$(nm)_$(rstr).csv")
            jl  = concentration(res, String(nm), rsym)
            report("$(nm)/$(rstr)", jl, ref)
        end
    end
    report("Water/enclosed", concentration(res, "Water", :enclosed),
           _readmat(refdir, "conc_Water_enclosed.csv"))

    if verbose
        @printf("\nworst over all fields: maxabs=%.3e  maxrel=%.3e\n", worst_abs, worst_rel)
        println(fail ? "RESULT: MISMATCH beyond tolerance." :
                       "RESULT: MATCH — Julia reproduces the MATLAB reference.")
    end
    return (match=!fail, worst_abs=worst_abs, worst_rel=worst_rel, flags_agree=flags_agree)
end

# Standalone entry point.
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: compare.jl <reference_dir>")
    r = golden_compare(ARGS[1]; verbose=true)
    exit((r.match && r.flags_agree) ? 0 : 1)
end
