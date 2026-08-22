# Golden-master comparison for the adaptive-CFL path: run the shared adaptive
# modelLund simulation in Julia and diff it against the MATLAB reference exported
# by export_adaptive_reference.m.
#
# Standalone (detailed diff table):
#   julia --project=julia julia/test/golden/compare_adaptive.jl <reference_dir>
# From the test suite: runtests.jl includes this and calls
#   golden_compare_adaptive(refdir; verbose=true).
#
# The run spec MUST match export_adaptive_reference.m exactly.

using SSF
using DelimitedFiles
using Printf

# --- shared run spec (keep in sync with export_adaptive_reference.m) ---------
const AD_N       = 20
const AD_NFRAMES = 5
#                          HET   PHO   POM   PAT    O2    IC    NH4   HPO4  DOM
const AD_INFLOW  = Float64[1e-2, 1e-3, 0.0,  1e-4, 1e-2, 5e-3, 4e-3, 1e-4, 1.0]
const AD_SIMTIME = 1e-5
const AD_INIT_DT = 1e-9
const AD_MAX_DT  = 1e-6
const AD_TOL     = 0.05
const AD_CFL     = 0.99

const GOLDEN_ADAPTIVE_REF = joinpath(@__DIR__, "reference_adaptive")

_admat(refdir, name) = readdlm(joinpath(refdir, name), ',', Float64)

function _adparse_meta(refdir)
    meta = Dict{String,String}()
    for line in eachline(joinpath(refdir, "meta.txt"))
        isempty(strip(line)) && continue
        k, v = split(line, '=', limit=2)
        meta[k] = v
    end
    return meta
end

"Run the shared adaptive golden-master simulation in Julia."
function golden_run_adaptive()
    filter = addgridpoints(SandFilter(), AD_N)
    model  = modelLund()
    return simulate(State(filter, model); inflow_concentrations=AD_INFLOW,
                    simulation_time=AD_SIMTIME, time_step=:adaptive,
                    adaptive_initial_dt=AD_INIT_DT, adaptive_max_dt=AD_MAX_DT,
                    adaptive_time_tolerance=AD_TOL, cfl_factor=AD_CFL,
                    n_frames=AD_NFRAMES, quiet=true)
end

function _addiffs(a, b; floor=1e-12)
    n = min(length(a), length(b))
    av, bv = a[1:n], b[1:n]
    absd = abs.(av .- bv)
    maxabs = maximum(absd; init=0.0)
    sig = abs.(bv) .> floor
    maxrel = any(sig) ? maximum(absd[sig] ./ abs.(bv[sig]); init=0.0) : 0.0
    return maxabs, maxrel
end

"""
    golden_compare_adaptive(refdir=GOLDEN_ADAPTIVE_REF; rtol=1e-4, atol=1e-7, verbose=false)

Compare the Julia adaptive run against the MATLAB reference: dt trajectory
(step_times), frame times/depths, per-region frame concentrations. A field
passes on abs OR rel within tolerance. The dt trajectory is compared over the
common-length prefix (step counts can differ by a few once dt reaches the CFL
plateau, where the bound depends on the sparse solve).
"""
function golden_compare_adaptive(refdir=GOLDEN_ADAPTIVE_REF; rtol=1e-4, atol=1e-7, verbose=false)
    meta = _adparse_meta(refdir)
    names = split(meta["names"], ',')
    res = golden_run_adaptive()
    flags_agree = String(strip(get(meta, "flag", "?"))) == res.flag

    ref_steps = Int(round(parse(Float64, meta["nsteps"])))
    jl_steps = length(res.simulation_data[:step_times]) - 1
    verbose && @printf("MATLAB flag=%s (%d steps) | Julia flag=%s (%d steps)\n\n",
                       get(meta, "flag", "?"), ref_steps, res.flag, jl_steps)

    worst_abs = 0.0; worst_rel = 0.0; fail = false
    function report(label, a, b)
        ma, mr = _addiffs(a, b)
        worst_abs = max(worst_abs, ma); worst_rel = max(worst_rel, mr)
        ok = ma <= atol || mr <= rtol
        fail |= !ok
        verbose && @printf("  %-26s maxabs=%.3e  maxrel=%.3e  %s\n", label, ma, mr, ok ? "ok" : "FAIL")
    end

    verbose && println("field                       differences (Julia vs MATLAB)")
    report("depths", depths(res), vec(_admat(refdir, "depths.csv")))
    report("frame times", times(res), vec(_admat(refdir, "times.csv")))
    # dt trajectory (compare the common prefix of the per-step time vector)
    report("step_times", res.simulation_data[:step_times], vec(_admat(refdir, "step_times.csv")))

    regions = [(:matrix, "matrix"), (:enclosed, "enclosed"), (:flowing, "flowing")]
    for nm in names, (rsym, rstr) in regions
        report("$(nm)/$(rstr)", concentration(res, String(nm), rsym),
               _admat(refdir, "conc_$(nm)_$(rstr).csv"))
    end
    report("Water/enclosed", concentration(res, "Water", :enclosed),
           _admat(refdir, "conc_Water_enclosed.csv"))

    if verbose
        @printf("\nworst: maxabs=%.3e maxrel=%.3e | step counts: MATLAB=%d Julia=%d\n",
                worst_abs, worst_rel, ref_steps, jl_steps)
        println(fail ? "RESULT: MISMATCH beyond tolerance." :
                       "RESULT: MATCH — Julia reproduces the MATLAB adaptive reference.")
    end
    return (match=!fail, worst_abs=worst_abs, worst_rel=worst_rel, flags_agree=flags_agree,
            ref_steps=ref_steps, jl_steps=jl_steps)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("usage: compare_adaptive.jl <reference_dir>")
    r = golden_compare_adaptive(ARGS[1]; verbose=true)
    exit((r.match && r.flags_agree) ? 0 : 1)
end
