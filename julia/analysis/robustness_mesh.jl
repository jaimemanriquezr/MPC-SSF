# Mesh-robustness of the OAT SENSITIVITY (not just the baseline). The mesh-
# convergence study showed the integral QoIs (Lmean) converge; this checks that
# their log-SENSITIVITY I_rms for the top parameters is also stable across meshes —
# the condition under which the coarse-mesh ranking is trustworthy despite the
# unresolved surface (δ, κ scales are sub-grid on any feasible mesh).
#
# Run: julia --project=julia julia/analysis/robustness_mesh.jl
using MPCSSF, Printf, Statistics
include(joinpath(@__DIR__, "log_oat_sensitivity.jl"))   # PARAMS, run_challenge, build_mature, TD, LN2

const TOP = ["attach_sand","sand_pathogen","influent_PAT","dispersivity","temperature","bacterivory"]
sel = [p for p in PARAMS if p.name in TOP]

function irms(ms, m0, p, nc, tpost, nframes)
    t0, L0, f0, tf0 = run_challenge(ms, m0, nothing, 0.0, nc, tpost, nframes; disturbance=:pulse)
    tp, Lp, fp, tfp = run_challenge(ms, m0, p, 2 * p.nominal,   nc, tpost, nframes; disturbance=:pulse)
    tm, Lm, fm, tfm = run_challenge(ms, m0, p, 0.5 * p.nominal, nc, tpost, nframes; disturbance=:pulse)
    tv = min(tf0, tfp, tfm); win = (t0 .>= TD) .& (t0 .<= tv + 1e-9)
    s = (Lp .- Lm) ./ (2LN2)
    return sqrt(mean(s[win] .^ 2)), "$fp/$fm"
end

function main(; meshes=(30, 60), tmature=2.0, tpost=1.0, nframes=40)
    R = Dict{Int,Dict{String,Float64}}()
    for nc in meshes
        m0 = pathogen_model(); ms = build_mature(m0, nc, tmature)
        R[nc] = Dict{String,Float64}()
        @printf("── ncells=%d ──\n", nc)
        for p in sel
            v, fl = irms(ms, m0, p, nc, tpost, nframes)
            R[nc][p.name] = v
            @printf("  %-14s I_rms=%.3e  (%s)\n", p.name, v, fl)
        end
    end
    println("\n── I_rms across meshes (ratio finest/coarsest) ──")
    lo, hi = minimum(meshes), maximum(meshes)
    for p in sel
        a = get(R[lo], p.name, NaN); b = get(R[hi], p.name, NaN)
        @printf("  %-14s  %d:%.3e  %d:%.3e   ratio=%.2f\n", p.name, lo, a, hi, b, b / a)
    end
end

main()
