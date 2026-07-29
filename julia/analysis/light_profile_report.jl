# Step 5 of the light–biofilm-profile study: collect the per-amplitude QoIs the
# sweep tasks wrote and report them against irradiance amplitude.
#
# The hypothesis under test (see .claude/plans/2026-07-29-light-biofilm-profile.md)
# is that raising irradiance puts MORE biofilm above z = 0. The photoinhibition
# reading predicts the opposite: the near-optimal light band sits where attenuation
# has brought I_eff down to ~1, so raising amplitude pushes that band DEEPER and
# leaves the supernatant pinned at minimum_light_factor. This script reports the
# direction rather than asserting it.
#
# Run: julia --project=. analysis/light_profile_report.jl [results_dir]

using Printf, DelimitedFiles, Statistics

function collect_qois(dir)
    files = filter(f -> startswith(f, "qoi_A") && endswith(f, ".csv"), readdir(dir))
    isempty(files) && error("no qoi_A*.csv in $dir — the sweep has not written yet")
    rows = NamedTuple[]
    for f in sort(files)
        d = readdlm(joinpath(dir, f), ',')
        hdr = string.(vec(d[1, :])); val = vec(d[2, :])
        g(name) = val[findfirst(==(name), hdr)]
        push!(rows, (amplitude = Float64(g("amplitude")),
                     flag       = string(g("flag")),
                     t_final    = Float64(g("t_final")),
                     z_peak     = Float64(g("z_peak_m")),
                     phib_peak  = Float64(g("phib_peak")),
                     mass_above = Float64(g("mass_above")),
                     mass_total = Float64(g("mass_total")),
                     wall       = Float64(g("wall_s"))))
    end
    return sort(rows; by = r -> r.amplitude)
end

function main(dir)
    rows = collect_qois(dir)
    incomplete = filter(r -> r.flag != "OK", rows)

    @printf("%-9s %-8s %-8s %11s %11s %12s %12s %8s\n",
            "amplitude", "flag", "t_final", "z_peak(mm)", "phib_peak",
            "M(z<0)", "M(total)", "frac>0")
    for r in rows
        frac = r.mass_total > 0 ? r.mass_above / r.mass_total : NaN
        @printf("%-9.4g %-8s %-8.3f %+11.2f %11.4g %12.4g %12.4g %8.3f\n",
                r.amplitude, r.flag, r.t_final, r.z_peak * 1e3, r.phib_peak,
                r.mass_above, r.mass_total, frac)
    end

    body = reduce(vcat, [permutedims(Any[
                r.amplitude, r.flag, r.t_final, r.z_peak, r.phib_peak,
                r.mass_above, r.mass_total,
                (r.mass_total > 0 ? r.mass_above / r.mass_total : NaN), r.wall])
            for r in rows])
    writedlm(joinpath(dir, "summary.csv"),
             vcat(["amplitude" "flag" "t_final" "z_peak_m" "phib_peak" "mass_above" "mass_total" "frac_above" "wall_s"],
                  body), ',')

    # Direction of the effect, stated only over runs that actually completed.
    ok = filter(r -> r.flag == "OK", rows)
    if length(ok) >= 2
        lo, hi = ok[1], ok[end]
        dz = (hi.z_peak - lo.z_peak) * 1e3
        dm = hi.mass_above - lo.mass_above
        println()
        @printf("amplitude %.4g -> %.4g:  peak moves %+.2f mm (%s),  M(z<0) %+.4g (%s)\n",
                lo.amplitude, hi.amplitude, dz, dz < 0 ? "deeper" : "shallower",
                dm, dm < 0 ? "falls" : "rises")
        println(dz < 0 && dm < 0 ?
                "Consistent with the photoinhibition reading: more light pushes biofilm deeper and empties the supernatant." :
                "NOT the photoinhibition prediction — re-examine the reasoning in the plan before reporting.")
    end

    if !isempty(incomplete)
        println()
        @printf("WARNING: %d of %d runs did not reach flag=OK; rows above are not comparable.\n",
                length(incomplete), length(rows))
        for r in incomplete
            @printf("  amplitude %-9.4g flag=%-10s stopped at t=%.4f of 7.0 d\n",
                    r.amplitude, r.flag, r.t_final)
        end
    end
    println("\nwrote ", joinpath(dir, "summary.csv"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results", "light_profile")
    main(dir)
end
