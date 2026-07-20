# Run chaining: resume a simulation from a previous run's final state, and
# stitch consecutive Results together.
#
# Port of the slow-sand-filtration @SDresults `get_filter_state` (resume) and
# `concatenate`/`plus` (append). `final_state(results)` rebuilds a [`State`](@ref)
# from a captured frame so `simulate(final_state(r); ...)` continues where `r`
# stopped; `concatenate(r1, r2)` (also `r1 + r2`) joins them into one Results.

"Mean density of a component group (mirrors simulate's assumption-1 averaging)."
_group_density(comps) = sum(c.density for c in comps) / length(comps)

"""
    final_state(results::Results; t=results.time_final) -> State

Rebuild a resumable [`State`](@ref) from the frame of `results` at time `t`
(default: the final frame). Concentrations, enclosed water, biofilm velocity and
time are taken from that frame, so `simulate(final_state(r); ...)` continues the
run from where `r` ended. Port of `@SDresults/get_filter_state.m`.
"""
function final_state(results::Results; t::Real=results.time_final)
    f = results.filter
    m = results.model
    P = particles(m); L = liquids(m)
    kP = length(P); kL = length(L)

    ft = times(results)
    idx = findfirst(>=(t), ft)
    idx === nothing && (idx = length(ft))

    cb = results.frames[:concentration_biofilm]   # N × nf × (2kP+kL)
    cf = results.frames[:concentration_flowing]   # N × nf × (kP+kL)
    cw = results.frames[:concentration_water]     # N × nf
    vb = results.frames[:velocity_biofilm]        # (N-1) × nf
    N = size(cb, 1)

    gc = GlobalConcentration(
        cb[:, idx, 1:kP],             # matrix
        cb[:, idx, kP+1:2kP],         # enclosed particles
        cf[:, idx, 1:kP],             # flowing particles
        cb[:, idx, 2kP+1:2kP+kL],     # enclosed liquids
        cf[:, idx, kP+1:kP+kL],       # flowing liquids
    )
    # water frame stores densityL * phiW (mean liquid density); invert it
    phiW = cw[:, idx] ./ _group_density(L)
    vel = Velocity(vb[:, idx], zeros(N + 1))   # flowing velocity is recomputed
    time = min(ft[idx], results.time_final)
    return State(f, m, Float64(time), gc, phiW, vel)
end

"Concatenate along the frame axis, dropping `b`'s first (duplicated) frame."
_catframes(a::AbstractArray{T,3}, b::AbstractArray{T,3}) where {T} = cat(a, b[:, 2:end, :]; dims=2)
_catframes(a::AbstractMatrix, b::AbstractMatrix) = hcat(a, b[:, 2:end])
_catframes(a::AbstractVector, b::AbstractVector) = vcat(a, b[2:end])

"""
    concatenate(r1::Results, r2::Results) -> Results
    r1 + r2

Append `r2` onto `r1`, returning a single Results spanning both. Requires the
same filter and model and `r1.time_final == r2.time_start` (the natural result
of `r2 = simulate(final_state(r1); ...)`). `r2`'s first frame — identical to
`r1`'s last — is dropped. Port of `@SDresults` `concatenate`/`plus`.
"""
function concatenate(r1::Results, r2::Results)
    r1.flag == "UNINITIATED" && return r2
    r2.flag == "UNINITIATED" && return r1
    (r1.filter === r2.filter && r1.model === r2.model) ||
        error("concatenate: filter/model differ between the two Results")
    isapprox(r1.time_final, r2.time_start; atol=0, rtol=1e-12) ||
        error("concatenate: r2.time_start ($(r2.time_start)) ≠ r1.time_final ($(r1.time_final))")

    out = Results(r1.filter, r1.model)
    out.time_start = r1.time_start
    out.time_final = r2.time_final
    out.flag = r2.flag
    for key in keys(r1.frames)
        haskey(r2.frames, key) || continue
        out.frames[key] = _catframes(r1.frames[key], r2.frames[key])
    end
    # carry r1's metadata, then splice the per-step time vector if both have it
    merge!(out.simulation_data, r1.simulation_data)
    if haskey(r1.simulation_data, :step_times) && haskey(r2.simulation_data, :step_times)
        out.simulation_data[:step_times] = _catframes(r1.simulation_data[:step_times],
                                                       r2.simulation_data[:step_times])
    end
    out.simulation_data[:time_final] = r2.time_final
    return out
end

Base.:+(r1::Results, r2::Results) = concatenate(r1, r2)
