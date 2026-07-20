module MPCSSFMakieExt

# Makie backend for MPCSSF plotting. Loads automatically when a Makie backend
# (CairoMakie / GLMakie) is imported alongside MPCSSF. Each function returns a
# Makie `Figure`; save static output with `CairoMakie.save("f.png", fig)`.
#
# Convention: depth z is on the y-axis, reversed so the filter surface (z=-height)
# is at the top and increasing depth runs downward.

using MPCSSF
using Makie

const _CMAP = :viridis

# select frame indices to draw (default: all captured frames)
_frames(r, frames) = frames === nothing ? (1:length(MPCSSF.times(r))) : frames

# add a colorbar mapping frame time -> color, over `ts`
function _time_colorbar!(fig, ts; col=2)
    Colorbar(fig[1, col]; limits=(minimum(ts), maximum(ts)),
             colormap=_CMAP, label="time (day)")
end

# draw one depth profile per selected frame, colored by time
function _profiles!(ax, x, z, ts, framesel)
    crange = (minimum(ts), maximum(ts))
    for j in framesel
        lines!(ax, x[:, j], z; color=ts[j], colorrange=crange, colormap=_CMAP)
    end
end

# ---------------------------------------------------------------------------
function MPCSSF.plot_concentration(r::Results, name::AbstractString, region::Symbol;
                                   frames=nothing, resolution=(640, 520))
    z = depths(r); ts = times(r); data = concentration(r, name, region)
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="concentration — $name ($region)",
              ylabel="depth z (m)", yreversed=true,
              title="$name concentration ($region)")
    _profiles!(ax, data, z, ts, _frames(r, frames))
    _time_colorbar!(fig, ts)
    return fig
end

function MPCSSF.plot_concentration_heatmap(r::Results, name::AbstractString, region::Symbol;
                                           resolution=(640, 520))
    z = depths(r); ts = times(r); data = concentration(r, name, region)
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="time (day)", ylabel="depth z (m)", yreversed=true,
              title="$name concentration ($region)")
    hm = heatmap!(ax, ts, z, permutedims(data); colormap=_CMAP)
    Colorbar(fig[1, 2], hm; label="concentration")
    return fig
end

function MPCSSF.plot_volume_fractions(r::Results; frame=length(times(r)), resolution=(640, 520))
    z = depths(r); vf = get_volume_fractions(r)
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="volume fraction", ylabel="depth z (m)", yreversed=true,
              title="volume fractions @ t = $(round(times(r)[frame], sigdigits=3)) day")
    lines!(ax, vf.biofilm[:, frame], z; label="biofilm φ_b", color=:seagreen)
    lines!(ax, vf.matrix[:, frame], z; label="matrix φ_M", color=:sienna)
    lines!(ax, vf.enclosed[:, frame], z; label="enclosed φ_e", color=:steelblue)
    axislegend(ax; position=:rb)
    return fig
end

function MPCSSF.plot_velocity(r::Results; frame=length(times(r)), resolution=(640, 520))
    # velocities live on interior cell faces (length N-1)
    zf = r.filter.grid.boundaries[2:end-1]
    vb = r.frames[:velocity_biofilm][:, frame]
    vfl = r.frames[:velocity_flowing][:, frame]
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="velocity (m/day)", ylabel="depth z (m)", yreversed=true,
              title="velocities @ t = $(round(times(r)[frame], sigdigits=3)) day")
    lines!(ax, vb, zf; label="biofilm v_b", color=:seagreen)
    lines!(ax, vfl, zf; label="flowing v_f", color=:steelblue)
    axislegend(ax; position=:rb)
    return fig
end

function MPCSSF.plot_cfl(r::Results; resolution=(640, 420))
    st = r.simulation_data[:step_times]
    dt = diff(st)
    fig = Figure(size=resolution)
    mode = get(r.simulation_data, :time_step, "?")
    ax = Axis(fig[1, 1]; xlabel="step", ylabel="dt (day)", yscale=log10,
              title="time step (mode: $mode)")
    lines!(ax, 1:length(dt), dt; color=:crimson)
    return fig
end

function MPCSSF.plot_reaction_rates(r::Results; frame=length(times(r)), resolution=(680, 520))
    z = depths(r)
    rr = reaction_rates(r)[:, frame, :]      # N × nRx
    names = reaction_names(r)
    fig = Figure(size=resolution)
    ax = Axis(fig[1, 1]; xlabel="reaction rate (1/day)", ylabel="depth z (m)", yreversed=true,
              title="biofilm reaction rates @ t = $(round(times(r)[frame], sigdigits=3)) day")
    for i in eachindex(names)
        lines!(ax, rr[:, i], z; label=names[i])
    end
    axislegend(ax; position=:rb)
    return fig
end

end # module
