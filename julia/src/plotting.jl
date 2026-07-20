# Plot-function generics. The real methods live in the Makie package extension
# (ext/MPCSSFMakieExt.jl), which loads automatically when a Makie backend is
# imported (`using CairoMakie` or `using GLMakie`). Each returns a Makie Figure.
#
# The broad `args...` fallbacks below are less specific than the extension's
# `(::Results, ...)` methods, so once a backend is loaded the real methods win;
# until then, calling one gives a clear hint instead of a raw MethodError.

_needs_makie(fn) = error(
    "$fn: MPCSSF plotting requires a Makie backend. Run `using CairoMakie` " *
    "(or GLMakie) to load the plotting extension, then call the plot function again.")

"""
    plot_concentration(results, name, region; frames=:, kwargs...)

Depth profiles of component `name` in `region` (`:matrix`/`:enclosed`/`:flowing`)
for the selected frames. Requires a Makie backend (`using CairoMakie`).
"""
function plot_concentration end
plot_concentration(args...; kwargs...) = _needs_makie("plot_concentration")

"""
    plot_concentration_heatmap(results, name, region; kwargs...)

Depth × time heatmap of component `name` in `region`. Requires a Makie backend.
"""
function plot_concentration_heatmap end
plot_concentration_heatmap(args...; kwargs...) = _needs_makie("plot_concentration_heatmap")

"""
    plot_volume_fractions(results; frame=end, kwargs...)

Biofilm / enclosed / matrix volume-fraction depth profiles. Requires a Makie backend.
"""
function plot_volume_fractions end
plot_volume_fractions(args...; kwargs...) = _needs_makie("plot_volume_fractions")

"""
    plot_velocity(results; frame=end, kwargs...)

Biofilm and flowing-phase velocity depth profiles. Requires a Makie backend.
"""
function plot_velocity end
plot_velocity(args...; kwargs...) = _needs_makie("plot_velocity")

"""
    plot_cfl(results; kwargs...)

Time-step `dt` versus step index (from `simulation_data[:step_times]`) — the
adaptive-CFL diagnostic. Requires a Makie backend.
"""
function plot_cfl end
plot_cfl(args...; kwargs...) = _needs_makie("plot_cfl")

"""
    plot_reaction_rates(results; frame=end, kwargs...)

Biofilm reaction-rate depth profiles per reaction (biological activity).
Requires a Makie backend.
"""
function plot_reaction_rates end
plot_reaction_rates(args...; kwargs...) = _needs_makie("plot_reaction_rates")
