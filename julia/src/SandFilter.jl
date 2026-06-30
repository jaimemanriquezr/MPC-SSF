# Port of src/@SandFilter/{SandFilter,addGridPoints,computePorosity,
#                            getCahnHilliardMatrices}.m
#
# Filter geometry, 1D staggered grid, porosity profile, light attenuation, and
# Cahn-Hilliard matrix assembly.
#
# MATLAB stored GridPoints as a struct with Boundaries + Centers and several
# Dependent properties (GridSize, GridZero, LightAttenuationEta{Water,Sand}).
# Here those become a Grid struct + accessor functions.
# TODO: port computePorosity, light-attenuation accessors, and
# getCahnHilliardMatrices; validate grid against addGridPoints.m.

"Staggered 1D grid: cell boundaries and centers (port of GridPoints struct)."
struct Grid
    boundaries::Vector{Float64}
    centers::Vector{Float64}
end

"""
    SandFilter(; height=1.0, depth=1.0, sand_porosity=0.4, sand_roughness=5e-3,
                 inflow_velocity=0.3*24, temperature=15+273,
                 light_irradiation=default_light, light_attenuation_water=0.32,
                 light_attenuation_sand=1500.0, grid=nothing)

Filter geometry and discretization. Defaults mirror SandFilter.m.
"""
Base.@kwdef mutable struct SandFilter
    height::Float64 = 1.0
    depth::Float64 = 1.0
    sand_porosity::Float64 = 0.4
    sand_roughness::Float64 = 5e-3
    inflow_velocity::Float64 = 0.3 * 24
    temperature::Float64 = 15 + 273
    light_irradiation::Function = default_light_irradiation
    light_attenuation_water::Float64 = 0.32
    light_attenuation_sand::Float64 = 1500.0
    grid::Union{Grid,Nothing} = nothing
end

"Default diel light forcing (port of the anonymous fn in SandFilter.m)."
default_light_irradiation(t) = 0.8 * max(sin(2π * (t - 13 / 48)) + 31 / 50, 0) / (1 + 31 / 50)

gridsize(f::SandFilter) = mean(diff(f.grid.boundaries))

"Port of addGridPoints.m — build the staggered grid with `n` interior cells."
function addgridpoints(f::SandFilter, n::Integer)
    error("addgridpoints not yet ported — see src/@SandFilter/addGridPoints.m")
end

"Port of computePorosity.m — porosity profile along depth z."
function computeporosity(f::SandFilter, z)
    error("computeporosity not yet ported — see src/@SandFilter/computePorosity.m")
end

"Port of getCahnHilliardMatrices.m — assemble (S, DD, D) sparse triplets."
function get_cahn_hilliard_matrices(f::SandFilter, model, upwinded::Bool=false)
    error("get_cahn_hilliard_matrices not yet ported — see getCahnHilliardMatrices.m")
end
