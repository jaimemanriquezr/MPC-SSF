# Port of src/@Model/{Model,computeReactionRates,isequal,plot}.m
#
# Assembles the ecological model: component lists, reactions, cohesion submodel,
# stoichiometric matrices, half-saturation constants, kinetic orders, and the
# global parameters (biofilm porosity beta, osmosis rate tau, water density).
# TODO: port computeReactionRates and the derived stoichiometric/half-saturation
# matrix assembly used by simulate.m; validate against presets/modelLund.m.

"""
    Model(; components=Component[], reactions=Reaction[],
            cohesion_submodel=nothing, biofilm_porosity=0.99,
            osmosis_rate=1e-5, density_water=1000.0)

Assembled ecological model. `components` mixes [`Particle`](@ref) and
[`Liquid`](@ref); use [`particles`](@ref)/[`liquids`](@ref) to filter.
"""
Base.@kwdef mutable struct Model
    components::Vector{Component} = Component[]
    reactions::Vector{Reaction} = Reaction[]
    cohesion_submodel::Union{CahnHilliardModel,Nothing} = nothing
    biofilm_porosity::Float64 = 0.99    # beta
    osmosis_rate::Float64 = 1e-5        # tau
    density_water::Float64 = 1000.0
end

particles(m::Model) = filter(c -> c isa Particle, m.components)
liquids(m::Model)   = filter(c -> c isa Liquid,   m.components)

"Port of computeReactionRates.m — temperature-corrected reaction-rate vector."
compute_reaction_rates(m::Model, temperature) =
    error("not yet ported — see src/@Model/computeReactionRates.m")

"Stoichiometric matrix restricted to particle components (port from simulate.m)."
stoichiometric_matrix_particles(m::Model) =
    error("not yet ported — see @Model + simulate.m assembly")

"Stoichiometric matrix restricted to liquid components."
stoichiometric_matrix_liquids(m::Model) =
    error("not yet ported — see @Model + simulate.m assembly")
