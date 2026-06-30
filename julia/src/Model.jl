# Port of src/@Model/{Model,computeReactionRates,isequal}.m
#
# Assembles the ecological model: component list, reactions, cohesion submodel,
# and global parameters (water density, biofilm porosity β, osmosis rate τ,
# detachment law). MATLAB exposes the derived quantities (particle/liquid lists,
# stoichiometric matrices, half-saturation constants, kinetic orders) as
# Dependent properties; here they are functions over the ecological-layer
# lookups, recomputed on demand (the values are cheap and not mutated).
#
# Defaults mirror Model.m: WaterDensity 998, BiofilmPorosity 0.99,
# OsmosisRate 1e-5, DetachmentFunction sqrt(|v|/qnom). The cohesion submodel has
# no usable default (CahnHilliardModel requires kappa/zeta_0), so it defaults to
# `nothing` and must be supplied (presets set it).

"""
    Model(; components=Component[], reactions=Reaction[],
            cohesion_submodel=nothing, water_density=998.0,
            biofilm_porosity=0.99, osmosis_rate=1e-5,
            detachment=(v -> zero(v)))
    Model(components, reactions=Reaction[]; kwargs...)

Assembled ecological model. `components` mixes [`Particle`](@ref) and
[`Liquid`](@ref); use [`particles`](@ref)/[`liquids`](@ref) to split them (order
preserved). `biofilm_porosity` is β, `osmosis_rate` is τ.
"""
Base.@kwdef struct Model
    components::Vector{Component} = Component[]
    reactions::Vector{Reaction} = Reaction[]
    cohesion_submodel::Union{CahnHilliardModel,Nothing} = nothing
    water_density::Float64 = 998.0
    biofilm_porosity::Float64 = 0.99    # β
    osmosis_rate::Float64 = 1e-5        # τ
    # Called by `simulate` as detachment(v) (one arg, matching simulate.m).
    # Default = no detachment; supply a closure, e.g. v -> sqrt.(abs.(v) ./ qnom).
    detachment::Function = (v -> zero(v))
end

# Ergonomic positional constructor: Model(components[, reactions]; kwargs...).
Model(components::AbstractVector{<:Component},
      reactions::AbstractVector{Reaction}=Reaction[]; kwargs...) =
    Model(; components=collect(Component, components),
            reactions=collect(Reaction, reactions), kwargs...)

# --- component views (port of get.Particles / get.Liquids) ------------------

"Particulate components, in declaration order."
particles(m::Model) = Particle[c for c in m.components if c isa Particle]

"Dissolved (liquid) components, in declaration order."
liquids(m::Model) = Liquid[c for c in m.components if c isa Liquid]

# --- derived matrices (ports of the Dependent properties) -------------------

"Reaction-rate vector at `temperature` (port of computeReactionRates.m)."
compute_reaction_rates(m::Model, temperature; kwargs...) =
    [compute_rate(r, temperature; kwargs...) for r in m.reactions]

"Stoichiometric coefficients over all components, (n_components × n_reactions)."
stoichiometric_coefficients(m::Model) =
    lookup_stoichiometric_coefficients(m.reactions, m.components)

"Stoichiometric matrix restricted to particles, (n_particles × n_reactions)."
stoichiometric_matrix_particles(m::Model) =
    lookup_stoichiometric_coefficients(m.reactions, particles(m))

"Stoichiometric matrix restricted to liquids, (n_liquids × n_reactions)."
stoichiometric_matrix_liquids(m::Model) =
    lookup_stoichiometric_coefficients(m.reactions, liquids(m))

"Half-saturation constants over all components, (n_components × n_reactions); absent ⇒ NaN."
half_saturation_constants(m::Model) =
    lookup_half_saturation_constants(m.reactions, m.components)

"Kinetic orders over all components, (n_components × n_reactions); absent ⇒ 0."
reaction_orders(m::Model) =
    lookup_order(m.reactions, m.components)

"Quotient/inhibition terms over all components (see [`lookup_quotients`](@ref))."
quotients(m::Model) =
    lookup_quotients(m.reactions, m.components)
