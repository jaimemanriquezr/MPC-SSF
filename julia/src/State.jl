# Port of src/@State/State.m
#
# State container: holds the filter, model, time, and the global concentration
# arrays split across regions (matrix, enclosed particles/liquids, flowing
# particles/liquids), enclosed water volume, and biofilm/flowing velocities.
# Dependent properties (GlobalConcentrationBiofilm/Flowing, VolumeFractions) in
# MATLAB become accessor functions here.
# TODO: port the constructor's initial-condition sizing and the volume-fraction
# accessors; validate against State.m.

"Region-split global concentrations (port of the GlobalConcentration struct)."
Base.@kwdef mutable struct GlobalConcentration
    matrix::Matrix{Float64}             # N x kP
    enclosed_particles::Matrix{Float64} # N x kP
    flowing_particles::Matrix{Float64}  # N x kP
    enclosed_liquids::Matrix{Float64}   # N x kL
    flowing_liquids::Matrix{Float64}    # N x kL
end

"Biofilm and flowing-phase velocities on cell faces."
Base.@kwdef mutable struct Velocity
    biofilm::Vector{Float64}   # length N-1
    flowing::Vector{Float64}   # length N+1
end

"""
    State(filter::SandFilter, model::Model; preset=:clean)

Simulation state. `preset=:clean` initialises an unclogged filter (port of the
`preset="clean"` path in State.m).
"""
Base.@kwdef mutable struct State
    filter::SandFilter
    model::Model
    time::Float64 = 0.0
    global_concentration::Union{GlobalConcentration,Nothing} = nothing
    enclosed_water_volume::Union{Vector{Float64},Nothing} = nothing
    velocity::Union{Velocity,Nothing} = nothing
end

function State(filter::SandFilter, model::Model; preset::Symbol=:clean)
    error("State constructor not yet ported — see src/@State/State.m")
end

"Concentrations in the biofilm region [matrix, enclosed P, enclosed L]."
global_concentration_biofilm(s::State) =
    error("not yet ported — see State.m get.GlobalConcentrationBiofilm")

"Concentrations in the flowing region [flowing P, flowing L]."
global_concentration_flowing(s::State) =
    error("not yet ported — see State.m get.GlobalConcentrationFlowing")

"Per-region volume fractions (port of get.VolumeFractions)."
volume_fractions(s::State) =
    error("not yet ported — see State.m get.VolumeFractions")
