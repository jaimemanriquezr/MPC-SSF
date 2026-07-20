# Port of src/@Results/{Results,getVolumeFractions}.m
#
# Post-processing container. simulate fills `frames` with raw 3D arrays
# (depth × frame × component-channel); MATLAB instead exposes a component×volume
# table, so we provide accessor functions that slice the right channel for a
# given (component, region) pair, plus get_volume_fractions. Plotting (the many
# plot_* methods) is deferred to a later phase.
#
# Frame channel layout (set by simulate):
#   concentration_biofilm  : [matrix particles (1:kP) | enclosed particles
#                             (kP+1:2kP) | enclosed liquids (2kP+1:2kP+kL)]
#   concentration_flowing  : [flowing particles (1:kP) | flowing liquids (kP+1:kP+kL)]
#   concentration_water    : enclosed water (depth × frame)

"""
    Results

Simulation output: per-frame concentration/velocity snapshots in `frames`
(a `Dict{Symbol,Any}`), run metadata, and a status `flag` ("OK", "CLOGGED",
"BIOFILM", "FLOWING", ...). Use [`concentration`](@ref),
[`get_volume_fractions`](@ref), [`times`](@ref), [`depths`](@ref) to read it.
"""
Base.@kwdef mutable struct Results
    filter::SandFilter
    model::Model
    time_start::Float64 = 0.0
    time_final::Float64 = 0.0
    flag::String = "UNINITIATED"
    frames::Dict{Symbol,Any} = Dict{Symbol,Any}()
    simulation_data::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Results(filter::SandFilter, model::Model) = Results(; filter, model)

"Frame times (days)."
times(r::Results) = r.frames[:time]

"Cell-center depths z (m)."
depths(r::Results) = r.filter.grid.centers

"""
    reaction_rates(r::Results) -> depth × frame × reaction array

Per-frame biofilm ecological reaction rates φ·μ·I·monod·product (captured by
`simulate`). Reaction order/names follow `r.model.reactions`; use
[`reaction_names`](@ref) for the labels.
"""
reaction_rates(r::Results) = r.frames[:reaction_rates]

"Names of the model's reactions, in column order of [`reaction_rates`](@ref)."
reaction_names(r::Results) = [rx.name for rx in r.model.reactions]

"""
    concentration(r::Results, name, region) -> depth × frame matrix

Concentration field for component `name` in `region ∈ (:matrix, :enclosed,
:flowing)`. Liquids have no matrix concentration (returns zeros); use
`name = "Water", region = :enclosed` for the enclosed water.
"""
function concentration(r::Results, name::AbstractString, region::Symbol)
    P = particles(r.model)
    L = liquids(r.model)
    kP = length(P)
    cb = r.frames[:concentration_biofilm]
    cf = r.frames[:concentration_flowing]

    pj = findfirst(p -> p.name == name, P)
    if pj !== nothing
        region === :matrix   && return cb[:, :, pj]
        region === :enclosed && return cb[:, :, kP+pj]
        region === :flowing  && return cf[:, :, pj]
    end
    lj = findfirst(l -> l.name == name, L)
    if lj !== nothing
        region === :matrix   && return zeros(size(cb, 1), size(cb, 2))
        region === :enclosed && return cb[:, :, 2kP+lj]
        region === :flowing  && return cf[:, :, kP+lj]
    end
    (name == "Water" && region === :enclosed) && return r.frames[:concentration_water]
    error("no concentration field for component \"$name\" in region :$region")
end

"""
    get_volume_fractions(r::Results) -> (biofilm, enclosed, matrix)

Per-frame volume fractions (port of getVolumeFractions.m): the matrix fraction
φ_M, the enclosed fraction φ_e (enclosed particles + liquids + water), and the
total biofilm fraction φ_b = φ_M + φ_e. Each is a depth × frame matrix.
"""
function get_volume_fractions(r::Results)
    P = particles(r.model)
    L = liquids(r.model)
    kP = length(P)
    cb = r.frames[:concentration_biofilm]
    nz, nt = size(cb, 1), size(cb, 2)

    phiM = zeros(nz, nt)
    phie = zeros(nz, nt)
    for (j, p) in enumerate(P)
        phiM .+= cb[:, :, j] ./ p.density
        phie .+= cb[:, :, kP+j] ./ p.density
    end
    for (j, l) in enumerate(L)
        phie .+= cb[:, :, 2kP+j] ./ l.density
    end
    phie .+= r.frames[:concentration_water] ./ r.model.water_density
    phib = phiM .+ phie
    return (biofilm=phib, enclosed=phie, matrix=phiM)
end
