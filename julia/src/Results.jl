# Port of src/@Results/{Results,getVolumeFractions,plotBiofilm,
#                        plotConcentration,plotConcentrations}.m
#
# Post-processing container: captured frames (time, concentrations, velocities),
# run metadata, and a status flag. Plotting is deferred to a later phase and
# will likely use Plots.jl or Makie.jl (decision pending).
# TODO: port getVolumeFractions and the plot_* methods.

"""
    Results

Simulation output: per-frame concentration/velocity snapshots plus run metadata
and a status `flag` ("OK", "CLOGGED", "BIOFILM", "FLOWING", ...).
"""
Base.@kwdef mutable struct Results
    filter::SandFilter
    model::Model
    time_start::Float64 = 0.0
    flag::String = "UNINITIATED"
    frames::Dict{Symbol,Any} = Dict{Symbol,Any}()
    simulation_data::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Results(filter::SandFilter, model::Model) = Results(; filter, model)

"Port of getVolumeFractions.m."
get_volume_fractions(r::Results) =
    error("not yet ported — see src/@Results/getVolumeFractions.m")
