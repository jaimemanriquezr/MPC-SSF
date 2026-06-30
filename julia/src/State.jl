# Port of src/@State/State.m
#
# State container: filter + model + the global concentration arrays split across
# regions (biofilm matrix, enclosed particles/liquids, flowing particles/liquids),
# the enclosed water volume, and the biofilm/flowing face velocities. The MATLAB
# constructor zero-allocates everything (a "clean" filter); there is no preset
# argument in the current src/ (the rollback dropped it).
#
# Dependent properties (GlobalConcentrationBiofilm/Flowing, VolumeFractions) become
# accessor functions. State is immutable, but its array fields are mutable so the
# solver can fill in / update them.

"Region-split global concentrations (port of the GlobalConcentration struct). Each is N × k."
struct GlobalConcentration
    matrix::Matrix{Float64}             # N × kP
    enclosed_particles::Matrix{Float64} # N × kP
    flowing_particles::Matrix{Float64}  # N × kP
    enclosed_liquids::Matrix{Float64}   # N × kL
    flowing_liquids::Matrix{Float64}    # N × kL
end

"Biofilm and flowing-phase velocities on cell faces."
struct Velocity
    biofilm::Vector{Float64}   # length N-1
    flowing::Vector{Float64}   # length N+1
end

"""
    State(filter::SandFilter, model::Model; time=0.0)

Simulation state for a discretized `filter`. All concentration/velocity arrays
are zero-initialized (a clean filter). The filter must already have a grid
(`addgridpoints`); otherwise an error is thrown.
"""
struct State
    filter::SandFilter
    model::Model
    time::Float64
    global_concentration::GlobalConcentration
    enclosed_water_volume::Vector{Float64}
    velocity::Velocity
end

function State(filter::SandFilter, model::Model; time::Real=0.0)
    filter.grid === nothing &&
        error("filter has no grid; call addgridpoints(filter, n) before building a State")
    N = length(filter.grid.centers)
    kP = length(particles(model))
    kL = length(liquids(model))

    gc = GlobalConcentration(zeros(N, kP), zeros(N, kP), zeros(N, kP),
                             zeros(N, kL), zeros(N, kL))
    vel = Velocity(zeros(N - 1), zeros(N + 1))
    return State(filter, model, Float64(time), gc, zeros(N), vel)
end

# --- dependent accessors -----------------------------------------------------

"Concentrations in the biofilm region: [matrix, enclosed particles, enclosed liquids] (N × 2kP+kL)."
global_concentration_biofilm(s::State) =
    hcat(s.global_concentration.matrix,
         s.global_concentration.enclosed_particles,
         s.global_concentration.enclosed_liquids)

"Concentrations in the flowing region: [flowing particles, flowing liquids] (N × kP+kL)."
global_concentration_flowing(s::State) =
    hcat(s.global_concentration.flowing_particles,
         s.global_concentration.flowing_liquids)

"""
    volume_fractions(s::State) -> NamedTuple

Per-region volume fractions: each region's concentration divided column-wise by
the corresponding component density (port of get.VolumeFractions). Returns
fields `matrix`, `enclosed_particles`, `enclosed_liquids`, `flowing_particles`,
`flowing_liquids`.
"""
function volume_fractions(s::State)
    pdens = [p.density for p in particles(s.model)]   # length kP
    ldens = [l.density for l in liquids(s.model)]     # length kL
    gc = s.global_concentration
    return (matrix             = gc.matrix ./ pdens',
            enclosed_particles = gc.enclosed_particles ./ pdens',
            enclosed_liquids   = gc.enclosed_liquids ./ ldens',
            flowing_particles  = gc.flowing_particles ./ pdens',
            flowing_liquids    = gc.flowing_liquids ./ ldens')
end
