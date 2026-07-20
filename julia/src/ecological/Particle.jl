# Port of src/ecological/Particle.m  (Particle < Component)
#
# MATLAB Particle adds AttachmentMatrix, AttachmentSand, Attenuation to the
# Component fields (Name, Density, Dispersivity, TransportRate).
#
# NOTE / TODO: the MATLAB Particle constructor declares an input named
# `Transport` (no "Rate") that is never assigned — almost certainly a leftover
# typo for `TransportRate`. We use `transport_rate` consistently here. Confirm
# against presets/modelLund.m that no preset relies on a field literally named
# `Transport`.

"""
    Particle(; name, density, dispersivity=0.0, transport_rate=0.0,
               attachment_matrix=0.0, attachment_sand=0.0, attenuation=0.0,
               sand_attachment_factor=1.0)

Particulate ecological component (e.g. heterotrophs `HET`, phototrophs `PHO`,
particulate organic matter `POM`, or a pathogen indicator `PAT`).

`sand_attachment_factor` scales this particle's attachment to *sand* (bare grains)
relative to its attachment onto existing *biofilm*, in the flowing phase. It is
the port of the MATLAB `ecological.sand_pathogen` parameter (applied there to the
pathogen only); the default `1.0` makes sand and biofilm attachment share the
same `attachment_sand` rate, reproducing the non-pathogen model exactly.

Units: `density` kg/m³, `dispersivity` m, `transport_rate`/`attachment_*` 1/day,
`attenuation` m²/kg, `sand_attachment_factor` dimensionless.
"""
Base.@kwdef struct Particle <: Component
    name::String
    density::Float64
    dispersivity::Float64 = 0.0
    transport_rate::Float64 = 0.0
    attachment_matrix::Float64 = 0.0
    attachment_sand::Float64 = 0.0
    attenuation::Float64 = 0.0
    sand_attachment_factor::Float64 = 1.0
end
