# Port of src/ecological/Particle.m
#
# Particulate ecological component. MATLAB fields (verify against Particle.m):
#   Name, Density, Dispersivity, AttachmentSand, TransportRate, Attenuation
# TODO: confirm field names/defaults against src/ecological/Particle.m and
# presets/modelLund.m before relying on these.

"""
    Particle(; name, density, dispersivity=0.0, attachment_sand=0.0,
               transport_rate=0.0, attenuation=0.0)

Particulate ecological component (e.g. HET, PHO, POM, pathogen).
"""
Base.@kwdef struct Particle <: Component
    name::String
    density::Float64
    dispersivity::Float64 = 0.0
    attachment_sand::Float64 = 0.0
    transport_rate::Float64 = 0.0
    attenuation::Float64 = 0.0
end
