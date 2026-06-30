# Port of src/ecological/Component.m
#
# Abstract supertype for the ecological components (particles and liquids).
# MATLAB used a Component base class with shared fields; in Julia we use an
# abstract type + shared fields on the concrete structs (see Particle, Liquid).

"""
    Component

Abstract supertype for ecological components ([`Particle`](@ref),
[`Liquid`](@ref)).
"""
abstract type Component end

name(c::Component) = c.name
density(c::Component) = c.density
