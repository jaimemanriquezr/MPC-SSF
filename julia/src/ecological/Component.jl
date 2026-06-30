# Port of src/ecological/Component.m
#
# MATLAB used a `Component` base class (matlab.mixin.Heterogeneous) with shared
# fields Name, Density, Dispersivity, TransportRate. Idiomatic Julia models this
# as an abstract supertype; the shared fields live on the concrete immutable
# structs (Particle, Liquid) so the compiler sees concrete field types.
#
# The MATLAB custom-display methods are presentation-only and intentionally not
# ported (a `Base.show` overload can be added later if desired).

"""
    Component

Abstract supertype for ecological components: [`Particle`](@ref) (particulate)
and [`Liquid`](@ref) (dissolved). Concrete subtypes carry, at minimum:

  - `name::String`
  - `density::Float64`       — kg/m³
  - `dispersivity::Float64`  — m
  - `transport_rate::Float64`— 1/day
"""
abstract type Component end
