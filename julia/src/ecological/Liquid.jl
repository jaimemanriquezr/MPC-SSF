# Port of src/ecological/Liquid.m  (Liquid < Component)
#
# MATLAB Liquid adds no fields beyond Component; it exists as a distinct type so
# the model can dispatch on particulate vs. dissolved components.

"""
    Liquid(; name, density, dispersivity=0.0, transport_rate=0.0)

Dissolved (liquid-phase) ecological component (e.g. dissolved organic matter
`DOM`, oxygen, substrate).

Units: `density` kg/m³, `dispersivity` m, `transport_rate` 1/day.
"""
Base.@kwdef struct Liquid <: Component
    name::String
    density::Float64
    dispersivity::Float64 = 0.0
    transport_rate::Float64 = 0.0
end
