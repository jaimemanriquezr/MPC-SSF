# Port of src/ecological/Liquid.m
#
# Dissolved/liquid ecological component. MATLAB fields (verify against Liquid.m):
#   Name, Density, Dispersivity, TransportRate
# TODO: confirm field names/defaults against src/ecological/Liquid.m.

"""
    Liquid(; name, density, dispersivity=0.0, transport_rate=0.0)

Dissolved (liquid-phase) ecological component (e.g. DOM, O2, substrate).
"""
Base.@kwdef struct Liquid <: Component
    name::String
    density::Float64
    dispersivity::Float64 = 0.0
    transport_rate::Float64 = 0.0
end
