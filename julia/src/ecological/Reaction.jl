# Port of src/ecological/@Reaction/*.m
#
# A biochemical reaction: nominal rate, stoichiometry, kinetic order, half-
# saturation (Monod) constants, and quotient (inhibition) terms. The MATLAB
# @Reaction folder splits this across:
#   Reaction.m, computeRate.m, lookupOrder.m, lookupHalfSaturationConstants.m,
#   lookupStoichiometricCoefficients.m, lookupQuotients.m
#
# In Julia we store stoichiometry/order/half-saturation as Dicts keyed by
# component name and provide lookup_* functions that resolve them against an
# ordered component list (matching how simulate.m assembles matrices).
#
# TODO: port computeRate (temperature correction: rate * Q10^(293 - T)) and the
# lookup_* helpers; validate against @Reaction outputs.

"""
    Reaction(; nominal_rate, stoichiometric_coefficients=Dict{String,Float64}(),
               order=Dict{String,Float64}(),
               half_saturation=Dict{String,Float64}(),
               temperature_factor=1.0,
               is_light_dependent=false, optimal_light_factor=1.0,
               minimum_light_factor=0.0)

A single ecological reaction. `stoichiometric_coefficients`, `order`, and
`half_saturation` are keyed by component name.
"""
Base.@kwdef struct Reaction
    nominal_rate::Float64
    stoichiometric_coefficients::Dict{String,Float64} = Dict{String,Float64}()
    order::Dict{String,Float64} = Dict{String,Float64}()
    half_saturation::Dict{String,Float64} = Dict{String,Float64}()
    temperature_factor::Float64 = 1.0          # Q10 / Arrhenius-style correction
    is_light_dependent::Bool = false
    optimal_light_factor::Float64 = 1.0
    minimum_light_factor::Float64 = 0.0
end

# --- ports of the @Reaction methods (stubs) ---------------------------------

"""
    compute_rate(rx::Reaction, temperature)

Temperature-corrected nominal rate: `nominal_rate * temperature_factor^(293 - T)`
(port of computeRate.m). TODO: confirm the exact correction form against MATLAB.
"""
function compute_rate(rx::Reaction, temperature)
    error("compute_rate not yet ported — see src/ecological/@Reaction/computeRate.m")
end

"Resolve stoichiometric coefficients into a vector aligned with `components`."
lookup_stoichiometric_coefficients(rx::Reaction, components) =
    error("not yet ported — see lookupStoichiometricCoefficients.m")

"Resolve kinetic orders into a vector aligned with `components`."
lookup_order(rx::Reaction, components) =
    error("not yet ported — see lookupOrder.m")

"Resolve half-saturation constants aligned with `components`."
lookup_half_saturation_constants(rx::Reaction, components) =
    error("not yet ported — see lookupHalfSaturationConstants.m")

"Resolve quotient (inhibition) terms; returns (K, denominator idx, numerator idx)."
lookup_quotients(rx::Reaction, components) =
    error("not yet ported — see lookupQuotients.m")
