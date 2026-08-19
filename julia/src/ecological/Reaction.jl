# Port of src/ecological/@Reaction/{Reaction,computeRate,lookupOrder,
#   lookupHalfSaturationConstants,lookupStoichiometricCoefficients,
#   lookupQuotients}.m
#
# A biochemical reaction: nominal rate + temperature correction, kinetic orders,
# Monod half-saturation constants, stoichiometry, biofilm/flowing efficiencies,
# and light dependence. The MATLAB dictionaries are keyed by component name (or
# by Component objects, which the constructor converts to names); we store them
# as `Dict{String,Float64}` and accept either key kind.

"""
    Reaction(; name="", nominal_rate=1.0, temperature_correction_factor=1.0,
               order=Dict(), half_saturation_constants=Dict(),
               stoichiometric_coefficients=Dict(),
               efficiency_biofilm=1.0, efficiency_flowing=1.0,
               is_light_dependent=false,
               minimum_light_factor=0.0, optimal_light_factor=0.0)

A single ecological reaction. The three dictionaries are keyed by component name;
keys may instead be [`Component`](@ref) objects (their `.name` is used).
`half_saturation_constants` may additionally contain quotient/inhibition keys of
the form `"numerator/denominator"` (see [`lookup_quotients`](@ref)).
"""
struct Reaction
    name::String
    nominal_rate::Float64
    temperature_correction_factor::Float64
    order::Dict{String,Float64}
    half_saturation_constants::Dict{String,Float64}
    stoichiometric_coefficients::Dict{String,Float64}
    efficiency_biofilm::Float64
    efficiency_flowing::Float64
    is_light_dependent::Bool
    minimum_light_factor::Float64
    optimal_light_factor::Float64
    # Dark-switch (Wolf2007 PHOBIA r6): when > 0, the reaction's light factor is
    # K/(K + I_local) with I_local the attenuated intensity in optimal-intensity
    # units — active in darkness, suppressed in light. Mutually exclusive with
    # is_light_dependent. 0.0 disables (default; goldens unchanged).
    light_inhibition::Float64
    # Complement mode: light factor = 1 − Steele(I) — exactly the complement of
    # the light-dependent growth factor (with its floor at 0), so the reaction
    # activates where and when photosynthesis idles. Mutually exclusive with
    # the other two light modes.
    is_light_complement::Bool
end

# Normalize a dict whose keys are either Strings or Components into String keys.
_namekeys(d::AbstractDict) =
    Dict{String,Float64}((k isa Component ? k.name : String(k)) => Float64(v)
                         for (k, v) in d)

function Reaction(; name::AbstractString="",
                  nominal_rate::Real=1.0,
                  temperature_correction_factor::Real=1.0,
                  order::AbstractDict=Dict{String,Float64}(),
                  half_saturation_constants::AbstractDict=Dict{String,Float64}(),
                  stoichiometric_coefficients::AbstractDict=Dict{String,Float64}(),
                  efficiency_biofilm::Real=1.0,
                  efficiency_flowing::Real=1.0,
                  is_light_dependent::Bool=false,
                  minimum_light_factor::Real=0.0,
                  optimal_light_factor::Real=0.0,
                  light_inhibition::Real=0.0,
                  is_light_complement::Bool=false)
    count((is_light_dependent, light_inhibition > 0, is_light_complement)) > 1 &&
        throw(ArgumentError("light modes (dependent/inhibited/complement) are mutually exclusive"))
    return Reaction(String(name), nominal_rate, temperature_correction_factor,
                    _namekeys(order), _namekeys(half_saturation_constants),
                    _namekeys(stoichiometric_coefficients),
                    efficiency_biofilm, efficiency_flowing,
                    is_light_dependent, minimum_light_factor, optimal_light_factor,
                    light_inhibition, is_light_complement)
end

"""
    compute_rate(rx::Reaction, temperature; scale=:celsius, nominal_temperature=20)

Temperature-corrected reaction rate:

    μ = μ₂₀ · θ^(T − T_nom)

with `μ₂₀ = rx.nominal_rate`, `θ = rx.temperature_correction_factor`. `θ`
(≈ 1.05–1.08) is calibrated against this temperature *difference* form
(θ^(T − 20 °C)), so the exponent is the difference, not the dimensionless ratio
`T/T_nom − 1`. Because it is a difference, the 273 K shift cancels: for
`scale = :celsius` the exponent is `temperature − nominal_temperature` in °C; for
`scale = :kelvin` pass both in kelvin (`nominal_temperature = 293`). Broadcast
over a vector of reactions: `compute_rate.(reactions, T)`.

NOTE: earlier ports (and MPC-SSF `computeRate.m`) used the dimensionless ratio
`θ^(T/T_nom − 1)`, which makes the response ~293× too weak unless θ is
recalibrated to θ^293; the authoritative slow-sand `run_biofilm.m` uses
`θ^(293 − T_K)`, which has the sign flipped. Both coincide with this form only at
the 20 °C reference.
"""
function compute_rate(rx::Reaction, temperature::Real;
                      scale::Symbol=:celsius, nominal_temperature::Real=20)
    μ20 = rx.nominal_rate
    θ = rx.temperature_correction_factor
    if scale === :celsius
        T = 273 + temperature
        Tnom = 273 + nominal_temperature
    elseif scale === :kelvin
        T = temperature
        Tnom = nominal_temperature
    else
        throw(ArgumentError("Invalid temperature scale $scale (use :celsius or :kelvin)"))
    end
    return μ20 * θ^(T - Tnom)
end

# --- (component × reaction) lookup matrices ---------------------------------
# Generic helper: build an (n_components × n_reactions) matrix by looking each
# component name up in `field(rx)` for every reaction, using `fallback` when the
# component is absent.
function _lookup_matrix(reactions::AbstractVector{Reaction},
                        components::AbstractVector{<:Component},
                        field::Function, fallback::Float64)
    names = [c.name for c in components]
    M = fill(fallback, length(names), length(reactions))
    for (j, rx) in pairs(reactions)
        d = field(rx)
        for (i, nm) in pairs(names)
            haskey(d, nm) && (M[i, j] = d[nm])
        end
    end
    return M
end

"Kinetic orders as an (n_components × n_reactions) matrix; absent ⇒ 0 (port of lookupOrder.m)."
lookup_order(reactions, components) =
    _lookup_matrix(reactions, components, rx -> rx.order, 0.0)

"Half-saturation constants, (n_components × n_reactions); absent ⇒ NaN (port of lookupHalfSaturationConstants.m)."
lookup_half_saturation_constants(reactions, components) =
    _lookup_matrix(reactions, components, rx -> rx.half_saturation_constants, NaN)

"Stoichiometric coefficients, (n_components × n_reactions); absent ⇒ 0 (port of lookupStoichiometricCoefficients.m)."
lookup_stoichiometric_coefficients(reactions, components) =
    _lookup_matrix(reactions, components, rx -> rx.stoichiometric_coefficients, 0.0)

"""
    lookup_quotients(reactions, components)
        -> (K, den_idx, num_idx)

Port of `lookupQuotients.m`. Scans each reaction's `half_saturation_constants`
for quotient/inhibition keys of the form `"numerator/denominator"`. Returns:

  - `K`        : (n_quotients × n_reactions) matrix; row `i` has its value placed
                 in the column of the reaction that owns quotient `i`, NaN else.
  - `den_idx`  : component index of each quotient's denominator.
  - `num_idx`  : component index of each quotient's numerator.

The MATLAB tuple order is `[K, denIdx, numIdx]`; this returns a NamedTuple with
those names for clarity.
"""
function lookup_quotients(reactions::AbstractVector{Reaction},
                          components::AbstractVector{<:Component})
    names = [c.name for c in components]
    nameindex(s) = something(findfirst(==(s), names), 0)

    rx_ids = Int[]
    values = Float64[]
    num_idx = Int[]
    den_idx = Int[]
    for (j, rx) in pairs(reactions)
        for (key, val) in rx.half_saturation_constants
            occursin("/", key) || continue
            num, den = split(key, "/"; limit=2)
            push!(rx_ids, j)
            push!(values, val)
            push!(num_idx, nameindex(String(num)))
            push!(den_idx, nameindex(String(den)))
        end
    end

    nq = length(rx_ids)
    K = fill(NaN, nq, length(reactions))
    for i in 1:nq
        K[i, rx_ids[i]] = values[i]
    end
    return (K=K, den_idx=den_idx, num_idx=num_idx)
end
