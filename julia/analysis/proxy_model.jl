# Fast PROXY of modelLund for quick evaluation runs — NOT physically faithful.
# modelLund is extremely stiff (dt~1e-7 to 1e-9): the CFL step is throttled by
# (a) osmosis 1/τ (τ=1e-7), (b) reaction rates, and (c) the tiny half-saturations
# (e.g. HPO4 K=1.4e-8), which make the liquid-consumption CFL weight X/(S+K) blow
# up as nutrients deplete. This proxy relaxes all three so runs finish in seconds,
# for wiring/plumbing/qualitative checks — do NOT use it for results.
#
#   modelLund_proxy(; rate_factor=0.1, tau=1e-3, halfsat_floor=1e-3, beta=0.90)
# rate_factor scales every reaction rate; tau is the osmosis time; halfsat_floor
# raises every (Monod) half-saturation to at least this value (kills the K→0 spike);
# beta is the biofilm porosity. modelLund uses beta=0.99, whose osmotic
# amplification β/(1−β)=99 makes φ_enclosed=99·φ_matrix — so when osmosis is
# SLOWED (large tau) the enclosed phase overshoots and the filter clogs (flag
# BIOFILM) within a day. Lowering beta shrinks that amplification so the relaxed
# proxy runs the full window without clogging. This is the trade you accept for
# speed: the proxy is NOT physically faithful.

using MPCSSF

_remake_rx(r::Reaction; kw...) = Reaction(; name=r.name, (; nominal_rate=r.nominal_rate,
    temperature_correction_factor=r.temperature_correction_factor, order=r.order,
    half_saturation_constants=r.half_saturation_constants,
    stoichiometric_coefficients=r.stoichiometric_coefficients, efficiency_biofilm=r.efficiency_biofilm,
    efficiency_flowing=r.efficiency_flowing, is_light_dependent=r.is_light_dependent,
    minimum_light_factor=r.minimum_light_factor, optimal_light_factor=r.optimal_light_factor, kw...)...)

# raise Monod half-sats to a floor (leave "num/den" quotient keys untouched)
_floor_K(d, floor) = Dict{String,Float64}(k => occursin("/", k) ? v : max(v, floor) for (k, v) in d)

"""
    modelLund_proxy(; rate_factor=0.1, tau=1e-3, halfsat_floor=1e-3, beta=0.90) -> Model

Fast, non-physical proxy of `modelLund`. Same components/stoichiometry, but
reaction rates ×`rate_factor`, osmosis time `tau`, half-saturations floored at
`halfsat_floor`, and biofilm porosity lowered to `beta` (from 0.99) so the
slowed osmosis does not clog. For quick evaluation / pipeline checks only.
"""
function modelLund_proxy(; rate_factor=0.1, tau=1e-3, halfsat_floor=1e-3, beta=0.90)
    m = modelLund()
    rxs = Reaction[_remake_rx(r; nominal_rate=rate_factor*r.nominal_rate,
                              half_saturation_constants=_floor_K(r.half_saturation_constants, halfsat_floor))
                   for r in m.reactions]
    Model(m.components, rxs; cohesion_submodel=m.cohesion_submodel, water_density=m.water_density,
          biofilm_porosity=beta, osmosis_rate=tau, detachment=m.detachment)
end
