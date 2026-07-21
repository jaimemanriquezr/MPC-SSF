# Fast PROXY of modelLund for quick evaluation runs — osmosis kept PHYSICAL.
#
# modelLund is cluster-stiff (dt~1e-7 to 1e-9). Three things throttle the step:
#   (a) osmosis 1/τ (τ=1e-7) — the dominant bound,
#   (b) reaction rates,
#   (c) the tiny half-saturations (e.g. HPO4 K=1.4e-8), which blow up the
#       liquid-consumption CFL weight X/(S+K) as nutrients deplete.
#
# (a) is handled WITHOUT touching the osmosis physics: pass `implicit_osmosis=true`
# to `simulate` (this repo's opt-in solver flag). It integrates the linear
# relaxation (β·φB−φe)/τ with backward Euler — unconditionally stable — so 1/τ
# drops out of the CFL while τ and β stay exactly at their modelLund values. That
# alone buys ~14× larger dt; the model here then relaxes (b) and (c) for the rest.
#
#   modelLund_proxy(; rate_factor=0.1, halfsat_floor=1e-3)   # osmosis UNMODIFIED
# rate_factor scales every reaction rate; halfsat_floor raises every Monod
# half-saturation to at least this value (kills the K→0 CFL spike). Osmosis time
# τ and biofilm porosity β are left at the physical modelLund values (1e-7, 0.99).
#
# Use `run_proxy(...)` to simulate with `implicit_osmosis=true` already set.
# NOT physically faithful in the reaction rates / half-sats — for wiring, plumbing
# and qualitative checks only. (The physical influent still clogs a coarse filter
# in ~0.2–0.5 d; that clogging is real, not a proxy artifact.)

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
    modelLund_proxy(; rate_factor=0.1, halfsat_floor=1e-3) -> Model

Fast, non-physical proxy of `modelLund` that keeps osmosis PHYSICAL (τ=1e-7,
β=0.99). Reaction rates are ×`rate_factor` and half-saturations floored at
`halfsat_floor`. To actually run fast, simulate it with `implicit_osmosis=true`
(see [`run_proxy`](@ref)), which removes the osmosis dt bound without changing τ.
For quick evaluation / pipeline checks only.
"""
function modelLund_proxy(; rate_factor=0.1, halfsat_floor=1e-3)
    m = modelLund()
    rxs = Reaction[_remake_rx(r; nominal_rate=rate_factor*r.nominal_rate,
                              half_saturation_constants=_floor_K(r.half_saturation_constants, halfsat_floor))
                   for r in m.reactions]
    # osmosis_rate (τ) and biofilm_porosity (β) kept at the modelLund values.
    Model(m.components, rxs; cohesion_submodel=m.cohesion_submodel, water_density=m.water_density,
          biofilm_porosity=m.biofilm_porosity, osmosis_rate=m.osmosis_rate, detachment=m.detachment)
end

"""
    run_proxy(state; simulation_time, kwargs...) -> Results

Convenience wrapper: `simulate` with `implicit_osmosis=true` and adaptive stepping
tuned for the proxy (large `adaptive_max_dt`). Extra kwargs pass through to
`simulate` and override the defaults.
"""
function run_proxy(state; simulation_time, inflow_concentrations,
                   cfl_factor=0.99, adaptive_initial_dt=1e-8, adaptive_max_dt=1e-2,
                   n_frames=200, quiet=false, kwargs...)
    simulate(state; inflow_concentrations, simulation_time, time_step=:adaptive,
             cfl_factor, adaptive_initial_dt, adaptive_max_dt, n_frames, quiet,
             implicit_osmosis=true, kwargs...)
end
