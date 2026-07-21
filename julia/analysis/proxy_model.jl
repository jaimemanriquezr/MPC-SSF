# Fast evaluation of modelLund — osmosis kept PHYSICAL, via the implicit_osmosis
# solver flag. Two tiers:
#
#   FAITHFUL (recommended):  run_proxy(State(f, modelLund()); ...)
#     Plain modelLund + implicit_osmosis + a capped dt. modelLund's step is bound
#     almost entirely by the osmosis 1/τ (=1e7 at τ=1e-7); integrating that linear
#     relaxation implicitly removes the bound while τ and β stay EXACTLY physical.
#     Measured: ~20× fewer steps than explicit, matching it to ~1e-4 relative at
#     t=0.15 (biofilm developed). Physics unchanged — this is the model, just faster.
#
#   AGGRESSIVE (lower fidelity):  run_proxy(State(f, modelLund_proxy()); ...)
#     Also floors the tiny half-saturations and scales reaction rates, which relaxes
#     the *reaction* CFL so dt grows further (to ~1e-3). This changes the reaction
#     physics AND, past dt~5e-6, lets fast-depleting species overshoot negative
#     (trips the BIOFILM guard). Use only for plumbing/qualitative checks where a
#     large dt matters more than correctness.
#
# WHY the dt cap: implicit osmosis only makes the *osmosis* term unconditionally
# stable. The reaction/transport terms are still explicit, and their CFL slightly
# under-bounds negativity once osmosis no longer limits dt — verified: explicit
# modelLund runs clean to t=0.55, the uncapped-implicit run goes negative at ~0.47.
# adaptive_max_dt=3e-6 sits safely below that (for ~30 cells) and still gives ~20×.

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

AGGRESSIVE proxy of `modelLund` (osmosis still physical: τ=1e-7, β=0.99). Reaction
rates ×`rate_factor` and half-saturations floored at `halfsat_floor` to relax the
reaction CFL for extra dt. Lower fidelity than plain `modelLund` and can trip the
negativity guard at very large dt — prefer plain `modelLund` with [`run_proxy`](@ref)
unless you specifically need dt beyond ~3e-6.
"""
function modelLund_proxy(; rate_factor=0.1, halfsat_floor=1e-3)
    m = modelLund()
    rxs = Reaction[_remake_rx(r; nominal_rate=rate_factor*r.nominal_rate,
                              half_saturation_constants=_floor_K(r.half_saturation_constants, halfsat_floor))
                   for r in m.reactions]
    Model(m.components, rxs; cohesion_submodel=m.cohesion_submodel, water_density=m.water_density,
          biofilm_porosity=m.biofilm_porosity, osmosis_rate=m.osmosis_rate, detachment=m.detachment)
end

"""
    run_proxy(state; simulation_time, inflow_concentrations, kwargs...) -> Results

Simulate with `implicit_osmosis=true` and a dt capped at `adaptive_max_dt` (default
3e-6, the faithful ceiling for ~30 cells). On plain `modelLund` this is ~20× faster
than explicit with physics exact. Extra kwargs pass through to / override `simulate`.
"""
function run_proxy(state; simulation_time, inflow_concentrations,
                   cfl_factor=0.99, adaptive_initial_dt=1e-8, adaptive_max_dt=3e-6,
                   n_frames=200, quiet=false, kwargs...)
    simulate(state; inflow_concentrations, simulation_time, time_step=:adaptive,
             cfl_factor, adaptive_initial_dt, adaptive_max_dt, n_frames, quiet,
             implicit_osmosis=true, kwargs...)
end
