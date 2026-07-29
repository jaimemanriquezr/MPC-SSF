# Simplified heterotroph/phototroph model for the light–biofilm-profile study.
#
# Reduced from `modelLund()` to isolate the light response: two particulates
# (HET, PHO) and three liquids (O2, DOM, IC). NH4, HPO4 and POM are dropped, and
# death routes straight to DOM instead of via particulate matter, so the carbon
# loop survives without the hydrolysis step.
#
# This is deliberately NOT a speed optimisation. Legacy timings are exactly linear
# in cell count (1002 / 2019 / 4042 s per simulated day at 20 / 50 / 100 cells),
# which shows the step size is pinned by `adaptive_max_dt`, not by a CFL or
# half-saturation bound. The point is interpretability: fewer reactions competing
# for attention when reading a light signal, and no HPO4 starvation artifact
# (influent HPO4 = 0 against K = 1.4e-8; see analysis/SESSION_LOG_2026-07.md:165).
#
# Rates, stoichiometric magnitudes and the light parameters are kept at their Lund
# values so the light behaviour stays comparable to the full model. Dropping the
# N and P terms makes the model non-conservative in those elements, which is
# consistent with not tracking them.
#
# See .claude/plans/2026-07-29-light-biofilm-profile.md.

using MPCSSF

"Component order of `INFLUENT_LIGHT`: HET, PHO, O2, IC, DOM."
const INFLUENT_LIGHT = [1.0e-3, 1.0e-3, 1.0e-2, 1.0e-2, 1.0e-3]

"""
    light_profile_model() -> Model

Two particles (HET, PHO), three liquids (O2, IC, DOM), four reactions: heterotroph
and phototroph growth and death. Phototroph growth is the only light-dependent
reaction, carrying Lund's `optimal_light_factor = 1.814e-2` and
`minimum_light_factor = 0.01`.
"""
function light_profile_model()
    density_particle = 1.117e3
    attenuation_particle = 0.094
    dispersivity_particle = 1.20e-2

    HET = Particle(name="HET", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                   attenuation=attenuation_particle)
    PHO = Particle(name="PHO", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                   attenuation=attenuation_particle)

    density_liquid = 0.998e3
    dispersivity_liquid = 1.20e-2
    O2  = Liquid(name="O2",  density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e1)
    IC  = Liquid(name="IC",  density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e1)
    DOM = Liquid(name="DOM", density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e1)

    # Lund's heterotroph growth minus the NH4 and HPO4 Monod terms and their
    # stoichiometry; the O2/IC/DOM coefficients are unchanged.
    heterotroph_growth = Reaction(name="Heterotroph growth", is_light_dependent=false,
        nominal_rate=1.81e-2, temperature_correction_factor=1.047,
        order=Dict("HET" => 1.0),
        half_saturation_constants=Dict("O2" => 3.00e-3, "DOM" => 2.00e-4),
        stoichiometric_coefficients=Dict("HET" => 1.0, "O2" => -1.2317,
                                         "IC" => 0.3848, "DOM" => -1.5873))

    # The only light-dependent reaction. Inorganic carbon is the sole remaining
    # nutrient limitation.
    phototroph_growth = Reaction(name="Phototroph growth", is_light_dependent=true,
        minimum_light_factor=0.01, optimal_light_factor=1.814e-2,
        nominal_rate=5.50, temperature_correction_factor=1.047,
        order=Dict("PHO" => 1.0),
        half_saturation_constants=Dict("IC" => 2.00e-5),
        stoichiometric_coefficients=Dict("PHO" => 1.0, "O2" => 0.9301, "IC" => -0.3600))

    # Death sends biomass to DOM directly, standing in for the dropped
    # POM + hydrolysis pathway; the carbon magnitudes are Lund's POM yields.
    heterotroph_death = Reaction(name="Heterotroph death", is_light_dependent=false,
        nominal_rate=2.00, temperature_correction_factor=1.066,
        order=Dict("HET" => 1.0),
        stoichiometric_coefficients=Dict("HET" => -1.0, "DOM" => 0.9123, "O2" => 0.0234))

    phototroph_death = Reaction(name="Phototroph death", is_light_dependent=false,
        nominal_rate=4.00e-1, temperature_correction_factor=1.080,
        order=Dict("PHO" => 1.0),
        stoichiometric_coefficients=Dict("PHO" => -1.0, "DOM" => 0.6316, "O2" => 0.2005))

    components = Component[HET, PHO, O2, IC, DOM]
    reactions  = Reaction[heterotroph_growth, phototroph_growth,
                          heterotroph_death, phototroph_death]

    # kappa/zeta_0 as modelLund; zeta_1 = 1/100 feeds the published potential
    # dpsi/du = u^2(u - 3*zeta_1/2) (see .claude/decisions/2026-07-28-published-cohesion-potential.md).
    cohesion = CahnHilliardModel(kappa=1.00e-6, zeta_0=1.00e6, zeta_1=1/100)

    return Model(components, reactions;
                 cohesion_submodel=cohesion,
                 biofilm_porosity=0.99,
                 osmosis_rate=1.00e-7,
                 detachment=(v -> sqrt.(v ./ 7.2)))
end

"""
    light_forcing(amplitude) -> Function

Diel forcing with the default shape, rescaled to `amplitude` at its peak. The
default `SandFilter` forcing is this with `amplitude = 0.8`; sweeping it is how
the study varies incidence, per
`.claude/plans/2026-07-29-light-biofilm-profile.md`.
"""
light_forcing(amplitude::Real) =
    t -> amplitude * max(sin(2π * (t - 13 / 48)) + 31 / 50, 0.0) / (1 + 31 / 50)

"""
    light_filter(; amplitude=0.8, ncells=500, delta=20e-3, temperature=15.0) -> SandFilter

Filter for the study. `delta` (`sand_roughness`) is held at 20 mm so the
near-optimal light band sits ~8 mm above the nominal sand surface, in water that
is ~36% solid, where bare-sand attachment is active.
"""
function light_filter(; amplitude::Real=0.8, ncells::Int=500,
                        delta::Real=20e-3, temperature::Real=15.0)
    f = SandFilter(temperature=temperature, sand_roughness=delta,
                   light_irradiation=light_forcing(amplitude))
    return addgridpoints(f, ncells)
end
