# Julia port of data/SimpleModel.mat (as used by examples/simulation_example.m).
#
# SimpleModel is a minimal 3-component / 3-reaction model used as the shared
# golden-master reference for MATLAB↔Julia validation because, unlike the full
# Lund preset, it runs to completion (flag "OK") on a clean filter with a modest
# inflow. The values below reproduce the stored .mat *with the two runtime
# overrides the example applies*: cohesion Kappa=1e-2 / Zeta0=1.0 and the
# one-arg detachment @(v) sqrt(abs(v)). The golden harness applies the identical
# overrides on the MATLAB side, so both codebases represent the same model.
#
# Components: Microorganism, POM (particulate), Nutrient (dissolved) — already in
# particles-then-liquids order (required by the reaction kernel).

"""
    simpleModel() -> Model

The minimal `SimpleModel` used for golden-master validation (port of
`data/SimpleModel.mat` with the example's cohesion/detachment overrides baked
in): microorganisms and particulate organic matter as particulates, a dissolved
nutrient, and three reactions (death, growth, hydrolysis).
"""
function simpleModel()
    Microorganism = Particle(name="Microorganism", density=1.1, dispersivity=0.0,
                             transport_rate=0.0, attachment_sand=1.0,
                             attachment_matrix=10.0, attenuation=0.0)
    POM = Particle(name="POM", density=1.1, dispersivity=0.0,
                   transport_rate=0.0, attachment_sand=1.0,
                   attachment_matrix=10.0, attenuation=0.0)
    Nutrient = Liquid(name="Nutrient", density=1.0, dispersivity=0.0, transport_rate=0.0)

    death = Reaction(name="Death", nominal_rate=1.0, temperature_correction_factor=1.0,
        order=Dict("Microorganism" => 1.0),
        stoichiometric_coefficients=Dict("Microorganism" => -1.0, "POM" => 1.0))

    growth = Reaction(name="Growth", nominal_rate=1.0, temperature_correction_factor=1.0,
        order=Dict("Microorganism" => 1.0),
        half_saturation_constants=Dict("Nutrient" => 1.0e-3),
        stoichiometric_coefficients=Dict("Microorganism" => 1.0, "Nutrient" => -1.0))

    hydrolysis = Reaction(name="Hydrolysis", nominal_rate=1.0, temperature_correction_factor=1.0,
        order=Dict("Microorganism" => 1.0),
        half_saturation_constants=Dict("Microorganism/POM" => 2.0e-5),
        stoichiometric_coefficients=Dict("POM" => -1.0, "Nutrient" => 1.0))

    components = Component[Microorganism, POM, Nutrient]
    reactions  = Reaction[death, growth, hydrolysis]

    # Example overrides: Kappa=1e-2, Zeta0=1.0 (Zeta1 default 0).
    cohesion = CahnHilliardModel(kappa=1.0e-2, zeta_0=1.0)

    return Model(components, reactions;
                 cohesion_submodel=cohesion,
                 water_density=1.0,
                 biofilm_porosity=0.99,
                 osmosis_rate=1.0e-5,
                 detachment=(v -> sqrt.(abs.(v))))
end
