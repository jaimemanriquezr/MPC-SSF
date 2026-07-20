# Port of the pathogen model stored in slow-sand-filtration/thesis_model.mat
# (the model exercised by @SDfilter/run_pathogen.m and example.m / srun_pathogen.m).
#
# It extends the Lund/Rosenqvist ecology with a pathogen indicator (PAT) and two
# extra reactions — pathogen inactivation (natural die-off) and bacterivory
# (predation of pathogens, promoted by heterotrophs). The values below are the
# verbatim `thesis_model.mat` contents (dumped via MATLAB); note this model uses
# its own unit convention (densities ~1, τ = 1e-3), distinct from `modelLund`.
#
# Two mechanisms specific to the pathogen model map onto general fields:
#   * `water_factor` (MATLAB `r_water = [0…0 1]·water_factor`) → per-reaction
#     `efficiency_flowing`: every reaction is inert in the flowing phase except
#     bacterivory, which acts there scaled by `water_factor`.
#   * `sand_pathogen` → PAT's `sand_attachment_factor`: pathogens attach to bare
#     sand at a fraction of their biofilm-attachment rate.
#
# The bacterivory Monod term is HET-based: run_pathogen computes
# `Xb(:,1)/(Xb(:,1)+kPred)` (HET is column 1), even though MATLAB parks `kPred`
# in the DOM row of the "Bacterivory" kinetic-parameter column. We therefore key
# its half-saturation on "HET", reproducing the solver's actual kinetics.

"""
    modelPathogen(; water_factor=1e-3, sand_pathogen=0.0, dark_respiration=1e-3) -> Model

The pathogen preset (port of `thesis_model.mat`): heterotrophs (HET), phototrophs
(PHO), particulate organic matter (POM) and a pathogen indicator (PAT) as
particulates; O2, IC, NH4, HPO4, DOM as liquids; seven reactions — the five Lund
reactions plus **Inactivation** (PAT die-off) and **Bacterivory** (HET-promoted
predation removing PAT).

`water_factor` scales bacterivory in the flowing phase (all other reactions are
inert there); `sand_pathogen` scales PAT's attachment to bare sand relative to
biofilm; `dark_respiration` is the phototroph baseline light factor. Defaults are
the base-case values from `example.m`.
"""
function modelPathogen(; water_factor::Real=1e-3,
                       sand_pathogen::Real=0.0,
                       dark_respiration::Real=1e-3)
    density_particle = 1.2
    dispersivity = 1.0e-3
    attenuation_particle = 0.094

    HET = Particle(name="HET", density=density_particle, dispersivity=dispersivity,
                   transport_rate=0.1, attachment_sand=1.0, attachment_matrix=1.0,
                   attenuation=attenuation_particle)
    PHO = Particle(name="PHO", density=density_particle, dispersivity=dispersivity,
                   transport_rate=0.1, attachment_sand=1.0, attachment_matrix=1.0,
                   attenuation=attenuation_particle)
    POM = Particle(name="POM", density=density_particle, dispersivity=dispersivity,
                   transport_rate=0.01, attachment_sand=0.0, attachment_matrix=0.0,
                   attenuation=attenuation_particle)
    PAT = Particle(name="PAT", density=density_particle, dispersivity=dispersivity,
                   transport_rate=1.0e-6, attachment_sand=2.0, attachment_matrix=2.0,
                   attenuation=attenuation_particle, sand_attachment_factor=sand_pathogen)

    density_liquid = 1.0
    O2   = Liquid(name="O2",   density=density_liquid, dispersivity=dispersivity, transport_rate=1.0)
    IC   = Liquid(name="IC",   density=density_liquid, dispersivity=dispersivity, transport_rate=1.0)
    NH4  = Liquid(name="NH4",  density=density_liquid, dispersivity=dispersivity, transport_rate=1.0)
    HPO4 = Liquid(name="HPO4", density=density_liquid, dispersivity=dispersivity, transport_rate=1.0)
    DOM  = Liquid(name="DOM",  density=density_liquid, dispersivity=dispersivity, transport_rate=0.1)

    heterotroph_growth = Reaction(name="Heterotroph growth", is_light_dependent=false,
        nominal_rate=0.02, temperature_correction_factor=1.047,
        order=Dict("HET" => 1.0),
        half_saturation_constants=Dict("O2" => 0.3, "NH4" => 0.4,
                                       "HPO4" => 2.0e-6, "DOM" => 0.02),
        stoichiometric_coefficients=Dict("HET" => 1.0, "O2" => -1.23168093,
                                         "IC" => 0.3847619048, "NH4" => -0.02476190476,
                                         "HPO4" => -0.01412698413, "DOM" => -1.587301587))

    phototroph_growth = Reaction(name="Phototroph growth", is_light_dependent=true,
        minimum_light_factor=dark_respiration, optimal_light_factor=1.08,
        nominal_rate=10.0, temperature_correction_factor=1.047,
        order=Dict("PHO" => 1.0),
        half_saturation_constants=Dict("IC" => 0.002, "NH4" => 1.2, "HPO4" => 0.02),
        stoichiometric_coefficients=Dict("PHO" => 1.0, "O2" => 0.9300460829,
                                         "IC" => -0.36, "NH4" => -0.06, "HPO4" => -0.01))

    heterotroph_death = Reaction(name="Heterotroph death", is_light_dependent=false,
        nominal_rate=2.0, temperature_correction_factor=1.066,
        order=Dict("HET" => 1.0),
        stoichiometric_coefficients=Dict("HET" => -1.0, "POM" => 0.9122807018,
                                         "O2" => 0.02336243835, "NH4" => 0.06526315789,
                                         "HPO4" => 0.02087719298))

    phototroph_death = Reaction(name="Phototroph death", is_light_dependent=false,
        nominal_rate=0.4, temperature_correction_factor=1.08,
        order=Dict("PHO" => 1.0),
        stoichiometric_coefficients=Dict("PHO" => -1.0, "POM" => 0.6315789474,
                                         "O2" => 0.2005093379, "NH4" => 0.02210526316,
                                         "HPO4" => 0.003684210526))

    hydrolysis = Reaction(name="Hydrolysis", is_light_dependent=false,
        nominal_rate=1.0, temperature_correction_factor=1.08,
        order=Dict("HET" => 1.0),
        half_saturation_constants=Dict("POM/HET" => 0.002),
        stoichiometric_coefficients=Dict("POM" => -1.0, "DOM" => 1.0))

    # Pathogen die-off: first-order in PAT, no Monod, inert in the flowing phase.
    inactivation = Reaction(name="Inactivation", is_light_dependent=false,
        nominal_rate=0.4, temperature_correction_factor=1.08,
        efficiency_flowing=0.0,
        order=Dict("PAT" => 1.0),
        stoichiometric_coefficients=Dict("PAT" => -1.0))

    # Bacterivory: first-order in PAT, HET-limited (Monod, K = kPred = 0.002),
    # removes PAT and returns biomass to HET. This is the only reaction active in
    # the flowing phase, scaled there by `water_factor`.
    bacterivory = Reaction(name="Bacterivory", is_light_dependent=false,
        nominal_rate=20.0, temperature_correction_factor=1.08,
        efficiency_flowing=float(water_factor),
        order=Dict("PAT" => 1.0),
        half_saturation_constants=Dict("HET" => 0.002),
        stoichiometric_coefficients=Dict("PAT" => -1.0, "HET" => 1.0))

    # All non-bacterivory reactions are inert in the flowing phase (r_water = 0);
    # rebuild them with efficiency_flowing = 0 (Reaction is immutable).
    mk_inert(rx) = Reaction(name=rx.name, nominal_rate=rx.nominal_rate,
        temperature_correction_factor=rx.temperature_correction_factor,
        order=rx.order, half_saturation_constants=rx.half_saturation_constants,
        stoichiometric_coefficients=rx.stoichiometric_coefficients,
        efficiency_biofilm=rx.efficiency_biofilm, efficiency_flowing=0.0,
        is_light_dependent=rx.is_light_dependent,
        minimum_light_factor=rx.minimum_light_factor,
        optimal_light_factor=rx.optimal_light_factor)

    components = Component[HET, PHO, POM, PAT, O2, IC, NH4, HPO4, DOM]
    reactions  = Reaction[mk_inert(heterotroph_growth), mk_inert(phototroph_growth),
                          mk_inert(heterotroph_death), mk_inert(phototroph_death),
                          mk_inert(hydrolysis), inactivation, bacterivory]

    cohesion = CahnHilliardModel(kappa=1.0e-6, zeta_0=1.0, zeta_1=1/100)

    return Model(components, reactions;
                 cohesion_submodel=cohesion,
                 water_density=density_liquid,
                 biofilm_porosity=0.99,
                 osmosis_rate=1.0e-3,
                 detachment=(v -> sqrt.(abs.(v) ./ 7.2)))
end
