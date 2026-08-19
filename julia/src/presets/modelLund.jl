# Port of src/presets/modelLund.m
#
# The "Lund" / "Rosenqvist" preset: a 9-component (4 particulate, 5 dissolved),
# 5-reaction ecological model with a Cahn-Hilliard cohesion submodel. Values are
# copied verbatim from modelLund.m so the two codebases share an identical model
# for golden-master validation.
#
# MATLAB builds the CohesionSubModel from the Kappa/Zeta0/Zeta1 arguments passed
# to Model(...); here we construct the CahnHilliardModel explicitly (its
# mobility/potential-gradient defaults already mirror CahnHilliardModel.m). The
# MATLAB DetachmentFunction is the one-arg `@(v) sqrt(v / 7.2)`.

"""
    modelLund(; phototroph_respiration=0.0) -> Model

The Lund/Rosenqvist preset model (port of `presets/modelLund.m`): heterotrophs
(HET), phototrophs (PHO), particulate organic matter (POM) and a pathogen
indicator (PAT) as particulates; O2, IC, NH4, HPO4, DOM as liquids; five
reactions (heterotroph/phototroph growth and death, hydrolysis). Used as the
shared reference model for MATLAB↔Julia validation.

`phototroph_respiration` > 0 enables the metabolic split (plan
2026-08-18-phototroph-respiration): a sixth reaction "Phototroph respiration" —
first-order in PHO, light-INdependent (maintenance respiration runs day and
night), O2-Monod-limited, stoichiometry the exact reverse of phototroph growth —
is appended with that nominal rate, and the growth reaction's dark floor
(minimum_light_factor) drops to 0, since photosynthesis is genuinely zero in
darkness once respiration carries the dark O2 cost. Literature value:
kra = 0.0020–0.0210 /h, avg 0.0115 /h = 0.276 /d (Campos2006 Table 3, Brown &
Barnwell 1987), θ_kra = 1.08. The default 0.0 reproduces the original 5-reaction
model bit-for-bit (golden suites unchanged).
"""
function modelLund(; phototroph_respiration::Real=0.0, pg_fraction::Real=0.2, pg_yield::Real=0.63, pg_excess::Bool=false)
    density_particle = 1.117e3
    attenuation_particle = 0.094
    dispersivity_particle = 1.20e-2

    HET = Particle(name="HET", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                   attenuation=attenuation_particle)
    PHO = Particle(name="PHO", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                   attenuation=attenuation_particle)
    POM = Particle(name="POM", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=0.0, attachment_matrix=0.0,
                   attenuation=attenuation_particle)
    PAT = Particle(name="PAT", density=density_particle, dispersivity=dispersivity_particle,
                   transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                   attenuation=attenuation_particle)

    density_liquid = 0.998e3
    dispersivity_liquid = 1.20e-2

    O2   = Liquid(name="O2",   density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e2)
    IC   = Liquid(name="IC",   density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e2)
    NH4  = Liquid(name="NH4",  density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e2)
    HPO4 = Liquid(name="HPO4", density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=6.00e2)
    DOM  = Liquid(name="DOM",  density=density_liquid, dispersivity=dispersivity_liquid, transport_rate=3.00e2)

    heterotroph_growth = Reaction(name="Heterotroph growth", is_light_dependent=false,
        nominal_rate=1.81e-2, temperature_correction_factor=1.047,
        order=Dict("HET" => 1.0),
        half_saturation_constants=Dict("O2" => 3.00e-3, "DOM" => 2.00e-4,
                                       "NH4" => 4.00e-3, "HPO4" => 1.40e-8),
        stoichiometric_coefficients=Dict("HET" => 1.0, "O2" => -1.2317, "IC" => 0.3848,
                                         "NH4" => -0.0248, "HPO4" => -0.0141, "DOM" => -1.5873))

    phototroph_growth = Reaction(name="Phototroph growth", is_light_dependent=true,
        minimum_light_factor=0.01, optimal_light_factor=1.814e-2,
        nominal_rate=5.50, temperature_correction_factor=1.047,
        order=Dict("PHO" => 1.0),
        half_saturation_constants=Dict("IC" => 2.00e-5, "NH4" => 1.20e-2, "HPO4" => 1.68e-4),
        stoichiometric_coefficients=Dict("PHO" => 1.0, "O2" => 0.9301, "IC" => -0.3600,
                                         "NH4" => -0.0600, "HPO4" => -0.0100))

    heterotroph_death = Reaction(name="Heterotroph death", is_light_dependent=false,
        nominal_rate=2.00, temperature_correction_factor=1.066,
        order=Dict("HET" => 1.0),
        stoichiometric_coefficients=Dict("HET" => -1.0, "POM" => 0.9123, "O2" => 0.0234,
                                         "NH4" => 0.0653, "HPO4" => 0.0209))

    phototroph_death = Reaction(name="Phototroph death", is_light_dependent=false,
        nominal_rate=4.00e-1, temperature_correction_factor=1.080,
        order=Dict("PHO" => 1.0),
        stoichiometric_coefficients=Dict("PHO" => -1.0, "POM" => 0.6316, "O2" => 0.2005,
                                         "NH4" => 0.0221, "HPO4" => 0.0037))

    hydrolysis = Reaction(name="Hydrolysis", is_light_dependent=false,
        nominal_rate=9.00e-2, temperature_correction_factor=1.080,
        order=Dict("HET" => 1.0),
        half_saturation_constants=Dict("POM/HET" => 2.00e-5),
        stoichiometric_coefficients=Dict("POM" => -1.0, "DOM" => 1.0))

    components = Component[HET, PHO, POM, PAT, O2, IC, NH4, HPO4, DOM]
    if phototroph_respiration > 0 && !pg_excess
        # Internally stored polyglucose (Wolf2007 PHOBIA): CH2O, COD 32/30,
        # carbon 0.4 kg C/kg, no N/P. Transport-identical to PHO (intracellular).
        # Appended as the 5th particle so Lund-structure indices survive.
        PG = Particle(name="PG", density=density_particle, dispersivity=dispersivity_particle,
                      transport_rate=5.47, attachment_sand=5.47e2, attachment_matrix=5.47e2,
                      attenuation=attenuation_particle)
        components = Component[HET, PHO, POM, PAT, PG, O2, IC, NH4, HPO4, DOM]
    end
    if phototroph_respiration > 0 && pg_excess
        # PG-in-excess variant (Jaime, 2026-08-19): the storage pool is assumed
        # never limiting and is NOT tracked. Growth keeps the original Lund row
        # (floor retired); respiration is Wolf2007 r6 WITHOUT the PG column,
        # normalized per unit (untracked) PG consumed, and its light factor is
        # the exact complement 1 − Steele(I) of the growth factor: biomass is
        # built in dark places, consuming NH4 and O2, releasing IC. Deliberately
        # mass-non-conservative toward the untracked pool.
        Y = float(pg_yield)
        phototroph_growth = Reaction(name="Phototroph growth", is_light_dependent=true,
            minimum_light_factor=0.0, optimal_light_factor=1.814e-2,
            nominal_rate=5.50, temperature_correction_factor=1.047,
            order=Dict("PHO" => 1.0),
            half_saturation_constants=Dict("IC" => 2.00e-5, "NH4" => 1.20e-2, "HPO4" => 1.68e-4),
            stoichiometric_coefficients=Dict("PHO" => 1.0, "O2" => 0.9301, "IC" => -0.3600,
                                             "NH4" => -0.0600, "HPO4" => -0.0100))
        phototroph_respiration_rx = Reaction(name="Phototroph respiration",
            nominal_rate=float(phototroph_respiration), temperature_correction_factor=1.08,
            is_light_complement=true,
            order=Dict("PHO" => 1.0),
            half_saturation_constants=Dict("O2" => 3.00e-3),
            stoichiometric_coefficients=Dict("PHO" => Y,
                                             "O2" => -(1.0667 - 0.9301Y),
                                             "IC" => 0.4 - 0.36Y,
                                             "NH4" => -0.0600Y, "HPO4" => -0.0100Y))
        reactions = Reaction[heterotroph_growth, phototroph_growth,
                             heterotroph_death, phototroph_death, hydrolysis,
                             phototroph_respiration_rx]
    elseif phototroph_respiration > 0
        # Photosynthesis with internal storage (Wolf2007): per unit PHO built,
        # the f-fraction goes additionally into PG (CH2O: +1.0667 O2, -0.4 C per
        # unit PG), coupled to the growth rate. Dark floor retired (respiration
        # carries the dark metabolism).
        f = float(pg_fraction)
        phototroph_growth = Reaction(name="Phototroph growth", is_light_dependent=true,
            minimum_light_factor=0.0, optimal_light_factor=1.814e-2,
            nominal_rate=5.50, temperature_correction_factor=1.047,
            order=Dict("PHO" => 1.0),
            half_saturation_constants=Dict("IC" => 2.00e-5, "NH4" => 1.20e-2, "HPO4" => 1.68e-4),
            stoichiometric_coefficients=Dict("PHO" => 1.0, "PG" => f,
                                             "O2" => 0.9301 + 1.0667f, "IC" => -(0.3600 + 0.4f),
                                             "NH4" => -0.0600, "HPO4" => -0.0100))
        # Wolf2007 (PHOBIA) r6, dark respiration = growth on stored polyglucose.
        # Per unit PG consumed with yield Y = pg_yield (Wolf leaves Y_PH/PG
        # unpinned; 0.63 = ASM heterotroph yield, ASSUMED): biomass produced,
        # NH4 consumed, O2 consumed, IC released. Columns are the difference of
        # the PG-synthesis and biomass rows, so COD/elements balance exactly.
        # Rate: 0.1*q_max (Tillmann & Rick 2001) = 0.55/d for mu_PHO = 5.5,
        # first-order in PHO, min-Monod over O2 and the PG/PHO quotient
        # (K_S,PH,PG = 0.005, Wolf Table VI), dark-only via K_inh/(K_inh + I)
        # with K_inh = 8e-5 normalized by I_opt = 1.814e-2.
        Y = float(pg_yield)
        phototroph_respiration_rx = Reaction(name="Phototroph respiration",
            nominal_rate=float(phototroph_respiration), temperature_correction_factor=1.08,
            light_inhibition=8e-5/1.814e-2,
            order=Dict("PHO" => 1.0),
            half_saturation_constants=Dict("O2" => 3.00e-3, "PG/PHO" => 0.005),
            stoichiometric_coefficients=Dict("PG" => -1.0, "PHO" => Y,
                                             "O2" => -(1.0667 - 0.9301Y),
                                             "IC" => 0.4 - 0.36Y,
                                             "NH4" => -0.0600Y, "HPO4" => -0.0100Y))
        reactions = Reaction[heterotroph_growth, phototroph_growth,
                             heterotroph_death, phototroph_death, hydrolysis,
                             phototroph_respiration_rx]
    else
        reactions = Reaction[heterotroph_growth, phototroph_growth,
                             heterotroph_death, phototroph_death, hydrolysis]
    end

    cohesion = CahnHilliardModel(kappa=1.00e-6, zeta_0=1.00e6, zeta_1=1/100)

    return Model(components, reactions;
                 cohesion_submodel=cohesion,
                 biofilm_porosity=0.99,
                 osmosis_rate=1.00e-7,
                 detachment=(v -> sqrt.(v ./ 7.2)))
end
