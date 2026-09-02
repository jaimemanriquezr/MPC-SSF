function model = modelLund(options)
    arguments
            options.WriteFile = false;
            % > 0 enables the phototroph metabolic split (plan
            % 2026-08-18-phototroph-respiration): appends a sixth reaction
            % "Phototroph respiration" (first-order in PHO, light-independent,
            % O2-Monod-limited, stoichiometry = reverse of phototroph growth)
            % at this nominal rate, and drops the growth reaction's dark floor
            % (MinimumLightFactor) to 0. Literature: kra = 0.0020-0.0210 /h,
            % avg 0.0115 /h = 0.276 /d (Campos2006 Table 3, Brown & Barnwell
            % 1987), theta_kra = 1.08. Default 0.0 reproduces the original
            % 5-reaction model exactly (goldens unchanged).
            options.PhototrophRespiration (1,1) {mustBeNumeric} = 0.0;
            % Wolf2007 storage fraction f (kg COD PG per kg COD PHO built) and
            % yield Y_PH/PG (ASSUMED 0.63, the ASM heterotroph yield; Wolf2007
            % Table VI leaves it unpinned).
            options.PGFraction (1,1) {mustBeNumeric} = 0.2;
            options.PGYield (1,1) {mustBeNumeric} = 0.63;
            % true = LightIrradiation(t) is ALREADY the normalized intensity
            % I_hat = I/I_opt, so the growth reaction's OptimalLightFactor
            % becomes 1.0. With the inherited curves (peak 0.8) the surface is
            % then sub-optimal all day instead of 44x optimal (units mismatch:
            % 1.814e-2 is Wolf2007's optimum in absolute PHOBIA units).
            options.NormalizedLight (1,1) logical = false;
            % Respiration light response is a Monod term K/(K + I_hat); K in
            % the same units as the (normalized) intensity. K = 1 halves
            % respiration at optimal light; Wolf's sharp dark switch ~ 4.4e-3.
            options.RespirationLightK (1,1) {mustBeNumeric} = 1.0;
            % true = PG-in-excess variant: the storage pool is assumed never
            % limiting and is NOT tracked (no PG component); respiration's
            % light factor is the complement 1 - Steele(I) of the growth
            % factor. Deliberately mass-non-conservative toward the pool.
            options.PGExcess (1,1) logical = false;
            % Which respiration process runs when PhototrophRespiration > 0:
            %   "reichert" - RWQM1 process (10), aerobic endogenous respiration of
            %                algae (Reichert et al. 2001): exact reverse of the
            %                growth row, first order in PHO, Monod on O2
            %                (K_O2,ALG = 0.2 g/m3), light-independent, no biomass
            %                produced, no storage pool. RWQM1 rate 0.1/d at 20 C,
            %                beta_ALG = 0.046/C (theta 1.047).
            %   "pgexcess" - the PG-in-excess variant above (legacy).
            %   "wolf"     - Wolf2007 r6 with a tracked PG pool (legacy).
            % "" resolves to the pre-2026-08-25 behaviour: PGExcess ? pgexcess : wolf.
            options.RespirationForm (1,1) string {mustBeMember(options.RespirationForm, ["", "reichert", "pgexcess", "wolf"])} = ""
    end
    respirationForm = options.RespirationForm;
    if respirationForm == ""
        if options.PGExcess, respirationForm = "pgexcess"; else, respirationForm = "wolf"; end
    end

    densityParticle = 1.117E+03;
    % Biofilm self-shading, m^2/kg. 52 = Tenore2021 Table 1 k_tot (210) converted from
    % DRY to WET mass, x f_dry = 0.25. Was 0.094 (Gallegos2000), which is light
    % attenuation by organic matter suspended in natural water, not by biofilm.
    %
    % Why the conversion is needed: Tenore's k_tot multiplies areal DRY biomass; this
    % coefficient multiplies the wet-mass concentration this model carries (rho_P = 1117
    % is a hydrated-cell density). Measured on chain_fld2x_lit_leg6 (60 d), the filter
    % holds 0.221 kg/m^2 wet = 0.055 kg/m^2 dry, against Tenore's ~0.05 kg/m^2 -- the
    % two systems carry the SAME areal loading, so his coefficient transfers once the
    % wet/dry basis is matched. At 52 the column reaches eta_P ~ 11.5, comparable to
    % Tenore's ~10.5. At the old 0.094 it was eta_P = 0.021, i.e. 2% attenuation over
    % the whole bed -- self-shading was effectively switched off.
    %
    % *** COUPLED TO f_dry (decision 2026-09-01, f_dry = 0.25, implement after Friday).
    % *** If f_dry is applied to rho_P and the influents, the state variable becomes DRY
    % *** mass and this MUST revert to Tenore's 210. Changing one without the other
    % *** double-counts the water. The two values are the same physics.
    attenuationParticle = 52;
    dispersivityParticle = 1.20E-02;
    HET =  Particle(Name="HET", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             Transport=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=5.47E+02,...
                             Attenuation=attenuationParticle);
    PHO =  Particle(Name="PHO", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             Transport=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=5.47E+02,...
                             Attenuation=attenuationParticle);
    POM =  Particle(Name="POM", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             Transport=5.47, AttachmentSand=0.0, AttachmentMatrix=0.0,...
                             Attenuation=attenuationParticle);
    PAT =  Particle(Name="PAT", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             Transport=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=5.47E+02,...
                             Attenuation=attenuationParticle);
                    
    densityLiquid = 0.998E+03;
    dispersivityLiquid = 1.20E-02;
    O2 =  Liquid(Name="O2", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            Transport=6.00E+02);
    IC =  Liquid(Name="IC", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            Transport=6.00E+02);
    NH4 =  Liquid(Name="NH4", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            Transport=6.00E+02);
    HPO4 =  Liquid(Name="HPO4", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            Transport=6.00E+02);
    DOM =  Liquid(Name="DOM", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            Transport=3.00E+02);


    heterotrophGrowth =  Reaction(Name="Heterotroph growth", IsLightDependent=false, ...
            NominalRate=4.80, ...   % Tenore2021 Table 1 mu_max,2 (f2 = heterotrophs)
            TemperatureCorrectionFactor=1.0725, ...
            Order=dictionary("HET", 1), ...
            HalfSaturationConstants=dictionary("O2", 2.00E-04, "DOM", 4.00E-03, "NH4", 1.00E-06, "HPO4", 2.00E-05), ...
            StoichiometricCoefficients=dictionary("HET", 1.0, ...
                                        "O2", -1.2317, "IC", 0.3848, "NH4", -0.0248, "HPO4", -0.0141, "DOM", -1.5873));

    phototrophGrowth =  Reaction(Name="Phototroph growth", IsLightDependent=true, ...
            MinimumLightFactor=0.0, ...
            OptimalLightFactor=ternaryOpt(options.NormalizedLight), ...
            NominalRate=2.00, ...
            TemperatureCorrectionFactor=1.047, ...
            Order=dictionary("PHO", 1), ...
            HalfSaturationConstants=dictionary("IC", 1.20E-03, "NH4", 2.00E-05, "HPO4", 2.00E-05), ...
            StoichiometricCoefficients=dictionary("PHO", 1.0, ...
                                        "O2", 0.9301, "IC", -0.3600, "NH4", -0.0600, "HPO4", -0.0100));

    heterotrophDeath =  Reaction(Name="Heterotroph death", IsLightDependent=false, ...
            NominalRate=0.40, ...
            TemperatureCorrectionFactor=1.0725, ...
            Order=dictionary("HET", 1), ...
            StoichiometricCoefficients=dictionary("HET", -1.0, "POM", 0.9123, ...
                                        "O2", 0.0234, "NH4", 0.0653, "HPO4", 0.0209));

    % Phototroph loss = Campos2006 k_ra, "algae loss due to the combined effects of
    % respiration and excretion" (Brown & Barnwell 1987): 0.048-0.504 /d, avg
    % 0.276 /d, theta 1.08 (Campos2006 Table 3). One loss process, as in Campos
    % and in Wolf2007 (b_ina,PH); no separate respiration reaction is used
    % (Jaime, 2026-08-26). Stoichiometry stays RWQM1 (11): the O2 cost of the dead
    % carbon is paid downstream by heterotroph growth on hydrolysed POM.
    phototrophDeath =  Reaction(Name="Phototroph death", IsLightDependent=false, ...
            NominalRate=0.276, ...
            TemperatureCorrectionFactor=1.080, ...
            Order=dictionary("PHO", 1), ...
            StoichiometricCoefficients=dictionary("PHO", -1.0, "POM", 0.6316, ...
                                        "O2", 0.2005, "NH4", 0.0221, "HPO4", 0.0037));

    hydrolysis =  Reaction(Name="Hydrolysis", IsLightDependent=false, ...
            NominalRate=3.00, ...
            TemperatureCorrectionFactor=1.0725, ...
            Order=dictionary("HET", 1), ...
            HalfSaturationConstants=dictionary("POM/HET", 0.1), ...
            StoichiometricCoefficients=dictionary("POM", -1.0, "DOM", 1.0));
                                    

    componentList = [HET; PHO; POM; PAT; 
                      O2; IC; NH4; HPO4; DOM];
    reactionList = [heterotrophGrowth; 
                    phototrophGrowth; 
                    heterotrophDeath; 
                    phototrophDeath; 
                    hydrolysis];
    if options.PhototrophRespiration > 0 && respirationForm == "pgexcess"
        % PG-in-excess variant (Jaime, 2026-08-19): growth keeps the original
        % Lund row (floor retired); respiration is Wolf2007 r6 WITHOUT the PG
        % column, light factor = 1 - Steele(I): biomass built in dark places,
        % consuming NH4 and O2, releasing IC.
        Y = options.PGYield;
        phototrophGrowth.MinimumLightFactor = 0.0;
        reactionList(2) = phototrophGrowth;
        % NH4/HPO4 half-saturations are DEPLETION PROTECTION (far below any
        % ambient level): without them the consumed pools cross zero once
        % drained (bottom-bed enclosed NH4 at fine grids).
        phototrophRespiration = Reaction(Name="Phototroph respiration", ...
                NominalRate=options.PhototrophRespiration, ...
                LightInhibition=options.RespirationLightK, ...
                TemperatureCorrectionFactor=1.08, ...
                Order=dictionary("PHO", 1), ...
                HalfSaturationConstants=dictionary(["O2", "NH4"], [3.00E-03, 1.0E-06]), ...
                StoichiometricCoefficients=dictionary( ...
                        ["PHO", "O2", "IC", "NH4"], ...
                        [Y, -(1.0667 - 0.9301*Y), 0.4 - 0.36*Y, -0.06*Y]));
        reactionList = [reactionList; phototrophRespiration];
    elseif options.PhototrophRespiration > 0 && respirationForm == "wolf"
        % Wolf2007 (PHOBIA): photosynthesis stores the f-fraction into an
        % internal polyglucose pool PG (CH2O: +1.0667 O2, -0.4 C per unit PG);
        % dark respiration (r6) grows PHO on PG, consuming NH4 and O2. PG is
        % transport-identical to PHO (intracellular) and appended as the 5th
        % particle so Lund-structure indices survive. Dark floor retired.
        f = options.PGFraction;
        Y = options.PGYield;
        PG = Particle(Name="PG", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                      Transport=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=5.47E+02, ...
                      Attenuation=attenuationParticle);
        componentList = [componentList(1:4); PG; componentList(5:end)];
        phototrophGrowth.MinimumLightFactor = 0.0;
        phototrophGrowth.StoichiometricCoefficients = dictionary( ...
                ["PHO", "PG", "O2", "IC", "NH4", "HPO4"], ...
                [1.0, f, 0.9301 + 1.0667*f, -(0.3600 + 0.4*f), -0.0600, -0.0100]);
        reactionList(2) = phototrophGrowth;
        % r6: rate 0.1*q_max = 0.55/d (Tillmann & Rick 2001), first-order in
        % PHO, min-Monod over O2 and the PG/PHO quotient (K_S,PH,PG = 0.005,
        % Wolf Table VI), dark-only via K_inh/(K_inh + I), K_inh = 8e-5
        % normalized by I_opt = 1.814e-2. Columns are the difference of the
        % PG-synthesis and biomass rows, so COD/elements balance exactly.
        phototrophRespiration = Reaction(Name="Phototroph respiration", ...
                NominalRate=options.PhototrophRespiration, ...
                LightInhibition=options.RespirationLightK, ...
                TemperatureCorrectionFactor=1.08, ...
                Order=dictionary("PHO", 1), ...
                HalfSaturationConstants=dictionary(["O2", "PG/PHO", "NH4"], ...
                        [3.00E-03, 0.005, 1.0E-06]), ...
                StoichiometricCoefficients=dictionary( ...
                        ["PG", "PHO", "O2", "IC", "NH4"], ...
                        [-1.0, Y, -(1.0667 - 0.9301*Y), 0.4 - 0.36*Y, -0.06*Y]));
        reactionList = [reactionList; phototrophRespiration];
    elseif options.PhototrophRespiration > 0
        % RWQM1 process (10): the growth row (9a) run backwards. With the algal
        % composition alpha_C 0.36, alpha_H 0.07, alpha_O 0.50, alpha_N 0.06,
        % alpha_P 0.01 the O2 coefficient 8a_C/3 + 8a_H - a_O - 12a_N/7 + 40a_P/31
        % is 0.9300, i.e. the same 0.9301 the growth row releases. Dark floor retired.
        phototrophGrowth.MinimumLightFactor = 0.0;
        reactionList(2) = phototrophGrowth;
        phototrophRespiration = Reaction(Name="Phototroph respiration", ...
                NominalRate=options.PhototrophRespiration, ...
                TemperatureCorrectionFactor=1.047, ...
                Order=dictionary("PHO", 1), ...
                HalfSaturationConstants=dictionary("O2", 2.00E-04), ...
                StoichiometricCoefficients=dictionary( ...
                        ["PHO", "O2", "IC", "NH4", "HPO4"], ...
                        [-1.0, -0.9301, 0.3600, 0.0600, 0.0100]));
        reactionList = [reactionList; phototrophRespiration];
    end
    model = Model(Components=componentList, ...
                 Kappa=1.00E-06, Zeta0=1.00E+06, Zeta1=1/100, ...
                 DetachmentFunction=@(v) sqrt(v / 7.2), ...
                 BiofilmPorosity=0.99, ...
                 OsmosisRate=1.00E-07, ...
                 Reactions=reactionList);

    if options.WriteFile
        save("./LundMPCModel.mat", "model");
    end
end

% MinimumLightFactor (dark-growth floor) is 0.0 since 2026-08-26. The published
% 0.01 let a covered filter photosynthesise at 1 % of optimum; it was retired
% inside the respiration branches only, so every driver that did not enable
% respiration kept it (season90 and the old cover data were contaminated).
% Growth, loss and hydrolysis rates and temperature factors (2026-08-25, second
% audit, .claude/decisions/2026-08-25-kinetics-wolf-reichert.md): rates from
% Reichert2001 Table 8 where Wolf2007 has none (mu_ALG 2.0, k_death,ALG 0.1,
% k_hyd 3.0) and Wolf2007 Table VI otherwise (b_ina,H 0.4, k_h 3.0); every theta
% is exp(beta) from Reichert (beta_ALG 0.046 -> 1.047, beta_H = beta_hyd 0.07 ->
% 1.0725), since PHOBIA has no temperature dependence. Old set for provenance:
%   mu_HET 1.81e-2 (th 1.047)  mu_PHO 5.5 (1.047)  d_HET 2.0 (1.066)
%   d_PHO 0.4 (1.08)  k_hyd 0.09 (1.08).  d_PHO was 0.10 (1.047, Reichert death
%   alone) between the second audit and the Campos merge of 2026-08-26.
% Half-saturation constants (2026-08-25): source-faithful values from Wolf2007
% Table VI and Reichert2001, replacing the published Table 2 set, which the
% 2026-08-20 audit (.claude/decisions/2026-08-20-kinetic-parameter-audit.md)
% showed to be adjacent-row slides of the cited tables. The old K_NH4,PHO =
% 1.2e-2 gave a Monod factor of 0.0017 at the influent, so phototrophs could
% not grow at any light level. Old set kept here for provenance:
%   HET: O2 3.0e-3, DOM 2.0e-4, NH4 4.0e-3, HPO4 1.4e-8
%   PHO: IC 2.0e-5, NH4 1.2e-2, HPO4 1.68e-4;  hydrolysis POM/HET 2.0e-5.
function opt = ternaryOpt(normalized)
if normalized, opt = 1.0; else, opt = 1.814E-02; end
end
