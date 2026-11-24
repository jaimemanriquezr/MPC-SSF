function model = modelLund(options)
% MODELLUND  The AUDITED MANRIQUEZ2026 WORKING SET -- not a general SSF model.
%
%   model = modelLund()
%   model = modelLund(Attenuation=52, ThetaSet="reichert", ...)
%
% READ THIS BEFORE REUSING ANY NUMBER BELOW. This file is not a neutral
% parameter library. It is the one specific parameterisation that Manriquez2026
% reports (tab:general-parameters, tab:rhs-parameters, tab:eco-parameters),
% reconstructed constant by constant from the audits of 2026-08 / 2026-09. The
% preset IS the working set used for the manuscript runs, so changing a default
% here changes the published model, and every regression golden in tests/ with it.
%
% PROVENANCE OF THE CONSTANTS. Most come from the literature and are cited in
% place beside the value. The ones that came from a RUN rather than a paper, or
% from a correction to the code, are:
%
%   Attenuation = 52 [m^2/kg]   Tenore2021 Table 1 k_tot = 210 converted from a
%       DRY to a WET mass basis with f_dry = 0.25. The conversion was VALIDATED
%       ON A RUN, not on a paper: the areal loading was measured on
%       `chain_fld2x_lit_leg6` (60 d) -- 0.221 kg/m^2 wet = 0.055 kg/m^2 dry
%       against Tenore's ~0.05 kg/m^2 -- which is what licenses transferring his
%       coefficient. *** COUPLED TO f_dry. If f_dry is
%       ever applied to rho_P and the influents, the state becomes DRY mass and
%       this MUST revert to 210. See the full note at the assignment below.
%
%   Zeta0 = 1, Zeta1 = 0.27, DetachmentFunction = 0.14*|v|/18
%       The 2026-09-05 correction. The preset had DRIFTED from the published
%       model on all three (Zeta0 1e6, Zeta1 1/100, detachment sqrt(v/7.2)) and
%       every production driver was silently overriding it back. Zeta0 and Zeta1
%       are the manuscript's tab:general-parameters values (own calibration, and
%       \citet{SSFnum}); the linear detachment law is tab:rhs-parameters' 0.14/d
%       at the nominal velocity 18 m/d. The old sqrt form additionally used the
%       INTAKE velocity 7.2 where 18 belongs and omitted abs(). Full note at the
%       Model() call at the bottom of this file.
%
%   ThetaSet = "reichert" (default since 2026-09-07, decision
%       2026-09-05-theta-set-all-reichert). "mixed" reproduces every run before
%       that date byte-identically. This is a CHOICE between literature sets,
%       not a measurement; see the option below.
%
% Reactions whose rate is a judgement value rather than a citation say so beside
% the value (e.g. the grazing rate in the EXPERIMENTAL block below).
%
% See also PATHOGENMODEL, MODEL, REACTION.
    arguments
            % > 0 enables the phototroph metabolic split (plan
            % 2026-08-18-phototroph-respiration): appends a sixth reaction
            % "Phototroph respiration" (first-order in PHO, light-independent,
            % O2-Monod-limited, stoichiometry = reverse of phototroph growth)
            % at this nominal rate, and drops the growth reaction's dark floor
            % (MinimumLightFactor) to 0. Literature: kra = 0.0020-0.0210 /h,
            % avg 0.0115 /h = 0.276 /d (Campos2006 Table 3, Brown & Barnwell
            % 1987), theta_kra = 1.08. Default 0.0 reproduces the original
            % 5-reaction model exactly (goldens unchanged).
            % Biofilm self-shading coefficient nu_P [m^2/kg]. Default is the preset
            % value; pass 0.094 to reproduce the pre-b1fcd03 (Gallegos2000) runs.
            % Multiplies AttachmentSand for the BIOMASS components HET, PHO and -- in the
            % "wolf" respiration branch only -- PG: attachment to
            % bare sand relative to attachment to biofilm. 1.0 = the current model (sand
            % and biofilm equally likely per unit area). Below 1, biomass accretes
            % preferentially on existing biofilm. POM is unaffected -- its AttachmentSand
            % is 0, so its factor is inert.
            options.SandBiomass (1,1) {mustBeNumeric} = 1.0;
            options.Attenuation (1,1) {mustBeNumeric} = 52;
            options.PhototrophRespiration (1,1) {mustBeNumeric} = 0.0;
            % Wolf2007 storage fraction f (kg COD PG per kg COD PHO built) and
            % yield Y_PH/PG (ASSUMED 0.63, the ASM heterotroph yield; Wolf2007
            % Table VI leaves it unpinned).
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
            % Which literature set the five base-model TemperatureCorrectionFactor
            % values (heterotroph/phototroph growth and death, hydrolysis) are drawn
            % from. "reichert" (DEFAULT since 2026-09-07) is all-Reichert2001
            % Table 8 p.28 (beta form, theta=exp(beta)). "mixed" is
            % the original per-process mix and reproduces every run before 2026-09-07
            % byte-identically; "campos" is all-Campos2006 Table 2 p.880 (theta form,
            % "assumed based on Chapra (1997)").
            options.ThetaSet (1,1) string {mustBeMember(options.ThetaSet, ["mixed", "reichert", "campos"])} = "reichert";
            % ================== EXPERIMENTAL -- NOT PART OF THE WORKING SET ==========
            % The five options below (ThetaDHET, ThetaDPHO, HETCardinal, Grazing,
            % GrazingTheta) exist only for the temperature-asymmetry screen, which asks whether
            % the model can produce more biomass in winter than in summer. None of them is
            % part of the published Manriquez2026 model, none has been accepted, and
            % EVERY ONE DEFAULTS TO OFF (NaN, NaN, [NaN NaN NaN], 0.0 respectively), so the
            % regression goldens are unaffected by their presence. Kept in the signature
            % deliberately (Jaime, 2026-09-22: "keep, mark as experimental") rather than
            % deleted; do not cite a number produced with any of them switched on without
            % saying which arm it came from.
            % -------------------------------------------------------------------------
            % Overrides the death theta only, for the temperature-asymmetry screen
            % NaN = the ThetaSet value.
            options.ThetaDHET (1,1) {mustBeNumeric} = NaN;
            options.ThetaDPHO (1,1) {mustBeNumeric} = NaN;
            % CTMI (Rosso 1993) on heterotroph GROWTH only, death stays exponential;
            % anchored at mu(20 C) = NominalRate by computeRate. [NaN NaN NaN] = off.
            options.HETCardinal (1,3) {mustBeNumeric} = [NaN NaN NaN];
            % Implicit-grazer predation closure on heterotrophs (plan
            % 2026-09-12-temperature-asymmetry-screen). 0.0 = off (default, goldens
            % unchanged); see the "Heterotroph grazing" reaction below for provenance.
            options.Grazing (1,1) {mustBeNumeric} = 0.0;
            options.GrazingTheta (1,1) {mustBeNumeric} = 1.08;
            % ================== end EXPERIMENTAL block ===============================
    end
    respirationForm = options.RespirationForm;
    if respirationForm == ""
        if options.PGExcess, respirationForm = "pgexcess"; else, respirationForm = "wolf"; end
    end

    switch options.ThetaSet
        case "mixed"
            thetaMuHet = 1.0725; thetaMuPho = 1.047; thetaDHet = 1.0725; thetaDPho = 1.080; thetaHyd = 1.0725;
        case "reichert"
            thetaMuHet = 1.0725; thetaMuPho = 1.047; thetaDHet = 1.0725; thetaDPho = 1.047; thetaHyd = 1.0725;
        case "campos"
            thetaMuHet = 1.047; thetaMuPho = 1.066; thetaDHet = 1.08; thetaDPho = 1.08; thetaHyd = 1.08;
    end
    % Death-theta overrides only, for the temperature-asymmetry screen (plan
    % 2026-09-12-temperature-asymmetry-screen). NaN = the ThetaSet value above.
    if ~isnan(options.ThetaDHET), thetaDHet = options.ThetaDHET; end
    if ~isnan(options.ThetaDPHO), thetaDPho = options.ThetaDPHO; end

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
    % *** COUPLED TO f_dry (proposed f_dry = 0.25, not applied).
    % *** If f_dry is applied to rho_P and the influents, the state variable becomes DRY
    % *** mass and this MUST revert to Tenore's 210. Changing one without the other
    % *** double-counts the water. The two values are the same physics.
    attenuationParticle = options.Attenuation;
    dispersivityParticle = 1.20E-02;
    HET =  mpcssf.Particle(Name="HET", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             TransportRate=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=2.735E+04,...
                             Attenuation=attenuationParticle);
    PHO =  mpcssf.Particle(Name="PHO", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             TransportRate=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=2.735E+04,...
                             Attenuation=attenuationParticle);
    POM =  mpcssf.Particle(Name="POM", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             TransportRate=5.47, AttachmentSand=0.0, AttachmentMatrix=0.0,...
                             Attenuation=attenuationParticle);
    PAT =  mpcssf.Particle(Name="PAT", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                             TransportRate=5.47, AttachmentSand=5.47E+02, AttachmentMatrix=2.735E+04,...
                             Attenuation=attenuationParticle);
                    
    densityLiquid = 0.998E+03;
    dispersivityLiquid = 1.20E-02;
    O2 =  mpcssf.Liquid(Name="O2", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            TransportRate=6.00E+02);
    IC =  mpcssf.Liquid(Name="IC", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            TransportRate=6.00E+02);
    NH4 =  mpcssf.Liquid(Name="NH4", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            TransportRate=6.00E+02);
    HPO4 =  mpcssf.Liquid(Name="HPO4", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            TransportRate=6.00E+02);
    DOM =  mpcssf.Liquid(Name="DOM", Density=densityLiquid, Dispersivity=dispersivityLiquid, ...
                            TransportRate=3.00E+02);


    % Particle is a VALUE class, so a loop variable would modify a copy. Assign
    % directly. POM is deliberately untouched: its AttachmentSand is 0, so its
    % factor is inert, and PAT is governed separately by SandPathogen.
    HET.AttachmentSand = options.SandBiomass*HET.AttachmentSand;
    PHO.AttachmentSand = options.SandBiomass*PHO.AttachmentSand;

    heterotrophGrowth =  mpcssf.Reaction(Name="Heterotroph growth", IsLightDependent=false, ...
            NominalRate=4.80, ...   % Tenore2021 Table 1 mu_max,2 (f2 = heterotrophs)
            TemperatureCorrectionFactor=thetaMuHet, ...
            Order=dictionary("HET", 1), ...
            HalfSaturationConstants=dictionary("O2", 2.00E-04, "DOM", 4.00E-03, "NH4", 1.00E-06, "HPO4", 2.00E-05), ...
            StoichiometricCoefficients=dictionary("HET", 1.0, ...
                                        "O2", -1.2317, "IC", 0.3848, "NH4", -0.0248, "HPO4", -0.0141, "DOM", -1.5873));
    % CTMI (Rosso 1993) on heterotroph GROWTH only, for the temperature-asymmetry
    % screen; death stays
    % exponential. mu(20 C) = NominalRate is anchored by computeRate itself.
    if all(isfinite(options.HETCardinal))
        heterotrophGrowth.TemperatureResponse = "cardinal";
        heterotrophGrowth.CardinalTemperatures = options.HETCardinal;
    end

    phototrophGrowth =  mpcssf.Reaction(Name="Phototroph growth", IsLightDependent=true, ...
            MinimumLightFactor=0.0, ...
            OptimalLightFactor=ternaryOpt(options.NormalizedLight), ...
            NominalRate=2.00, ...
            TemperatureCorrectionFactor=thetaMuPho, ...
            Order=dictionary("PHO", 1), ...
            HalfSaturationConstants=dictionary("IC", 1.20E-03, "NH4", 2.00E-05, "HPO4", 2.00E-05), ...
            StoichiometricCoefficients=dictionary("PHO", 1.0, ...
                                        "O2", 0.9301, "IC", -0.3600, "NH4", -0.0600, "HPO4", -0.0100));

    heterotrophDeath =  mpcssf.Reaction(Name="Heterotroph death", IsLightDependent=false, ...
            NominalRate=0.40, ...
            TemperatureCorrectionFactor=thetaDHet, ...
            Order=dictionary("HET", 1), ...
            StoichiometricCoefficients=dictionary("HET", -1.0, "POM", 0.9123, ...
                                        "O2", 0.0234, "NH4", 0.0653, "HPO4", 0.0209));

    % Phototroph loss = Campos2006 k_ra, "algae loss due to the combined effects of
    % respiration and excretion" (Brown & Barnwell 1987): 0.048-0.504 /d, avg
    % 0.276 /d, theta 1.08 (Campos2006 Table 3). One loss process, as in Campos
    % and in Wolf2007 (b_ina,PH); no separate respiration reaction is used
    % (Jaime, 2026-08-26). Stoichiometry stays RWQM1 (11): the O2 cost of the dead
    % carbon is paid downstream by heterotroph growth on hydrolysed POM.
    % mixed-set default theta_d_PHO = 1.080 = Campos2006 Table 2 p.880, theta_kra.
    % VERIFIED in the source 2026-09-04, and it pairs with NominalRate 0.276 =
    % Campos Table 3 p.881 kra (0.0115/h x 24). Briefly changed to 1.047 on
    % 2026-09-04 on a mistaken miscitation report; reverted. Now also settable
    % via ThetaSet.
    phototrophDeath =  mpcssf.Reaction(Name="Phototroph death", IsLightDependent=false, ...
            NominalRate=0.276, ...
            TemperatureCorrectionFactor=thetaDPho, ...
            Order=dictionary("PHO", 1), ...
            StoichiometricCoefficients=dictionary("PHO", -1.0, "POM", 0.6316, ...
                                        "O2", 0.2005, "NH4", 0.0221, "HPO4", 0.0037));

    hydrolysis =  mpcssf.Reaction(Name="Hydrolysis", IsLightDependent=false, ...
            NominalRate=3.00, ...
            TemperatureCorrectionFactor=thetaHyd, ...
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
        Y = 0.63;      % was options.PGYield, removed 2026-09-22: no caller ever set it
        phototrophGrowth.MinimumLightFactor = 0.0;
        reactionList(2) = phototrophGrowth;
        % NH4/HPO4 half-saturations are DEPLETION PROTECTION (far below any
        % ambient level): without them the consumed pools cross zero once
        % drained (bottom-bed enclosed NH4 at fine grids).
        phototrophRespiration = mpcssf.Reaction(Name="Phototroph respiration", ...
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
        %
        % "Transport-identical to PHO" is taken literally: PG carries PHO's
        % AttachmentMatrix (2.735e4, since 2026-09-05) and PHO's SandBiomass
        % scaling of AttachmentSand. It was missed on 2026-09-05 and left at
        % 5.47e2 with no scaling, which was an oversight, not a choice
        % (Jaime, 2026-09-07). Only the "wolf" branch builds PG and no
        % production driver uses it, so no recorded run is affected.
        f = 0.2;       % was options.PGFraction, removed 2026-09-22: no caller ever set it
        Y = 0.63;      % was options.PGYield, removed 2026-09-22: no caller ever set it
        PG = mpcssf.Particle(Name="PG", Density=densityParticle, Dispersivity=dispersivityParticle, ...
                      TransportRate=5.47, AttachmentSand=options.SandBiomass*5.47E+02, ...
                      AttachmentMatrix=2.735E+04, ...
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
        phototrophRespiration = mpcssf.Reaction(Name="Phototroph respiration", ...
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
        phototrophRespiration = mpcssf.Reaction(Name="Phototroph respiration", ...
                NominalRate=options.PhototrophRespiration, ...
                TemperatureCorrectionFactor=1.047, ...
                Order=dictionary("PHO", 1), ...
                HalfSaturationConstants=dictionary("O2", 2.00E-04), ...
                StoichiometricCoefficients=dictionary( ...
                        ["PHO", "O2", "IC", "NH4", "HPO4"], ...
                        [-1.0, -0.9301, 0.3600, 0.0600, 0.0100]));
        reactionList = [reactionList; phototrophRespiration];
    end
    % Cohesion and detachment: the AUDITED WORKING SET of Manriquez2026
    % tab:general-parameters and tab:rhs-parameters, as used for the manuscript
    % runs. Corrected 2026-09-05 --
    % the preset had drifted from the published model on all three:
    %
    %   Zeta0   1.00E+06 -> 1       tab:general-parameters, own calibration.
    %                               1e6 made the cohesive velocity dominate the
    %                               CFL bound so hard that a fixed dt = 1e-6 went
    %                               unphysical within a few steps on a developed
    %                               biofilm.
    %   Zeta1   1/100 -> 0.27       tab:general-parameters, \citet{SSFnum}.
    %   b^det   sqrt(v/7.2)         tab:rhs-parameters: nominal rate 0.14 /d and
    %           -> 0.14*|v|/18      nominal velocity 18 m/d, in the LINEAR form
    %                               the working set runs (the sqrt law agrees
    %                               only at v = 18). The old form also carried
    %                               the INTAKE velocity 7.2 in the slot the
    %                               nominal detachment velocity 18 belongs in,
    %                               and omitted abs(), so a negative face
    %                               velocity produced a complex rate.
    %
    % Every production driver already rebuilt the Model with exactly these three
    % values, so this changes nothing they compute -- it makes the preset itself
    % agree with the model that was published, instead of requiring three silent
    % overrides to get there.
    if options.Grazing > 0
        % Implicit-grazer predation closure on heterotrophs (PARTNER-2026-09-07 §C(b)):
        % protozoa are not a state variable; the grazer is taken proportional to its
        % prey, so the loss is HET*HET/(K+HET) -- the same shape as the marker
        % bacterivory of pathogenModel.m. K = 2e-3 is the marker's K_pred (uncited).
        % Egestion: 1 - E_h = 0.4 of grazed carbon returns as POM (Campos2006 Eq. 27
        % p.877; E_h = 0.6, Table 3 p.881); the assimilated 0.6 is mineralised to IC
        % (protozoa at quasi-steady state, k_dp + k_rp). Unlike detachment this is a
        % TRUE sink. GrazingTheta 1.08 = Campos2006 Table 2 p.880 theta_Cga; 1.0 =
        % temperature-independent (screen arm). Rate 0.8 /d in the screen is a
        % judgement value, not a citation. Default 0 = off, goldens unchanged.
        % O2 half-saturation 2e-4 (the same depletion protection as heterotrophGrowth)
        % added because this reaction's rate is HET-only-limited and would otherwise
        % drive O2 negative under its -0.6*1.0667 stoichiometric draw.
        heterotrophGrazing = mpcssf.Reaction(Name="Heterotroph grazing", IsLightDependent=false, ...
                NominalRate=options.Grazing, ...
                TemperatureCorrectionFactor=options.GrazingTheta, ...
                Order=dictionary("HET", 1), ...
                HalfSaturationConstants=dictionary("HET", 2.00E-03, "O2", 2.00E-04), ...
                StoichiometricCoefficients=dictionary("HET", -1.0, "POM", 0.4*0.9123, ...
                                            "IC", 0.6*0.3600, "O2", -0.6*1.0667, ...
                                            "NH4", 0.0653, "HPO4", 0.0209));
        reactionList = [reactionList; heterotrophGrazing];
    end
    model = mpcssf.Model(Components=componentList, ...
                 Kappa=1.00E-06, Zeta0=1, Zeta1=0.27, ...
                 DetachmentFunction=@(v) 0.14*abs(v)/18, ...
                 BiofilmPorosity=0.99, ...
                 OsmosisRate=1.00E-07, ...
                 Reactions=reactionList);

end

% MinimumLightFactor (dark-growth floor) is 0.0 since 2026-08-26. The published
% 0.01 let a covered filter photosynthesise at 1 % of optimum; it was retired
% inside the respiration branches only, so every driver that did not enable
% respiration kept it (season90 and the old cover data were contaminated).
% Growth, loss and hydrolysis rates and temperature factors (2026-08-25, second
% audit): rates from
% Reichert2001 Table 8 where Wolf2007 has none (mu_ALG 2.0, k_death,ALG 0.1,
% k_hyd 3.0) and Wolf2007 Table VI otherwise (b_ina,H 0.4, k_h 3.0); every theta
% is exp(beta) from Reichert (beta_ALG 0.046 -> 1.047, beta_H = beta_hyd 0.07 ->
% 1.0725), since PHOBIA has no temperature dependence. Old set for provenance:
%   mu_HET 1.81e-2 (th 1.047)  mu_PHO 5.5 (1.047)  d_HET 2.0 (1.066)
%   d_PHO 0.4 (1.08)  k_hyd 0.09 (1.08).  d_PHO was 0.10 (1.047, Reichert death
%   alone) between the second audit and the Campos merge of 2026-08-26.
% Half-saturation constants (2026-08-25): source-faithful values from Wolf2007
% Table VI and Reichert2001, replacing the published Table 2 set, which the
% 2026-08-20 audit showed to be adjacent-row slides of the cited tables. The old K_NH4,PHO =
% 1.2e-2 gave a Monod factor of 0.0017 at the influent, so phototrophs could
% not grow at any light level. Old set kept here for provenance:
%   HET: O2 3.0e-3, DOM 2.0e-4, NH4 4.0e-3, HPO4 1.4e-8
%   PHO: IC 2.0e-5, NH4 1.2e-2, HPO4 1.68e-4;  hydrolysis POM/HET 2.0e-5.
function opt = ternaryOpt(normalized)
if normalized, opt = 1.0; else, opt = 1.814E-02; end
end
