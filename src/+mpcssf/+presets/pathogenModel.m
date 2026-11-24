function model = pathogenModel(options)
% PATHOGENMODEL  Physical pathogen model: MATLAB port of the Julia campaign
% model `pathogen_model()` of the SSF.jl study scripts.
%
%   model = pathogenModel()
%   model = pathogenModel(WaterFactor=1e-3, SandPathogen=0.0)
%
% = modelLund (rho_P = 1117, physical units) + the three marker reactions of
% Manriquez2026 Table B.4: aerobic growth r7 (mu_PAT = 0.2/d), inactivation r8
% (d_PAT = 0.02/d), bacterivory r9 (p_PAT = 8.0/d, HET-promoted, scaled in the
% flowing phase by WaterFactor = E^PAT_{f,p} = 0.001). Every other reaction runs
% at efficacy 1 in all three phases, as ecomodel.tex:99 states -- see
% LundFlowingInert below, whose default flipped on 2026-09-07. PAT gets Transport x40
% (srun_pathogen convention) and AttachmentSand *= SandPathogen. Cohesion
% Zeta0 and the detachment law default to the Julia campaign values (1e2 and
% 1.4e-5*sqrt(|v|/18)) but are now OPTIONS -- see Zeta0/DetachForm/DetachScale
% below -- so the preset can be run on the working set (zeta0 = 1, linear
% detachment) without rebuilding the Model by hand.
%
% AUDIT 2026-08-27.
% The Lund base is INHERITED, not duplicated: every Lund constant comes from
% modelLund() below and is verified field-by-field by tests/diagnostics/testPathogen.m.
% The marker reactions' Monod set is now the audited heterotroph-growth set
% (half-saturations corrected 2026-08-25), because the
% manuscript defines r7 with the heterotroph stoichiometry and gives no
% PAT-specific half-saturations of its own (tab:eco-parameters lists mu/d/p for
% PAT and nothing else). HPO4 was ADDED: the manuscript's pathogen section
% defines M^PAT as min over {O2, NH4, HPO4, DOM}, and the solver combines Monod
% terms with min (evaluateReactions in @State/simulate.m), so the omission was a departure
% from the manuscript's own equation.
%
% CAVEATS. (1) tab:eco-parameters lists no K_pred; the thesis_model value 2e-3
% is kept. (2) It lists no theta for any PAT process either; 1.047 (growth) and
% 1.08 (inactivation, bacterivory) are kept, matching the Lund growth/death
% thetas of the published table. (3) pathogen.tex:60 sets the flowing-phase
% bacterivory efficiency E^PAT_{f,p} = 0.001, which is exactly WaterFactor's
% default here; the "0.01" this comment used to cite was a misreading of the tex,
% which has said 0.001 since April (corrected 2026-09-07). Bacterivory is
% therefore the manuscript's OWN documented exception to "all efficacies are 1"
% (Jaime, 2026-09-07) and the only reaction still damped in the flowing
% suspension. Distinct from modelPathogen (thesis units) -- do NOT mix.
%
% THIS IS THE AUDITED MANRIQUEZ2026 WORKING SET, NOT A GENERAL PATHOGEN MODEL.
% Its Lund base is modelLund() -- see that file's header for the working-set
% warning and for the run-derived constants (Attenuation = 52 from
% `chain_fld2x_lit_leg6`, COUPLED TO f_dry; the 2026-09-05 Zeta0/Zeta1/detachment
% correction). Everything that warning says about modelLund applies here through
% inheritance. On top of it, TWO marker constants below are conventions inherited
% from the Julia campaign rather than results in the manuscript, and both change
% what any sensitivity conclusion means:
%
%   TransportRate x40 (line ~119). The PAT transport rate is multiplied by a bare
%       40 after modelLund builds the component. This is the "srun_pathogen
%       convention" carried over from the Julia campaign driver; the factor has no
%       citation in Manriquez2026 and is not derived from a measured diffusivity.
%       It makes the marker 40x more mobile than the biomass particles it is
%       otherwise identical to.
%
%   SandPathogen = 0.0 (default, line ~120). This multiplies the marker's
%       AttachmentSand, so BY DEFAULT THE MARKER NEVER ATTACHES TO BARE SAND --
%       it can only be captured by existing biofilm (AttachmentMatrix). Every
%       sensitivity conclusion this model has produced about attachment is
%       conditional on that default: with SandPathogen = 0 the clean-bed removal
%       pathway is switched off by construction, so "attachment does not matter
%       early" is a statement about this configuration, not about the physics.
%       Pass SandPathogen > 0 to re-enable it.
arguments
    options.WaterFactor (1,1) {mustBeNumeric} = 1e-3;
    options.SandPathogen (1,1) {mustBeNumeric} = 0.0;
    % Pathogen-type knobs. Both default to 1,
    % i.e. the audited model above is unchanged unless they are passed.
    %   MatrixPathogen    multiplies PAT AttachmentMatrix. 0 = a marker that is
    %                     never captured by the biofilm matrix and can only be
    %                     removed after f -> e transport.
    %   TransportPathogen multiplies PAT TransportRate AFTER the x40 convention
    %                     below, so 1 = x40, 10 = x400, 0.1 = x4.
    options.MatrixPathogen (1,1) {mustBeNumeric} = 1.0;
    options.TransportPathogen (1,1) {mustBeNumeric} = 1.0;
    options.GrowthRate (1,1) {mustBeNumeric} = 0.2;
    options.InactivationRate (1,1) {mustBeNumeric} = 0.02;
    options.BacterivoryRate (1,1) {mustBeNumeric} = 8.0;
    options.KPred (1,1) {mustBeNumeric} = 2e-3;
    % Cohesion strength. Default 1e2 = the Julia campaign / manuscript value
    % (results.tex:79); the E1-E6 working set runs at 1.
    options.Zeta0 (1,1) {mustBeNumeric} = 1e2;
    % Detachment law k_det(v_f) in 1/d:
    %   "campaign" - 1.4e-5*sqrt(|v|/18)   (the Julia campaign value; default,
    %                so every existing caller is unchanged)
    %   "sqrt"     - 0.14*sqrt(|v|/18)     (manuscript Table rhs-parameters)
    %   "linear"   - 0.14*|v|/18           (the E3 working set)
    % All three are multiplied by DetachScale. "sqrt"/"linear" match the
    % detachment laws of the manuscript runs exactly.
    options.DetachForm (1,1) string ...
        {mustBeMember(options.DetachForm, ["campaign", "sqrt", "linear"])} = "campaign";
    options.DetachScale (1,1) {mustBeNumeric} = 1.0;
    % false (DEFAULT since 2026-09-07) = the inherited Lund reactions keep
    % modelLund's own flowing efficiencies, all 1.0. This is THE model:
    % ecomodel.tex:99 states that all efficacy constants are 1 for every
    % star in {b, e, f}, and since 2026-09-07 the code follows the manuscript.
    %
    % true = every inherited Lund reaction is made inert in the flowing
    % suspension (EfficiencyFlowing = 0), the Julia campaign convention that
    % EVERY run before 2026-09-07 used. Pass it true to reproduce those runs, and
    % for pulses on snapshots grown before 2026-09-07: those were themselves grown
    % flowing-inert -- a DIFFERENT host model.
    options.LundFlowingInert (1,1) logical = false;
    % Forwarded to modelLund (0.0 = off; see modelLund.m for provenance).
    options.PhototrophRespiration (1,1) {mustBeNumeric} = 0.0;
    options.PGExcess (1,1) logical = false;
    options.NormalizedLight (1,1) logical = false;
    options.RespirationLightK (1,1) {mustBeNumeric} = 1.0;
    % nu_P [m^2/kg]; NaN = modelLund's preset (52). 0.094 = pre-b1fcd03.
    % Multiplies HET/PHO AttachmentSand; NaN = modelLund's preset (1.0).
    options.SandBiomass (1,1) {mustBeNumeric} = NaN;
    options.Attenuation (1,1) {mustBeNumeric} = NaN;
    options.RespirationForm (1,1) string = "";
    % Forwarded to modelLund. "reichert" = the preset default since 2026-09-07.
    % "mixed"
    % reproduces every run before that date byte-identically.
    options.ThetaSet (1,1) string {mustBeMember(options.ThetaSet, ["mixed", "reichert", "campos"])} = "reichert";
    % Forwarded to modelLund (temperature-asymmetry options).
    % Defaults are inert (NaN/NaN/[NaN NaN NaN]/0.0/1.08); see modelLund.m for
    % provenance.
    options.ThetaDHET (1,1) {mustBeNumeric} = NaN;
    options.ThetaDPHO (1,1) {mustBeNumeric} = NaN;
    options.HETCardinal (1,3) {mustBeNumeric} = [NaN NaN NaN];
    options.Grazing (1,1) {mustBeNumeric} = 0.0;
    options.GrazingTheta (1,1) {mustBeNumeric} = 1.08;
end

lundArgs = {"PhototrophRespiration", options.PhototrophRespiration, ...
            "PGExcess", options.PGExcess, ...
            "NormalizedLight", options.NormalizedLight, ...
            "RespirationLightK", options.RespirationLightK, ...
            "RespirationForm", options.RespirationForm, ...
            "ThetaSet", options.ThetaSet, ...
            "ThetaDHET", options.ThetaDHET, "ThetaDPHO", options.ThetaDPHO, ...
            "HETCardinal", options.HETCardinal, ...
            "Grazing", options.Grazing, "GrazingTheta", options.GrazingTheta};
if ~isnan(options.Attenuation)
    lundArgs = [lundArgs, {"Attenuation", options.Attenuation}];
end
if ~isnan(options.SandBiomass)
    lundArgs = [lundArgs, {"SandBiomass", options.SandBiomass}];
end
lund = mpcssf.presets.modelLund(lundArgs{:});
components = lund.Components;
for i = 1:length(components)
    if isa(components(i), "mpcssf.Particle") && components(i).Name == "PAT"
        % x40: the "srun_pathogen convention" inherited from the Julia campaign
        % driver. Uncited in Manriquez2026 and not derived from a diffusivity --
        % it simply makes the marker 40x more mobile than HET/PHO/POM.
        components(i).TransportRate = options.TransportPathogen*40*components(i).TransportRate;
        % SandPathogen defaults to 0, i.e. NO attachment to bare sand: the marker
        % is captured only by existing biofilm. Every attachment sensitivity
        % result from this preset is conditional on that default (see the header).
        components(i).AttachmentSand = options.SandPathogen*components(i).AttachmentSand;
        components(i).AttachmentMatrix = options.MatrixPathogen*components(i).AttachmentMatrix;
    end
end

reactions = lund.Reactions;
if options.LundFlowingInert
    for i = 1:length(reactions)
        reactions(i).EfficiencyFlowing = 0.0;   % Lund reactions inert in the flowing phase
    end
end

% Monod set = the AUDITED heterotroph-growth set (2026-08-25), since r7 carries
% the heterotroph stoichiometry and the manuscript gives no PAT-specific K.
%   O2   2.0e-4  Reichert2001 RWQM1 K_O2,H = 0.2 g/m3   (was 3.0e-3)
%   NH4  1.0e-6  Wolf2007 K_S,H,NH3 ~ 0 ("N never limits"); depletion-protective
%                stand-in                                (was 4.0e-3)
%   HPO4 2.0e-5  Reichert2001 K_HPO4,H = 0.02 g P/m3     (ABSENT before)
%   DOM  4.0e-3  Wolf2007 K_S,H,SS                       (was 2.0e-4)
% Stoichiometry = the FULL heterotroph liquid column of pathogen.tex:61-66,
% identical to modelLund's "Heterotroph growth" (modelLund.m:136-137) with HET
% replaced by PAT. Before 2026-09-07 only O2 and DOM were carried, so the marker
% Monod-limited on NH4 and HPO4 without ever consuming them; the code now
% follows the manuscript.
markerGrowth = mpcssf.Reaction(Name="MarkerGrowth", ...
    NominalRate=options.GrowthRate, TemperatureCorrectionFactor=1.047, ...
    EfficiencyFlowing=1.0, Order=dictionary("PAT", 1.0), ...
    HalfSaturationConstants=dictionary(["O2", "NH4", "HPO4", "DOM"], ...
                                       [2.0e-4, 1.0e-6, 2.0e-5, 4.0e-3]), ...
    StoichiometricCoefficients=dictionary(["PAT", "O2", "IC", "NH4", "HPO4", "DOM"], ...
                                          [1.0, -1.2317, 0.3848, -0.0248, -0.0141, -1.5873]));

% POM +1 is the r8 row of pathogen.tex:57. Before 2026-09-07 the code carried
% PAT -1 alone, i.e. inactivation destroyed mass; the marker's carbon now becomes
% particulate organic matter, as the manuscript's stoichiometric matrix says.
inactivation = mpcssf.Reaction(Name="Inactivation", ...
    NominalRate=options.InactivationRate, TemperatureCorrectionFactor=1.08, ...
    EfficiencyFlowing=1.0, Order=dictionary("PAT", 1.0), ...
    StoichiometricCoefficients=dictionary(["PAT", "POM"], [-1.0, 1.0]));

bacterivory = mpcssf.Reaction(Name="Bacterivory", ...
    NominalRate=options.BacterivoryRate, TemperatureCorrectionFactor=1.08, ...
    EfficiencyFlowing=options.WaterFactor, Order=dictionary("PAT", 1.0), ...
    HalfSaturationConstants=dictionary("HET", options.KPred), ...
    StoichiometricCoefficients=dictionary(["PAT", "HET"], [-1.0, 1.0]));

model = mpcssf.Model(Components=components, ...
    Reactions=[reactions; markerGrowth; inactivation; bacterivory], ...
    Kappa=lund.CohesionSubModel.Kappa, Zeta0=options.Zeta0, Zeta1=lund.CohesionSubModel.Zeta1, ...
    DetachmentFunction=detachmentLaw(options.DetachForm, options.DetachScale), ...
    WaterDensity=lund.WaterDensity, BiofilmPorosity=lund.BiofilmPorosity, ...
    OsmosisRate=lund.OsmosisRate);
end

function f = detachmentLaw(form, scale)
% k_det(v_f) in 1/d. "sqrt"/"linear" are the manuscript runs' two forms verbatim;
% "campaign" is the historical pathogenModel law (1e4x smaller prefactor).
switch form
    case "campaign", f = @(v) scale*1.4e-5*sqrt(abs(v)/18);
    case "sqrt",     f = @(v) scale*0.14*sqrt(abs(v)/18);
    case "linear",   f = @(v) scale*0.14*abs(v)/18;
end
end
