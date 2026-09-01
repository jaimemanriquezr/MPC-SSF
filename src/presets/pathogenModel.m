function model = pathogenModel(options)
% PATHOGENMODEL  Physical pathogen model: MATLAB port of the Julia campaign
% model `pathogen_model()` in julia/analysis/pathogen_repro.jl.
%
%   model = pathogenModel()
%   model = pathogenModel(WaterFactor=1e-3, SandPathogen=0.0)
%
% = modelLund (rho_P = 1117, physical units) + the three marker reactions of
% Manriquez2026 Table B.4: aerobic growth r7 (mu_PAT = 0.2/d), inactivation r8
% (d_PAT = 0.02/d), bacterivory r9 (p_PAT = 8.0/d, HET-promoted, the only
% reaction active in the flowing phase, scaled there by WaterFactor). The Lund
% reactions are made inert in the flowing phase. PAT gets Transport x40
% (srun_pathogen convention) and SandAttachmentFactor = SandPathogen. Cohesion
% Zeta0 and the detachment law default to the Julia campaign values (1e2 and
% 1.4e-5*sqrt(|v|/18)) but are now OPTIONS -- see Zeta0/DetachForm/DetachScale
% below -- so the preset can be run on the working set (zeta0 = 1, linear
% detachment) without rebuilding the Model by hand.
%
% AUDIT 2026-08-27 (.claude/decisions/2026-08-27-pathogen-model-audit.md).
% The Lund base is INHERITED, not duplicated: every Lund constant comes from
% modelLund() below and is verified field-by-field by analysis/testPathogen.m.
% The marker reactions' Monod set is now the audited heterotroph-growth set of
% .claude/decisions/2026-08-25-half-saturations-corrected.md, because the
% manuscript defines r7 with the heterotroph stoichiometry and gives no
% PAT-specific half-saturations of its own (tab:eco-parameters lists mu/d/p for
% PAT and nothing else). HPO4 was ADDED: manuscripts/AWR-SSF/pathogen.tex:33
% defines M^PAT as min over {O2, NH4, HPO4, DOM}, and the solver combines Monod
% terms with min (src/@State/simulate.m:1129), so the omission was a departure
% from the manuscript's own equation.
%
% CAVEATS. (1) tab:eco-parameters lists no K_pred; the thesis_model value 2e-3
% is kept. (2) It lists no theta for any PAT process either; 1.047 (growth) and
% 1.08 (inactivation, bacterivory) are kept, matching the Lund growth/death
% thetas of the published table. (3) pathogen.tex:60 sets the flowing-phase
% bacterivory efficiency E^PAT_{f,p} = 0.01, but WaterFactor defaults to 1e-3
% here (the Julia campaign value); NOT changed, because it would move every
% number in the completed OAT/Sobol campaign. Pass WaterFactor=1e-2 for the
% manuscript value. Distinct from modelPathogen (thesis units) -- do NOT mix.
arguments
    options.WaterFactor (1,1) {mustBeNumeric} = 1e-3;
    options.SandPathogen (1,1) {mustBeNumeric} = 0.0;
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
    % All three are multiplied by DetachScale. "sqrt"/"linear" match
    % probeChain.m's detachmentFunction exactly.
    options.DetachForm (1,1) string ...
        {mustBeMember(options.DetachForm, ["campaign", "sqrt", "linear"])} = "campaign";
    options.DetachScale (1,1) {mustBeNumeric} = 1.0;
    % true (default) = every inherited Lund reaction is made inert in the
    % flowing suspension (EfficiencyFlowing = 0), the Julia campaign convention.
    % KEEP IT TRUE for pulses on the E1-E7 chain snapshots: probeChain.m:77
    % builds its model from pathogenModel(), so those snapshots were themselves
    % grown with the Lund reactions flowing-inert. false restores modelLund's
    % own flowing efficiencies (all 1.0) -- a DIFFERENT host model.
    options.LundFlowingInert (1,1) logical = true;
    % Forwarded to modelLund (0.0 = off; see modelLund.m for provenance).
    options.PhototrophRespiration (1,1) {mustBeNumeric} = 0.0;
    options.PGExcess (1,1) logical = false;
    options.NormalizedLight (1,1) logical = false;
    options.RespirationLightK (1,1) {mustBeNumeric} = 1.0;
    options.RespirationForm (1,1) string = "";
end

lund = modelLund(PhototrophRespiration=options.PhototrophRespiration, ...
                 PGExcess=options.PGExcess, ...
                 NormalizedLight=options.NormalizedLight, ...
                 RespirationLightK=options.RespirationLightK, ...
                 RespirationForm=options.RespirationForm);
components = lund.Components;
for i = 1:length(components)
    if isa(components(i), "Particle") && components(i).Name == "PAT"
        components(i).TransportRate = 40*components(i).TransportRate;
        components(i).SandAttachmentFactor = options.SandPathogen;
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
markerGrowth = Reaction(Name="MarkerGrowth", ...
    NominalRate=options.GrowthRate, TemperatureCorrectionFactor=1.047, ...
    EfficiencyFlowing=0.0, Order=dictionary("PAT", 1.0), ...
    HalfSaturationConstants=dictionary(["O2", "NH4", "HPO4", "DOM"], ...
                                       [2.0e-4, 1.0e-6, 2.0e-5, 4.0e-3]), ...
    StoichiometricCoefficients=dictionary(["PAT", "O2", "DOM"], [1.0, -1.2317, -1.5873]));

inactivation = Reaction(Name="Inactivation", ...
    NominalRate=options.InactivationRate, TemperatureCorrectionFactor=1.08, ...
    EfficiencyFlowing=0.0, Order=dictionary("PAT", 1.0), ...
    StoichiometricCoefficients=dictionary("PAT", -1.0));

bacterivory = Reaction(Name="Bacterivory", ...
    NominalRate=options.BacterivoryRate, TemperatureCorrectionFactor=1.08, ...
    EfficiencyFlowing=options.WaterFactor, Order=dictionary("PAT", 1.0), ...
    HalfSaturationConstants=dictionary("HET", options.KPred), ...
    StoichiometricCoefficients=dictionary(["PAT", "HET"], [-1.0, 1.0]));

model = Model(Components=components, ...
    Reactions=[reactions; markerGrowth; inactivation; bacterivory], ...
    Kappa=lund.CohesionSubModel.Kappa, Zeta0=options.Zeta0, Zeta1=lund.CohesionSubModel.Zeta1, ...
    DetachmentFunction=detachmentLaw(options.DetachForm, options.DetachScale), ...
    WaterDensity=lund.WaterDensity, BiofilmPorosity=lund.BiofilmPorosity, ...
    OsmosisRate=lund.OsmosisRate);
end

function f = detachmentLaw(form, scale)
% k_det(v_f) in 1/d. "sqrt"/"linear" are probeChain.m's two forms verbatim;
% "campaign" is the historical pathogenModel law (1e4x smaller prefactor).
switch form
    case "campaign", f = @(v) scale*1.4e-5*sqrt(abs(v)/18);
    case "sqrt",     f = @(v) scale*0.14*sqrt(abs(v)/18);
    case "linear",   f = @(v) scale*0.14*abs(v)/18;
end
end
