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
% Zeta0 is overridden to 1e2 and the detachment law to 1.4e-5*sqrt(|v|/18),
% both matching the Julia campaign model.
%
% CAVEATS carried over verbatim from pathogen_repro.jl: (1) Table B.4 lists no
% K_pred; the thesis_model value 2e-3 is kept. (2) The r7 stoichiometry is
% modelled like heterotroph growth (O2/NH4/DOM Monod, O2 consumed); confirm
% against Section 4. Distinct from modelPathogen (thesis units) -- do NOT mix.
arguments
    options.WaterFactor (1,1) {mustBeNumeric} = 1e-3;
    options.SandPathogen (1,1) {mustBeNumeric} = 0.0;
    options.GrowthRate (1,1) {mustBeNumeric} = 0.2;
    options.InactivationRate (1,1) {mustBeNumeric} = 0.02;
    options.BacterivoryRate (1,1) {mustBeNumeric} = 8.0;
    options.KPred (1,1) {mustBeNumeric} = 2e-3;
    % Forwarded to modelLund (0.0 = off; see modelLund.m for provenance).
    options.PhototrophRespiration (1,1) {mustBeNumeric} = 0.0;
    options.PGExcess (1,1) logical = false;
    options.NormalizedLight (1,1) logical = false;
    options.RespirationLightK (1,1) {mustBeNumeric} = 1.0;
end

lund = modelLund(PhototrophRespiration=options.PhototrophRespiration, ...
                 PGExcess=options.PGExcess, ...
                 NormalizedLight=options.NormalizedLight, ...
                 RespirationLightK=options.RespirationLightK);
components = lund.Components;
for i = 1:length(components)
    if isa(components(i), "Particle") && components(i).Name == "PAT"
        components(i).TransportRate = 40*components(i).TransportRate;
        components(i).SandAttachmentFactor = options.SandPathogen;
    end
end

reactions = lund.Reactions;
for i = 1:length(reactions)
    reactions(i).EfficiencyFlowing = 0.0;   % Lund reactions inert in the flowing phase
end

markerGrowth = Reaction(Name="MarkerGrowth", ...
    NominalRate=options.GrowthRate, TemperatureCorrectionFactor=1.047, ...
    EfficiencyFlowing=0.0, Order=dictionary("PAT", 1.0), ...
    HalfSaturationConstants=dictionary(["O2", "NH4", "DOM"], [3.0e-3, 4.0e-3, 2.0e-4]), ...
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
    Kappa=lund.CohesionSubModel.Kappa, Zeta0=1e2, Zeta1=lund.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 1.4e-5*sqrt(abs(v)/18), ...
    WaterDensity=lund.WaterDensity, BiofilmPorosity=lund.BiofilmPorosity, ...
    OsmosisRate=lund.OsmosisRate);
end
