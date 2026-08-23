function results = simulate(obj, options, parameters)
arguments
    obj  State
    options struct = struct.empty;
    parameters.InflowConcentrations = []
    parameters.SimulationTime (1,1) {mustBeNumeric} = 1.0;
    parameters.TimeStep (1,1) = 1E-5;
    parameters.FrameNumber (1,1) {mustBeNumeric} = 200;
    parameters.CloggingFraction (1,1) {mustBeNumeric} = .99;
    parameters.IsUpwinded (1,1) = false;
    parameters.Quiet = false;

    % Adaptive (CFL) time-stepping: set TimeStep="adaptive" to enable. Ported
    % from slow-sand-filtration @SDfilter/run_biofilm.m (the authoritative
    % reference). Assumes the Lund/Rosenqvist model structure (5 reactions in
    % the order growth/growth/death/death/hydrolysis; particles HET,PHO,POM,PAT;
    % liquids incl. the two growth half-saturations) -- see modelLund.
    parameters.CFLFactor (1,1) {mustBeNumeric} = .99;
    parameters.AdaptiveVelocityFactor (1,1) {mustBeNumeric} = 0;
    parameters.AdaptiveTimeTolerance (1,1) {mustBeNumeric} = 1e-3;
    parameters.AdaptiveInitialDt (1,1) {mustBeNumeric} = 5e-7;
    parameters.AdaptiveMaxDt (1,1) {mustBeNumeric} = 5e-7;

    % Backward-Euler integration of the osmosis relaxation (port of Julia
    % simulate.jl `implicit_osmosis`; see implicit-osmosis.md). The relaxation is
    % (beta*phiB - phiE)/tau = A_osm - kosm*phiW with kosm = (1-beta)/tau, a stiff
    % linear decay in phiW. Damping the rate by 1/(1+dt*kosm) is exactly backward
    % Euler on that term, unconditionally stable, and lets the CFL bound drop
    % 1/tau so dt is no longer osmosis-bound (~20x speedup on modelLund).
    parameters.ImplicitOsmosis (1,1) = false;

    % Liebig-limitation diagnostic. When true, every frame additionally records
    % which substrate is currently binding the min over Monod terms, the value
    % of that min, and the two light arrays, for the biofilm and flowing
    % regions. Purely observational: the reaction rates are bit-identical with
    % it on or off, and nothing is allocated when it is off.
    parameters.RecordLimitation (1,1) logical = false;
    % Phase 0 of the Bailo investigation (.claude/plans/2026-08-23-bailo-scheme.md):
    % attribute the adaptive CFL bound to region and term, and price the
    % counterfactual in which the cohesive part of v_b were treated implicitly.
    % Exact counters accumulated every step; no per-step arrays.
    parameters.RecordCflBudget (1,1) logical = false;
    % Clamp the Cahn-Hilliard mobility to its positive part,
    % lambda = max(zeta_0*M(u), 0). For the REAL solve this is a no-op: u comes
    % from Solver B, which is guarded against negativity, so u is in [0,1] and
    % zeta_0*u(1-u) >= 0 already -- the goldens re-export bit-identically with
    % this on, which is the proof. It matters for the free-running parallel state
    % (RecordCflBudget), where u is NOT re-seeded: there u goes negative almost
    % immediately, mobility follows it negative, diffusion becomes
    % anti-diffusion, and the solution blows up. NOTE this is a pointwise clamp
    % of the centred mobility, NOT Bailo's two-point upwind M(x,y) =
    % zeta_0*(x)^+ (1-y)^+ -- it removes negative mobility but carries no
    % bound-preservation proof.
    parameters.PositiveMobility (1,1) logical = false;
    % Treat the flowing-phase dispersive flux implicitly. Phase 0 measured
    % dispersion at 94% of the binding region's sum at N=500 with the matrix
    % region's margin at 0.0688, so this is worth at most ~14.5x on dt; after it,
    % advection binds at dz/q = 2.78e-04 d. Components do not couple through
    % dispersion, so this is one tridiagonal solve per component. The dispersion
    % coefficient (|v_f|, phi_b) and phiFlowing are LAGGED, which keeps the
    % operator linear -- the lag has its own accuracy ceiling, which is exactly
    % what needs measuring. See PLANS.md, "Implicit Solver B".
    parameters.ImplicitDispersion (1,1) logical = false;
end
if ~isempty(options)
    for field = string(fieldnames(options)).'
        if isfield(parameters, field)
            parameters.(field) = options.(field);
        end
    end
end
disp("Loading parameters...")
filter = obj.SandFilter;
model = obj.Model;

temperature = filter.Temperature;
%================= LOADING OPTIONS ONTO VARIABLES ==================%
simulationTime = parameters.SimulationTime;
timeStep = parameters.TimeStep;
numFrames = parameters.FrameNumber;

[quotientK, quotientDenIdx, quotientNumIdx] = model.Reactions.lookupQuotients(model.Components);
halfSaturationK = model.HalfSaturationConstants;
listK = permute([halfSaturationK; quotientK], [3 1 2]);
listOrder = zeros(size(listK));
listOrder(1, 1:length(model.Components), :) = model.Order;

newlist = cell(2, size(listK, 3));
for i = 1:size(listK, 3)
newlist{1, i} = ~isnan(listK(:, :, i));
    if any(newlist{1, i})
        newlist{2, i} = listK(:, newlist{1, i}, i);
    end
end
listK = newlist;
[listOrder, ~] = ind2sub(size(listOrder, [2 3]), find(listOrder));

if isempty(parameters.InflowConcentrations)
    inflowConcentrations = zeros(1, length(model.Components));
elseif isa(parameters.InflowConcentrations, "dictionary")
    if isa(parameters.InflowConcentrations.keys, "string")
        inflowConcentrations = lookup(parameters.InflowConcentrations, [model.Components.Name], FallbackValue=0.0);
    else
        inflowConcentrations = lookup(parameters.InflowConcentrations, model.Components, FallbackValue=0.0);
    end
else
    inflowConcentrations = parameters.InflowConcentrations;
end
if ~isa(inflowConcentrations, 'function_handle')
    inflowConcentrations = inflowConcentrations(:).';
end

%================= I. SAND FILTER PARAMETERS ====================%
depthCenters = filter.GridPoints.Centers;
porosityCenters = computePorosity(filter, depthCenters);

depthBoundaries = filter.GridPoints.Boundaries;
porosityBoundaries = computePorosity(filter, depthBoundaries);

dz = filter.GridSize;
n0 = filter.GridZero;

[S, DD, D] = getCahnHilliardMatrices(filter, model, parameters.IsUpwinded);
S = sparse(S.Rows, S.Columns, S.Values, 2*n0, 2*n0);
DD = sparse(DD.Rows, DD.Columns, DD.Values, 2*n0, 2*n0);

Id = speye(2*n0);
CH0 = Id - DD;

etaWater = filter.LightAttenuationEtaWater;
etaSand = filter.LightAttenuationEtaSand;

volumetricFlow = filter.InflowVelocity;
volumeAvgVelocity = volumetricFlow./porosityBoundaries;
%========================================================%

%================= II. MODEL PARAMETERS ====================%
kP = length(model.Particles);
kL = length(model.Liquids);
densityP = [model.Particles.Density];
densityL = [model.Liquids.Density];

attachmentRates = [model.Particles.AttachmentSand];
transportParticleRates = [model.Particles.TransportRate];
transportLiquidRates = [model.Liquids.TransportRate];
alpha = [model.Particles.Dispersivity, model.Liquids.Dispersivity];

beta = model.BiofilmPorosity;
tau = model.OsmosisRate;

sigmaParticles = model.StoichiometricMatrixParticles;
sigmaLiquids = model.StoichiometricMatrixLiquids;

dpsi_fun = model.CohesionSubModel.PotentialGradient;
zeta_0 = model.CohesionSubModel.Zeta0;
mobility = model.CohesionSubModel.MobilityFunction;

% assumption 1: all components (of the same type) have the same density
densityP = mean(densityP);
densityL = mean(densityL);
%========================================================%

%=================== III. INITIAL CONDITIONS ======================%
startingConditions = obj;
timeStart = obj.Time;

globalBiofilm = startingConditions.GlobalConcentrationBiofilm;
globalFlowing = startingConditions.GlobalConcentrationFlowing;
phiW = startingConditions.EnclosedWaterVolume;
velBiofilm = startingConditions.Velocity.Biofilm;

% compute initial vf based on initial vb and phib(Cb,phiWe)
phiBiofilm = phiW + sum(globalBiofilm(:,1:kP),2)/densityP ...
            + sum(globalBiofilm(:,kP+1:kP+kP),2)/densityP ...
            + sum(globalBiofilm(:,kP+kP+1:kP+kP+kL),2)/densityL;
phiBiofilmBoundaries = .5*(phiBiofilm(2:end) + phiBiofilm(1:end-1));
velFlowing = [volumeAvgVelocity(1); (volumeAvgVelocity(2:end-1) - velBiofilm.*phiBiofilmBoundaries)./(1 - phiBiofilmBoundaries); volumeAvgVelocity(end)];
%=========================================================================%

%=================== IV. PRE-ALLOCATION ======================%
% ---------- GENERAL VARIABLES --------------%
timeFrames = zeros(numFrames,1);
concFramesBiofilm = zeros(length(depthCenters),numFrames, kP + kP + kL);
concFramesWater = zeros(length(depthCenters),numFrames, 1);
concFramesFlowing = zeros(length(depthCenters),numFrames, kP + kL);
[velFramesBiofilm, velFramesFlowing] = deal(zeros(length(depthCenters)-1,numFrames,1));

% ---------- LIMITATION DIAGNOSTIC (opt-in) --------------%
% Allocated only when requested. uint8 for the argmin, single for the values:
% ~13 MB for a 100-cell/2160-frame/6-reaction run.
if parameters.RecordCflBudget
    cflBudget = struct( ...
        "Steps", 0, "CapBound", 0, ...          % dt pinned to AdaptiveMaxDt, not the CFL
        "RegionWins", zeros(1,5), ...           % argmax over [matrix, encP, encL, flowP, flowL]
        "TermSums", zeros(1,4), ...             % share of w_v/dz, w_a/dz^2, w_b, w_s in that region
        "XSum", 0, "XMin", Inf, "XMax", 0, ...  % dt_CFL(advective v_b only) / dt_CFL(actual)
        "VbCohSum", 0, "VbCohMax", 0, ...       % cohesive share of max|v_b|
        "RegionShareSum", zeros(1,5), ...       % mean total(region)/max(total): the MARGIN
        "RegionShareMax", zeros(1,5), ...       % worst-case approach of each region
        "MatrixShareNoCoh", 0, ...              % same for region 1 with cohesion removed
        "WsP50Sum", 0, "WsP90Sum", 0, "WsP99Sum", 0, ...   % spread of w_s ACROSS CELLS
        "WsCellMaxSum", 0, "WsBoundSum", 0, ...  % per-cell max vs the bound actually used
        "WsConcSum", 0, "WsConcMax", 0, "WsConcN", 0);  % max/median: few cells setting dt?
    phibCHFrames = nan(filter.GridZero, numFrames);
    % Free-running CH state: seeded once from phi_b, then advanced by Solver A
    % alone and NEVER re-seeded from Solver B. phibCHFrames above is re-seeded
    % every step, so it can only ever show a one-step (O(dt)) difference; this
    % one accumulates, which is what the component-elimination design actually
    % depends on -- there Solver A's phi_b becomes authoritative for the whole run.
    phibParFrames = nan(filter.GridZero, numFrames);
    uPar = [];
    parDiag = struct("tFirstNeg", NaN, "tFirstAbove1", NaN, "tFirstNaN", NaN, ...
                     "minSeen", Inf, "maxSeen", -Inf, "dead", false);
end

if parameters.RecordLimitation
    nRx = numel(model.Reactions);   % == size(muRates, 2), but muRates is built later
    limFramesBiofilm  = zeros(length(depthCenters), numFrames, nRx, "uint8");
    limFramesFlowing  = zeros(length(depthCenters), numFrames, nRx, "uint8");
    monodFramesBiofilm = nan(length(depthCenters), numFrames, nRx, "single");
    monodFramesFlowing = nan(length(depthCenters), numFrames, nRx, "single");
    lightFramesAtten  = nan(length(depthCenters), numFrames, "single");
    lightFramesFactor = nan(length(depthCenters), numFrames, nRx, "single");
end
%========================================================%

%=================== V. OUTPUT RESULTS ======================%
results =  Results(filter, model);
results.TimeStart = timeStart;
results.SimulationData.time_final_intended = timeStart + simulationTime;
results.SimulationData.time = timeStart;
results.Flag = "OK";

timeSnap = linspace(timeStart,timeStart + simulationTime,numFrames);
counter = 1;
timeFrames(counter) = timeStart;
concFramesBiofilm(:,counter,:) = globalBiofilm;
concFramesWater(:,counter,:) = densityL*phiW;
concFramesFlowing(:,counter,:) = globalFlowing;

velFramesBiofilm(:,counter,:) = velBiofilm;
velFramesFlowing(:,counter,:) = velFlowing(2:end-1);

counter = counter + 1;
%=========================================================================%

%=============== VI. TIME INTEGRATION ====================================%
% Time-stepping mode: a numeric TimeStep is a fixed dt; TimeStep="adaptive"
% recomputes dt each step from a per-step CFL bound (ported from
% slow-sand-filtration @SDfilter/run_biofilm.m).
if isnumeric(timeStep)
    adaptivity = "fixed";
    dt = timeStep;
    results.SimulationData.TimeStep = dt;
else
    adaptivity = "adaptive";
    dt = parameters.AdaptiveInitialDt;
    results.SimulationData.TimeStep = "adaptive";
end
cflFactor = parameters.CFLFactor;
adaptiveVelocityFactor = parameters.AdaptiveVelocityFactor;
adaptiveTimeTolerance = parameters.AdaptiveTimeTolerance;
adaptiveMaxDt = parameters.AdaptiveMaxDt;

% CFL constants (Lund-structured model assumption; see the arguments block).
alphaP = alpha(1);
alphaL = alpha(kP + 1);
% liquid half-saturation constants for the two growth reactions; absent -> inf
% so the corresponding 1/(S+K) contribution vanishes.
halfSatLiquids = model.HalfSaturationConstants(kP+1:kP+kL, :);
% Guarded for non-Lund models with < 2 reactions (constants only feed the
% adaptive CFL bound; fixed-step runs never read them).
[K_HetGrowth, K_PhoGrowth] = deal(inf(kL, 1));
if size(halfSatLiquids, 2) >= 1, K_HetGrowth = halfSatLiquids(:, 1); end
if size(halfSatLiquids, 2) >= 2, K_PhoGrowth = halfSatLiquids(:, 2); end
K_HetGrowth(isnan(K_HetGrowth)) = inf;
K_PhoGrowth(isnan(K_PhoGrowth)) = inf;
K_CFL = [K_HetGrowth, K_PhoGrowth];
% hydrolysis quotient (POM/HET) half-saturation constant
if isempty(quotientK)
    K_Hyd = inf;
else
    K_Hyd = max(quotientK(:));
end

t = timeStart;
% Use the Model method (arrayfun over reactions); calling
% model.Reactions.computeRate(...) on the reaction array returns only the first
% reaction's rate, which is wrong whenever rates differ between reactions.
muRates = reshape(model.computeReactionRates(temperature), 1, []);
lightOptimal = max([model.Reactions.OptimalLightFactor]);
if isempty(lightOptimal) || lightOptimal <= 0
    lightOptimal = 1.0;   % inhibited-only models: keep the normalization finite
end
lightInhibition = [model.Reactions.LightInhibition];
inhibitionDependency = lightInhibition > 0;
complementDependency = [model.Reactions.IsLightComplement];
if any((inhibitionDependency + [model.Reactions.IsLightDependent] + complementDependency) > 1)
    error("light modes (dependent/inhibited/complement) are mutually exclusive");
end
attenuationParticles = [model.Particles.Attenuation];

minimumLight = [model.Reactions.MinimumLightFactor];
lightDependency = [model.Reactions.IsLightDependent];
minimumLight = minimumLight(lightDependency);
lightFactor = ones(length(filter.GridPoints.Centers), length(model.Reactions));

% Per-reaction phase efficiencies (1 x nRx, default 1.0). EfficiencyFlowing is
% the general form of the pathogen model's water_factor (every reaction inert in
% the flowing phase except bacterivory, scaled there by water_factor);
% EfficiencyBiofilm scales the biofilm/enclosed phases. Mirrors
% julia/src/simulate.jl:225-226.
efficiencyBiofilm = reshape([model.Reactions.EfficiencyBiofilm], 1, []);
efficiencyFlowing = reshape([model.Reactions.EfficiencyFlowing], 1, []);
% Per-particle bare-sand attachment scaling (1 x kP, default 1.0); the pathogen
% model's sand_pathogen. Mirrors julia/src/simulate.jl:228.
sandFactors = reshape([model.Particles.SandAttachmentFactor], 1, kP);

if isnumeric(inflowConcentrations)
    globalConcInflow = [inflowConcentrations(:)].';
end
disp("Starting simulation...")
while t < timeStart + simulationTime
    if isa(inflowConcentrations, 'function_handle')
        globalConcInflow = inflowConcentrations(t);
    end

    % ============= I. MAIN: computing state variables ================%
    %% GLOBAL CONCENTRATIONS
    % MATRIX
    globalMatrix = globalBiofilm(:, 1:kP);
    % ENCLOSED REGION
    globalEnclosedP = globalBiofilm(:, kP+1:kP+kP);
    globalEnclosedL = globalBiofilm(:, kP+kP+1:kP+kP+kL);
    % FLOWING REGION
    globalFlowingP = globalFlowing(:, 1:kP);
    globalFlowingL = globalFlowing(:, kP+1:kP+kL);

    %% VOLUME FRACTIONS
    phiMatrix = sum(globalMatrix,2)/densityP;
    phiEnclosed = phiW + sum(globalEnclosedP,2)/densityP + sum(globalEnclosedL,2)/densityL;
    phiBiofilm = phiMatrix + phiEnclosed;
    phiFlowing = 1 - phiBiofilm;

    % %============= CHECK FOR filter CLOGGING =========================%
    problemCellClogging = find(phiBiofilm > parameters.CloggingFraction, 1);
    if ~isempty(problemCellClogging)
        warning('Clogged! %.0f%% biofilm!\nT = %f\nCELL = %i', 100*max(phiBiofilm), t, problemCellClogging - n0)
        results.Flag = "CLOGGED";
        results.SimulationData.error.description = "Biofilm volume fraction has surpassed clogging value.";
        results.SimulationData.error.problem_cells = problemCellClogging;
        results.SimulationData.error.time = t;
        break;
    end
    % %================================================================%

    %% LOCAL CONCENTRATIONS
    localBiofilmX = globalMatrix./(phiBiofilm + realmin);
    % localBiofilmX(isnan(localBiofilmX)) = 0;
    localBiofilmS = globalEnclosedL./(phiBiofilm + realmin);
    % localBiofilmS(isnan(localBiofilmS)) = 0;
    localBiofilmXS = [localBiofilmX, localBiofilmS];
    localBiofilmQuotients = localBiofilmXS(:, quotientNumIdx) ./ (localBiofilmXS(:, quotientDenIdx) + realmin);
    localBiofilm = [localBiofilmXS, localBiofilmQuotients];

    localEnclosedX = globalEnclosedP./(phiEnclosed + realmin);
    % localEnclosedX(isnan(localEnclosedX)) = 0;
    localEnclosedS = globalEnclosedL./(phiEnclosed + realmin);
    % localEnclosedS(isnan(localEnclosedS)) = 0;
    localEnclosedXS = [localEnclosedX, localEnclosedS];
    localEnclosedQuotients = localEnclosedXS(:, quotientNumIdx) ./ (localEnclosedXS(:, quotientDenIdx) + realmin);
    localEnclosed = [localEnclosedXS, localEnclosedQuotients];

    localFlowingX = globalFlowingP./phiFlowing;
    localFlowingS = globalFlowingL./phiFlowing;
    localFlowingXS = [localFlowingX, localFlowingS];
    localFlowingQuotients = localFlowingXS(:, quotientNumIdx) ./ (localFlowingXS(:, quotientDenIdx) + realmin);
    localFlowing = [localFlowingXS, localFlowingQuotients];

    %% ATTACHMENT RATES (likelihood based on available volume)
    attachmentEnclosedFactor = phiMatrix./phiBiofilm;
    attachmentEnclosedFactor(isnan(attachmentEnclosedFactor)) = 0;
    attachmentFlowingFactor = (1 - porosityCenters) + porosityCenters.*phiBiofilm;

    %================= III. MAIN: compute reaction terms ===============================================%
    %% REACTION TERMS
    % LIGHT IRRADIATION
    globalParticles = globalMatrix + globalEnclosedP + globalFlowingP;
    etaParticles = cumsum(sum(attenuationParticles .* globalParticles, 2))*dz;
    eta = etaWater + etaSand + etaParticles;

    light = filter.LightIrradiation(t);
    lightAttenuated = light * exp(-eta) / lightOptimal;
    lightEffective = lightAttenuated.*exp(1 - lightAttenuated);
    % Dark-respiration FLOOR (authoritative slow-sand-filtration form
    % I = max(fdark, I_eff*e^{1-I_eff}); @SDfilter/run_biofilm.m, run_pathogen.m).
    % Was the additive (minimumLight + lightEffective + |...|)/2 = max(0, min+eff),
    % which double-counts the baseline at high light. MinimumLightFactor stands in
    % for the global dark_respiration. Equal to the old form when either is 0.
    if any(lightDependency)
        lightFactor(:, lightDependency) = max(minimumLight, lightEffective);
    end
    % Dark-switch reactions (Wolf2007 r6): K/(K + I_local), same normalization
    % as lightAttenuated.
    for jInh = find(inhibitionDependency)
        K = lightInhibition(jInh);
        lightFactor(:, jInh) = K./(K + lightAttenuated);
    end
    % Complement reactions: 1 - Steele(I) (Steele <= 1, factor stays in [0,1])
    % -- on in darkness, zero at optimal light.
    for jCmp = find(complementDependency)
        lightFactor(:, jCmp) = 1 - lightEffective;
    end

    % Phase efficiencies scale each reaction per region (defaults 1.0, so this
    % reduces to the unscaled rates for every pre-existing preset).
    if parameters.RecordLimitation
        [rxB, monodB, limB] = evaluateReactions(localBiofilm, listK, phiBiofilm, muRates, lightFactor, listOrder);
        [rxF, monodF, limF] = evaluateReactions(localFlowing, listK, phiFlowing, muRates, lightFactor, listOrder);
        ecoRxBiofilm = efficiencyBiofilm.*rxB;
        ecoRxFlowing = efficiencyFlowing.*rxF;
    else
        ecoRxBiofilm = efficiencyBiofilm.*evaluateReactions(localBiofilm, listK, phiBiofilm, muRates, lightFactor, listOrder);
        ecoRxFlowing = efficiencyFlowing.*evaluateReactions(localFlowing, listK, phiFlowing, muRates, lightFactor, listOrder);
    end
    ecoRxEnclosed = efficiencyBiofilm.*evaluateReactions(localEnclosed, listK, phiEnclosed, muRates, lightFactor, listOrder);

    ecoRxM = ecoRxBiofilm*sigmaParticles';
    ecoRxPe = ecoRxEnclosed*sigmaParticles';
    ecoRxLe = (ecoRxBiofilm + ecoRxEnclosed)*sigmaLiquids';
    ecoRxPf = ecoRxFlowing*sigmaParticles';
    ecoRxLf = ecoRxFlowing*sigmaLiquids';

    attE = attachmentEnclosedFactor.*globalEnclosedP.*attachmentRates;
    % Flowing attachment splits into a bare-sand term (scaled per particle by
    % SandAttachmentFactor) and a biofilm term. With every factor 1.0 this is
    % identically attachmentFlowingFactor.*globalFlowingP.*attachmentRates.
    % Mirrors julia/src/simulate.jl:385-386. (N x 1) .* (1 x kP) -> N x kP.
    attachmentFlowingFactorP = (1 - porosityCenters).*sandFactors + porosityCenters.*phiBiofilm;
    attF = attachmentFlowingFactorP.*globalFlowingP.*attachmentRates;

    velFlowingCenters = .5*(velFlowing(2:end) + velFlowing(1:end-1));
    detM = model.DetachmentFunction(velFlowingCenters).*globalMatrix;

    discreteFickP = localFlowingX - localEnclosedX;
    transP = (phiEnclosed/beta).*discreteFickP.*transportParticleRates;
    discreteFickL = localFlowingS - localEnclosedS;
    transL = (phiEnclosed/beta).*discreteFickL.*transportLiquidRates;

    reactionsMatrix = ecoRxM + attE + attF - detM;
    reactionsEnclosedParticles = ecoRxPe - attE + transP;
    reactionsEnclosedLiquids = ecoRxLe + transL;
    reactionsFlowingParticles = ecoRxPf - attF + detM - transP;
    reactionsFlowingLiquids = ecoRxLf - transL;

    rhsBiofilm = [reactionsMatrix, reactionsEnclosedParticles, reactionsEnclosedLiquids];
    rhsFlowing = [reactionsFlowingParticles, reactionsFlowingLiquids];
    rhsEnclosedWater = (beta*phiBiofilm - phiEnclosed)/tau;
    kosm = (1 - beta)/tau;
    if parameters.ImplicitOsmosis
        rhsEnclosedWaterA = rhsEnclosedWater./(1 + dt*kosm);
    else
        rhsEnclosedWaterA = rhsEnclosedWater;
    end

    %=================== IV. SOLVER A: compute biofilm velocity ======================%
    rhsBiofilmVolume = sum(reactionsMatrix,2)/densityP ...
                        + sum(reactionsEnclosedParticles,2)/densityP ...
                        + sum(reactionsEnclosedLiquids,2)/densityL ...
                        + rhsEnclosedWaterA;
    rhsBiofilmVolume = rhsBiofilmVolume(1:n0);

    u = phiBiofilm(1:n0);
    uB = .5*(u(2:end) + u(1:end-1));

    lambda = zeta_0*mobility(uB);
    if parameters.PositiveMobility, lambda = max(lambda, 0); end
    lhsCH = CH0 - dt*(S + sparse(D.Rows, D.Columns, D.Values.*(lambda).', 2*n0, 2*n0));
    rhsCH = [u + dt*rhsBiofilmVolume(1:n0); dpsi_fun(u)];

    % solving linear system
    xCH = lhsCH \ rhsCH;
    % uCH = xCH(1:n0);
    muCH = xCH(n0+1:end);
    if parameters.RecordCflBudget
        % (a) Solver A's own phi_b for THIS step, normally discarded. Re-seeded
        % from Solver B every step, so this measures local (one-step) consistency.
        uCHrec = xCH(1:n0);

        % (b) The same operator run in PARALLEL, seeded once and never re-seeded.
        % Its own state sets its own mobility, so it is a genuinely independent
        % CH evolution. It shares the reaction source rhsBiofilmVolume, which is
        % deliberate: the source depends on all components, not just phi_b, so it
        % cannot be recomputed from uPar -- holding it common isolates the CH
        % transport discretisation as the sole origin of any drift.
        % Diagnostic only: nothing here feeds muCH, v_b, or Solver B.
        if isempty(uPar), uPar = u; end
        uBpar = .5*(uPar(2:end) + uPar(1:end-1));
        lambdaPar = zeta_0*mobility(uBpar);
        if parameters.PositiveMobility, lambdaPar = max(lambdaPar, 0); end
        lhsPar = CH0 - dt*(S + sparse(D.Rows, D.Columns, D.Values.*(lambdaPar).', 2*n0, 2*n0));
        if ~parDiag.dead
            xPar = lhsPar \ [uPar + dt*rhsBiofilmVolume(1:n0); dpsi_fun(uPar)];
            uPar = xPar(1:n0);
            % Record WHERE it first leaves the physical range, not just how far
            % it ends up: the failure mode is the point, since a mobility
            % zeta_0*u(1-u) that goes negative turns diffusion into
            % anti-diffusion and the blow-up is then self-reinforcing.
            parDiag.minSeen = min(parDiag.minSeen, min(uPar));
            parDiag.maxSeen = max(parDiag.maxSeen, max(uPar));
            if isnan(parDiag.tFirstNeg)    && any(uPar < 0),  parDiag.tFirstNeg    = t; end
            if isnan(parDiag.tFirstAbove1) && any(uPar > 1),  parDiag.tFirstAbove1 = t; end
            if any(~isfinite(uPar))
                parDiag.tFirstNaN = t;
                parDiag.dead = true;   % stop wasting a solve per step on NaNs
            end
        end
    end

    % computing vb
    velBiofilm(1:n0-1) = volumeAvgVelocity(2:n0) - zeta_0*(1 - uB).*diff(muCH)/dz;
    %================================================================%
    %====================== V. SOLVER B: compute concentrations =================================%
    phiBiofilmBoundaries = .5*(phiBiofilm(2:end) + phiBiofilm(1:end-1));
    velFlowing = [volumeAvgVelocity(1);
        (volumeAvgVelocity(2:end-1) - velBiofilm.*phiBiofilmBoundaries)./(1 - phiBiofilmBoundaries);
        volumeAvgVelocity(end)];

    %===================== TIME ADAPTIVITY (CFL) ==========================%
    % Per-step CFL bound, ported from slow-sand-filtration
    % @SDfilter/run_biofilm.m. dt grows by at most (1 + tolerance) per step
    % toward the CFL limit, capped at AdaptiveMaxDt. Weights: w_v advection,
    % w_a dispersion, w_b exchange (attach/detach/transfer/osmosis), w_s
    % ecological source terms; one entry per region (matrix, enclosed P,
    % enclosed L, flowing P, flowing L).
    if adaptivity == "adaptive"
        vbmax = (1 + adaptiveVelocityFactor)*max(abs(velBiofilm));
        vfmax = (1 + adaptiveVelocityFactor)*max(abs(velFlowing));
        maxPhib = max(phiBiofilm);
        phie_f_max = max(phiEnclosed./phiFlowing);
        det_vf = model.DetachmentFunction(velFlowingCenters);

        % liquid-consumption bound from the two growth reactions (HET, PHO)
        Xb_r = permute(localBiofilmX(:, [1 2]), [3 2 1]);   % 1 x 2 x N
        Xe_r = permute(localEnclosedX(:, [1 2]), [3 2 1]);
        Xf_r = permute(localFlowingX(:, [1 2]), [3 2 1]);
        sigmaLiquidsGrowth = sigmaLiquids(:, 1:2);          % kL x 2
        muGrowth = muRates(1:2);                            % 1 x 2
        L_b = -sum(sigmaLiquidsGrowth.*(muGrowth.*Xb_r)./(permute(localBiofilmS, [2 3 1]) + K_CFL), 2);
        L_e = -sum(sigmaLiquidsGrowth.*(muGrowth.*Xe_r)./(permute(localEnclosedS, [2 3 1]) + K_CFL), 2);
        L_f = -sum(sigmaLiquidsGrowth.*(muGrowth.*Xf_r)./(permute(localFlowingS, [2 3 1]) + K_CFL), 2);

        % hydrolysis quotient monod (POM/HET) per region
        Xi_b = localBiofilmX(:,3)./(localBiofilmX(:,3) + K_Hyd*localBiofilmX(:,1)); Xi_b(isnan(Xi_b)) = 0;
        Xi_e = localEnclosedX(:,3)./(localEnclosedX(:,3) + K_Hyd*localEnclosedX(:,1)); Xi_e(isnan(Xi_e)) = 0;
        Xi_f = localFlowingX(:,3)./(localFlowingX(:,3) + K_Hyd*localFlowingX(:,1)); Xi_f(isnan(Xi_f)) = 0;

        % death-reaction source bound (HET death rx 3, PHO death rx 4)
        ws_t0 = max(abs(sigmaParticles(1,3))*muRates(3), abs(sigmaParticles(2,4))*muRates(4));

        w_v = [2*vbmax*[1 1 1], 2*vfmax*[1 1]];
        if parameters.ImplicitDispersion
            w_a = zeros(1,5);      % solved implicitly, so out of the CFL bound
        else
            w_a = [0 0 0, ...
                   2*vfmax*alphaP*(1 + 1/(1 - maxPhib)), ...
                   2*vfmax*alphaL*(1 + 1/(1 - maxPhib))];
        end
        w_b = [max(det_vf), ...
               max(attachmentRates)*max(attachmentEnclosedFactor) + max(transportParticleRates)/beta, ...
               max([max(transportLiquidRates)/beta, ...
                    ~parameters.ImplicitOsmosis/tau]), ...
               max(attachmentRates)*max(attachmentFlowingFactor) + max(transportParticleRates)/beta*phie_f_max, ...
               max(transportLiquidRates)/beta*phie_f_max];
        w_s = [max(ws_t0, abs(sigmaParticles(3,5))*muRates(5)*max(Xi_b)), ...
               max(ws_t0, abs(sigmaParticles(3,5))*muRates(5)*max(Xi_e)), ...
               max(abs(L_b), [], 'all') + max(abs(L_e), [], 'all'), ...
               max(ws_t0, abs(sigmaParticles(3,5))*muRates(5)*max(Xi_f)), ...
               max(abs(L_f), [], 'all')];

        dt_CFL = cflFactor/max(w_v/dz + w_a/dz^2 + w_b + w_s);
        dt = min([dt_CFL, (1 + adaptiveTimeTolerance)*dt, adaptiveMaxDt]);

        if parameters.RecordCflBudget
            % 0a/0b: which region attains the max, and each term's share of it.
            totals = w_v/dz + w_a/dz^2 + w_b + w_s;
            [~, kReg] = max(totals);
            terms = [w_v(kReg)/dz, w_a(kReg)/dz^2, w_b(kReg), w_s(kReg)];
            % 0c: counterfactual with the cohesive part of v_b made implicit.
            % velBiofilm(1:n0-1) = volumeAvgVelocity(2:n0) - zeta_0*(1-uB).*diff(muCH)/dz,
            % so removing cohesion leaves the advective part alone; entries n0..end
            % (the bed, where v_b == 0) are unchanged.
            vbAdv = velBiofilm; vbAdv(1:n0-1) = volumeAvgVelocity(2:n0);
            vbmaxAdv = (1 + adaptiveVelocityFactor)*max(abs(vbAdv));
            w_v_adv = [2*vbmaxAdv*[1 1 1], 2*vfmax*[1 1]];
            dt_CFL_adv = cflFactor/max(w_v_adv/dz + w_a/dz^2 + w_b + w_s);
            X = dt_CFL_adv/dt_CFL;

            cflBudget.Steps = cflBudget.Steps + 1;
            cflBudget.CapBound = cflBudget.CapBound + (dt_CFL >= adaptiveMaxDt);
            cflBudget.RegionWins(kReg) = cflBudget.RegionWins(kReg) + 1;
            cflBudget.TermSums = cflBudget.TermSums + terms/sum(terms);
            cflBudget.XSum = cflBudget.XSum + X;
            cflBudget.XMin = min(cflBudget.XMin, X);
            cflBudget.XMax = max(cflBudget.XMax, X);
            % How close did each region come to winning? Reporting only the
            % argmax hides whether the matrix region lost by 1% or by 100x --
            % and that margin is what predicts finer meshes, since region 1 has
            % no dispersion term and its cohesive part scales as zeta_0*kappa/dz^4
            % while region 4 scales as dz^-2.
            shares = totals/max(totals);
            cflBudget.RegionShareSum = cflBudget.RegionShareSum + shares;
            cflBudget.RegionShareMax = max(cflBudget.RegionShareMax, shares);
            totalsAdv = w_v_adv/dz + w_a/dz^2 + w_b + w_s;
            cflBudget.MatrixShareNoCoh = cflBudget.MatrixShareNoCoh ...
                + totalsAdv(1)/max(totalsAdv);

            % How concentrated is w_s across cells? The bound uses max over ALL
            % cells and liquids, so one nearly-depleted cell can set dt for the
            % whole domain -- L ~ 1/(S + K) with K small. If max >> median the
            % constraint is an artefact of a few cells and could be attacked far
            % more cheaply than by making reactions implicit.
            LbC = reshape(max(abs(L_b), [], 1), [], 1);   % per cell, max over liquids
            LeC = reshape(max(abs(L_e), [], 1), [], 1);
            wsCell = LbC + LeC;                            % region-3 ecology, per cell
            sw = sort(wsCell); nsw = numel(sw);
            p50 = sw(max(1, round(0.50*nsw)));
            p90 = sw(max(1, round(0.90*nsw)));
            p99 = sw(max(1, round(0.99*nsw)));
            wsCellMax = sw(end);
            cflBudget.WsP50Sum = cflBudget.WsP50Sum + p50;
            cflBudget.WsP90Sum = cflBudget.WsP90Sum + p90;
            cflBudget.WsP99Sum = cflBudget.WsP99Sum + p99;
            cflBudget.WsCellMaxSum = cflBudget.WsCellMaxSum + wsCellMax;
            % w_s(3) sums two INDEPENDENT maxima (over L_b and over L_e), which
            % may be attained in different cells, so it is >= the per-cell max.
            % That gap is extra conservatism costing nothing to remove.
            cflBudget.WsBoundSum = cflBudget.WsBoundSum + w_s(3);
            % Only meaningful once the median is nonzero: at startup every cell
            % is zero and max/p50 would be Inf, poisoning the running mean.
            if p50 > 0
                conc = wsCellMax/p50;
                cflBudget.WsConcSum = cflBudget.WsConcSum + conc;
                cflBudget.WsConcMax = max(cflBudget.WsConcMax, conc);
                cflBudget.WsConcN   = cflBudget.WsConcN + 1;
            end

            cohShare = 1 - vbmaxAdv/max(vbmax, eps);
            cflBudget.VbCohSum = cflBudget.VbCohSum + cohShare;
            cflBudget.VbCohMax = max(cflBudget.VbCohMax, cohShare);
        end
    end
    %======================================================================%

    %%%% FLUX COMPUTING %%%%
    % BIOFILM
    fluxBiofilm = [
        0*globalBiofilm(1, :);
        globalBiofilm(1:end-1, :).*max(0, velBiofilm) ...
        + globalBiofilm(2:end, :).*min(0, velBiofilm);
        0*globalBiofilm(end, :)
        ];
    fluxBiofilmIn = porosityBoundaries(1:end-1).*fluxBiofilm(1:end-1, :);
    fluxBiofilmOut = porosityBoundaries(2:end).*fluxBiofilm(2:end, :);

    % FLOWING SUSPENSION
    localFlowingGradient = diff(localFlowingXS(1:end-1, :));
    dispersionStrength = abs(velFlowing(2:end-2)).*(1 - phiBiofilmBoundaries(1:end-1));
    dispersionFlux = [
        0*globalFlowing(1, :);
        (dispersionStrength.*localFlowingGradient.*alpha) / dz;
        0*globalFlowing(end-1, :);
        0*globalFlowing(end, :)
        ];

    convectionFlux = [
        volumeAvgVelocity(1)*globalConcInflow;
        globalFlowing(1:end-1, :).*max(0, velFlowing(2:end-1)) ...
        + globalFlowing(2:end, :).*min(0, velFlowing(2:end-1));
        volumeAvgVelocity(end)*globalFlowing(end, :)
        ];

    if parameters.ImplicitDispersion
        fluxFlowing = convectionFlux;          % dispersion applied after the update
    else
        fluxFlowing = convectionFlux - dispersionFlux;
    end
    fluxFlowingIn = porosityBoundaries(1:end-1).*fluxFlowing(1:end-1,:);
    fluxFlowingOut = porosityBoundaries(2:end).*fluxFlowing(2:end,:);

    % ENCLOSED WATER COMPONENT
    fluxWater = [
        0*phiW(1, :);
        phiW(1:end-1, :).*max(0, velBiofilm) ...
        + phiW(2:end, :).*min(0, velBiofilm);
        0*phiW(end, :)
        ];
    fluxWaterIn = porosityBoundaries(1:end-1).*fluxWater(1:end-1, :);
    fluxWaterOut = porosityBoundaries(2:end).*fluxWater(2:end, :);

    %======================= VI. MAIN: update cell values =================================%
    % ORDERING (time splitting, biofilm resolved first). The biofilm block and the
    % enclosed water carry NO dispersion -- w_a is zero for those regions and
    % fluxBiofilm/fluxWater are pure upwind advection -- and every source and flux
    % for them is evaluated at time n. So both advance from time-n data alone, and
    % phi_f^{n+1} = 1 - phiBiofilm^{n+1} is known BEFORE the flowing dispersion
    % solve. Level A of the sketch in PLANS.md: the flowing solve then uses
    % phi_f^{n+1} rather than the lagged phi_f^n, at no extra cost and still as 9
    % independent tridiagonals.
    globalBiofilm = globalBiofilm + (dt/dz)*(fluxBiofilmIn - fluxBiofilmOut)./porosityCenters + dt*rhsBiofilm;
    if parameters.ImplicitOsmosis
        rhsEnclosedWaterB = rhsEnclosedWater./(1 + dt*kosm);
    else
        rhsEnclosedWaterB = rhsEnclosedWater;
    end
    phiW = phiW + (dt/dz)*(fluxWaterIn - fluxWaterOut)./porosityCenters + dt*rhsEnclosedWaterB;

    globalFlowing = globalFlowing + (dt/dz)*(fluxFlowingIn - fluxFlowingOut)./porosityCenters + dt*rhsFlowing;
    if parameters.ImplicitDispersion
        % phi_f at n+1, from the blocks just advanced (mirrors :337-340).
        phiBiofilmNew = sum(globalBiofilm(:, 1:kP), 2)/densityP + phiW ...
            + sum(globalBiofilm(:, kP+1:kP+kP), 2)/densityP ...
            + sum(globalBiofilm(:, kP+kP+1:kP+kP+kL), 2)/densityL;
        phiFlowingNew = 1 - phiBiofilmNew;
        % dispersionStrength = |v_f| * (1 - phiBiofilm at faces), and that second
        % factor IS phi_f at the faces, so it moves to n+1 too. Only |v_f| stays at
        % time n, since it needs Solver A.
        phiFlowingFaceNew = .5*(phiFlowingNew(2:end) + phiFlowingNew(1:end-1));
        dispersionStrengthNew = abs(velFlowing(2:end-2)).*phiFlowingFaceNew(1:end-1);
        globalFlowing = solveImplicitDispersion(globalFlowing, dispersionStrengthNew, ...
            alpha, phiFlowingNew, porosityBoundaries, porosityCenters, dz, dt);
    end

    %========== CHECK IF CONCENTRATIONS ARE NEGATIVE ============%
    % Underflow floor: a sink acting on an exactly-zero pool (r6 consuming PG
    % before any is stored; POM hydrolysis) leaves O(realmin) negatives from
    % the Monod regularizer that would trip the strict guard below. Clamp only
    % denormal-scale negatives; genuine instabilities overshoot far beyond -1e-20.
    globalBiofilm(globalBiofilm < 0 & globalBiofilm > -1e-20) = 0;
    globalFlowing(globalFlowing < 0 & globalFlowing > -1e-20) = 0;
    problemCellBiofilm = find(globalBiofilm(:) < 0 | isnan(globalBiofilm(:)));
    if ~isempty(problemCellBiofilm)
        [pRow, pCol] = ind2sub(size(globalBiofilm), problemCellBiofilm(1));
        fprintf('Unphysical concentration in biofilm. \nTIME = %e\nCELL = %d (z = %g), STATE COLUMN = %d, VALUE = %e\n', ...
            t, pRow, depthCenters(pRow), pCol, globalBiofilm(pRow, pCol))

        results.Flag = "BIOFILM";
        results.SimulationData.error.description = "Concentration in biofilm volume has reached unphysical values.";
        results.SimulationData.error.problem_cells = problemCellBiofilm;
        results.SimulationData.error.time = t;
        break;
    end

    problemCellFlowing = mod(find(globalFlowing(:) < 0 | isnan(globalFlowing(:))), size(globalFlowing, 1));
    if ~isempty(problemCellFlowing)
        fprintf('Unphysical concentration in flowing suspension. \nTIME = %e\n',t)
        fprintf("CELL = %i\n", problemCellFlowing(1) - n0)
        fprintf("HEIGHT = %i\n", depthCenters(problemCellFlowing(1)));
        results.Flag = "BIOFILM";
        results.SimulationData.error.description = "Concentration in flowing volume has reached unphysical values.";
        results.SimulationData.error.problem_cells = problemCellFlowing;
        results.SimulationData.error.time = t;
        break;
    end
    %===========================================================%

    % NEXT STEP
    t = t + dt;
    results.SimulationData.time(end+1) = t;

    %============ SAVING FRAMES ================================%
    if t >= timeSnap(counter)
        timeFrames(counter) = t;

        concFramesBiofilm(:, counter, :) = globalBiofilm;
        concFramesFlowing(:, counter, :) = globalFlowing;
        concFramesWater(:, counter, :) = densityL*phiW;

        velFramesBiofilm(:, counter, :) = velBiofilm;
        velFramesFlowing(:, counter, :) = velFlowing(2:end-1);

        if parameters.RecordCflBudget
            phibCHFrames(:, counter) = uCHrec;
            phibParFrames(:, counter) = uPar;
        end

        if parameters.RecordLimitation
            limFramesBiofilm(:, counter, :)   = limB;
            limFramesFlowing(:, counter, :)   = limF;
            monodFramesBiofilm(:, counter, :) = single(monodB);
            monodFramesFlowing(:, counter, :) = single(monodF);
            lightFramesAtten(:, counter)      = single(lightAttenuated);
            lightFramesFactor(:, counter, :)  = single(lightFactor);
        end
        counter = counter + 1;
        if ~parameters.Quiet
            fprintf("t = %.4e\n", t);
        end
    end
    %==========================================================%
end
disp("Simulation ended.")
disp("Saving results...")

%====================== VII. OUTPUTS ===================%
results.Frames.Time = timeFrames;

concentrations_cell = cell(kP + kL + 1, 3);
for j = 1:kP
    concentrations_cell{j, 1} = concFramesBiofilm(:, :, j);
    concentrations_cell{j, 2} = concFramesBiofilm(:, :, kP + j);
    concentrations_cell{j, 3} = concFramesFlowing(:, :, j);
end

for j = 1:kL
    concentrations_cell{kP + j, 1} = 0;
    concentrations_cell{kP + j, 2} = concFramesBiofilm(:, :, kP + kP + j);
    concentrations_cell{kP + j, 3} = concFramesFlowing(:, :, kP + j);
end
concentrations_cell{kP + kL + 1, 1} = 0;
concentrations_cell{kP + kL + 1, 2} = concFramesWater;
concentrations_cell{kP + kL + 1, 3} = Inf;

results.Frames.Concentrations = cell2table(concentrations_cell,  ...
    'VariableNames', {'Matrix',  'Enclosed',  'Flowing'},  ...
    'rowNames', [results.Model.Particles.Name,  results.Model.Liquids.Name,  'Water'],  ...
    'DimensionNames', {'Component', 'Volume'});
results.Frames.Velocity.Biofilm = velFramesBiofilm;
results.Frames.Velocity.Flowing = velFramesFlowing;

if parameters.RecordCflBudget
    % SimulationData is the existing home for run diagnostics (see :193-221);
    % adding a Results property would change the class for every consumer.
    results.SimulationData.CflBudget = cflBudget;
    results.SimulationData.PhibCH = phibCHFrames;
    results.SimulationData.PhibPar = phibParFrames;
    results.SimulationData.PhibParDiag = parDiag;
end

if parameters.RecordLimitation
    results.Frames.Limitation.Biofilm = limFramesBiofilm;
    results.Frames.Limitation.Flowing = limFramesFlowing;
    results.Frames.Limitation.MonodBiofilm = monodFramesBiofilm;
    results.Frames.Limitation.MonodFlowing = monodFramesFlowing;
    results.Frames.Limitation.LightAttenuated = lightFramesAtten;
    results.Frames.Limitation.LightFactor = lightFramesFactor;
    results.Frames.Limitation.Names = limitationNames(model, quotientNumIdx, quotientDenIdx);
    results.Frames.Limitation.ReactionNames = [model.Reactions.Name];
end

results.TimeFinal = t;
disp("Results saved.")
end

function [rx, monod, lim] = evaluateReactions(local, K, phi, mu, I, orders)
    % monod = permute(min(min((local + realmin)./(K + local + realmin), [], 2), 1), [1 3 2]);
    % `lim` records which column of `local` supplied the min for each cell and
    % reaction (0 = the reaction has no Monod terms, so monod stays 1). The
    % argmin is already computed by min(); returning it costs nothing and `rx`
    % is unchanged. Column indices run over [Particles, Liquids, Quotients] --
    % see limitationNames() for the matching labels.
    monod = ones(size(phi, 1), size(mu, 2));
    lim = zeros(size(phi, 1), size(mu, 2), "uint8");
    for i = 1:size(monod, 2)
        ii = K{1, i};
        if any(ii)
            k = K{2, i};
            mon_term = (local(:, ii) + realmin)./(k + local(:, ii) + realmin);
            [monod(:, i), j] = min(mon_term, [], 2);
            gidx = uint8(find(ii));
            lim(:, i) = gidx(j);
        end
    end
    %product = permute(prod(local.^orders, 2, "omitnan"), [1 3 2]);
    product = local(:, orders);
    rx = phi.*mu.*I.*monod.*product;
end

function names = limitationNames(model, numIdx, denIdx)
    % Labels for the columns of `local`, i.e. the values `lim` takes.
    componentNames = [model.Components.Name];
    names = componentNames;
    for q = 1:numel(numIdx)
        names(end+1) = componentNames(numIdx(q)) + "/" + componentNames(denIdx(q)); %#ok<AGROW>
    end
end


function g = solveImplicitDispersion(g, S, alpha, phiFlowing, poroB, poroC, dz, dt)
% Backward-Euler the flowing-phase dispersive flux.
%
% Explicit form being replaced (see the FLUX COMPUTING block):
%   dispFlux_f = S_{f-1} * alpha_j * (c_f - c_{f-1}) / dz,   c = g ./ phiFlowing
%   g <- g + (dt/dz)*(poroB_l*flux_l - poroB_r*flux_r)/poroC,  flux = conv - disp
% so the dispersive contribution to the update is -DISP(g) with
%   DISP(g) = (dt/dz)*(poroB_l*dispFlux_l - poroB_r*dispFlux_r)/poroC,
% and taking it at n+1 gives (I + DISP)g^{n+1} = <explicit part>, already in g.
%
% Dispersion acts on the LOCAL concentration g./phiFlowing, not on g, so the
% diag(1./phiFlowing) belongs inside the operator. phiFlowing and S are lagged,
% which is what keeps this linear.
%
% Components couple through phiFlowing, which is 1 - phiBiofilm and so sums all 13
% biofilm-block components plus enclosed water (:337-340). They decouple here only
% because phiFlowing is supplied as a KNOWN vector -- at n+1 under the time
% splitting, since the biofilm block carries no dispersion and is advanced first.
% Given that, only alpha differs between components, so the operator is assembled
% once and rescaled. Tridiagonal.
%
% A fully implicit variant (needed only if the biofilm block also goes implicit)
% would couple the flowing block to 14 further unknowns per cell -- but only
% through the SCALAR phi_f, so that coupling is rank one and is best handled by
% carrying phi_f as one auxiliary unknown (10 per cell, not 23). See PLANS.md.
nC = size(g, 1);
if nC < 3, return, end

% cells -> faces: face f carries -/+ S(f-1)/dz at cells f-1, f, for f = 2..nC-1.
% Faces 1, nC, nC+1 carry no dispersive flux, matching the zero rows of
% dispersionFlux in the explicit branch.
f = (2:nC-1)';
Gr = sparse([f; f], [f-1; f], [-S(:)/dz; S(:)/dz], nC+1, nC);

% faces -> cells: left face i with weight poroB(i), right face i+1 with poroB(i+1)
i = (1:nC)';
SelL = sparse(i, i,   poroB(1:nC),   nC, nC+1);
SelR = sparse(i, i+1, poroB(2:nC+1), nC, nC+1);

Base = spdiags((dt/dz)./poroC(:), 0, nC, nC) * (SelL - SelR) * Gr * ...
       spdiags(1./phiFlowing(:), 0, nC, nC);

Id = speye(nC);
for j = 1:size(g, 2)
    g(:, j) = (Id + alpha(j)*Base) \ g(:, j);
end
end
