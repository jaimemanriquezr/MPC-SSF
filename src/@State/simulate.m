function results = simulate(obj, parameters)
arguments
    obj  State
    parameters.InflowConcentrations = []
    parameters.SimulationTime (1,1) {mustBeNumeric} = 1.0;
    parameters.TimeStep (1,1) = 1E-5;
    parameters.FrameNumber (1,1) {mustBeNumeric} = 200;
    parameters.CloggingFraction (1,1) {mustBeNumeric} = .99;
    parameters.IsUpwinded (1,1) = false;
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
else
    inflowConcentrations = parameters.InflowConcentrations;
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

attachmentRates = [model.Particles.AttachmentRate];
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
dt = timeStep;
results.SimulationData.TimeStep = dt;

t = timeStart;
muRates = reshape(model.Reactions.computeRate(temperature), 1, []);
lightOptimal = max([model.Reactions.OptimalLightFactor]);
attenuationParticles = [model.Particles.Attenuation];

minimumLight = [model.Reactions.MinimumLightFactor];
lightDependency = [model.Reactions.IsLightDependent];
minimumLight = minimumLight(lightDependency);
lightFactor = ones(length(filter.GridPoints.Centers), length(model.Reactions));

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
    lightFactor(:, lightDependency) = (minimumLight + lightEffective + abs(minimumLight + lightEffective))/2;

    ecoRxBiofilm = evaluateReactions(localBiofilm, listK, phiBiofilm, muRates, lightFactor, listOrder);
    ecoRxEnclosed = evaluateReactions(localEnclosed, listK, phiEnclosed, muRates, lightFactor, listOrder);
    ecoRxFlowing = evaluateReactions(localFlowing, listK, phiFlowing, muRates, lightFactor, listOrder);

    ecoRxM = ecoRxBiofilm*sigmaParticles';
    ecoRxPe = ecoRxEnclosed*sigmaParticles';
    ecoRxLe = (ecoRxBiofilm + ecoRxEnclosed)*sigmaLiquids';
    ecoRxPf = ecoRxFlowing*sigmaParticles';
    ecoRxLf = ecoRxFlowing*sigmaLiquids';

    attE = attachmentEnclosedFactor.*globalEnclosedP.*attachmentRates;
    attF = attachmentFlowingFactor.*globalFlowingP.*attachmentRates;

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

    %=================== IV. SOLVER A: compute biofilm velocity ======================%
    rhsBiofilmVolume = sum(reactionsMatrix,2)/densityP ...
                        + sum(reactionsEnclosedParticles,2)/densityP ...
                        + sum(reactionsEnclosedLiquids,2)/densityL ...
                        + rhsEnclosedWater;
    rhsBiofilmVolume = rhsBiofilmVolume(1:n0);

    u = phiBiofilm(1:n0);
    uB = .5*(u(2:end) + u(1:end-1));

    lambda = zeta_0*mobility(uB);
    lhsCH = CH0 - dt*(S + sparse(D.Rows, D.Columns, D.Values.*(lambda).', 2*n0, 2*n0));
    rhsCH = [u + dt*rhsBiofilmVolume(1:n0); dpsi_fun(u)];

    % solving linear system
    xCH = lhsCH \ rhsCH;
    % uCH = xCH(1:n0);
    muCH = xCH(n0+1:end);

    % computing vb
    velBiofilm(1:n0-1) = volumeAvgVelocity(2:n0) - zeta_0*(1 - uB).*diff(muCH)/dz;
    %================================================================%
    %====================== V. SOLVER B: compute concentrations =================================%
    phiBiofilmBoundaries = .5*(phiBiofilm(2:end) + phiBiofilm(1:end-1));
    velFlowing = [volumeAvgVelocity(1);
        (volumeAvgVelocity(2:end-1) - velBiofilm.*phiBiofilmBoundaries)./(1 - phiBiofilmBoundaries);
        volumeAvgVelocity(end)];

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

    fluxFlowing = convectionFlux - dispersionFlux;
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
    globalBiofilm = globalBiofilm + (dt/dz)*(fluxBiofilmIn - fluxBiofilmOut)./porosityCenters + dt*rhsBiofilm;
    globalFlowing = globalFlowing + (dt/dz)*(fluxFlowingIn - fluxFlowingOut)./porosityCenters + dt*rhsFlowing;
    phiW = phiW + (dt/dz)*(fluxWaterIn - fluxWaterOut)./porosityCenters + dt*rhsEnclosedWater;

    %========== CHECK IF CONCENTRATIONS ARE NEGATIVE ============%
    problemCellBiofilm = mod(find(globalBiofilm(:) < 0 | isnan(globalBiofilm(:))), size(globalBiofilm, 1));
    if ~isempty(problemCellBiofilm)
        fprintf('Unphysical concentration in biofilm. \nTIME = %e\n',t)
        fprintf("CELL = %i\n", problemCellBiofilm(1) - n0)
        fprintf("HEIGHT = %i\n", depthCenters(problemCellBiofilm(1)));
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
        counter = counter + 1;
        fprintf("t = %.4e\n", t);
    end
    %==========================================================%
end
disp("Simulation ended.")
disp("Saving results...")

%====================== VII. OUTPUTS ===================%
results.Frames.time = timeFrames;

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
    'VariableNames', {'matrix',  'enclosed',  'flowing'},  ...
    'rowNames', [results.Model.Particles.Name,  results.Model.Liquids.Name,  'WATER'],  ...
    'DimensionNames', {'Component', 'Volume'});
results.Frames.Velocity.Biofilm = velFramesBiofilm;
results.Frames.Velocity.Flowing = velFramesFlowing;

results.TimeFinal = t;
disp("Results saved.")
end

function rx = evaluateReactions(local, K, phi, mu, I, orders)
    % monod = permute(min(min((local + realmin)./(K + local + realmin), [], 2), 1), [1 3 2]);
    monod = ones(size(phi, 1), size(mu, 2));
    for i = 1:size(monod, 2)
        ii = K{1, i};
        if any(ii)
            k = K{2, i};
            mon_term = (local(:, ii) + realmin)./(k + local(:, ii) + realmin);
            monod(:, i) = min(mon_term, [], 2);
        end
    end
    %product = permute(prod(local.^orders, 2, "omitnan"), [1 3 2]);
    product = local(:, orders);
    rx = phi.*mu.*I.*monod.*product;
end

