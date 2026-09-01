function out = plotBiologicalActivity(opts)
% PLOTBIOLOGICALACTIVITY  Depth-resolved and column-integrated reaction rates.
%
% "Biological activity" measured from the model's OWN reaction terms rather than
% from a standing stock: for every reaction j the volumetric extent rate
%
%   r_j(z) = phi * mu_j(T) * I_j(z) * min_k Monod_k * X_driver_local        [kg/(m^3 d)]
%
% summed over the three phases the solver evaluates it in (biofilm matrix,
% enclosed, flowing). The assembly below is copied from src/@State/simulate.m --
% see the cited line ranges -- because a re-derivation is exactly how a figure
% that looks right stops meaning what it says. simulate.m:501-505 anticipates
% this: the Monod and limitation factors are pure functions of the local
% concentrations, so evaluateReactions can be re-run on a stored frame.
%
% Sources mirrored (src/@State/simulate.m):
%   :70-84    listK / listOrder assembly from HalfSaturationConstants + quotients
%   :103-137  porosity, dz, eta_water/eta_sand, inflow velocity
%   :374-401  muRates, lightOptimal, light-mode flags, attenuationParticles,
%             lightFactor init, per-reaction phase efficiencies
%   :415-462  global -> phi -> local concentrations in each phase
%   :472-497  eta, light attenuation, Steele curve, lightFactor per mode
%   :506-514  ecoRx{Biofilm,Flowing,Enclosed} and the stoichiometric projections
%   :966,984  the Solver B update: the reaction source enters as dt*rhs on a
%             quantity whose conserved form is eps*g, hence the eps weighting of
%             every column integral here
%   :1115     local function evaluateReactions (reproduced as evalRx below)
%
% QUASI-STEADY DIEL AVERAGE. snap.Time is a whole number of days, where the chain
% light curve is identically zero, so evaluating at the snapshot instant returns
% zero phototroph growth. The state is frozen and the rates are evaluated over 48
% times spanning one day; what is plotted is the DAILY MEAN. Legitimate only
% because the state is quasi-steady by >= 30 d (E9: 0.001 %/10 d profile change).
%
% Reactions 6-8 (MarkerGrowth, Inactivation, Bacterivory) are driven by PAT,
% which is identically zero in the E4/E9 influent, so they are identically zero
% and are omitted from the figure.

arguments
    opts.LitTag  (1,1) string = "fld2x_lit"
    opts.DarkTag (1,1) string = "fld2x_dark"
    opts.Legs    (1,:) double = 1:11
    opts.NTimes  (1,1) double = 48        % diel sampling points over one day
    opts.OutDir  (1,1) string = ""
    opts.Split   (1,1) logical = false    % two files instead of one 2-panel figure
end
here = fileparts(mfilename("fullpath")); W = fileparts(here);
addpath(genpath(fullfile(W, "src"))); addpath(here); addpath(fullfile(here, "probes"));
dataDir = fullfile(here, "probes", "data", "chain");
if opts.OutDir == ""
    opts.OutDir = fullfile(here, "results", "figures", "recreation");
end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

rxShow  = ["Heterotroph growth", "Phototroph growth", "Heterotroph death", ...
           "Phototroph death", "Hydrolysis"];
rxLabel = ["HET growth", "PHO growth", "HET death", "PHO death", "Hydrolysis"];

arms = struct("tag", {opts.LitTag, opts.DarkTag}, "light", {1, 0}, ...
              "name", {"lit", "dark"});

out = struct();
for a = 1:2
    [f, m] = buildWorkingSet(arms(a).light);
    nL = numel(opts.Legs);
    tLeg = nan(1, nL);
    prof = [];  colInt = [];
    for q = 1:nL
        fn = fullfile(dataDir, sprintf("chain_%s_leg%d.mat", arms(a).tag, opts.Legs(q)));
        if ~isfile(fn), error("missing %s", fn); end
        D = load(fn, "snap");
        R = dielMeanRates(f, m, D.snap, opts.NTimes);
        tLeg(q) = D.snap.Time;
        if isempty(prof)
            prof   = nan(numel(R.z), numel(R.names), nL);
            colInt = nan(numel(R.names), nL);
            zAxis  = R.z;  epsC = R.eps;  dz = R.dz;
            rxNames = R.names;
        end
        prof(:, :, q) = R.rateTotal;                       % kg/(m^3 d), all phases
        colInt(:, q)  = (epsC(:).' * R.rateTotal).' * dz;  % kg/(m^2 d)
        if opts.Legs(q) == max(opts.Legs)
            out.(arms(a).name).check = R;                  % keep the last frame whole
        end
    end
    out.(arms(a).name).t      = tLeg;
    out.(arms(a).name).prof   = prof;
    out.(arms(a).name).colInt = colInt;
    out.(arms(a).name).z      = zAxis;
    out.(arms(a).name).eps    = epsC;
    out.(arms(a).name).dz     = dz;
    out.(arms(a).name).names  = rxNames;
end

% ------------------------- VALIDATION ------------------------------------- %
% Reaction-side O2 sink at the final leg, lit, against the advective deficit.
% The solver advances eps*g (simulate.m:966,984: fluxes carry porosityBoundaries
% and the divergence is divided by porosityCenters), so the column integral of a
% source is int eps*rhs dz. The inflow flux is q*C_in with q = filter.InflowVelocity
% (porosityBoundaries(1) = 1 cancels volumeAvgVelocity(1) = q/1, simulate.m:137).
Rc = out.lit.check;
o2Sink = -sum(Rc.eps(:).*Rc.o2Net)*Rc.dz;                  % kg O2/(m^2 d), > 0 = sink
qFilt  = Rc.q;
o2Adv  = qFilt*(Rc.o2In - Rc.o2Out);
relErr = (o2Sink - o2Adv)/o2Adv;
fprintf("\n== O2 mass balance, %s, t = %g d ==\n", opts.LitTag, out.lit.t(end));
fprintf("  reaction sink   %8.3f g O2/(m^2 d)\n", 1e3*o2Sink);
fprintf("  advective q(Cin-Cout) %8.3f g O2/(m^2 d)  (q = %g m/d, Cin = %.3e, Cout = %.3e kg/m^3)\n", ...
    1e3*o2Adv, qFilt, Rc.o2In, Rc.o2Out);
fprintf("  relative difference  %+6.1f %%\n", 100*relErr);
out.validation = struct("O2SinkGm2d", 1e3*o2Sink, "O2AdvGm2d", 1e3*o2Adv, ...
    "RelErr", relErr, "q", qFilt, "Cin", Rc.o2In, "Cout", Rc.o2Out);

jPho = find(out.dark.names == "Phototroph growth");
darkPho = max(abs(out.dark.prof(:, jPho, :)), [], "all");
fprintf("== dark arm phototroph growth: max |rate| = %.3e kg/(m^3 d) (must be 0)\n", darkPho);
jPhoL = find(out.lit.names == "Phototroph growth");
pL = out.lit.prof(:, jPhoL, end);  zL = out.lit.z;
lit1 = zL(pL > 1e-3*max(pL));
if isempty(lit1), lit1 = NaN(1,2); end
fprintf("== lit arm phototroph growth confined to z in [%.4f, %.4f] m\n", min(lit1), max(lit1));
out.validation.DarkPhoMax = darkPho;
out.validation.LitPhoZone = [min(lit1) max(lit1)];

% Phase split of the column-integrated activity. EfficiencyFlowing = 0 for every
% reaction except Bacterivory (which is PAT-driven and so identically zero), so
% the flowing phase is inert here BY CONSTRUCTION -- see simulate.m:394-396.
wB = sum(Rc.eps(:).*Rc.rateBiofilm, 1)*Rc.dz;
wE = sum(Rc.eps(:).*Rc.rateEnclosed, 1)*Rc.dz;
wF = sum(Rc.eps(:).*Rc.rateFlowing, 1)*Rc.dz;
tot = sum(wB + wE + wF);
fprintf("== phase split of total column activity: biofilm %.6f, enclosed %.3e, flowing %.3e\n", ...
    sum(wB)/tot, sum(wE)/tot, sum(wF)/tot);
out.validation.PhaseShare = [sum(wB) sum(wE) sum(wF)]/tot;

% -------------------------- FIGURE ---------------------------------------- %
set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.6, ...
    "defaultFigureVisible", "off");
jShow = arrayfun(@(s) find(out.lit.names == s), rxShow);
cmap  = lines(numel(jShow));

tagStr = strrep(opts.LitTag, "_", "\_") + " / " + strrep(opts.DarkTag, "_", "\_");
if opts.Split
    figA = figure(Position=[60 60 720 660], Visible="off");
    axA = axes(figA); %#ok<LAXES>
    hA = drawProfiles(axA, out, jShow, rxLabel, cmap);
    legend(axA, hA, Location="southwest", FontSize=8, Box="on");
    saveBoth(figA, fullfile(opts.OutDir, "rec_BiologicalActivity_profiles"));
    figB = figure(Position=[60 60 720 660], Visible="off");
    axB = axes(figB); %#ok<LAXES>
    hB = drawSeries(axB, out, jShow, rxLabel, cmap);
    legend(axB, hB, Location="southeast", FontSize=8, Box="on", NumColumns=2);
    saveBoth(figB, fullfile(opts.OutDir, "rec_BiologicalActivity_timeseries"));
else
    fig = figure(Position=[60 60 1250 700], Visible="off");
    tl = tiledlayout(fig, 1, 2, TileSpacing="compact", Padding="compact");
    hA = drawProfiles(nexttile(tl), out, jShow, rxLabel, cmap);
    drawSeries(nexttile(tl), out, jShow, rxLabel, cmap);
    lgd = legend(hA, FontSize=9, Box="on", NumColumns=5, Orientation="horizontal");
    lgd.Layout.Tile = "south";
    title(tl, "Biological activity from the reaction terms (" + tagStr + ...
        ", working set, T = 19 \circC)", Interpreter="tex");
    subtitle(tl, sprintf("daily mean over %d times with the state frozen at the snapshot " + ...
        "(quasi-steady); solid = lit, dashed = dark; PAT-driven reactions omitted (PAT \\equiv 0)", ...
        opts.NTimes), Interpreter="tex", FontSize=9);
    saveBoth(fig, fullfile(opts.OutDir, "rec_BiologicalActivity"));
end

% ------------------------- SUMMARY TABLE ---------------------------------- %
fprintf("\n%-20s %14s %14s %9s\n", "reaction", "lit [g/m2/d]", "dark [g/m2/d]", "lit/dark");
for k = 1:numel(jShow)
    L = 1e3*out.lit.colInt(jShow(k), end);
    Dk = 1e3*out.dark.colInt(jShow(k), end);
    fprintf("%-20s %14.4f %14.4f %9.3f\n", rxLabel(k), L, Dk, L/Dk);
end
end

% ========================================================================== %
function [f, m] = buildWorkingSet(lightScale)
% E4/E9 working set, copied from analysis/probes/probeChain.m:74-120 with the
% E9 option values (EXPERIMENTS.md:619): Zeta0 = 1, Zeta1 = 0.27, KDOM = 3e-4,
% TransferScale = 10, KHPO4 = 1e-6, DetachForm = "linear", EtaSand = 1500,
% Delta = 5e-3, Temperature = 19, LightForm = "chain", MuPHO/MuHET/DPHO preset.
lightFn = @(t) lightScale*0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50);
f = SandFilter(Temperature=19, LightAttenuationCoeffSand=1500, SandRoughness=5e-3, ...
    LightIrradiation=lightFn);
f = f.addGridPoints(500);
mp = pathogenModel(PhototrophRespiration=0, RespirationForm="reichert", NormalizedLight=true);
rx = mp.Reactions;
rx([rx.Name] == "Phototroph growth").MinimumLightFactor = 0.0;
iH = find([rx.Name] == "Heterotroph growth");
H = rx(iH).HalfSaturationConstants; H("DOM") = 3e-4; rx(iH).HalfSaturationConstants = H;
iP = find([rx.Name] == "Phototroph growth");
P = rx(iP).HalfSaturationConstants; P("HPO4") = 1e-6; rx(iP).HalfSaturationConstants = P;
mp.Reactions = rx;
cs = mp.Components;                                   % TransferScale = 10
for q = 1:numel(cs)
    if isa(cs(q), "Liquid"), cs(q).TransportRate = 10*cs(q).TransportRate; end
end
mp.Components = cs;
m = Model(mp.Components, mp.Reactions, Kappa=1e-6, Zeta0=1, Zeta1=0.27, ...
    DetachmentFunction=@(v) 0.14*abs(v)/18, ...   % probeChain.m:325, DetachForm "linear"
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, ...
    OsmosisRate=mp.OsmosisRate);
end

% ========================================================================== %
function R = dielMeanRates(filter, model, snap, nT)
% Daily-mean reaction extent rates for a frozen state. Every block is the
% corresponding block of simulate.m; the line numbers are in the header.

% --- listK / listOrder (simulate.m:70-84) --------------------------------- %
[quotientK, quotientDenIdx, quotientNumIdx] = model.Reactions.lookupQuotients(model.Components);
halfSaturationK = model.HalfSaturationConstants;
listK = permute([halfSaturationK; quotientK], [3 1 2]);
listOrder = zeros(size(listK));
listOrder(1, 1:length(model.Components), :) = model.Order;
newlist = cell(2, size(listK, 3));
for i = 1:size(listK, 3)
    newlist{1, i} = ~isnan(listK(:, :, i));
    if any(newlist{1, i}), newlist{2, i} = listK(:, newlist{1, i}, i); end
end
listK = newlist;
[listOrder, ~] = ind2sub(size(listOrder, [2 3]), find(listOrder));

% --- filter / model constants (simulate.m:103-137, 141-159, 374-401) ------ %
depthCenters = filter.GridPoints.Centers;
porosityCenters = computePorosity(filter, depthCenters);
dz = filter.GridSize;
etaWater = filter.LightAttenuationEtaWater;
etaSand  = filter.LightAttenuationEtaSand;
kP = length(model.Particles);  kL = length(model.Liquids);
densityP = mean([model.Particles.Density]);
densityL = mean([model.Liquids.Density]);
sigmaLiquids = model.StoichiometricMatrixLiquids;
muRates = reshape(model.computeReactionRates(filter.Temperature), 1, []);
lightOptimal = max([model.Reactions.OptimalLightFactor]);
if isempty(lightOptimal) || lightOptimal <= 0, lightOptimal = 1.0; end
lightInhibition = [model.Reactions.LightInhibition];
inhibitionDependency = lightInhibition > 0;
complementDependency = [model.Reactions.IsLightComplement];
attenuationParticles = [model.Particles.Attenuation];
minimumLight = [model.Reactions.MinimumLightFactor];
lightDependency = [model.Reactions.IsLightDependent];
minimumLight = minimumLight(lightDependency);
efficiencyBiofilm = reshape([model.Reactions.EfficiencyBiofilm], 1, []);
efficiencyFlowing = reshape([model.Reactions.EfficiencyFlowing], 1, []);

% --- state -> local concentrations (simulate.m:415-462) ------------------- %
globalMatrix    = snap.Matrix;
globalEnclosedP = snap.EnclosedParticles;
globalEnclosedL = snap.EnclosedLiquids;
globalFlowingP  = snap.FlowingParticles;
globalFlowingL  = snap.FlowingLiquids;
phiW = snap.EnclosedWaterVolume;

phiMatrix   = sum(globalMatrix, 2)/densityP;
phiEnclosed = phiW + sum(globalEnclosedP, 2)/densityP + sum(globalEnclosedL, 2)/densityL;
phiBiofilm  = phiMatrix + phiEnclosed;
phiFlowing  = 1 - phiBiofilm;

localBiofilmXS = [globalMatrix./(phiBiofilm + realmin), globalEnclosedL./(phiBiofilm + realmin)];
localBiofilm = [localBiofilmXS, ...
    localBiofilmXS(:, quotientNumIdx)./(localBiofilmXS(:, quotientDenIdx) + realmin)];
localEnclosedXS = [globalEnclosedP./(phiEnclosed + realmin), globalEnclosedL./(phiEnclosed + realmin)];
localEnclosed = [localEnclosedXS, ...
    localEnclosedXS(:, quotientNumIdx)./(localEnclosedXS(:, quotientDenIdx) + realmin)];
localFlowingXS = [globalFlowingP./phiFlowing, globalFlowingL./phiFlowing];
localFlowing = [localFlowingXS, ...
    localFlowingXS(:, quotientNumIdx)./(localFlowingXS(:, quotientDenIdx) + realmin)];

% --- attenuation (simulate.m:472-474): state-dependent, so computed once --- %
globalParticles = globalMatrix + globalEnclosedP + globalFlowingP;
etaParticles = cumsum(sum(attenuationParticles .* globalParticles, 2))*dz;
eta = etaWater + etaSand + etaParticles;

% --- diel sweep with the state frozen ------------------------------------- %
nRx = length(model.Reactions);
nC  = numel(depthCenters);
tGrid = snap.Time + (0:nT-1)/nT;
accB = zeros(nC, nRx); accE = accB; accF = accB;
noonB = []; noonE = []; noonF = [];
for it = 1:nT
    lightFactor = ones(nC, nRx);
    light = filter.LightIrradiation(tGrid(it));
    lightAttenuated = light * exp(-eta) / lightOptimal;          % :477
    lightEffective  = lightAttenuated.*exp(1 - lightAttenuated); % :478
    if any(lightDependency)
        lightFactor(:, lightDependency) = max(minimumLight, lightEffective);   % :485
    end
    for jInh = find(inhibitionDependency)                                      % :489-492
        K = lightInhibition(jInh);
        lightFactor(:, jInh) = K./(K + lightAttenuated);
    end
    for jCmp = find(complementDependency)                                      % :495-497
        lightFactor(:, jCmp) = 1 - lightEffective;
    end
    rxB = efficiencyBiofilm.*evalRx(localBiofilm,  listK, phiBiofilm,  muRates, lightFactor, listOrder);
    rxF = efficiencyFlowing.*evalRx(localFlowing,  listK, phiFlowing,  muRates, lightFactor, listOrder);
    rxE = efficiencyBiofilm.*evalRx(localEnclosed, listK, phiEnclosed, muRates, lightFactor, listOrder);
    accB = accB + rxB/nT;  accE = accE + rxE/nT;  accF = accF + rxF/nT;
    if abs(mod(tGrid(it), 1) - 0.5) < 0.5/nT
        noonB = rxB; noonE = rxE; noonF = rxF;
    end
end

R.z = depthCenters;  R.eps = porosityCenters;  R.dz = dz;
R.names = [model.Reactions.Name];
R.rateBiofilm = accB;  R.rateEnclosed = accE;  R.rateFlowing = accF;
% The three phase rates enter the mass balance additively: the matrix block gets
% rxB*sigmaP', the enclosed block rxE*sigmaP', the liquids (rxB + rxE)*sigmaL'
% and the flowing block rxF*sigma' (simulate.m:510-514). So the TOTAL extent of
% reaction j per unit pore volume is the plain sum of the three.
R.rateTotal = accB + accE + accF;
if ~isempty(noonB), R.rateNoon = noonB + noonE + noonF; else, R.rateNoon = nan(nC, nRx); end
% Net O2 source (kg O2/(m^3 pore d)); negative where O2 is consumed.
R.o2Net = R.rateTotal*sigmaLiquids(1, :).';
R.q     = filter.InflowVelocity;
R.o2In  = 9.10e-3;                                    % E4/E9 influent, EXPERIMENTS.md:619
R.o2Out = localFlowingXS(end, kP + 1);                % local flowing O2 in the last cell
R.kP = kP; R.kL = kL;
end

% ========================================================================== %
function rx = evalRx(local, K, phi, mu, I, orders)
% Verbatim from the local function evaluateReactions in src/@State/simulate.m:1115
% (limitation-argmin output dropped, which does not touch rx).
monod = ones(size(phi, 1), size(mu, 2));
for i = 1:size(monod, 2)
    ii = K{1, i};
    if any(ii)
        k = K{2, i};
        monod(:, i) = min((local(:, ii) + realmin)./(k + local(:, ii) + realmin), [], 2);
    end
end
rx = phi.*mu.*I.*monod.*local(:, orders);
end

% ========================================================================== %
function h = drawProfiles(ax, out, jShow, labels, cmap)
hold(ax, "on");
h = gobjects(0);
zi = out.lit.z >= -0.02 & out.lit.z <= 0.30;
vmax = 0;
for k = 1:numel(jShow)
    v = 1e3*out.lit.prof(:, jShow(k), end);
    h(end+1) = plot(ax, v, out.lit.z, "-", Color=cmap(k,:), ...
        DisplayName=labels(k)); %#ok<AGROW>
    vmax = max(vmax, max(v(zi)));
end
for k = 1:numel(jShow)
    d = 1e3*out.dark.prof(:, jShow(k), end);
    if max(abs(d)) == 0, continue, end
    plot(ax, d, out.dark.z, "--", Color=cmap(k,:), HandleVisibility="off");
    vmax = max(vmax, max(d(zi)));
end
set(ax, YDir="reverse", XScale="log");
yline(ax, 0, "k:", HandleVisibility="off");
ylim(ax, [-0.02 0.30]);
xlim(ax, [1e-3, 10^ceil(log10(1.2*vmax))]);
xlabel(ax, "daily-mean rate [g m^{-3} d^{-1}]");
ylabel(ax, "depth z [m]");
title(ax, sprintf("(a) rate profiles at t = %g d (solid lit, dashed dark)", out.lit.t(end)));
grid(ax, "on"); box(ax, "on");
end

% ========================================================================== %
function h = drawSeries(ax, out, jShow, labels, cmap)
hold(ax, "on");
h = gobjects(0);
for k = 1:numel(jShow)
    h(end+1) = plot(ax, out.lit.t, 1e3*out.lit.colInt(jShow(k), :), "-o", ...
        Color=cmap(k,:), MarkerSize=4, MarkerFaceColor=cmap(k,:), ...
        DisplayName=labels(k)); %#ok<AGROW>
end
for k = 1:numel(jShow)
    d = 1e3*out.dark.colInt(jShow(k), :);
    if max(abs(d)) == 0, continue, end
    plot(ax, out.dark.t, d, "--s", Color=cmap(k,:), MarkerSize=4, ...
        HandleVisibility="off");
end
set(ax, YScale="log");
xlim(ax, [0 110]);
xlabel(ax, "t [d]");
ylabel(ax, "column-integrated daily-mean rate [g m^{-2} d^{-1}]");
title(ax, "(b) ripening of column activity (solid lit, dashed dark)");
grid(ax, "on"); box(ax, "on");
end

% ========================================================================== %
function saveBoth(fig, stem)
exportgraphics(fig, stem + ".png", Resolution=200);
exportgraphics(fig, stem + ".pdf", ContentType="vector");
fprintf("wrote %s.{png,pdf}\n", stem);
end
