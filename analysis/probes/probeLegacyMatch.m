function res = probeLegacyMatch(season, opts)
% PROBELEGACYMATCH  Reproduce the PUBLISHED seasons run inside MPC-SSF.
%
%   probeLegacyMatch("summer")   probeLegacyMatch("winter")
%
% WHY ------------------------------------------------------------------------
% Manuscript Fig. 7 (fig:seasons-results) shows summer > winter. Every MPC-SSF
% arm run so far -- including probeSeason90's "manuscript" variant -- gives
% winter > summer. Before treating that as a regression it has to be established
% that MPC-SSF was ever running the published configuration. It was not.
%
% The published figure was produced by the LEGACY implementation:
%   code/1d/slow-sand-filtration/srun_biofilm.m
%     parameters_file = "SDparameters.mat"          (rho 1117/998, Table B.2/B.4)
%     model.detachment = @(v) 0.14*sqrt(abs(v)/(7.2/0.4))       -> qnom = 18
%     model.cahn_hilliard.zeta_0 = 1e2
%     filter = SDfilter().add_cells(500);  filter.temperature = 293
%     options.total_experiment_length = 90   -> batch_seasons
%   and batch_simulations/batch_seasons.m sets 19 C / 3 C with the same two
%   light curves manuscriptExperiments uses.
%
% Divergences of probeSeason90("*","manuscript") from that, largest first:
%   1. It adds a SIXTH reaction -- PG-excess phototroph respiration at 0.55/d --
%      which the published model does not contain (SDparameters.mat has exactly
%      five: Het_Growth, Pho_Growth, Het_Death, Pho_Death, Hydrolysis).
%   2. Passing PhototrophRespiration>0 makes modelLund zero MinimumLightFactor
%      (modelLund.m:119). The published run has dark_respiration = 0.01.
%   3. I_opt: the legacy SI model carries 1.08. Manuscript Table B.4 prints
%      1.814e-2, which is what modelLund implements -- and which puts the surface
%      at ~44x optimal, deep in Steele photoinhibition. The NormalizedLight
%      convention (I_opt = 1.0) lands back near the legacy value by a different
%      route. THE TABLE VALUE APPEARS TO BE A TRANSCRIPTION ERROR.
%   4. Osmosis rate: legacy tau = 1e-5, modelLund tail OsmosisRate = 1e-7 (100x).
%      Manuscript Table B.2 prints 1e5, matching neither.
%   5. 500 cells vs 100.
%
% This probe reconstructs (1)-(4) and leaves resolution as an option, so the
% seasonal ordering can be attributed. Items 1 and 2 are the prime suspects: the
% kinetic audit found the PG-excess form to be "a net biomass source fed by an
% infinite untracked pool, so darkness manufactures biomass", and measured winter
% carrying 3.6x summer's phototroph mass under it.
%
% See also PROBESEASON90, MANUSCRIPTEXPERIMENTS.

arguments
    season (1,1) string {mustBeMember(season, ["summer","winter"])}
    opts.NCells (1,1) double = 100      % 500 = published; 100 for turnaround
    opts.Tsim   (1,1) double = 90.0
    opts.MaxDt  (1,1) double = 3e-6
    opts.Variant (1,1) string {mustBeMember(opts.Variant, ...
        ["published","plusResp","publishedIopt1814"])} = "published"
    % Winter photoperiod. "manuscript" is the curve batch_seasons.m and
    % manuscriptExperiments both use; "realistic" corrects it -- see below.
    opts.WinterLight (1,1) string {mustBeMember(opts.WinterLight, ...
        ["manuscript","realistic"])} = "manuscript"
    opts.Force  (1,1) logical = false
end

here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

tag = sprintf("legmatch_%s_%s_n%d", season, opts.Variant, opts.NCells);
if opts.WinterLight == "realistic", tag = tag + "_wreal"; end
out = fullfile(S, tag + ".mat");
if isfile(out) && ~opts.Force
    fprintf("LEGMATCH skip %s\n", tag); res = load(out).res; return
end

%% ---- model: the published five-reaction Lund ecology -------------------
% modelLund, NOT pathogenModel: SDparameters.mat has no PAT reactions, and
% pathogenModel would add MarkerGrowth/Inactivation/Bacterivory plus its own
% Zeta0 and detachment overrides.
switch opts.Variant
    case "published"
        m = modelLund(NormalizedLight=false);          % 5 reactions, dark floor 0.01
        iOpt = 1.08;                                   % legacy SDparameters.mat
    case "publishedIopt1814"
        m = modelLund(NormalizedLight=false);
        iOpt = 1.814e-2;                               % manuscript Table B.4 as printed
    case "plusResp"
        % Published + the post-publication respiration, to isolate its effect.
        m = modelLund(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=false);
        iOpt = 1.08;
end

rx = m.Reactions;
for k = 1:numel(rx)
    if rx(k).OptimalLightFactor > 0, rx(k).OptimalLightFactor = iOpt; end
end

% Published cohesion/transport overrides from srun_biofilm.m and SDparameters.mat.
m = Model(m.Components, rx, ...
    Kappa = 1e-7, ...                                  % SDparameters + srun
    Zeta0 = 1e2, ...                                   % srun_biofilm override
    Zeta1 = 1e-2, ...
    DetachmentFunction = @(v) 0.14*sqrt(abs(v)/18), ...% srun: 7.2/0.4 = 18
    WaterDensity = m.WaterDensity, ...
    BiofilmPorosity = 0.99, ...
    OsmosisRate = 1e-5);                               % legacy tau

%% ---- filter ------------------------------------------------------------
if season == "summer"
    f = SandFilter(Temperature=19);
    light = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
else
    f = SandFilter(Temperature=3);
    % The manuscript winter curve gives 13.54 h of daylight and a peak of 0.60,
    % i.e. 0.75x summer's peak and 0.63x its daily integral. At Lund (55.7 N),
    % where the "Scandinavian winter" scenario is set, 21 December has 6 h 55 m
    % of daylight and a noon solar elevation of 10.8 deg against 57.8 deg on
    % 21 June -- a true daylight ratio of 0.395 and a true peak-irradiance ratio
    % of sin(10.8)/sin(57.8) = 0.221. So the manuscript curve understates the
    % seasonal light contrast by roughly sixfold in integrated terms.
    %
    % Offset 0.81 fixes it and is self-consistent on BOTH measures
    % independently: 6.89 h of daylight (actual 6.92) and a peak ratio of 0.237
    % (physical 0.221). Summer needs no correction -- 16.92 h against an actual
    % 17.5 h.
    %
    % This matters for the seasonal ordering: with winter light barely reduced,
    % the 19 -> 3 C temperature drop dominates, and since standing stock is
    % supply / theta-suppressed decay, slower winter decay wins.
    if opts.WinterLight == "realistic"
        light = @(t) max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.81, 0);
    else
        light = @(t) max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
    end
end
f = f.addGridPoints(opts.NCells);
f.LightIrradiation = light;

% SDparameters.mat inflow_concentrations, 8 components (no PAT).
infl = dictionary( ...
    ["HET","PHO","POM","PAT","O2","IC","NH4","HPO4","DOM"], ...
    [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]);

% Report the light curve as actually instantiated, not as labelled: integrate it
% so the printed daylight hours can be checked against the intent.
tt = linspace(0,1,20001); Ic = arrayfun(f.LightIrradiation, tt);
fprintf("LEGMATCH start %s variant=%s wlight=%s Iopt=%g nRx=%d n=%d t=%gd T=%g | peak=%.2f daylight=%.2fh\n", ...
    tag, opts.Variant, opts.WinterLight, iOpt, numel(m.Reactions), opts.NCells, opts.Tsim, ...
    f.Temperature, max(Ic), 24*mean(Ic > 0));
tRun = tic;

r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Tsim, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=max(round(2*opts.Tsim), 12), ImplicitOsmosis=true, Quiet=true, ...
    RecordLimitation=true);

%% ---- the plotted quantity ----------------------------------------------
C = r.Frames.Concentrations;
z = f.GridPoints.Centers(:); dz = f.GridSize;
poros = computePorosity(f, z);
dL = mean([m.Liquids.Density]); dP = mean([m.Particles.Density]);
phiB = C{"Water","Enclosed"}{1}/dL;
for nm = [m.Particles.Name], phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [m.Liquids.Name],   phiB = phiB + C{nm,"Enclosed"}{1}/dL; end

res = struct("tag", tag, "season", season, "variant", opts.Variant, ...
    "winterLight", opts.WinterLight, "temperatureC", f.Temperature, ...
    "lightPeak", max(Ic), "daylightHours", 24*mean(Ic > 0), ...
    "iOpt", iOpt, "nReactions", numel(m.Reactions), "ncells", opts.NCells, ...
    "flag", string(r.Flag), "tFinal", r.TimeFinal, "wallMin", toc(tRun)/60, ...
    "phib_max", max(phiB(:,end)), "phib_sup", max(phiB(z<0,end)), ...
    "phib_int_eps", sum(poros.*phiB(:,end))*dz, ...
    "pho_mass_eps", sum(poros.*(C{"PHO","Matrix"}{1}(:,end) + C{"PHO","Enclosed"}{1}(:,end)))*dz, ...
    "o2_out", mean(C{"O2","Flowing"}{1}(end,max(1,end-47):end))*1000);

profile = [z, poros, phiB(:,end)];
save(out, "res", "profile", "-v7");
writematrix(profile, fullfile(here, "data", "legmatch_profile_" + tag + ".csv"));

fprintf("LEGMATCH done %s flag=%s wall=%.1fmin | phib_max=%.4f phib_sup=%.4f phib_int=%.4f PHO=%.4f O2out=%.2f\n", ...
    tag, res.flag, res.wallMin, res.phib_max, res.phib_sup, res.phib_int_eps, res.pho_mass_eps, res.o2_out);
end
