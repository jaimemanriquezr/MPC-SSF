function plotSolverADrift(opts)
% PLOTSOLVERADRIFT  Free-running Solver A against the authoritative scheme.
%
% The authoritative phi_b is the sum over biofilm components produced by
% Solver B. The free-running phi_b is Solver A advanced on its own, seeded once
% and never re-seeded -- what the component-elimination design would rely on.
%
% Both are plotted because the summary statistics hide the shape: the drift
% SATURATES into a band rather than growing, and the band oscillates on what
% looks like the diel light period. A single "max drift" number shows neither.

arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
    opts.OutDir (1,1) string = ""
end
if opts.Days < 3
    error("plotSolverADrift:inactiveRegime", "Days = %g is below the 3 d minimum.", opts.Days);
end
here = fileparts(mfilename("fullpath")); W = fileparts(here);
addpath(genpath(fullfile(W,"src")));
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

schemes = ["shin" "bailo"]; R = cell(1,2);
for k = 1:2
    t0 = tic;
    R{k} = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
        TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, FrameNumber=72, ...
        ImplicitOsmosis=true, RecordCflBudget=true, CohesionScheme=schemes(k), Quiet=true);
    fprintf("  %-6s %s  %.0f s\n", schemes(k), R{k}.Flag, toc(t0));
end

n0 = f.GridZero; z = f.GridPoints.Centers(1:n0);
set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.5, ...
    "defaultFigureVisible", "off");
cols = lines(2); ls = ["-" "--"];
fig = figure(Position=[70 70 1200 760], Visible="off");
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

D = cell(1,2); phiA = cell(1,2); phiS = cell(1,2); tt = cell(1,2); phiC = cell(1,2);
for k = 1:2
    phiA{k} = R{k}.SimulationData.PhibPar;              % free-running Solver A
    phiS{k} = biofilmFrac(R{k}); phiS{k} = phiS{k}(1:n0, :);   % authoritative sum
    phiC{k} = R{k}.SimulationData.PhibCH;               % re-seeded, one-step
    % PhibCH is written from frame 1; PhibPar only once the free state is seeded.
    % Masking each by its own NaNs leaves arrays of different width, so use the
    % free state's mask throughout -- it is the narrower one.
    ok = ~isnan(phiA{k}(1,:));
    tt{k} = R{k}.Frames.Time(ok).';
    D{k} = max(abs(phiA{k}(:,ok) - phiS{k}(:,ok)), [], 1);
    phiA{k} = phiA{k}(:,ok); phiS{k} = phiS{k}(:,ok); phiC{k} = phiC{k}(:,ok);
end
% Cache, so re-plotting does not cost another full pair of simulations.
cache = struct("z", z, "tt", {tt}, "phiA", {phiA}, "phiS", {phiS}, "phiC", {phiC}, ...
               "D", {D}, "schemes", schemes, "NCells", opts.NCells, "Kappa", opts.Kappa);
save(fullfile(opts.OutDir, sprintf("solverA_drift_n%d.mat", opts.NCells)), "cache", "-v7");
den = max(cellfun(@(p) max(p(:)), phiS));

nexttile; hold on                                        % drift vs time
for k = 1:2, plot(tt{k}, D{k}, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)); end
xlabel("t [d]"); ylabel("max_z |\phi_b^{free} - \phi_b^{sum}|");
% Envelope, because the diel oscillation hides the trend: reading sampled values
% I first called this "saturating", but the PEAKS rise monotonically
% (shin 0.0122 -> 0.0137 -> 0.0149 -> 0.0169 over 3 d).
for k = 1:2
    D = D{k};  li = find(D(2:end-1) > D(1:end-2) & D(2:end-1) > D(3:end)) + 1;  pk = D(li);  % no Signal Toolbox
    if numel(pk) >= 2
        plot(tt{k}(li), pk, "o-", Color=cols(k,:)*0.6, LineWidth=1.0, ...
             MarkerSize=4, DisplayName=schemes(k)+" peaks");
    end
end
title("drift vs time -- peaks CREEP UPWARD, not saturating"); legend(Location="southeast"); grid on; box on

nexttile; hold on                                        % same, relative, log
for k = 1:2, semilogy(tt{k}, 100*D{k}/den, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)); end
set(gca, YScale="log"); xlabel("t [d]"); ylabel("% of max \phi_b");
title("relative drift (log)"); legend(Location="southeast"); grid on; box on

nexttile; hold on                                        % bounds of the free state
for k = 1:2
    plot(tt{k}, min(phiA{k}, [], 1), Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)+" min");
    plot(tt{k}, max(phiA{k}, [], 1), Color=cols(k,:), LineStyle=":", DisplayName=schemes(k)+" max");
end
yline(0, "k-", HandleVisibility="off");
xlabel("t [d]"); ylabel("\phi_b^{free}"); title("bounds of the free state");
legend(Location="east", FontSize=7); grid on; box on

nexttile; hold on                                        % final profiles
for k = 1:2
    plot(phiS{k}(:,end), z, Color=cols(k,:), LineWidth=2.6, DisplayName=schemes(k)+" authoritative");
    plot(phiA{k}(:,end), z, Color=cols(k,:), LineStyle="--", DisplayName=schemes(k)+" free");
end
set(gca, YDir="reverse"); xlabel("\phi_b"); ylabel("z [m]");
title(sprintf("profiles at t = %.0f d", opts.Days)); legend(Location="southeast", FontSize=7);
grid on; box on

nexttile; hold on                                        % pointwise drift profile
for k = 1:2
    plot(phiA{k}(:,end) - phiS{k}(:,end), z, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k));
end
set(gca, YDir="reverse"); xline(0, "k:", HandleVisibility="off");
xlabel("\phi_b^{free} - \phi_b^{sum}"); ylabel("z [m]");
title("where the drift lives"); legend(Location="southeast"); grid on; box on

nexttile                                                 % one-step vs accumulated
oneStep = nan(1,2);
for k = 1:2
    oneStep(k) = max(abs(phiC{k} - phiS{k}), [], "all");
end
b = bar(categorical(schemes, schemes), [oneStep(:), cellfun(@max, D)'], "grouped");
b(1).DisplayName = "one-step (re-seeded)"; b(2).DisplayName = "accumulated (free)";
set(gca, YScale="log"); ylabel("max |difference|");
title("local vs accumulated"); legend(Location="northwest", FontSize=8); grid on; box on

title(tl, sprintf("Free-running Solver A vs the authoritative scheme (N = %d, \\kappa = %g, %.0f d)", ...
    opts.NCells, opts.Kappa, opts.Days));
out = fullfile(opts.OutDir, sprintf("fig8_solverA_drift_n%d.png", opts.NCells));
exportgraphics(fig, out, Resolution=200);
fprintf("wrote %s\n", out);
for k = 1:2
    fprintf("  %-6s one-step %.3e | accumulated max %.3e (%.2f%%) | free range [%.3e, %.3e]\n", ...
        schemes(k), oneStep(k), max(D{k}), 100*max(D{k})/den, min(phiA{k}(:)), max(phiA{k}(:)));
end
close all
end

function p = biofilmFrac(r)
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
p = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], p = p + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   p = p + C{nm,"Enclosed"}{1}/dL; end
end
