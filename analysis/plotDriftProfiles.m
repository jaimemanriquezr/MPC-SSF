function plotDriftProfiles(opts)
% PLOTDRIFTPROFILES  phi_b(z) at successive times for the free-running Solver A
% against the authoritative Solver B sum. Uses the cached frames written by
% plotSolverADrift, so it costs seconds rather than re-simulating.

arguments
    opts.Cache (1,1) string = ""
    opts.OutDir (1,1) string = ""
    opts.Tag (1,1) string = ""
end
here = fileparts(mfilename("fullpath"));
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if opts.Cache == "", opts.Cache = fullfile(opts.OutDir, "solverA_drift_n100.mat"); end
c = load(opts.Cache).cache;
z = c.z; schemes = c.schemes;
set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.4, ...
    "defaultFigureVisible", "off");

% Pick ~6 times spread over the run
nT = numel(c.tt{1}); idx = unique(round(linspace(1, nT, 6)));
cmap = parula(numel(idx));

fig = figure(Position=[60 60 1250 800], Visible="off");
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

for k = 1:2                                     % row 1: profiles, free vs authoritative
    nexttile; hold on
    for j = 1:numel(idx)
        plot(c.phiS{k}(:,idx(j)), z, Color=cmap(j,:), LineWidth=2.2, ...
             DisplayName=sprintf("t=%.2f auth", c.tt{k}(idx(j))));
        plot(c.phiA{k}(:,idx(j)), z, Color=cmap(j,:), LineStyle="--", ...
             HandleVisibility="off");
    end
    set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off"); ylim([min(z) 0.05]);
    xlabel("\phi_b"); ylabel("z [m]");
    title(schemes(k) + ": solid = authoritative, dashed = free");
    if k == 1, legend(Location="southeast", FontSize=7); end
    grid on; box on
end

nexttile; hold on                                % relative growth of the difference
for k = 1:2
    dd = max(abs(c.phiA{k} - c.phiS{k}), [], 1);
    plot(c.tt{k}, dd, DisplayName=schemes(k));
end
xlabel("t [d]"); ylabel("max_z |difference|"); title("difference vs time");
legend(Location="southeast"); grid on; box on

for k = 1:2                                      % row 2: difference profiles
    nexttile; hold on
    for j = 1:numel(idx)
        plot(c.phiA{k}(:,idx(j)) - c.phiS{k}(:,idx(j)), z, Color=cmap(j,:), ...
             DisplayName=sprintf("t = %.2f d", c.tt{k}(idx(j))));
    end
    set(gca, YDir="reverse"); xline(0, "k:", HandleVisibility="off");
    yline(0, "k:", HandleVisibility="off"); ylim([min(z) 0.05]);
    xlabel("\phi_b^{free} - \phi_b^{sum}"); ylabel("z [m]");
    title(schemes(k) + ": difference profile");
    if k == 1, legend(Location="southwest", FontSize=7); end
    grid on; box on
end

nexttile; hold on                                % where the max difference sits
for k = 1:2
    [~, im] = max(abs(c.phiA{k} - c.phiS{k}), [], 1);
    plot(c.tt{k}, z(im), "o-", MarkerSize=3, DisplayName=schemes(k));
end
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off");
xlabel("t [d]"); ylabel("z of max |difference| [m]");
title("where the drift sits"); legend(Location="best"); grid on; box on

title(tl, sprintf("Free-running vs authoritative \\phi_b profiles (N = %d, \\kappa = %g)", ...
    c.NCells, c.Kappa));
out = fullfile(opts.OutDir, "fig10_drift_profiles" + opts.Tag + ".png");
exportgraphics(fig, out, Resolution=200); fprintf("wrote %s\n", out);

fprintf("\n%8s", "t [d]");
for k = 1:2, fprintf(" | %-28s", schemes(k) + ": max|diff| @ z"); end
fprintf("\n");
for j = 1:numel(idx)
    fprintf("%8.2f", c.tt{1}(idx(j)));
    for k = 1:2
        [mv, im] = max(abs(c.phiA{k}(:,idx(j)) - c.phiS{k}(:,idx(j))));
        fprintf(" | %12.4e @ %+7.4f m ", mv, z(im));
    end
    fprintf("\n");
end
close all
end
