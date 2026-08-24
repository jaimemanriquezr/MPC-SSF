function plotDriftFinal(opts)
% PLOTDRIFTFINAL  The drift picture at a SINGLE time (default: the last frame).
%
% fig10 overlays six times, which shows the front migrating but crowds the
% detail. This is the same comparison at one instant, with the vertical range
% cropped to where phi_b is actually nonzero.

arguments
    opts.Cache (1,1) string = ""
    opts.OutDir (1,1) string = ""
    opts.Time (1,1) double = NaN      % NaN -> last frame
end
here = fileparts(mfilename("fullpath"));
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if opts.Cache == "", opts.Cache = fullfile(opts.OutDir, "solverA_drift_n100.mat"); end
c = load(opts.Cache).cache;
z = c.z; schemes = c.schemes; dz = mean(diff(z));

j = zeros(1,2); tsel = zeros(1,2);
for k = 1:2
    if isnan(opts.Time), j(k) = numel(c.tt{k});
    else, [~, j(k)] = min(abs(c.tt{k} - opts.Time)); end
    tsel(k) = c.tt{k}(j(k));
end

% crop to where there is biofilm, else 90% of the panel is empty
act = false(size(z));
for k = 1:2, act = act | c.phiS{k}(:,j(k)) > 1e-4; end
zlo = min(z(act)) - 0.02; zhi = max(z(act)) + 0.01;

set(groot, "defaultAxesFontSize", 12, "defaultLineLineWidth", 1.8, ...
    "defaultFigureVisible", "off");
cols = lines(2);
fig = figure(Position=[60 60 1250 480], Visible="off");
tl = tiledlayout(fig, 1, 3, TileSpacing="compact", Padding="compact");

nexttile; hold on                                  % profiles
for k = 1:2
    plot(c.phiS{k}(:,j(k)), z, "-",  Color=cols(k,:), LineWidth=2.6, ...
         DisplayName=schemes(k)+" authoritative");
    plot(c.phiA{k}(:,j(k)), z, "--", Color=cols(k,:), LineWidth=1.6, ...
         DisplayName=schemes(k)+" free-running");
end
set(gca, YDir="reverse"); ylim([zlo zhi]); yline(0, "k:", HandleVisibility="off");
xlabel("\phi_b [-]"); ylabel("z [m]");
title(sprintf("\\phi_b at t = %.2f d", tsel(1)));
legend(Location="southeast", FontSize=9); grid on; box on

nexttile; hold on                                  % difference
for k = 1:2
    d = c.phiA{k}(:,j(k)) - c.phiS{k}(:,j(k));
    plot(d, z, Color=cols(k,:), DisplayName=schemes(k));
    [~, im] = max(abs(d));
    plot(d(im), z(im), "o", Color=cols(k,:), MarkerFaceColor=cols(k,:), ...
         MarkerSize=7, HandleVisibility="off");
    text(d(im), z(im), sprintf("  %.2e @ %.3f m", d(im), z(im)), FontSize=8, ...
         Color=cols(k,:)*0.7);
end
set(gca, YDir="reverse"); ylim([zlo zhi]);
xline(0, "k:", HandleVisibility="off"); yline(0, "k:", HandleVisibility="off");
xlabel("\phi_b^{free} - \phi_b^{sum}"); ylabel("z [m]");
title("difference (dipole = displaced front)");
legend(Location="southeast", FontSize=9); grid on; box on

nexttile; hold on                                  % cumulative mass from the top
for k = 1:2
    cf = cumsum(c.phiA{k}(:,j(k)))*dz;  ca = cumsum(c.phiS{k}(:,j(k)))*dz;
    plot(ca, z, "-",  Color=cols(k,:), LineWidth=2.6, DisplayName=schemes(k)+" authoritative");
    plot(cf, z, "--", Color=cols(k,:), LineWidth=1.6, DisplayName=schemes(k)+" free-running");
    fprintf("%-6s t=%.2f d | integral free %.5e  auth %.5e  rel diff %.3e | " + ...
            "pointwise max %.4e (%.2f%% of max phi_b)\n", schemes(k), tsel(k), ...
            cf(end), ca(end), abs(cf(end)-ca(end))/ca(end), ...
            max(abs(c.phiA{k}(:,j(k)) - c.phiS{k}(:,j(k)))), ...
            100*max(abs(c.phiA{k}(:,j(k)) - c.phiS{k}(:,j(k))))/max(c.phiS{k}(:,j(k))));
end
set(gca, YDir="reverse"); ylim([zlo zhi]); yline(0, "k:", HandleVisibility="off");
xlabel("cumulative \int\phi_b dz from top [m]"); ylabel("z [m]");
title("cumulative mass -- curves converge");
legend(Location="southeast", FontSize=9); grid on; box on

title(tl, sprintf("Free-running vs authoritative at t = %.2f d (N = %d, \\kappa = %g)", ...
    tsel(1), c.NCells, c.Kappa));
out = fullfile(opts.OutDir, "fig11_drift_final.png");
exportgraphics(fig, out, Resolution=200); fprintf("wrote %s\n", out);
close all
end
