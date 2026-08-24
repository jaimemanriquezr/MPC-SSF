function plotPhibProfiles(opts)
% PLOTPHIBPROFILES  phi_b(z) at successive times, from the E3 production run.
%
% Six views of the same data, because a single full-column plot hides everything:
% phi_b is ~0 over 90% of the domain and all the structure sits in a few
% centimetres around z = 0.

arguments
    opts.Root (1,1) string = ""
    opts.OutDir (1,1) string = ""
end
here = fileparts(mfilename("fullpath"));
% E3 writes to OutRoot relative to the REPO ROOT (the sbatch passes
% 'results/manuscript_k1e6_20d'), not under analysis/ where manuscriptExperiments
% puts its own default.
if opts.Root == "", opts.Root = fullfile(fileparts(here), "results", "manuscript_k1e6_20d"); end
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

S    = readmatrix(fullfile(opts.Root, "E3_longterm", "phib_snapshots_summer.csv"));
days = readmatrix(fullfile(opts.Root, "E3_longterm", "phib_snapshot_days.csv"));
mS   = readmatrix(fullfile(opts.Root, "E3_longterm", "total_biomass_summer.csv"));
W    = readmatrix(fullfile(opts.Root, "E1_seasons", "phib_final_winter.csv"));
mW   = readmatrix(fullfile(opts.Root, "E1_seasons", "total_biomass_winter.csv"));

z = S(:,1); P = S(:,2:end); zw = W(:,1); pw = W(:,2);
set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.6, ...
    "defaultFigureVisible", "off");
cmap = parula(numel(days));

fig = figure(Position=[60 60 1250 800], Visible="off");
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

nexttile; hold on                                   % full column
for k = 1:numel(days)
    plot(P(:,k), z, Color=cmap(k,:), DisplayName=sprintf("day %g", days(k)));
end
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off");
xlabel("\phi_b [-]"); ylabel("z [m]"); title("full column");
legend(Location="southeast", FontSize=8); grid on; box on

nexttile; hold on                                   % schmutzdecke zoom
for k = 1:numel(days), plot(P(:,k), z, Color=cmap(k,:)); end
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off");
ylim([-0.06 0.10]); xlabel("\phi_b [-]"); ylabel("z [m]");
title("schmutzdecke (z \in [-6, 10] cm)"); grid on; box on

nexttile; hold on                                   % log x, to see the tail
for k = 1:numel(days), semilogx(max(P(:,k), 1e-8), z, Color=cmap(k,:)); end
set(gca, XScale="log", YDir="reverse"); yline(0, "k:", HandleVisibility="off");
xlim([1e-6 1]); xlabel("\phi_b [-] (log)"); ylabel("z [m]");
title("log scale -- penetration depth"); grid on; box on

nexttile; hold on                                   % summer vs winter, final
plot(P(:,end), z, "r-", LineWidth=2.4, DisplayName=sprintf("summer, day %g", days(end)));
plot(pw, zw, "b--", LineWidth=2.4, DisplayName="winter, day 20");
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off"); ylim([-0.06 0.30]);
xlabel("\phi_b [-]"); ylabel("z [m]"); title("summer vs winter at 20 d");
legend(Location="southeast", FontSize=8); grid on; box on

nexttile; hold on                                   % biomass
plot(mS(:,1), mS(:,2), "r-", DisplayName="summer");
plot(mW(:,1), mW(:,2), "b--", DisplayName="winter");
xlabel("t [d]"); ylabel("total biomass [kg/m^2]"); title("total biomass");
legend(Location="northwest", FontSize=8); grid on; box on

nexttile; hold on                                   % peak and its depth
[pk, ipk] = max(P, [], 1);
yyaxis left;  plot(days, pk, "o-"); ylabel("peak \phi_b");
yyaxis right; plot(days, z(ipk), "s--"); ylabel("depth of peak [m]"); set(gca, YDir="reverse");
xlabel("t [d]"); title("peak \phi_b and where it sits"); grid on; box on

title(tl, "E3 production run: \phi_b profiles (\kappa = 10^{-6}, N = 500, 20 d)");
out = fullfile(opts.OutDir, "fig9_phib_profiles_E3.png");
exportgraphics(fig, out, Resolution=200); fprintf("wrote %s\n", out);

fprintf("\n%6s %12s %12s %14s %12s\n", "day", "peak phi_b", "depth [m]", "phi_b at z=0", "mass [kg/m2]");
iz0 = find(z >= 0, 1);
for k = 1:numel(days)
    mk = interp1(mS(:,1), mS(:,2), days(k), "linear", "extrap");
    fprintf("%6g %12.5f %12.4f %14.5f %12.5f\n", days(k), pk(k), z(ipk(k)), P(iz0,k), mk);
end
fprintf("%6s %12.5f %12.4f %14.5f %12.5f\n", "wint", max(pw), zw(find(pw==max(pw),1)), ...
    pw(find(zw>=0,1)), mW(end,2));
close all
end
