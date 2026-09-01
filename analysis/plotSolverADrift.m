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
%
% This is now the PLOT half only: each scheme is run by probeSolverADrift, which
% writes analysis/probes/data/solverA_drift_<scheme>_n<N>.mat. That split is what
% lets the schemes run as a Slurm array (slurm/solverA_drift.sbatch) instead of
% in sequence. Any scheme whose cache is missing is run locally on demand, so
% the old one-call usage still works -- but on cosmos, submit the array.

arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 20
    opts.Schemes (1,:) string = ["shin", "matched", "bailo"]
    opts.OutDir (1,1) string = ""
end
here = fileparts(mfilename("fullpath")); W = fileparts(here);
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(here,"probes"));
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end
dataDir = fullfile(here, "probes", "data");

schemes = opts.Schemes;  nS = numel(schemes);
R = cell(1, nS);
for k = 1:nS
    fn = fullfile(dataDir, sprintf("solverA_drift_%s_n%d.mat", schemes(k), opts.NCells));
    if isfile(fn)
        L = load(fn, "rec");  R{k} = L.rec;
        fprintf("  loaded %s (%s, %.0f d)\n", schemes(k), R{k}.flag, R{k}.Days);
    else
        fprintf("  no cache for %s -- running locally\n", schemes(k));
        R{k} = probeSolverADrift(schemes(k), NCells=opts.NCells, ...
            Kappa=opts.Kappa, Days=opts.Days);
    end
end

z = R{1}.z;
tt   = cellfun(@(r) r.tt,   R, UniformOutput=false);
phiA = cellfun(@(r) r.phiA, R, UniformOutput=false);
phiS = cellfun(@(r) r.phiS, R, UniformOutput=false);
phiC = cellfun(@(r) r.phiC, R, UniformOutput=false);
D    = cellfun(@(r) r.D,    R, UniformOutput=false);
den  = max(cellfun(@(p) max(p(:)), phiS));

set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.5, ...
    "defaultFigureVisible", "off");
cols = lines(nS);
lsAll = ["-", "--", ":", "-."];  ls = lsAll(mod(0:nS-1, numel(lsAll)) + 1);
fig = figure(Position=[70 70 1200 760], Visible="off");
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

nexttile; hold on                                        % drift vs time
for k = 1:nS, plot(tt{k}, D{k}, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)); end
xlabel("t [d]"); ylabel("max_z |\phi_b^{free} - \phi_b^{sum}|");
% Envelope, because the diel oscillation hides the trend: reading sampled values
% I first called this "saturating", but the PEAKS rise monotonically
% (shin 0.0122 -> 0.0137 -> 0.0149 -> 0.0169 over 3 d).
for k = 1:nS
    % Dk, NOT D: `D = D{k}` clobbered the cell array on the first iteration and
    % made the second one throw. The peaks overlay has never rendered.
    Dk = D{k};
    li = find(Dk(2:end-1) > Dk(1:end-2) & Dk(2:end-1) > Dk(3:end)) + 1;  % no Signal Toolbox
    pk = Dk(li);
    if numel(pk) >= 2
        plot(tt{k}(li), pk, "o-", Color=cols(k,:)*0.6, LineWidth=1.0, ...
             MarkerSize=4, DisplayName=schemes(k)+" peaks");
    end
end
title("drift vs time -- peaks CREEP UPWARD, not saturating"); legend(Location="southeast"); grid on; box on

nexttile; hold on                                        % same, relative, log
for k = 1:nS, semilogy(tt{k}, 100*D{k}/den, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)); end
set(gca, YScale="log"); xlabel("t [d]"); ylabel("% of max \phi_b");
title("relative drift (log)"); legend(Location="southeast"); grid on; box on

nexttile; hold on                                        % bounds of the free state
for k = 1:nS
    plot(tt{k}, min(phiA{k}, [], 1), Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k)+" min");
    plot(tt{k}, max(phiA{k}, [], 1), Color=cols(k,:), LineStyle=":", DisplayName=schemes(k)+" max");
end
yline(0, "k-", HandleVisibility="off");
xlabel("t [d]"); ylabel("\phi_b^{free}"); title("bounds of the free state");
legend(Location="east", FontSize=7); grid on; box on

nexttile; hold on                                        % final profiles
for k = 1:nS
    plot(phiS{k}(:,end), z, Color=cols(k,:), LineWidth=2.6, DisplayName=schemes(k)+" authoritative");
    plot(phiA{k}(:,end), z, Color=cols(k,:), LineStyle="--", DisplayName=schemes(k)+" free");
end
set(gca, YDir="reverse"); xlabel("\phi_b"); ylabel("z [m]");
title(sprintf("profiles at t = %.0f d", opts.Days)); legend(Location="southeast", FontSize=7);
grid on; box on

nexttile; hold on                                        % pointwise drift profile
for k = 1:nS
    plot(phiA{k}(:,end) - phiS{k}(:,end), z, Color=cols(k,:), LineStyle=ls(k), DisplayName=schemes(k));
end
set(gca, YDir="reverse"); xline(0, "k:", HandleVisibility="off");
xlabel("\phi_b^{free} - \phi_b^{sum}"); ylabel("z [m]");
title("where the drift lives"); legend(Location="southeast"); grid on; box on

nexttile                                                 % one-step vs accumulated
oneStep = nan(1,nS);
for k = 1:nS
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
fprintf("\n%-8s %11s %11s %9s %12s %12s %8s\n", ...
    "scheme", "one-step", "accum max", "% of max", "free min", "free max", "wall s");
for k = 1:nS
    fprintf("%-8s %11.3e %11.3e %8.2f%% %12.3e %12.3e %8.0f\n", ...
        schemes(k), oneStep(k), max(D{k}), 100*max(D{k})/den, ...
        min(phiA{k}(:)), max(phiA{k}(:)), R{k}.wall);
end
close all
end
