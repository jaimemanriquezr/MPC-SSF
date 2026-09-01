function plotPulseOutflow(opts)
% PLOTPULSEOUTFLOW  Reproduce Manriquez2026 fig:pulse-outflow on the working set.
%
%   plotPulseOutflow
%
% Two panels, saved separately, matching pathogen.tex:80-98:
%   (a) fig:pat-low-inact  - marker concentration at the inflow and the outflow
%       for two AMPLITUDES of a pulse influx from day 30 to 32 (1x and 100x
%       c_ref). Source runs: pulse_fig_p1x, pulse_fig_p100x.
%   (b) fig:pat-high-inact - the 1x pulse for three RATE SETS:
%         Experiment 1  d_PAT = 0.02, p_PAT = 8.0   (tab:eco-parameters)
%         Experiment 2  d_PAT = 2.0                 (caption)
%         Experiment 3  d_PAT = 2e-6, p_PAT = 2.0   (caption)
%       Source runs: pulse_fig_p1x (= Experiment 1), pulse_fig_exp2, pulse_fig_exp3.
%
% CONCENTRATIONS, not log removals. The log-removal ceiling of EXPERIMENTS.md
% P1 caveat 4 (L defensible only to ~2.9) does not apply to this figure, but the
% same outflow term carries an ~8.4e-5 relative error on the exported mass, so
% the deep tail of the outflow curve is at the edge of what the solver resolves.
% A floor line is drawn at that level for exactly that reason.
%
% The published panels are semilog-y (the tail spans many decades and the
% caption is about tail LENGTH), so the y-axis is log here too. Inflow is drawn
% as the dashed rectangular pulse it is.
arguments
    opts.Dir (1,1) string = "";      % default: analysis/probes/data/pulse
    opts.OutDir (1,1) string = "";   % default: analysis/results/figures/repro
    opts.Floor (1,1) double = 8.39e-5;  % relative resolution of the outflow term
end
here = fileparts(mfilename("fullpath")); W = fileparts(here);
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis", "probes"));
D = opts.Dir;    if D == "",    D = fullfile(here, "probes", "data", "pulse"); end
O = opts.OutDir; if O == "",    O = fullfile(here, "results", "figures", "repro"); end
if ~isfolder(O), mkdir(O); end

% ---- panel (a): amplitudes ------------------------------------------------
A = struct( ...
    "tag",   ["fig_p1x", "fig_p100x"], ...
    "label", ["Pulse ($1\times c_{\mathrm{ref}}$)", "BigPulse ($100\times c_{\mathrm{ref}}$)"]);
makePanel(D, O, A.tag, A.label, "fig_pulse-outflow_low-inact", ...
    "Marker concentration, in- and outflow: pulse amplitude", opts.Floor);

% ---- panel (b): inactivation / predation rates ----------------------------
B = struct( ...
    "tag",   ["fig_p1x", "fig_exp2", "fig_exp3"], ...
    "label", ["Experiment 1: $d^{\mathrm{PAT}}_{20}=0.02$, $p^{\mathrm{PAT}}_{20}=8$", ...
              "Experiment 2: $d^{\mathrm{PAT}}_{20}=2$", ...
              "Experiment 3: $d^{\mathrm{PAT}}_{20}=2\!\times\!10^{-6}$, $p^{\mathrm{PAT}}_{20}=2$"]);
makePanel(D, O, B.tag, B.label, "fig_pulse-outflow_high-inact", ...
    "Marker concentration, in- and outflow: inactivation and predation", opts.Floor);
end

% =============================================================================
function makePanel(D, O, tags, labels, base, ttl, floorRel)
% `labels` is a local copy of the argument and is relabelled below for any run
% that aborted, so the legend states the clog rather than hiding it.
recs = cell(1, numel(tags));
for k = 1:numel(tags)
    fn  = fullfile(D, "pulse_" + tags(k) + ".mat");
    fnA = fullfile(D, "pulse_" + tags(k) + "_ABORTED.mat");
    if isfile(fn)
        S = load(fn, "rec"); recs{k} = S.rec;
    elseif isfile(fnA)
        % A run that CLOGGED still has a valid trajectory up to the abort; plot
        % it, trimmed at the abort, and say so in the legend. Frames past the
        % abort are zeros (probePulse writes the whole frame buffer).
        S = load(fnA, "rec", "aborted");
        r = S.rec;
        keep = r.tRel > 0 & r.tRel <= (S.aborted.tFinal - r.tSnapshot) + 1e-12;
        keep(1) = true;
        r.tRel = r.tRel(keep); r.effPAT = r.effPAT(keep); r.L = r.L(keep);
        recs{k} = r;
        labels(k) = labels(k) + sprintf(" -- CLOGGED at day %.2f", S.aborted.tFinal);
        warning("plotPulseOutflow:aborted", ...
            "%s aborted (%s at t = %g d); curve is truncated at the abort", ...
            tags(k), S.aborted.flag, S.aborted.tFinal);
    else
        warning("plotPulseOutflow:missing", "%s not found -- panel %s skipped", fn, base);
        return
    end
end

f = figure("Visible", "off", "Position", [100 100 760 520]);
ax = axes(f); hold(ax, "on"); set(ax, "YScale", "log");
co = [0.00 0.45 0.70; 0.84 0.37 0.00; 0.00 0.62 0.45];

for k = 1:numel(recs)
    r = recs{k};
    day = r.tRel(:) + r.tSnapshot;               % absolute day
    % Inflow: the rectangular pulse. NaN outside the window rather than 0 --
    % on a log axis a zero baseline drags the y-range down by 300 decades and
    % the whole curve collapses into the top gridline.
    cIn = nan(size(day));
    inPulse = r.tRel(:) >= r.PulseT0 & r.tRel(:) < r.PulseT1;
    cIn(inPulse) = r.cPeak;
    % outflow
    cOut = max(r.effPAT(:), realmin);
    c = co(mod(k-1, size(co,1)) + 1, :);
    % draw the inflow only once per distinct amplitude
    if k == 1 || r.cPeak ~= recs{1}.cPeak
        plot(ax, day, cIn, "--", "Color", c, "LineWidth", 1.4, ...
             "DisplayName", "inflow, " + labels(k));
    end
    plot(ax, day, cOut, "-", "Color", c, "LineWidth", 1.8, ...
         "DisplayName", "outflow, " + labels(k));
end

% Resolution floor of the outflow term (EXPERIMENTS.md P1 caveat 4): the
% cumulative export carries ~floorRel relative error, so an outflow value this
% far below the pulse peak is not resolved by the solver.
peak = max(cellfun(@(r) r.cPeak, recs));
yfloor = floorRel*peak;
% Fix the range explicitly: from one decade below the resolution floor to just
% above the largest inflow. Everything below the floor is solver noise and
% showing it (the curves run to 1e-300) makes the figure unreadable.
ylim(ax, [yfloor/10, 3*peak]);
yl = yline(ax, yfloor, ":", sprintf("outflow-term resolution (%.1e of peak)", floorRel));
yl.Color = [0.4 0.4 0.4];
yl.LabelHorizontalAlignment = "right";
yl.LabelVerticalAlignment = "bottom";
yl.FontSize = 8;
yl.Annotation.LegendInformation.IconDisplayStyle = "off";

% pulse window shading
r1 = recs{1};
xp = [r1.PulseT0 r1.PulseT1] + r1.tSnapshot;
yy = ylim(ax);
p = patch(ax, [xp(1) xp(2) xp(2) xp(1)], [yy(1) yy(1) yy(2) yy(2)], [0.9 0.9 0.9], ...
      "EdgeColor", "none", "FaceAlpha", 0.35);
uistack(p, "bottom");
p.Annotation.LegendInformation.IconDisplayStyle = "off";
ylim(ax, yy);

grid(ax, "on"); box(ax, "on");
xlabel(ax, "time (d)");
ylabel(ax, "c^{PAT} (kg m^{-3})");
title(ax, ttl);
lg = legend(ax, "Location", "southwest", "Interpreter", "latex");
lg.FontSize = 8;
xlim(ax, [r1.tSnapshot, r1.tRel(end) + r1.tSnapshot]);

exportgraphics(f, fullfile(O, base + ".png"), "Resolution", 200);
exportgraphics(f, fullfile(O, base + ".pdf"), "ContentType", "vector");
close(f);
fprintf("wrote %s.{png,pdf}\n", fullfile(O, base));
end
