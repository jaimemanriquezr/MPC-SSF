function recreateFigsManriquez2026(what, opts)
% RECREATEFIGSMANRIQUEZ2026  Redraw the Manriquez2026 results figures, panel for
% panel, with the working-set data of this repo.
%
%   recreateFigsManriquez2026()             % every figure whose data exists
%   recreateFigsManriquez2026("seasons")    % one figure family
%
% The goal is VISUAL one-to-one comparability with the published PDFs in
% ../../manuscripts/AWR-SSF/figures/: same panel structure, axis orientation,
% ranges, labels, units, legends, line styles and zoom windows. The DATA is the
% audited working set (E4/E7), so the curves themselves are not expected to
% agree with the published ones.
%
% Figure keys (each writes PNG + PDF into analysis/results/figures/recreation/):
%   light       results_light-seasons-time.pdf
%   seasons     results_biofilm_seasons{,_zoomed}.pdf
%   roofed      results_biofilm_covered{,_zoomed,-ultrazoomed_zoomed}.pdf
%   month       results_evolution_summer{,_zoomed,_mass}.pdf
%   scrape      results_evolution_Scraped_GP{0,4,8,12}_zoomed.pdf, _GP_mass.pdf
%   outflow1d   Outflow1D{O2,IC,NH4,HPO4,DOM}.pdf
%   outflow2d   Outflow2D{O2,IC,NH4,HPO4,DOM}.pdf
%   patpulse    March-11/PATPulse_{HighLow,FastSlow}.pdf
%   filtration  March-11/FiltrationRate_HETPHO.pdf
%   pat2d       March-11/PAT{Flowing,Matrix,Enclosed}{Reference,LowInactivation}.pdf
%
% Two further keys have NO published counterpart -- they are manuscript
% CANDIDATES (see .claude/CANDIDATE-FIGURES.md), drawn from the same chains:
%   composition rec_composition_depth.pdf   community structure vs depth, four
%               scenarios at the common horizon, paired with phi_b(z)
%   tenore      rec_composition_time.pdf    Tenore2021-style: depth-integrated
%               biomass fractions f_i(t) over the bed, paired with areal biomass
%
% Every figure whose input data is missing prints a "SKIP <key>: ..." line and
% is left out; nothing here starts a simulation. Deviations are recorded in
% analysis/results/figures/recreation/README.md.
arguments
    what (1,1) string = "all"
    opts.LitTag    (1,1) string = "fld2x_lit"      % summer / uncovered arm
    opts.WinterTag (1,1) string = "fld2x_winter"   % 3 C winter arm
    opts.CovTag    (1,1) string = "fld2x_cov01"    % 1 % light (covered)
    opts.DarkTag   (1,1) string = "fld2x_dark"     % 0 % light (covered fallback)
end
here    = fileparts(mfilename("fullpath"));
repo    = fileparts(here);
addpath(genpath(fullfile(repo, "src")));
dataDir = fullfile(here, "probes", "data", "chain");
pulseDir= fullfile(here, "probes", "data", "pulse");
scrapeDir = fullfile(here, "results", "figures", "repro", "data");
outDir  = fullfile(here, "results", "figures", "recreation");
if ~isfolder(outDir), mkdir(outDir); end

keys = ["light" "seasons" "roofed" "month" "scrape" "outflow1d" "outflow2d" ...
        "patpulse" "filtration" "pat2d" "composition" "tenore"];
if what == "all"
    todo = keys;
elseif ismember(what, keys)
    todo = what;
else
    error("recreateFigsManriquez2026:what", "what must be ""all"" or one of: %s", strjoin(keys, ", "));
end

for k = todo
    switch k
        case "light",      figLight(outDir);
        case "seasons",    figSeasons(dataDir, outDir, opts);
        case "roofed",     figRoofed(dataDir, outDir, opts);
        case "month",      figMonth(dataDir, outDir, opts);
        case "scrape",     figScrape(scrapeDir, outDir);
        case "outflow1d",  figOutflow1D(dataDir, outDir, opts);
        case "outflow2d",  figOutflow2D(dataDir, outDir, opts);
        case "patpulse",   figPatPulse(pulseDir, outDir);
        case "filtration", figFiltration(pulseDir, outDir);
        case "pat2d",      figPat2D(pulseDir, outDir);
        case "composition",figComposition(dataDir, outDir, opts);
        case "tenore",     figTenore(dataDir, outDir, opts);
    end
end
fprintf("recreation figures in %s\n", outDir);
end

% ======================================================================== %
% fig:seasons-light -- results_light-seasons-time.pdf
% ======================================================================== %
function figLight(outDir)
t  = linspace(0, 1, 2000);
Ls = max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);          % published summer
Lw = max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);          % published winter
Lp = 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50);  % working-set chain

[f, ax] = squareAxes();
xlabel(ax, "Time [d]", Interpreter="latex", FontSize=20);
ylabel(ax, "Incoming daylight irradiation [--]", Interpreter="latex", FontSize=20);
plot(ax, t, Ls, "k-",  LineWidth=1.0, DisplayName="Summer");
plot(ax, t, Lw, "k--", LineWidth=3.0, DisplayName="Winter");
plot(ax, t, Lp, "k:",  LineWidth=2.0, DisplayName="Working set");
xlim(ax, [0 1]); ylim(ax, [0 1]);
boxedLegend(ax, "northeast");
saveBoth(f, outDir, "rec_light-seasons-time");
end

% ======================================================================== %
% fig:seasons-results -- results_biofilm_seasons{,_zoomed}.pdf
% ======================================================================== %
function figSeasons(dataDir, outDir, opts)
if isempty(chainLegs(dataDir, opts.LitTag))
    fprintf("SKIP seasons: no chain_%s_leg*.mat\n", opts.LitTag); return
end
S = loadChain(dataDir, opts.LitTag, MaxT=90);
[phiS, tS] = profileAt(S, 90);
haveW = ~isempty(chainLegs(dataDir, opts.WinterTag));
if haveW
    W = loadChain(dataDir, opts.WinterTag, MaxT=90);
    [phiW, tW] = profileAt(W, 90);
else
    fprintf("PARTIAL seasons: no chain_%s_leg*.mat -- Winter curve omitted\n", opts.WinterTag);
end

windows = {[-0.2 1.0], [-0.02 0.04]};
names   = ["rec_biofilm_seasons", "rec_biofilm_seasons_zoomed"];
locs    = ["southeast", "northeast"];
for q = 1:2
    [f, ax] = profileAxes();
    plot(ax, phiS, S.z, "k-",  LineWidth=1.0, DisplayName="Summer");
    if haveW, plot(ax, phiW, W.z, "k--", LineWidth=3.0, DisplayName="Winter"); end
    ylim(ax, windows{q});
    if q == 2, yticks(ax, -0.02:0.01:0.04); end
    boxedLegend(ax, locs(q));
    saveBoth(f, outDir, names(q));
end
if haveW
    fprintf("  seasons: summer t = %.2f d, winter t = %.2f d\n", tS, tW);
else
    fprintf("  seasons: summer t = %.2f d\n", tS);
end
end

% ======================================================================== %
% fig:roofed-results -- results_biofilm_covered{,_zoomed,-ultrazoomed_zoomed}
% ======================================================================== %
function figRoofed(dataDir, outDir, opts)
if isempty(chainLegs(dataDir, opts.LitTag))
    fprintf("SKIP roofed: no chain_%s_leg*.mat\n", opts.LitTag); return
end
covTag = opts.CovTag;
if isempty(chainLegs(dataDir, covTag))
    if isempty(chainLegs(dataDir, opts.DarkTag))
        fprintf("SKIP roofed: neither chain_%s_leg*.mat nor chain_%s_leg*.mat\n", ...
            opts.CovTag, opts.DarkTag); return
    end
    covTag = opts.DarkTag;
    fprintf("PARTIAL roofed: %s missing -- ""Covered"" drawn from %s (0 %% light)\n", ...
        opts.CovTag, covTag);
end
U = loadChain(dataDir, opts.LitTag, MaxT=30);
C = loadChain(dataDir, covTag,      MaxT=30);
[phiU, tU] = profileAt(U, 30);
[phiC, tC] = profileAt(C, 30);
fprintf("  roofed: uncovered t = %.2f d, covered t = %.2f d (%s)\n", tU, tC, covTag);

% panel 1 (whole filter) and panel 2 (zoom), published windows
windows = {[-0.2 1.0], [-0.025 0.01]};
names   = ["rec_biofilm_covered", "rec_biofilm_covered_zoomed"];
locs    = ["southeast", "northeast"];
for q = 1:2
    [f, ax] = profileAxes();
    plot(ax, phiC, C.z, "k-",  LineWidth=1.0, DisplayName="Covered");
    plot(ax, phiU, U.z, "k--", LineWidth=3.0, DisplayName="Uncovered");
    ylim(ax, windows{q});
    if q == 2, yticks(ax, -0.025:0.005:0.01); end
    boxedLegend(ax, locs(q));
    saveBoth(f, outDir, names(q));
end

% panel 3: ultra zoom. Published window is z in [0.10 0.11], phi in
% [0.39 0.405] -- 1 cm of depth against a 0.015-wide phi window. The phi window
% is re-centred on our own curves (see README) so that the two lines are inside
% the axes; the width (0.015) and the depth window are the published ones.
[f, ax] = profileAxes();
zw = C.z >= 0.10 & C.z <= 0.11;
vals = [phiC(zw); phiU(zw)];
c = (min(vals) + max(vals))/2;
plot(ax, phiC, C.z, "k-",  LineWidth=1.0, DisplayName="Covered");
plot(ax, phiU, U.z, "k--", LineWidth=3.0, DisplayName="Uncovered");
ylim(ax, [0.10 0.11]); yticks(ax, 0.10:0.002:0.11);
xlim(ax, [c - 0.0075, c + 0.0075]);
xticks(ax, round((c - 0.0075):0.005:(c + 0.0075), 4));
boxedLegend(ax, "southeast");
saveBoth(f, outDir, "rec_biofilm_covered-ultrazoomed_zoomed");
end

% ======================================================================== %
% fig:month-change -- results_evolution_summer{,_zoomed,_mass}.pdf
% ======================================================================== %
function figMonth(dataDir, outDir, opts)
if isempty(chainLegs(dataDir, opts.LitTag))
    fprintf("SKIP month: no chain_%s_leg*.mat\n", opts.LitTag); return
end
days = [1 2 5 10 20];
mk   = ["o" "+" "*" "." "x"];
A = loadChain(dataDir, opts.LitTag, MaxT=max(days));
windows = {[-0.2 1.0], [-0.02 0.04]};
names   = ["rec_evolution_summer", "rec_evolution_summer_zoomed"];
locs    = ["southeast", "northeast"];
for q = 1:2
    [f, ax] = profileAxes();
    for j = 1:numel(days)
        p = profileAt(A, days(j));
        plot(ax, p, A.z, "k-", Marker=mk(j), MarkerSize=6, LineWidth=0.8, ...
            MarkerIndices=1:12:numel(A.z), DisplayName=sprintf("$t = %g$", days(j)));
    end
    ylim(ax, windows{q});
    if q == 2, yticks(ax, -0.02:0.01:0.04); end
    boxedLegend(ax, locs(q));
    saveBoth(f, outDir, names(q));
end

% mass panel: 0 .. 90 d, single thick black curve
M = loadChain(dataDir, opts.LitTag, MaxT=90, NeedMass=true);
if isempty(M.mass)
    fprintf("PARTIAL month: no results_py in the legs -- mass panel skipped\n"); return
end
[f, ax] = squareAxes();
xlabel(ax, "Time [d]", Interpreter="latex", FontSize=20);
ylabel(ax, "Biofilm mass [kg/m$^2$]", Interpreter="latex", FontSize=20);
plot(ax, M.ts, M.mass, "k-", LineWidth=3.0);
xlim(ax, [0 90]); xticks(ax, 0:10:90);
ylim(ax, [0 max(1.05*max(M.mass), eps)]);
legend(ax, "off");
saveBoth(f, outDir, "rec_evolution_summer_mass");
fprintf("  month: mass %.3f -> %.3f kg/m^2 over %.1f d\n", M.mass(1), M.mass(end), M.ts(end));
end

% ======================================================================== %
% fig:month-scrape -- results_evolution_Scraped_GP{0,4,8,12}_zoomed, _GP_mass
% ======================================================================== %
function figScrape(scrapeDir, outDir)
g = dir(fullfile(scrapeDir, "scrape_GP*.mat"));
if isempty(g)
    fprintf("SKIP scrape: no scrape_GP*.mat in %s\n", scrapeDir); return
end
S = struct("cm", {}, "sc", {});
for k = 1:numel(g)
    D = load(fullfile(g(k).folder, g(k).name), "sc");
    S(k).sc = D.sc; S(k).cm = round(D.sc.depth*100);
end
[~, o] = sort([S.cm]); S = S(o);

days = [30 31 32 35 40 50];
mk   = ["o" "+" "*" "." "x" "square"];
tMax = 30 + max(arrayfun(@(s) max(s.sc.ts), S));
if tMax < 50
    days = days(days <= tMax); mk = mk(1:numel(days));
    fprintf("PARTIAL scrape: regrowth data stop at t = %.0f d -- t = 50 d curve omitted\n", tMax);
end
for k = 1:numel(S)
    sc = S(k).sc;
    [f, ax] = profileAxes();
    for j = 1:numel(days)
        [~, i] = min(abs(sc.ts + 30 - days(j)));
        plot(ax, sc.phi(:, i), sc.z, "k-", Marker=mk(j), MarkerSize=6, LineWidth=0.8, ...
            MarkerIndices=1:12:numel(sc.z), DisplayName=sprintf("$t = %g$", days(j)));
    end
    ylim(ax, [-0.1 0.2]); yticks(ax, -0.1:0.05:0.2);
    boxedLegend(ax, "northeast");
    saveBoth(f, outDir, sprintf("rec_evolution_Scraped_GP%d_zoomed", S(k).cm));
end

% mass panel: one curve per scraping depth
lw = [1.0 1.5 3.0 3.0]; ls = ["-" "--" "--" ":"];
[f, ax] = squareAxes();
xlabel(ax, "Time [d]", Interpreter="latex", FontSize=20);
ylabel(ax, "Biofilm mass [kg/m$^2$]", Interpreter="latex", FontSize=20);
mAll = [];
for k = 1:numel(S)
    j = min(k, numel(ls));
    plot(ax, S(k).sc.ts + 30, S(k).sc.mass, "k" + ls(j), LineWidth=lw(j), ...
        DisplayName=sprintf("$z_S = %d$ [cm]", S(k).cm));
    mAll = [mAll, S(k).sc.mass]; %#ok<AGROW>
end
xlim(ax, [30 max(60, tMax)]);
ylim(ax, [min(mAll) 1.02*max(mAll)]);
boxedLegend(ax, "southeast");
saveBoth(f, outDir, "rec_evolution_Scraped_GP_mass");
end

% ======================================================================== %
% fig:1d-outflow-liquids -- Outflow1D{O2,IC,NH4,HPO4,DOM}.pdf, one file each
% ======================================================================== %
function figOutflow1D(dataDir, outDir, opts)
if isempty(chainLegs(dataDir, opts.LitTag))
    fprintf("SKIP outflow1d: no chain_%s_leg*.mat\n", opts.LitTag); return
end
A = loadChain(dataDir, opts.LitTag, MaxT=30);
keep = A.ts <= 30;
for nm = ["O2" "IC" "NH4" "HPO4" "DOM"]
    j = find(A.effNames == nm);
    c = A.eff(j, keep)*1e3;                 % kg/m3 = g/L -> mg/L
    [f, ax] = squareAxes(Grid=false);
    xlabel(ax, "Time [days]", Interpreter="latex", FontSize=20);
    ylabel(ax, "Concentration [mg/L]", Interpreter="latex", FontSize=20);
    plot(ax, A.ts(keep), c, "k-", LineWidth=1.2);
    xlim(ax, [0 30]); xticks(ax, 0:5:30);
    ylim(ax, [0 max(1.05*max(c), eps)]);
    legend(ax, "off");
    saveBoth(f, outDir, "rec_Outflow1D" + nm);
end
end

% ======================================================================== %
% fig:2d-plots -- Outflow2D{O2,IC,NH4,HPO4,DOM}.pdf, one file each
% ======================================================================== %
function figOutflow2D(dataDir, outDir, opts)
if isempty(chainLegs(dataDir, opts.LitTag))
    fprintf("SKIP outflow2d: no chain_%s_leg*.mat\n", opts.LitTag); return
end
A = loadChain(dataDir, opts.LitTag, MaxT=30, NeedFields=true);
if isempty(fieldnames(A.fields))
    fprintf("SKIP outflow2d: no results_py fields in the legs\n"); return
end
keep = A.fieldT <= 30;
[T, Z] = meshgrid(A.fieldT(keep), A.z);
for nm = ["O2" "IC" "NH4" "HPO4" "DOM"]
    C = A.fields.(nm)(:, keep);             % kg/m3 = g/L, as published
    f = figure("Visible", "off", "Position", [0 0 900 560], Color="w");
    ax = axes(f);
    surf(ax, T, Z, C, LineStyle="none");
    colormap(ax, "jet");
    cb = colorbar(ax); cb.Label.String = "Concentration [g/L]";
    cb.Label.Interpreter = "latex"; cb.Label.FontSize = 16;
    cb.TickLabelInterpreter = "latex";
    xlabel(ax, "Time [days]", Interpreter="latex", FontSize=20);
    ylabel(ax, "Depth [m]", Interpreter="latex", FontSize=20);
    % YDir reverse puts depth -1 on the far left and +1 at the axis vertex, as
    % in the published Outflow2D panels.
    set(ax, TickLabelInterpreter="latex", FontSize=14, YDir="reverse");
    xlim(ax, [0 30]); ylim(ax, [-1 1]);
    zlim(ax, [0 max(max(C(:)), eps)*1.02]);
    % RASTERISED, not vector like the rest (saveBoth). These are surf plots:
    % vector export writes one path per face, which made the PDFs large and slow
    % to render in the compiled manuscript. 300 dpi at this canvas size is well
    % above print need.
    stem = fullfile(outDir, "rec_Outflow2D" + nm);
    exportgraphics(f, stem + ".png", Resolution=300);
    exportgraphics(f, stem + ".pdf", ContentType="image", Resolution=300);
    close(f);
    fprintf("  rec_Outflow2D%s.{png,pdf} (raster, 300 dpi)\n", nm);
end
end

% ======================================================================== %
% fig:pat-low-inact / fig:pat-high-inact -- PATPulse_{HighLow,FastSlow}.pdf
% ======================================================================== %
function figPatPulse(pulseDir, outDir)
[p1, n1] = loadPulse(pulseDir, "p1x");
[p2, n2] = loadPulse(pulseDir, "p1em3");
[e2, m2] = loadPulse(pulseDir, "exp2");
[e3, m3] = loadPulse(pulseDir, "exp3");

% --- panel (a): PATPulse_HighLow
if isempty(p1)
    fprintf("SKIP patpulse (HighLow): no pulse_fig_p1x*.mat\n");
else
    [f, ax] = pulseAxes();
    plotPulseIn (ax, p1, "k-", "o", true,  "Pulse 1: IN");
    if ~isempty(p2)
        plotPulseIn(ax, p2, "k-", "v", true,  "Pulse 2: IN");
    end
    plotPulseOut(ax, p1, "k-", "o", false, "Pulse 1: OUT");
    if ~isempty(p2)
        plotPulseOut(ax, p2, "k-", "v", false, "Pulse 2: OUT");
    else
        fprintf("PARTIAL patpulse (HighLow): no pulse_fig_p1em3*.mat -- Pulse 2 omitted\n");
    end
    boxedLegend(ax, "northeast");
    saveBoth(f, outDir, "rec_PATPulse_HighLow");
    fprintf("  patpulse HighLow: %s%s\n", n1, ternary(isempty(p2), "", " + " + n2));
end

% --- panel (b): PATPulse_FastSlow
if isempty(p1)
    fprintf("SKIP patpulse (FastSlow): no pulse_fig_p1x*.mat\n"); return
end
[f, ax] = pulseAxes();
plotPulseIn (ax, p1, "k-", "o",      true,  "Pulse: IN");
plotPulseOut(ax, p1, "k-", "o",      false, "Experiment 1");
if ~isempty(e2), plotPulseOut(ax, e2, "k-", "v",      false, "Experiment 2");
else, fprintf("PARTIAL patpulse (FastSlow): no pulse_fig_exp2*.mat\n"); end
if ~isempty(e3), plotPulseOut(ax, e3, "k-", "square", false, "Experiment 3");
else, fprintf("PARTIAL patpulse (FastSlow): no pulse_fig_exp3*.mat\n"); end
boxedLegend(ax, "northeast");
saveBoth(f, outDir, "rec_PATPulse_FastSlow");
fprintf("  patpulse FastSlow: %s / %s / %s\n", n1, m2, m3);
end

function [f, ax] = pulseAxes()
[f, ax] = squareAxes(Grid=false);
xlabel(ax, "Time [days]", Interpreter="latex", FontSize=20);
ylabel(ax, "Concentration [g/L]", Interpreter="latex", FontSize=20);
set(ax, YScale="log");
xlim(ax, [30 37]); xticks(ax, 30:37);
ylim(ax, [1e-12 1e0]); yticks(ax, 10.^(-12:2:0));
end

function plotPulseIn(ax, r, sty, mk, filled, name)
t = pulseTime(r);
c = zeros(size(t));
c(t >= r.tSnapshot + r.PulseT0 & t <= r.tSnapshot + r.PulseT1) = r.cPeak;
c(c == 0) = NaN;
h = plot(ax, t, c, sty, LineWidth=1.2, Marker=mk, MarkerSize=7, ...
    MarkerIndices=markerIdx(c), DisplayName=name);
if filled, h.MarkerFaceColor = "k"; else, h.MarkerFaceColor = "w"; end
end

function plotPulseOut(ax, r, sty, mk, filled, name)
t = pulseTime(r);
c = r.effPAT(:).';
c(c <= 0) = NaN;
h = plot(ax, t, c, sty, LineWidth=1.2, Marker=mk, MarkerSize=7, ...
    MarkerIndices=markerIdx(c), DisplayName=name);
if filled, h.MarkerFaceColor = "k"; else, h.MarkerFaceColor = "w"; end
end

function idx = markerIdx(c)
ok = find(isfinite(c));
if isempty(ok), idx = 1; return, end
idx = unique(round(linspace(ok(1), ok(end), min(10, numel(ok)))));
end

function t = pulseTime(r)
if isfield(r, "tRel"), t = r.tSnapshot + r.tRel(:).';
else,                  t = r.tSnapshot + r.ts(:).'; end
end

% ======================================================================== %
% fig:pat-filtration -- March-11/FiltrationRate_HETPHO.pdf
% ======================================================================== %
function figFiltration(pulseDir, outDir)
[p1, n1] = loadPulse(pulseDir, "p1x");
[p2, n2] = loadPulse(pulseDir, "p1em3");
if isempty(p1)
    fprintf("SKIP filtration: no pulse_fig_p1x*.mat\n"); return
end
red = [0.85 0 0]; green = [0 0.5 0];
[f, ax] = squareAxes(Grid=false);
xlabel(ax, "", Interpreter="latex");     % published panel carries no axis labels
ylabel(ax, "", Interpreter="latex");
sets = {p1, p2}; nmSet = ["Pulse 1", "Pulse 2"];
sty = ["-", "--"]; lw = [1.2 3.0];
nm = ["HET" "PHO"]; col = {red, green};
allL = [];
for q = 1:2                                   % species outer: legend order as published
    for s = 1:2
        r = sets{s};
        if isempty(r), continue, end
        t = pulseTime(r);
        j   = find(r.effNames == nm(q));         % effluent rows: effNames order
        cin = r.Influent(influentIndex(nm(q)));  % influent vector: probeChain order
        L = log10(cin ./ max(r.effluent(j, :), realmin));
        allL = [allL, L]; %#ok<AGROW>
        plot(ax, t, L, sty(s), Color=col{q}, LineWidth=lw(s), ...
            DisplayName=sprintf("%s: %s", nm(q), nmSet(s)));
    end
end
if isempty(p2)
    fprintf("PARTIAL filtration: no pulse_fig_p1em3*.mat -- Pulse 2 curves omitted\n");
end
xlim(ax, [30 37]); xticks(ax, 30:37);
% Published y window is [4 7] (three decades of log removal). Our removals are
% lower, so the window is slid to cover them while keeping the published span of
% three decades and the 0.5 tick spacing (see README).
lo = floor(2*min(allL))/2; hi = max(lo + 3, ceil(2*max(allL))/2);
ylim(ax, [lo hi]); yticks(ax, lo:0.5:hi);
boxedLegend(ax, "northeast");
saveBoth(f, outDir, "rec_FiltrationRate_HETPHO");
fprintf("  filtration: log removal range %.2f .. %.2f\n", min(allL), max(allL));
fprintf("  filtration: %s%s\n", n1, ternary(isempty(p2), "", " + " + n2));
end

% ======================================================================== %
% fig:2d-pat-plots-ref -- March-11/PAT{Flowing,Matrix,Enclosed}Reference.pdf
% ======================================================================== %
function figPat2D(pulseDir, outDir)
fn = pulseFile(pulseDir, "p1x");
if fn == ""
    fprintf("SKIP pat2d: no pulse_fig_p1x*.mat\n"); return
end
w = whos("-file", fn);
if ~ismember("results", string({w.name}))
    fprintf("SKIP pat2d: %s has no ""results"" object\n", fn); return
end
D = load(fn, "results", "rec");
if ~isa(D.results, "Results")
    fprintf("SKIP pat2d: ""results"" in %s did not load as a Results object\n", fn); return
end
z = D.results.SandFilter.GridPoints.Centers(:);
t = D.results.Frames.Time(:).';
% Frames.Time may be relative to the snapshot or already absolute -- shift only
% if it is clearly relative (starts at ~0).
if isfield(D.rec, "tSnapshot") && t(1) < D.rec.tSnapshot
    t = t + D.rec.tSnapshot;
end
for vol = ["Flowing" "Matrix" "Enclosed"]
    C = D.results.Frames.Concentrations{"PAT", vol}{:};
    if isscalar(C) || isempty(C)
        fprintf("SKIP pat2d %s: empty concentration table entry\n", vol); continue
    end
    [T, Z] = meshgrid(t, z);
    f = figure("Visible", "off", "Position", [0 0 900 560], Color="w");
    ax = axes(f);
    surf(ax, T, Z, max(C, realmin), LineStyle="none");
    view(ax, 2); colormap(ax, "jet");
    set(ax, ColorScale="log", YDir="reverse", TickLabelInterpreter="latex", FontSize=14);
    cb = colorbar(ax); cb.Label.String = "Concentration [g/L]";
    cb.Label.Interpreter = "latex"; cb.Label.FontSize = 16;
    cb.TickLabelInterpreter = "latex";
    xlabel(ax, "Time [days]", Interpreter="latex", FontSize=20);
    ylabel(ax, "Depth [m]", Interpreter="latex", FontSize=20);
    xlim(ax, [t(1) t(end)]); ylim(ax, [-1 1]);
    saveBoth(f, outDir, "rec_PAT" + vol + "Reference");
end

% Low-inactivation arm (InactivationRate = 2e-6, vs. NaN/nominal for p1x):
% same three panels, sourced from pulse_fig_exp3_d37.mat.
fn2 = pulseFile(pulseDir, "exp3");
if fn2 == ""
    fprintf("SKIP pat2d-lowinactivation: no pulse_fig_exp3*.mat\n"); return
end
w2 = whos("-file", fn2);
if ~ismember("results", string({w2.name}))
    fprintf("SKIP pat2d-lowinactivation: %s has no ""results"" object\n", fn2); return
end
D2 = load(fn2, "results", "rec");
if ~isa(D2.results, "Results")
    fprintf("SKIP pat2d-lowinactivation: ""results"" in %s did not load as a Results object\n", fn2); return
end
z2 = D2.results.SandFilter.GridPoints.Centers(:);
t2 = D2.results.Frames.Time(:).';
if isfield(D2.rec, "tSnapshot") && t2(1) < D2.rec.tSnapshot
    t2 = t2 + D2.rec.tSnapshot;
end
for vol = ["Flowing" "Matrix" "Enclosed"]
    C = D2.results.Frames.Concentrations{"PAT", vol}{:};
    if isscalar(C) || isempty(C)
        fprintf("SKIP pat2d-lowinactivation %s: empty concentration table entry\n", vol); continue
    end
    [T, Z] = meshgrid(t2, z2);
    f = figure("Visible", "off", "Position", [0 0 900 560], Color="w");
    ax = axes(f);
    surf(ax, T, Z, max(C, realmin), LineStyle="none");
    view(ax, 2); colormap(ax, "jet");
    set(ax, ColorScale="log", YDir="reverse", TickLabelInterpreter="latex", FontSize=14);
    cb = colorbar(ax); cb.Label.String = "Concentration [g/L]";
    cb.Label.Interpreter = "latex"; cb.Label.FontSize = 16;
    cb.TickLabelInterpreter = "latex";
    xlabel(ax, "Time [days]", Interpreter="latex", FontSize=20);
    ylabel(ax, "Depth [m]", Interpreter="latex", FontSize=20);
    xlim(ax, [t2(1) t2(end)]); ylim(ax, [-1 1]);
    saveBoth(f, outDir, "rec_PAT" + vol + "LowInactivation");
end
end

% ======================================================================== %
% CANDIDATE (no published counterpart) -- community structure against depth
% rec_composition_depth.pdf
% ======================================================================== %
function figComposition(dataDir, outDir, opts)
[sc, T] = scenarioSet(dataDir, opts);
if isempty(sc), fprintf("SKIP composition: no fld2x chains on disk\n"); return, end
if isnan(T),    fprintf("SKIP composition: no common horizon\n"); return, end

n = numel(sc);
f = figure("Visible", "off", "Position", [0 0 420*n 760], Color="w");
tl = tiledlayout(f, 2, n, TileSpacing="compact", Padding="compact");
for i = 1:n
    fn = legForTime(dataDir, sc(i).tag, T);
    w  = whos("-file", fn);
    if ~ismember("results", string({w.name}))
        fprintf("SKIP composition: %s has no ""results"" object\n", fn); close(f); return
    end
    D = load(fn, "results");
    ax = nexttile(tl, i);
    set(ax, TickLabelInterpreter="latex", LineWidth=1.0);
    plotBiofilmComposition(D.results, T, AxisHandle=ax, DepthLimits=[0 1], ...
        Phase="Biofilm", FontSize=15, Grid=false);
    title(ax, sc(i).name, Interpreter="latex", FontSize=17);
    xlabel(ax, "Depth $z$ [m]", Interpreter="latex", FontSize=17);
    if i == 1
        ylabel(ax, "Relative volume fraction of biofilm [--]", Interpreter="latex", FontSize=17);
    else
        ylabel(ax, ""); yticklabels(ax, []);
    end
    lg = legend(ax);
    if i < n
        set(lg, Visible="off");
    else
        set(lg, Interpreter="latex", FontSize=14, Box="on", EdgeColor="k", ...
            Orientation="horizontal");
        lg.Layout.Tile = "south";
    end
end
% Bottom row: the ABSOLUTE measure. Relative fractions alone are misleading --
% phi_b peaks near 3e-3 here, so a composition plot describes a nearly empty
% filter. Never show the top row without this one.
phiMax = 0;
for i = 1:n
    A = loadChain(dataDir, sc(i).tag, MaxT=T);
    [phi, tGot] = profileAt(A, T);
    sc(i).phi = phi; sc(i).z = A.z; sc(i).tGot = tGot;
    phiMax = max(phiMax, max(phi(A.z >= 0)));
end
for i = 1:n
    ax = nexttile(tl, n + i);
    set(ax, NextPlot="add", TickLabelInterpreter="latex", FontSize=15, LineWidth=1.0);
    box(ax, "on"); grid(ax, "on");
    plot(ax, sc(i).z, max(sc(i).phi, 1e-6), "k-", LineWidth=1.4);
    % phi_b spans four decades between the schmutzdecke and the bottom of the
    % bed; on a linear axis everything but the z = 0 spike collapses onto zero.
    set(ax, YScale="log");
    xlim(ax, [0 1]); ylim(ax, [1e-5 2*phiMax]);
    xlabel(ax, "Depth $z$ [m]", Interpreter="latex", FontSize=17);
    if i == 1
        ylabel(ax, "Biofilm volume fraction $\phi_{\rm b}$ [--]", Interpreter="latex", FontSize=17);
    else
        yticklabels(ax, []);
    end
    fprintf("  composition %-14s t=%.2f d  max phi_b(bed)=%.3e\n", ...
        sc(i).tag, sc(i).tGot, max(sc(i).phi(sc(i).z >= 0)));
end
title(tl, sprintf("Biofilm community structure at $t = %g$ d (chains \\texttt{fld2x\\_*})", T), ...
    Interpreter="latex", FontSize=19);
saveBoth(f, outDir, "rec_composition_depth");
end

% ======================================================================== %
% CANDIDATE (no published counterpart) -- Tenore2021-style f_i(t)
% rec_composition_time.pdf
% ======================================================================== %
function figTenore(dataDir, outDir, opts)
% Tenore et al. (2021) present phototroph-heterotroph competition as biomass
% FRACTIONS f_i evolving in time (their Figs. 2, 3, 6). The analogue here is the
% bed-integrated composition: what fraction of the sessile biomass in 0 <= z <= 1 m
% is HET / PHO / POM / PAT at each time, per scenario, with the absolute areal
% biomass drawn on top so the reader sees how much biomass the fractions describe.
[sc, T] = scenarioSet(dataDir, opts);
if isempty(sc), fprintf("SKIP tenore: no fld2x chains on disk\n"); return, end
if isnan(T),    fprintf("SKIP tenore: no common horizon\n"); return, end

n = numel(sc);
for i = 1:n
    [sc(i).t, sc(i).areal, sc(i).names] = arealSeries(dataDir, sc(i).tag, T);
end
mMax = max(arrayfun(@(s) max(sum(s.areal, 2)), sc));

f  = figure("Visible", "off", "Position", [0 0 470*ceil(n/2), 760], Color="w");
tl = tiledlayout(f, 2, ceil(n/2), TileSpacing="compact", Padding="compact");
for i = 1:n
    ax = nexttile(tl, i);
    set(ax, NextPlot="add", TickLabelInterpreter="latex", FontSize=15, LineWidth=1.0);
    box(ax, "on");
    tot = sum(sc(i).areal, 2);
    frac = cumsum(sc(i).areal, 2) ./ max(tot, realmin);
    kp = numel(sc(i).names);
    % Cumulative bands drawn largest-first, exactly as plotBiofilmComposition
    % does, so the two candidate figures share colours and legend order.
    for q = 0:kp-1
        area(ax, sc(i).t, frac(:, end-q), DisplayName=sc(i).names(end-q));
    end
    xlim(ax, [0 T]); ylim(ax, [0 1]);
    yyaxis(ax, "right");
    plot(ax, sc(i).t, tot, "k-", LineWidth=2.0, DisplayName="areal biomass");
    ylim(ax, [0 1.05*mMax]);
    set(ax, YColor="k");
    if mod(i, ceil(n/2)) == 0
        ylabel(ax, "Areal biomass [kg\,m$^{-2}$] (wet)", Interpreter="latex", FontSize=17);
    else
        yticklabels(ax, []);
    end
    yyaxis(ax, "left"); set(ax, YColor="k");
    if mod(i, ceil(n/2)) == 1
        ylabel(ax, "Fraction of sessile biomass [--]", Interpreter="latex", FontSize=17);
    else
        yticklabels(ax, []);
    end
    xlabel(ax, "Time [d]", Interpreter="latex", FontSize=17);
    title(ax, sc(i).name, Interpreter="latex", FontSize=17);
    if i == n
        % The bands fill the axes, so an in-axes legend always covers data --
        % park it in the layout's own tile instead.
        lg = legend(ax, "show");
        set(lg, Interpreter="latex", FontSize=14, Box="on", EdgeColor="k", ...
            Orientation="horizontal");
        lg.Layout.Tile = "south";
    end
    fprintf("  tenore %-14s t=%.1f d  areal=%.4f kg/m^2  f=[%s]\n", sc(i).tag, ...
        sc(i).t(end), tot(end), ...
        strjoin(compose("%s %.3f", sc(i).names(:), (sc(i).areal(end, :)./tot(end)).'), ", "));
end
title(tl, "Bed-integrated biofilm composition over time (chains \texttt{fld2x\_*})", ...
    Interpreter="latex", FontSize=19);
saveBoth(f, outDir, "rec_composition_time");
end

% ======================================================================== %
% candidate-figure helpers
% ======================================================================== %
function [sc, T] = scenarioSet(dataDir, opts)
% The four fld2x arms, in the order the manuscript discusses them, restricted to
% the ones actually on disk, plus the LATEST horizon that ALL of them reach.
cand = struct( ...
    "tag",  {opts.LitTag, opts.WinterTag, opts.CovTag, opts.DarkTag}, ...
    "name", {"Summer (19$^\circ$C, uncovered)", "Winter (3$^\circ$C)", ...
             "Covered (1\% light)", "Dark (0\% light)"});
sc = cand([]); T = NaN; ends = [];
for c = cand
    g = chainLegs(dataDir, c.tag);
    if isempty(g)
        fprintf("PARTIAL: no chain_%s_leg*.mat -- arm omitted\n", c.tag); continue
    end
    D = load(fullfile(g(end).folder, g(end).name), "rec");
    sc(end+1) = c; %#ok<AGROW>
    ends(end+1) = D.rec.ts(end); %#ok<AGROW>
end
if ~isempty(ends), T = min(ends); end
end

function fn = legForTime(dataDir, tag, T)
% The leg whose record contains time T (the last one if T is past the chain).
g = chainLegs(dataDir, tag);
fn = string(fullfile(g(end).folder, g(end).name));
for k = 1:numel(g)
    p = fullfile(g(k).folder, g(k).name);
    D = load(p, "rec");
    if D.rec.ts(1) <= T && T <= D.rec.ts(end), fn = string(p); return, end
end
end

function [t, areal, names] = arealSeries(dataDir, tag, T)
% Bed-integrated (z >= 0) sessile areal biomass per particulate, kg/m^2 WET,
% stitched across legs. Sessile = Matrix + Enclosed; Flowing is excluded.
% Same leg bookkeeping as loadChain -- do not diverge.
g = chainLegs(dataDir, tag);
t = []; areal = []; names = strings(0); prevEnd = -inf;
for k = 1:numel(g)
    D = load(fullfile(g(k).folder, g(k).name), "results_py");
    py = D.results_py;
    if isempty(names), names = string(strsplit(py.particleNames, '|')); end
    tk = py.time(:).';
    if k > 1 && abs(tk(1) - 2*prevEnd) < abs(tk(1) - prevEnd), tk = tk - prevEnd; end
    keep = 1:numel(tk); if k > 1, keep = 2:numel(tk); end
    bed = py.z(:) >= 0;
    w   = py.porosity(:) .* bed * py.dz;
    m   = squeeze(sum(w .* (py.particles(:, :, :, 1) + py.particles(:, :, :, 2)), 1));
    t     = [t, tk(keep)];               %#ok<AGROW>
    areal = [areal; m(keep, :)];         %#ok<AGROW>
    prevEnd = tk(end);
    if tk(end) >= T, break, end
end
keepT = t <= T + eps(T);
t = t(keepT).'; areal = areal(keepT, :);
end

% ======================================================================== %
% shared style helpers
% ======================================================================== %
function [f, ax] = squareAxes(o)
arguments, o.Grid (1,1) logical = true; end
f  = figure("Visible", "off", "Position", [0 0 640 620], Color="w");
ax = axes(f, "NextPlot", "add");
axis(ax, "square"); box(ax, "on");
if o.Grid, grid(ax, "on"); end
set(ax, TickLabelInterpreter="latex", FontSize=16, LineWidth=1.0);
end

function [f, ax] = profileAxes()
[f, ax] = squareAxes();
set(ax, YDir="reverse");
xlabel(ax, "Biofilm volume fraction $\phi_{\rm b}$", Interpreter="latex", FontSize=20);
ylabel(ax, "Depth $z$ [m]", Interpreter="latex", FontSize=20);
xlim(ax, [0 0.6]); xticks(ax, 0:0.1:0.6);
end

function boxedLegend(ax, loc)
lg = legend(ax, "show");
set(lg, Interpreter="latex", FontSize=15, Location=loc, Box="on", ...
    EdgeColor="k", Color="w");
end

function saveBoth(f, outDir, name)
exportgraphics(f, fullfile(outDir, name + ".png"), Resolution=200);
exportgraphics(f, fullfile(outDir, name + ".pdf"), ContentType="vector");
close(f);
fprintf("  %s.{png,pdf}\n", name);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

% ======================================================================== %
% data loading
% ======================================================================== %
function g = chainLegs(dataDir, tag)
g = dir(fullfile(dataDir, sprintf("chain_%s_leg*.mat", tag)));
if isempty(g), return, end
g = g(~contains(string({g.name}), "ABORTED"));
if isempty(g), return, end
n = arrayfun(@(x) str2double(regexp(string(x.name), "leg(\d+)\.mat", "tokens", "once")), g);
[~, o] = sort(n); g = g(o);
end

function A = loadChain(dataDir, tag, o)
% Stitch the legs of one chain into a single time series. Mirrors the leg
% bookkeeping of analysis/reproduceManriquez2026.m (do not diverge).
arguments
    dataDir (1,1) string
    tag (1,1) string
    o.MaxT (1,1) double = inf
    o.NeedMass (1,1) logical = false
    o.NeedFields (1,1) logical = false
end
g = chainLegs(dataDir, tag);
needPy = o.NeedMass || o.NeedFields;
liq = ["O2" "IC" "NH4" "HPO4" "DOM"];
A = struct("tag", tag, "z", [], "eps", [], "dz", [], "delta", [], "ts", [], ...
    "phi", [], "eff", [], "effNames", strings(0), "influent", [], ...
    "mass", [], "fields", struct(), "fieldT", []);
prevEnd = -inf;
for k = 1:numel(g)
    if needPy
        D = load(fullfile(g(k).folder, g(k).name), "rec", "results_py");
    else
        D = load(fullfile(g(k).folder, g(k).name), "rec");
    end
    rec = D.rec;
    t = rec.ts(:).';
    if k > 1 && abs(t(1) - 2*prevEnd) < abs(t(1) - prevEnd), t = t - prevEnd; end
    keep = 1:numel(t); if k > 1, keep = 2:numel(t); end
    A.z = rec.z(:); A.eps = rec.eps(:); A.dz = rec.dz; A.delta = rec.delta;
    A.effNames = rec.effNames; A.influent = rec.influent;
    A.ts  = [A.ts,  t(keep)];
    A.phi = [A.phi, rec.phiT(:, keep)];
    A.eff = [A.eff, rec.effluent(:, keep)];
    if needPy && isfield(D, "results_py")
        py = D.results_py;
        if o.NeedMass
            pN = string(strsplit(py.particleNames, '|'));
            tot = 0;
            for nm = ["HET" "PHO" "POM"]
                j = pN == nm;
                tot = tot + py.particles(:, :, j, 1) + py.particles(:, :, j, 2);
            end
            A.mass = [A.mass, sum(py.porosity(:).*tot(:, keep), 1)*py.dz];
        end
        if o.NeedFields
            lN = string(strsplit(py.liquidNames, '|'));
            for nm = liq
                if ~isfield(A.fields, nm), A.fields.(nm) = []; end
                A.fields.(nm) = [A.fields.(nm), py.liquids(:, keep, lN == nm, 2)];
            end
            A.fieldT = [A.fieldT, t(keep)];
        end
    end
    prevEnd = t(end);
    if t(end) >= o.MaxT, break, end
end
end

function j = influentIndex(nm)
% TWO DIFFERENT ORDERINGS live in the same record -- do not mix them up:
%   rec.effluent rows / rec.effNames : ["O2" "IC" "NH4" "HPO4" "DOM" "HET" "PHO" "POM" "PAT"]
%   rec.influent / rec.Influent      : probeChain's InflowConcentrations argument
%                                      order, ["HET" "PHO" "POM" "PAT" "O2" "IC"
%                                      "NH4" "HPO4" "DOM"]
% (Same convention as liquidIndex in analysis/reproduceManriquez2026.m.) Reading
% the influent with the effNames index gives HET = 6.23e-3 instead of 3.0e-4 and
% inflates the HET log removal by ~1.3 decades.
j = find(["HET" "PHO" "POM" "PAT" "O2" "IC" "NH4" "HPO4" "DOM"] == nm);
end

function [p, tGot] = profileAt(A, tWant)
[~, k] = min(abs(A.ts - tWant));
p = A.phi(:, k); tGot = A.ts(k);
end

function fn = pulseFile(pulseDir, stem)
% Prefer the day-30..37 reruns (pulse_fig_<stem>_d37.mat) over the older
% day-30..35 files (pulse_fig_<stem>.mat).
fn = "";
for cand = ["pulse_fig_" + stem + "_d37.mat", "pulse_fig_" + stem + ".mat"]
    p = fullfile(pulseDir, cand);
    if isfile(p), fn = string(p); return, end
end
end

function [r, name] = loadPulse(pulseDir, stem)
r = []; name = "";
fn = pulseFile(pulseDir, stem);
if fn == "", return, end
D = load(fn, "rec"); r = D.rec; [~, b, e] = fileparts(fn); name = string(b) + e;
end
