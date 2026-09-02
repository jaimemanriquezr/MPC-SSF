function plotSupernatantHeightChains(opts)
% PLOTSUPERNATANTHEIGHTCHAINS  Supernatant biofilm height against time, one line per
% scenario, for the fld2x chain family.
%
%   plotSupernatantHeightChains()
%   plotSupernatantHeightChains(Reference=0.25)
%
% The measure is the equivalent height of Results/getSupernatantHeight:
%
%       h(t) = ( \int_{z<0} eps phi_b dz ) / phi_ref
%
% i.e. the thickness the biofilm standing above the sand surface would have if
% compacted to volume fraction phi_ref. A threshold height was tried first and
% rejected -- at N = 500 it can only return 0, 2 or 4 mm and separates nothing.
% See .claude/plans/2026-09-02-supernatant-height.md.
%
% Writes rec_supernatant_height.{png,pdf} into analysis/results/figures/recreation/.
% Nothing here starts a simulation; a scenario whose chain is absent is skipped.
%
% The legs are stitched exactly as analysis/recreateFigsManriquez2026.m/loadChain
% does (absolute frame times, drop the repeated first frame of every leg after the
% first). Do not diverge from it.
arguments
    opts.Reference (1,1) double {mustBePositive} = 0.3;   % phi_ref, see above
    opts.SurfaceDepth (1,1) double = 0;                   % z of the sand surface
    opts.Tags (1,:) string = ["fld2x_lit" "fld2x_winter" "fld2x_cov01" "fld2x_dark"];
    opts.Names (1,:) string = ["Summer (19 C, lit)" "Winter (3 C)" ...
                               "Covered, 1\% light" "Covered, dark"];
    opts.Verify (1,1) logical = true;   % cross-check rec.phiT against the @Results method
end
here    = fileparts(mfilename("fullpath"));
repo    = fileparts(here);
addpath(genpath(fullfile(repo, "src")));
dataDir = fullfile(here, "probes", "data", "chain");
outDir  = fullfile(here, "results", "figures", "recreation");
if ~isfolder(outDir), mkdir(outDir); end

sty = ["-" "--" "-." ":"];
lw  = [1.6 2.6 1.6 2.2];

f  = figure("Visible", "off", "Position", [0 0 760 560], Color="w");
ax = axes(f, "NextPlot", "add"); box(ax, "on"); grid(ax, "on");
set(ax, TickLabelInterpreter="latex", FontSize=16, LineWidth=1.0);
xlabel(ax, "Time [d]", Interpreter="latex", FontSize=20);
ylabel(ax, "Supernatant biofilm height $h$ [mm]", Interpreter="latex", FontSize=20);

tMax = 0; drew = 0;
for q = 1:numel(opts.Tags)
    tag = opts.Tags(q);
    g = chainLegs(dataDir, tag);
    if isempty(g)
        fprintf("SKIP %s: no chain_%s_leg*.mat\n", tag, tag); continue
    end
    A = loadChainPhi(g);
    sup = A.z < opts.SurfaceDepth;
    h = sum(A.eps(sup) .* A.phi(sup, :), 1) * A.dz / opts.Reference;
    j = min(q, numel(sty));
    plot(ax, A.ts, 1e3*h, "k" + sty(j), LineWidth=lw(j), ...
        DisplayName=sprintf("%s (%d legs, $t \\le %g$ d)", opts.Names(q), numel(g), A.ts(end)));
    tMax = max(tMax, A.ts(end)); drew = drew + 1;
    fprintf("  %-14s legs %2d  t = %6.2f d  h = %6.3f mm (final)  peak %6.3f mm  h/dz = %.2f\n", ...
        tag, numel(g), A.ts(end), 1e3*h(end), 1e3*max(h), h(end)/A.dz);
    for T = [10 30 60 90]
        if A.ts(end) >= T - 0.05
            [~, k] = min(abs(A.ts - T));
            fprintf("      t = %3d d : h = %6.3f mm\n", T, 1e3*h(k));
        end
    end
    if opts.Verify
        verifyAgainstResults(g(end), opts, 1e3*h(end));
    end
end
if drew == 0
    fprintf("nothing to plot\n"); close(f); return
end

xlim(ax, [0 tMax]);
% The chain tag family and the reference fraction are part of the measure, so they
% are stated on the figure rather than only in the caption.
title(ax, {sprintf("Equivalent height above the sand surface, $\\phi_{\\rm ref} = %g$", opts.Reference), ...
           sprintf("chains \\texttt{chain\\_fld2x\\_*}, horizon $t \\le %g$ d, $N = 500$", tMax)}, ...
    Interpreter="latex", FontSize=14);
lg = legend(ax, "show");
set(lg, Interpreter="latex", FontSize=13, Location="southeast", Box="on", ...
    EdgeColor="k", Color="w");
saveBoth(f, outDir, "rec_supernatant_height");
end

% ======================================================================== %
function verifyAgainstResults(leg, opts, hFromRec)
% Recompute the same number from the saved Results object through the @Results
% method, so the figure and the class method are known to agree.
D = load(fullfile(leg.folder, leg.name), "results");
if ~isfield(D, "results") || ~isa(D.results, "Results")
    fprintf("      verify: no Results object in %s\n", leg.name); return
end
hM = getSupernatantHeight(D.results, ...
    Reference=opts.Reference, SurfaceDepth=opts.SurfaceDepth);
fprintf("      verify: @Results method %.4f mm vs rec path %.4f mm (%s)\n", ...
    1e3*hM(end), hFromRec, leg.name);
end

function g = chainLegs(dataDir, tag)
g = dir(fullfile(dataDir, sprintf("chain_%s_leg*.mat", tag)));
if isempty(g), return, end
g = g(~contains(string({g.name}), "ABORTED"));
if isempty(g), return, end
n = arrayfun(@(x) str2double(regexp(string(x.name), "leg(\d+)\.mat", "tokens", "once")), g);
[~, o] = sort(n); g = g(o);
end

function A = loadChainPhi(g)
A = struct("z", [], "eps", [], "dz", [], "ts", [], "phi", []);
prevEnd = -inf;
for k = 1:numel(g)
    D = load(fullfile(g(k).folder, g(k).name), "rec");
    rec = D.rec;
    t = rec.ts(:).';
    if k > 1 && abs(t(1) - 2*prevEnd) < abs(t(1) - prevEnd), t = t - prevEnd; end
    keep = 1:numel(t); if k > 1, keep = 2:numel(t); end
    A.z = rec.z(:); A.eps = rec.eps(:); A.dz = rec.dz;
    A.ts  = [A.ts,  t(keep)];
    A.phi = [A.phi, rec.phiT(:, keep)];
    prevEnd = t(end);
end
end

function saveBoth(f, outDir, name)
exportgraphics(f, fullfile(outDir, name + ".png"), Resolution=200);
exportgraphics(f, fullfile(outDir, name + ".pdf"), ContentType="vector");
close(f);
fprintf("  %s.{png,pdf}\n", name);
end
