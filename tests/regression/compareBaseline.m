function out = compareBaseline(names)
% COMPAREBASELINE  How FAR does the current code differ from the references?
%
%   out = compareBaseline()          % all configs
%   out = compareBaseline("theta_*")
%
% checkBaseline answers "identical?" -- the right question for a refactor that
% must change nothing. This answers "by how much?", which is the question you
% have when a change is DELIBERATELY reassociating arithmetic: rewriting
%     ((1-eps) + eps*phi_b) * g * b   as   ((1-eps)*b + eps*phi_b*b) * g
% is exact in real arithmetic and inexact in floating point, so the gate goes
% red and the only thing worth knowing is whether the drift is at machine
% epsilon (rounding, accept and re-capture) or above it (a real change, stop).
%
% Reports, per config, the worst RELATIVE difference over every stored array,
% and where it occurs. Judge against eps = 2.2e-16: a few ulps accumulated over
% ~17000 steps is rounding; 1e-9 is not.
%
% This does NOT re-capture anything. Accepting a drift is a deliberate act:
% run captureBaseline yourself, and say in the journal why.

arguments
    names (1,:) string = "*"
end

here = fileparts(mfilename("fullpath"));
addpath(here);
dataDir = fullfile(here, "data");

cfgs = regressionConfigs();
sel = selectConfigs(cfgs, names);
out = struct("name", {}, "maxRel", {}, "where", {}, "maxAbs", {});

fprintf("%-26s %12s  %-28s %12s\n", "config", "max rel", "worst block", "max abs");
fprintf("%s\n", repmat('-', 1, 84));
for k = sel
    cfg = cfgs(k);
    f = fullfile(dataDir, "baseline_" + cfg.name + ".mat");
    if ~isfile(f), fprintf("%-26s   NO REFERENCE\n", cfg.name); continue, end
    ref = load(f).d;
    got = runConfig(cfg);

    worstRel = 0; worstAbs = 0; where = "-";
    if isequaln(ref.Labels, got.Labels)
        for b = 1:numel(ref.Blocks)
            A = ref.Blocks{b}; B = got.Blocks{b};
            if ~isequal(size(A), size(B)), where = "SIZE CHANGED"; worstRel = Inf; break, end
            dA = abs(A - B);
            scale = max(abs(A), abs(B));
            % Relative where there is something to be relative to; absolute
            % differences on an exactly-zero entry are reported separately.
            rel = dA ./ scale;
            rel(scale == 0) = 0;
            [r, i] = max(rel(:));
            if r > worstRel
                worstRel = r;
                [c_, fr] = ind2sub(size(A), i);
                where = sprintf("%s/%s c%d f%d", ref.Labels(b,1), ref.Labels(b,2), c_, fr);
            end
            worstAbs = max(worstAbs, max(dA(:)));
        end
    else
        where = "LAYOUT CHANGED"; worstRel = Inf;
    end

    fprintf("%-26s %12.3e  %-28s %12.3e\n", cfg.name, worstRel, where, worstAbs);
    out(end+1) = struct("name", string(cfg.name), "maxRel", worstRel, ...
        "where", string(where), "maxAbs", worstAbs); %#ok<AGROW>
end

if ~isempty(out)
    m = max([out.maxRel]);
    fprintf("%s\n", repmat('-', 1, 84));
    fprintf("worst relative difference anywhere: %.3e   (eps = %.3e, ratio %.1f ulp)\n", ...
        m, eps, m/eps);
end
end
