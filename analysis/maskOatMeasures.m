function T = maskOatMeasures(resultsDir, opts)
% MASKOATMEASURES  Resolution-masked OAT measures.
%
% The raw campaign measures integrate the sensitivity s(t) over the whole
% post-disturbance window, including times where L(t) exceeds the solver's
% log-removal resolution (the O(dt) closure bound of
% .claude/decisions/2026-09-01-pat-export-closure-resolved.md: ~2.9 log at the
% operating MaxDt). Pre-breakthrough and deep-tail values of L there are
% numerical noise, and parameters that move the front arrival (dispersivity
% above all) inherit inflated I_rms/I_max from them.
%
% This recomputes I_rms/I_max/D_min from each curves_<param>.csv keeping only
% times where baseline AND both perturbed arms lie inside the defensible band:
%   max(L0, L0+dLp, L0+dLm) <= LCap.
% Writes measures_masked.csv beside the per-parameter dirs and returns the
% combined ranking (with kept-fraction per parameter, so a low-coverage row is
% visibly resting on little data).
arguments
    resultsDir (1,1) string = "analysis/results/oat"
    opts.LCap (1,1) double = 2.9
    opts.OutFile (1,1) string = ""   % default: <resultsDir>/measures_masked.csv
    % Scope to one scenario's dirs when several share resultsDir
    % (log_oat_pulse_*, log_oat_startup_*, ...): the param-dedupe would
    % otherwise silently keep whichever scenario sorts first.
    opts.Pattern (1,1) string = "log_oat_*"
end
d = dir(fullfile(resultsDir, opts.Pattern, "measures.csv"));
assert(~isempty(d), "maskOatMeasures:noMeasures", "no per-parameter measures under %s", resultsDir);
rows = {}; compat = {};
for i = 1:numel(d)
    M = readtable(fullfile(d(i).folder, d(i).name), TextType="string");
    p = M.param(1);
    curves = readmatrix(fullfile(d(i).folder, "curves_" + p + ".csv"));   % t, dLp, dLm, s
    L0 = readmatrix(fullfile(d(i).folder, "L0.csv"));                     % t, L0
    t = curves(:, 1); dLp = curves(:, 2); dLm = curves(:, 3); s = curves(:, 4);
    L0v = interp1(L0(:, 1), L0(:, 2), t, "linear", "extrap");
    mask = max([L0v, L0v + dLp, L0v + dLm], [], 2) <= opts.LCap;
    kept = nnz(mask)/numel(mask);
    if any(mask)
        Irms = sqrt(mean(s(mask).^2));
        Imax = max(abs(s(mask)));
        Dmin = max(min(L0v(mask) + dLp(mask)) - min(L0v(mask)), ...
                   min(L0v(mask) + dLm(mask)) - min(L0v(mask)));
        % denominator of the log-OAT sensitivity, recovered pointwise from
        % s = (dLp - dLm)/denom (constant per parameter by construction)
        nz = mask & abs(s) > 0;
        denom = median((dLp(nz) - dLm(nz))./s(nz));
        Imin = (min(L0v(mask) + dLp(mask)) - min(L0v(mask) + dLm(mask)))/denom;
        dp = dLp(mask); dm = dLm(mask);
        asym = sqrt(sum((dp + dm).^2))/(sqrt(sum((dp - dm).^2)) + 1e-30);
    else
        Irms = NaN; Imax = NaN; Dmin = NaN; Imin = NaN; asym = NaN;
    end
    rows(end+1, :) = {p, M.block(1), Irms, Imax, Dmin, kept, ...
        M.I_rms(1), M.clog_driver(1), M.source(1)}; %#ok<AGROW>
    compat(end+1, :) = {p, M.block(1), Irms, Imax, Dmin, Imin, asym, ...
        M.flag_plus(1), M.flag_minus(1), M.clog_driver(1), M.source(1)}; %#ok<AGROW>
end
T = cell2table(rows, VariableNames=["param", "block", "I_rms_masked", ...
    "I_max_masked", "D_min_masked", "kept_fraction", "I_rms_raw", ...
    "clog_driver", "source"]);
T = sortrows(T, "I_rms_masked", "descend");
outFile = opts.OutFile;
if outFile == "", outFile = fullfile(resultsDir, "measures_masked.csv"); end
writetable(T, outFile);
% a measures.csv-compatible copy, so plotOatTornado/makeOatRankingTable can
% run unchanged on the masked ranking: <resultsDir>/masked/measures.csv
C = cell2table(compat, VariableNames=["param", "block", "I_rms", "I_max", ...
    "D_min", "I_min", "asymmetry", "flag_plus", "flag_minus", "clog_driver", "source"]);
maskedDir = fullfile(resultsDir, "masked");
if ~isfolder(maskedDir), mkdir(maskedDir); end
writetable(sortrows(C, "I_rms", "descend"), fullfile(maskedDir, "measures.csv"));
fprintf("wrote %s (%d parameters, LCap=%.2f) + masked/measures.csv\n", ...
    outFile, height(T), opts.LCap);
disp(T(1:min(10, height(T)), ["param", "I_rms_masked", "I_rms_raw", "kept_fraction"]));
end
