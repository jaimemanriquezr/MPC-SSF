function ok = checkBaseline(names)
% CHECKBASELINE  The gate. Re-run every config and demand BIT-IDENTICAL output.
%
%   ok = checkBaseline()            % all configs
%   ok = checkBaseline("theta_*")   % a subset, by name or wildcard
%
% Comparison is isequaln, not a tolerance. This is deliberate: the cleanup it
% guards (the 2026-09 simulation cleanup) changes NO physics and
% NO numerics, so every stored double must come back with the same bits. A
% tolerance would let a reordered floating-point sum through, and a reordered
% sum is exactly the kind of accident a 1285-line refactor produces.
%
% Prints, for each failure, the first block that differs and where -- component,
% phase, frame and cell -- so a red gate points at the defect rather than just
% announcing one.

arguments
    names (1,:) string = "*"
end

here = fileparts(mfilename("fullpath"));
addpath(here);
dataDir = fullfile(here, "data");

cfgs = regressionConfigs();
sel = selectConfigs(cfgs, names);

nPass = 0; nFail = 0; nMissing = 0;
failed = strings(0);

fprintf("Checking %d configuration(s) against tests/regression/data.\n\n", numel(sel));
for k = sel
    cfg = cfgs(k);
    f = fullfile(dataDir, "baseline_" + cfg.name + ".mat");
    fprintf("%-26s ", cfg.name);
    if ~isfile(f)
        fprintf("NO REFERENCE (run captureBaseline)\n");
        nMissing = nMissing + 1;
        continue
    end
    ref = load(f).d;
    got = runConfig(cfg);
    diffs = compareDigests(ref, got);
    if isempty(diffs)
        fprintf("ok    (%s, %d steps, entry=%s)\n", got.Flag, got.NSteps, got.entry);
        nPass = nPass + 1;
    else
        fprintf("FAIL\n");
        for m = 1:numel(diffs)
            fprintf("        %s\n", diffs(m));
        end
        nFail = nFail + 1;
        failed(end+1) = string(cfg.name); %#ok<AGROW>
    end
end

fprintf("\n%d passed, %d failed, %d without a reference.\n", nPass, nFail, nMissing);
if nFail > 0
    fprintf("Failing: %s\n", strjoin(failed, ", "));
end
ok = (nFail == 0) && (nMissing == 0);
end

function diffs = compareDigests(ref, got)
% Scalars and metadata first -- a changed Flag or step count explains every
% array difference below it, so report it first and don't bury it.
diffs = strings(0);

    function chk(field, label)
        if ~isequaln(ref.(field), got.(field))
            diffs(end+1) = sprintf("%s: %s -> %s", label, ...
                fmt(ref.(field)), fmt(got.(field))); %#ok<AGROW>
        end
    end

chk("Flag", "Flag");
chk("NSteps", "step count");
chk("TimeFinal", "TimeFinal");
chk("Solver", "SolverOptions");

if ~isequaln(ref.Time, got.Time)
    diffs(end+1) = "frame times differ";
end
if ~isequaln(ref.TimeGrid, got.TimeGrid)
    n = min(numel(ref.TimeGrid), numel(got.TimeGrid));
    i = find(ref.TimeGrid(1:n) ~= got.TimeGrid(1:n), 1);
    if isempty(i)
        diffs(end+1) = "dt history: same prefix, different length";
    else
        diffs(end+1) = sprintf("dt history first differs at step %d: %.17g -> %.17g", ...
            i, ref.TimeGrid(i), got.TimeGrid(i));
    end
end

% Then the concentration blocks, named.
if ~isequaln(ref.Labels, got.Labels)
    diffs(end+1) = "component/phase layout changed";
    return
end
for b = 1:numel(ref.Blocks)
    A = ref.Blocks{b}; B = got.Blocks{b};
    if isequaln(A, B), continue, end
    lbl = ref.Labels(b, 1) + "/" + ref.Labels(b, 2);
    if ~isequal(size(A), size(B))
        diffs(end+1) = sprintf("%s: size %s -> %s", lbl, ...
            mat2str(size(A)), mat2str(size(B))); %#ok<AGROW>
        continue
    end
    i = find(~(A == B | (isnan(A) & isnan(B))), 1);
    [cell_, frame] = ind2sub(size(A), i);
    diffs(end+1) = sprintf("%s: first differs at cell %d, frame %d: %.17g -> %.17g", ...
        lbl, cell_, frame, A(i), B(i)); %#ok<AGROW>
end

if ~isequaln(ref.VelocityBiofilm, got.VelocityBiofilm)
    diffs(end+1) = "biofilm velocity differs";
end
if ~isequaln(ref.VelocityFlowing, got.VelocityFlowing)
    diffs(end+1) = "flowing velocity differs";
end
end

function s = fmt(v)
if isstring(v) || ischar(v)
    s = string(v);
elseif isscalar(v)
    s = sprintf("%.17g", v);
else
    s = "<" + strjoin(string(size(v)), "x") + ">";
end
end
