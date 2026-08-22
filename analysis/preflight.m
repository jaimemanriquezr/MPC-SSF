function ok = preflight(targets)
% PREFLIGHT  Cheap checks to run before submitting anything to the cluster.
%
%   ok = preflight()                       % everything under analysis/ and src/
%   ok = preflight("analysis/probes")      % a subtree
%
% Exists because a defect that costs nothing to catch here cost a 3-hour,
% six-task cosmos array (job 3530219) everything but its exit code: every task
% completed its 90-day run and saved valid results, then died on its final
% summary line and was recorded FAILED. The cause was invisible to checkcode.
%
% Checks:
%   1. checkcode (MATLAB's own lint) on every .m file
%   2. bracket-array format strings -- fprintf(["a" "b"], ...) builds a STRING
%      ARRAY, which fprintf/sprintf/error/warning reject at RUNTIME with
%      "Invalid file identifier". checkcode does not flag it. Use "a" + "b".
%   3. hardcoded absolute paths -- the legacy probes embed a worktree path that
%      no longer resolves, so they only run under an sbatch that supplies its
%      own addpath. New code should derive paths from mfilename("fullpath").
%
% Returns true if all checks pass. Prints every finding.

arguments
    targets (1,:) string = ["analysis", "src"]
end

here = fileparts(mfilename("fullpath"));
W = fileparts(here);
self = string(mfilename("fullpath")) + ".m";   % this file matches its own patterns

files = strings(0);
for t = targets
    d = dir(fullfile(W, t, "**", "*.m"));
    for k = 1:numel(d)
        if contains(d(k).folder, ".git"), continue, end
        f = string(fullfile(d(k).folder, d(k).name));
        if f == self, continue, end            % skip self: the pattern literals below
        files(end+1) = f; %#ok<AGROW>
    end
end
fprintf("preflight: %d files under %s\n", numel(files), strjoin(targets, ", "));

% ---- 1. checkcode -------------------------------------------------------
nLint = 0;
for f = files
    info = checkcode(f, "-struct");
    for k = 1:numel(info)
        fprintf("  LINT  %s:%d  %s\n", relPath(f, W), info(k).line, info(k).message);
        nLint = nLint + 1;
    end
end

% ---- 2. bracket-array format strings ------------------------------------
nFmt = 0;
pat = "(fprintf|sprintf|error|warning)\(\[";
for f = files
    L = readlines(f);
    hits = find(~cellfun(@isempty, regexp(L, pat, "once")));
    for k = hits(:)'
        fprintf("  FORMAT %s:%d  bracket-array format -- concatenate with + , not with []\n", ...
            relPath(f, W), k);
        nFmt = nFmt + 1;
    end
end

% ---- 3. hardcoded absolute paths ----------------------------------------
nPath = 0;
for f = files
    L = readlines(f);
    hits = find(contains(L, "/Users/") | contains(L, "/home/"));
    for k = hits(:)'
        if startsWith(strtrim(L(k)), "%"), continue, end     % comments are fine
        fprintf("  PATH  %s:%d  hardcoded absolute path\n", relPath(f, W), k);
        nPath = nPath + 1;
    end
end

fprintf("preflight: %d lint, %d format, %d path\n", nLint, nFmt, nPath);
ok = (nLint + nFmt + nPath) == 0;
if ok, fprintf("preflight: PASS\n"); else, fprintf("preflight: FAIL\n"); end
end

function p = relPath(f, W)
p = erase(f, string(W) + string(filesep));
end
