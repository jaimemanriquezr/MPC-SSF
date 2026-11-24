function captureBaseline(names)
% CAPTUREBASELINE  Record the reference output of every regression config.
%
%   captureBaseline()            % all configs
%   captureBaseline("theta_*")   % a subset, by name or wildcard
%
% Writes tests/regression/data/baseline_<name>.mat, one per config. Run it ONCE on the
% commit the cleanup starts from; after that, checkBaseline is the gate and
% this is only re-run for a config whose output changed DELIBERATELY -- in
% which case the journal says which one and why.
%
% Re-capturing to make a red gate go green is the one way to make this whole
% harness worthless. Don't.

arguments
    names (1,:) string = "*"
end

here = fileparts(mfilename("fullpath"));
addpath(here);
outDir = fullfile(here, "data");
if ~isfolder(outDir), mkdir(outDir); end

cfgs = regressionConfigs();
sel = selectConfigs(cfgs, names);

fprintf("Capturing %d configuration(s).\n\n", numel(sel));
for k = sel
    cfg = cfgs(k);
    fprintf("[%2d/%2d] %-26s ", find(sel == k), numel(sel), cfg.name);
    tic;
    d = runConfig(cfg); %#ok<NASGU>
    f = fullfile(outDir, "baseline_" + cfg.name + ".mat");
    save(f, "d", "-v7");
    fprintf("%-6s %6.1f s  ->  %s\n", d.Flag, toc, "baseline_" + cfg.name + ".mat");
end
fprintf("\nDone. %d reference file(s) in %s\n", numel(sel), outDir);
end
