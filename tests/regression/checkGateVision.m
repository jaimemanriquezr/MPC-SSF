function blind = checkGateVision()
% CHECKGATEVISION  Does the gate actually SEE the branches it claims to pin?
%
%   blind = checkGateVision()
%
% For every config, revert each option it sets beyond the default production
% path and demand that the output MOVES. A config whose own switch can be
% reverted with no change is blind: it would stay green through a change to the
% very branch it exists to protect.
%
% WHY THIS EXISTS. On 2026-09-07 the suite was green 19/19 while every config
% ran at night -- the phototroph growth, the light attenuation and the whole PG
% pool were unexercised, and respiration_wolf had been re-captured with PG
% identically zero.
% A green gate is not evidence that the gate can see anything. Run this after
% any captureBaseline, and whenever a config is added.
%
% WHAT IT DOES NOT COVER. Configs that differ from the production path by their
% MODEL rather than by a simulate option (pathogen_sandfactor, respiration_*,
% theta_*, dark_window) report "(none: it IS the default)" -- their
% distinguishing feature is not a name-value pair. For those, distinctness
% WITHIN a family is the evidence: the three theta_* baselines differ from each
% other, so ThetaSet demonstrably moves the output. Run-length knobs
% (SimulationTime, FrameNumber) are skipped: reverting SimulationTime to its 1 d
% default turns a 1.5 s config into a 26 min one, and TimeFinal plus the whole
% time grid are in the digest anyway.
%
% KNOWN BLIND, 2026-09-22 (recorded, not fixed):
%   diag_trackcfl/TrackCFL, diag_posmob/TrackCFL -- TrackCFL writes diagnostic
%     fields that runConfig's digest does not capture. These configs still pin
%     the shared time loop, which is what they are needed for.
%   diag_posmob/PositiveMobility -- the clamp (simulateDiagnostic.m:608,
%     lambda = max(lambda,0)) is real but never fires at the seeded state,
%     because lambda is already non-negative there. The branch is unprotected.
% --- gate 1: pairwise distinctness of the stored references -------------------
% Free, no simulation. Two configs with identical output means at least one pins
% nothing. This is also what covers the configs whose distinguishing feature is
% their MODEL rather than a simulate option: the three theta_* differ from one
% another, so ThetaSet demonstrably moves the output.
here = fileparts(mfilename("fullpath"));
d = dir(fullfile(here, "data", "baseline_*.mat"));
dig = cell(1, numel(d)); nm = strings(1, numel(d));
for k = 1:numel(d)
    S = load(fullfile(d(k).folder, d(k).name));
    f = fieldnames(S); dig{k} = S.(f{1});
    nm(k) = extractBetween(string(d(k).name), "baseline_", ".mat");
end
dupes = 0;
for i = 1:numel(d)
    for j = i+1:numel(d)
        if isequaln(dig{i}, dig{j})
            fprintf("IDENTICAL BASELINES: %s == %s\n", nm(i), nm(j)); dupes = dupes + 1;
        end
    end
end
fprintf("%d baselines, %d identical pair(s)\n\n", numel(d), dupes);

% --- gate 2: revert each config's own switch and demand the output move -------
cfgs = regressionConfigs();
ir = find(strcmp(string({cfgs.name}), "shin_neumann_adaptive"), 1);
assert(~isempty(ir), "no reference config");
b = cfgs(ir).build(); refArgs = b{2};
refOpts = string(refArgs(1:2:end));

fprintf("%-24s %-22s %-8s %s\n", "config", "option reverted", "entry", "verdict");
blind = strings(0,1);
for k = 1:numel(cfgs)
    c = cfgs(k);
    bb = c.build(); args = bb{2};
    names = string(args(1:2:end));
    % options this config adds, or sets differently from the production default
    extra = strings(0,1);
    for i = 1:numel(names)
        j = find(refOpts == names(i), 1);
        if isempty(j)
            extra(end+1) = names(i); %#ok<AGROW>
        elseif ~isequaln(args{2*i}, refArgs{2*j})
            extra(end+1) = names(i); %#ok<AGROW>
        end
    end
    % Run-length knobs are excluded: reverting SimulationTime to its 1 d default
    % turns a 1.5 s config into a 26 min one, and its sensitivity is trivial
    % anyway -- TimeFinal and the whole time grid are in the digest.
    extra = extra(~ismember(extra, ["SimulationTime", "FrameNumber"]));
    base = runConfig(c);
    if isempty(extra) && c.entry == "simulate"
        fprintf("%-24s %-22s %-8s %s\n", c.name, "(none: it IS the default)", c.entry, "n/a");
        continue
    end
    for i = 1:numel(extra)
        keep = true(1, numel(args));
        % ALL occurrences: COMMON may set an option that the config then overrides,
        % and MATLAB lets the later pair win -- removing only the first changes nothing.
        for p = find(names == extra(i))
            keep(2*p-1:2*p) = false;
        end
        c2 = c; c2.build = @() {bb{1}, args(keep)};
        try
            d2 = runConfig(c2);
            moved = ~isequaln(rmfield(base,{'name','about'}), rmfield(d2,{'name','about'}));
        catch err
            moved = true;   % reverting it makes the run impossible: also "seen"
            fprintf("%-24s %-22s %-8s ERRORS (%s)\n", c.name, extra(i), c.entry, err.identifier);
            continue
        end
        fprintf("%-24s %-22s %-8s %s\n", c.name, extra(i), c.entry, string(moved));
        if ~moved, blind(end+1) = c.name + "/" + extra(i); end %#ok<AGROW>
    end
end
fprintf("\nBLIND (reverting the option changed nothing): %d\n", numel(blind));
if ~isempty(blind), fprintf("  %s\n", blind); end
end
