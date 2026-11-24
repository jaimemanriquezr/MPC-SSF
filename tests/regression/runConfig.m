function d = runConfig(cfg)
% RUNCONFIG  Run one regression configuration and reduce it to a digest.
%
%   d = runConfig(cfg)   % cfg is one element of regressionConfigs()
%
% The digest holds the run's full numeric output rather than a hash. A hash
% answers "did something change?"; the arrays answer "which component, which
% phase, which frame" -- and that is the question you actually have at 2 a.m.
% after a refactor. At 30 cells and 6 frames the whole suite is ~1 MB.
%
% `entry` selects the function under test. Before Step 1 of the cleanup plan,
% `simulateDiagnostic` does not exist and every config runs through `simulate`;
% after it, the diag_* configs run through the fork and must reproduce the
% references captured here. Resolving it by `exist` (rather than hard-coding)
% is what lets the same reference files span the fork.

entry = cfg.entry;
if entry == "simulateDiagnostic" && isempty(which("simulateDiagnostic"))
    entry = "simulate";     % pre-fork: the diagnostics still live in simulate
end

built = cfg.build();
state = built{1};
args  = built{2};

fn = str2func(entry);
r = fn(state, args{:});

d = struct();
d.name      = string(cfg.name);
d.about     = string(cfg.about);
d.entry     = string(entry);
d.Flag      = string(r.Flag);
d.TimeStart = r.TimeStart;
d.TimeFinal = r.TimeFinal;
d.Time      = r.Frames.Time;
d.Solver    = string(r.SolverOptions.summary());

% Concentrations is a table of cells, one row per component and Water, three
% phase columns. Flatten to a labelled numeric array so isequaln can compare it
% and a mismatch can be named. Empty-phase entries are stored as the scalars
% simulate writes there (0 for the absent matrix phase, Inf for flowing water).
C = r.Frames.Concentrations;
rows = string(C.Properties.RowNames);
cols = string(C.Properties.VariableNames);
d.Labels = strings(0, 2);
d.Blocks = {};
for i = 1:numel(rows)
    for j = 1:numel(cols)
        d.Labels(end+1, :) = [rows(i), cols(j)]; %#ok<AGROW>
        d.Blocks{end+1} = C{rows(i), cols(j)}{1}; %#ok<AGROW>
    end
end

d.VelocityBiofilm = r.Frames.Velocity.Biofilm;
d.VelocityFlowing = r.Frames.Velocity.Flowing;

% Step count and the adaptive dt history: a refactor that leaves the frames
% intact but changes how many steps produced them has still changed the solver,
% and this is the cheapest place to notice.
d.NSteps   = numel(r.SimulationData.time);
d.TimeGrid = r.SimulationData.time(:);
end
