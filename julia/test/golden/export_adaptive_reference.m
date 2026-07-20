function export_adaptive_reference(outdir)
% EXPORT_ADAPTIVE_REFERENCE  Run the shared adaptive-CFL golden-master
% simulation in MATLAB and dump results as CSV for the Julia comparer
% (compare_adaptive.jl).
%
%   export_adaptive_reference(OUTDIR) writes reference data into OUTDIR.
%
% Requires the adaptive-CFL simulate.m (branch matlab-claude). The run spec here
% MUST match compare_adaptive.jl:
%   filter  = SandFilter() defaults, addGridPoints(N)
%   model   = modelLund()
%   inflow  = INFLOW (ordered by model.Components, particles-then-liquids)
%   simulate(TimeStep="adaptive", ...) with the options below
%
% Run from the repo root with the matlab-claude MATLAB source checked out, e.g.:
%   git checkout matlab-claude -- src/@State/simulate.m src/@Model/Model.m
%   matlab -batch "initpath; export_adaptive_reference('/path/to/out')"
%   git checkout HEAD -- src/@State/simulate.m src/@Model/Model.m

    if nargin < 1 || isempty(outdir)
        error("export_adaptive_reference:noOutdir", "Provide an output directory.");
    end
    if ~exist(outdir, "dir"); mkdir(outdir); end

    % --- shared run spec (keep in sync with compare_adaptive.jl) ---------
    N       = 20;
    NFRAMES = 5;
    %                HET   PHO   POM  PAT   O2    IC    NH4   HPO4  DOM
    INFLOW = [1e-2, 1e-3, 0.0, 1e-4, 1e-2, 5e-3, 4e-3, 1e-4, 1.0];
    SIMTIME  = 1e-5;
    INIT_DT  = 1e-9;
    MAX_DT   = 1e-6;
    TOL      = 0.05;
    CFL      = 0.99;

    % --- build + run -----------------------------------------------------
    filter = SandFilter();
    filter = filter.addGridPoints(N);
    model  = modelLund();
    state  = State(filter, model);
    results = state.simulate(InflowConcentrations=INFLOW, SimulationTime=SIMTIME, ...
                 TimeStep="adaptive", AdaptiveInitialDt=INIT_DT, AdaptiveMaxDt=MAX_DT, ...
                 AdaptiveTimeTolerance=TOL, CFLFactor=CFL, FrameNumber=NFRAMES, Quiet=true);

    % --- metadata --------------------------------------------------------
    names = [string([model.Particles.Name]), string([model.Liquids.Name])];
    stepTimes = results.SimulationData.time(:);
    writematrix(results.SandFilter.GridPoints.Centers(:), fullfile(outdir, "depths.csv"));
    writematrix(results.Frames.Time(:),                    fullfile(outdir, "times.csv"));
    writematrix(stepTimes,                                 fullfile(outdir, "step_times.csv"));

    meta = fopen(fullfile(outdir, "meta.txt"), "w");
    fprintf(meta, "kP=%d\n", numel(model.Particles));
    fprintf(meta, "kL=%d\n", numel(model.Liquids));
    fprintf(meta, "N=%d\n",  numel(results.SandFilter.GridPoints.Centers));
    fprintf(meta, "nframes=%d\n", NFRAMES);
    fprintf(meta, "nsteps=%d\n", numel(stepTimes) - 1);
    fprintf(meta, "flag=%s\n", string(results.Flag));
    fprintf(meta, "names=%s\n", strjoin(names, ","));
    fclose(meta);

    % --- concentration fields per component/region -----------------------
    C = results.Frames.Concentrations;
    regions = ["Matrix", "Enclosed", "Flowing"];
    for i = 1:numel(names)
        for r = regions
            cell_ij = C{names(i), r};
            val = cell_ij{1};
            if isscalar(val)
                val = zeros(numel(results.SandFilter.GridPoints.Centers), NFRAMES);
            end
            writematrix(val, fullfile(outdir, sprintf("conc_%s_%s.csv", names(i), lower(r))));
        end
    end
    water = C{"Water", "Enclosed"};
    writematrix(water{1}, fullfile(outdir, "conc_Water_enclosed.csv"));

    writematrix(results.Frames.Velocity.Biofilm, fullfile(outdir, "vel_biofilm.csv"));
    writematrix(results.Frames.Velocity.Flowing, fullfile(outdir, "vel_flowing.csv"));

    fprintf("Adaptive reference exported to %s (flag=%s, %d steps)\n", ...
        outdir, string(results.Flag), numel(stepTimes) - 1);
end
