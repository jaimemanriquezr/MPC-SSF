function export_reference(outdir)
% EXPORT_REFERENCE  Run the shared golden-master simulation in MATLAB and dump
% the results as plain CSV files for the Julia comparer (compare.jl).
%
%   export_reference(OUTDIR) writes reference data into OUTDIR.
%
% The run spec here MUST be kept identical to compare.jl and to the Julia
% simpleModel() preset:
%   filter  = SandFilter() defaults, addGridPoints(N)
%   model   = SimpleModel.mat with example overrides (Kappa=1e-2, Zeta0=1.0,
%             detachment @(v) sqrt(abs(v)))
%   inflow  = Microorganism=1e-2, Nutrient=1.0
%   simulate(SimulationTime=TSIM, TimeStep=DT, FrameNumber=NFRAMES, Quiet=true)
%
% Run from the repo root, e.g.:
%   matlab -batch "initpath; export_reference('/path/to/out')"

    if nargin < 1 || isempty(outdir)
        error("export_reference:noOutdir", "Provide an output directory.");
    end
    if ~exist(outdir, "dir"); mkdir(outdir); end

    % --- shared run spec (keep in sync with compare.jl) ------------------
    N       = 20;
    TSIM    = 1e-3;
    DT      = 1e-5;
    NFRAMES = 5;

    % --- build + run -----------------------------------------------------
    filter = SandFilter();
    filter = filter.addGridPoints(N);
    model  = load("./data/SimpleModel.mat").model;
    model.CohesionSubModel.Kappa = 1e-2;
    model.CohesionSubModel.Zeta0 = 1.0;
    model.DetachmentFunction = @(v) sqrt(abs(v));
    inflow = dictionary("Microorganism", 1e-2, "Nutrient", 1.0);

    state = State(filter, model);
    results = state.simulate(InflowConcentrations=inflow, SimulationTime=TSIM, ...
                             TimeStep=DT, FrameNumber=NFRAMES, Quiet=true);

    % --- dump metadata ---------------------------------------------------
    names  = [string([model.Particles.Name]), string([model.Liquids.Name])];
    writematrix(results.SandFilter.GridPoints.Centers(:), fullfile(outdir, "depths.csv"));
    writematrix(results.Frames.Time(:),                    fullfile(outdir, "times.csv"));

    meta = fopen(fullfile(outdir, "meta.txt"), "w");
    fprintf(meta, "kP=%d\n", numel(model.Particles));
    fprintf(meta, "kL=%d\n", numel(model.Liquids));
    fprintf(meta, "N=%d\n",  numel(results.SandFilter.GridPoints.Centers));
    fprintf(meta, "nframes=%d\n", NFRAMES);
    fprintf(meta, "flag=%s\n", string(results.Flag));
    fprintf(meta, "names=%s\n", strjoin(names, ","));
    fclose(meta);

    % --- dump concentration fields (depth x frame) per component/region --
    C = results.Frames.Concentrations;   % table, rows=component, cols=Matrix/Enclosed/Flowing
    regions = ["Matrix", "Enclosed", "Flowing"];
    for i = 1:numel(names)
        for r = regions
            cell_ij = C{names(i), r};
            val = cell_ij{1};
            if isscalar(val)   % liquids Matrix = 0 (stored as scalar 0)
                val = zeros(numel(results.SandFilter.GridPoints.Centers), NFRAMES);
            end
            fname = sprintf("conc_%s_%s.csv", names(i), lower(r));
            writematrix(val, fullfile(outdir, fname));
        end
    end

    % enclosed water
    water = C{"Water", "Enclosed"};
    writematrix(water{1}, fullfile(outdir, "conc_Water_enclosed.csv"));

    % velocities (depth-1 x frame)
    writematrix(results.Frames.Velocity.Biofilm, fullfile(outdir, "vel_biofilm.csv"));
    writematrix(results.Frames.Velocity.Flowing, fullfile(outdir, "vel_flowing.csv"));

    fprintf("Reference exported to %s (flag=%s)\n", outdir, string(results.Flag));
end
