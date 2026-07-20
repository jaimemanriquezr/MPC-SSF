function export_pathogen_reference(outdir)
% EXPORT_PATHOGEN_REFERENCE  Run the pathogen golden-master simulation using the
% AUTHORITATIVE slow-sand-filtration code (@SDfilter/run_pathogen.m + the
% thesis_model.mat pathogen model) and dump the frames as CSV for compare_pathogen.jl.
%
% Unlike export_reference.m (which drives MPC-SSF's own class-based simulate),
% the pathogen model exists only in ../slow-sand-filtration. This script must be
% run with THAT repo's classes on the path. From the MPC-SSF repo root:
%
%   matlab -batch "run('julia/test/golden/run_export_pathogen.m')"
%
% (run_export_pathogen.m cds into slow-sand-filtration, calls initpath, then this.)
%
% The run spec is chosen to remove every cross-codebase confounder so the two
% solvers are numerically comparable (see compare_pathogen.jl for the mirror):
%   * grid:  SDfilter.add_cells(N) == Julia addgridpoints(N) (identical dz/centers)
%   * temp:  293 K -> theta^(293-T)=theta^0=1, i.e. mu = nominal (Julia uses 20 C,
%            which gives theta^(293/293-1)=1 as well)
%   * light: dark_respiration = 0 makes MATLAB max(fdark, I*exp(1-I)) equal Julia
%            max(0, min_light + I*exp(1-I)) since the light term is >= 0
%   * detachment: set explicitly to sqrt(|v|/7.2) to match modelPathogen()

    if nargin < 1 || isempty(outdir)
        error("export_pathogen_reference:noOutdir", "Provide an output directory.");
    end
    if ~exist(outdir, "dir"); mkdir(outdir); end

    % --- shared run spec (keep in sync with compare_pathogen.jl) ---------
    N       = 20;
    TSIM    = 1e-5;
    DT      = 1e-6;
    NFRAMES = 5;

    % --- filter (293 K removes the temperature correction) ---------------
    filter = SDfilter(temperature = 293);
    filter = filter.add_cells(N);

    % --- pathogen model, confounders removed -----------------------------
    model = load("thesis_model.mat").model;
    model.ecological.water_factor   = 1e-3;
    model.ecological.sand_pathogen  = 0.1;   % exercise differential PAT->sand attachment
    model.ecological.dark_respiration = 0.0;
    model.detachment = @(v) sqrt(abs(v) / 7.2);

    % inflow order: [HET PHO POM PAT | O2 IC NH4 HPO4 DOM]
    inflow = [1e-3, 1e-3, 0, 1e-5, 1e-2, 1e-2, 1e-5, 0, 1e-4] / 10;

    % Seed a uniform mature biofilm so the biofilm-phase reactions (death,
    % hydrolysis, pathogen inactivation, bacterivory) actually fire — a clean
    % filter stays clean (growth is order-1 in existing biomass). Layout of the
    % biofilm block: [matrix P (1:4) | enclosed P (5:8) | enclosed L (9:13)].
    Nc = numel(filter.mesh.cell_centers);
    ic = struct();
    ic.biofilm = zeros(Nc, 13);
    ic.biofilm(:, 1:4)  = repmat([0.05, 0.05, 0.02, 0.01], Nc, 1);  % HET PHO POM PAT
    ic.biofilm(:, 9:13) = 1e-3;                                     % enclosed liquids
    ic.enclosed_water = 0.1 * ones(Nc, 1);

    results = run_pathogen(filter, model, ...
        time_step = DT, time_run = TSIM, n_frames = NFRAMES, ...
        inflow_concentrations = inflow, initial_conditions = ic);

    % --- metadata --------------------------------------------------------
    pnames = string(model.particles.Name(:)).';
    lnames = string(model.liquids.Name(:)).';
    names  = [pnames, lnames];
    centers = filter.mesh.cell_centers;

    writematrix(centers(:),               fullfile(outdir, "depths.csv"));
    writematrix(results.frames.time(:),   fullfile(outdir, "times.csv"));

    meta = fopen(fullfile(outdir, "meta.txt"), "w");
    fprintf(meta, "kP=%d\n", numel(pnames));
    fprintf(meta, "kL=%d\n", numel(lnames));
    fprintf(meta, "N=%d\n",  numel(centers));
    fprintf(meta, "nframes=%d\n", NFRAMES);
    fprintf(meta, "flag=%s\n", string(results.flag));
    fprintf(meta, "names=%s\n", strjoin(names, ","));
    fclose(meta);

    % --- concentration fields (depth x frame) per component/region -------
    C = results.frames.concentrations;   % rows=component, cols matrix/enclosed/flowing
    regions = ["matrix", "enclosed", "flowing"];
    for i = 1:numel(names)
        for r = regions
            val = C{names(i), r}{1};
            if isscalar(val)   % liquids have no matrix region (stored as scalar 0)
                val = zeros(numel(centers), NFRAMES);
            end
            fname = sprintf("conc_%s_%s.csv", names(i), r);
            writematrix(val, fullfile(outdir, fname));
        end
    end

    water = C{"WATER", "enclosed"}{1};
    writematrix(water, fullfile(outdir, "conc_Water_enclosed.csv"));

    writematrix(results.frames.velocity.biofilm, fullfile(outdir, "vel_biofilm.csv"));
    writematrix(results.frames.velocity.flowing, fullfile(outdir, "vel_flowing.csv"));

    fprintf("Pathogen reference exported to %s (flag=%s)\n", outdir, string(results.flag));
end
