function out = compare_pathogen(refdir, options)
% COMPARE_PATHOGEN  Golden-master comparison for the PATHOGEN model, MATLAB side.
%
% MATLAB mirror of compare_pathogen.jl: run MPC-SSF's own class-based
% simulate on modelPathogen() and diff it field-by-field against the committed
% reference in REFDIR, which was exported by export_pathogen_reference.m from the
% AUTHORITATIVE slow-sand-filtration code (@SDfilter/run_pathogen.m +
% thesis_model.mat). The reference is codebase-agnostic CSV, so no MATLAB legacy
% checkout is needed here.
%
%   compare_pathogen                      % refdir = ./reference_pathogen
%   compare_pathogen('/abs/path/to/ref')
%   out = compare_pathogen(refdir, Verbose=true, Atol=1e-7, Rtol=1e-6)
%
% OUT is a struct with fields Match, WorstAbs, WorstRel, FlagsAgree, and Table
% (the per-field diff table).
%
% The run spec MUST match export_pathogen_reference.m / compare_pathogen.jl
% exactly. Confounders are removed so the two codebases are numerically
% comparable:
%   * grid:  addGridPoints(20) == SDfilter.add_cells(20) (identical dz/centers)
%   * temp:  20 degC here (SandFilter.Temperature is CELSIUS and computeRate uses
%            theta^(T - 20)) == 293 K there; both give mu = nominal
%   * light: DarkRespiration = 0 and the default diel forcing, which is dark over
%            this short run, so the light factor is trivially equal on both sides
%   * detachment: sqrt(|v|/7.2) on both sides (set by modelPathogen)
%
% To exercise the pathogen physics the run seeds a uniform MATURE BIOFILM (a clean
% filter stays clean, since growth is order-1 in existing biomass). Combined with
% SandPathogen = 0.1 and WaterFactor = 1e-3 this covers both pathogen-specific
% solver knobs, Particle.SandAttachmentFactor and Reaction.EfficiencyFlowing.

    arguments
        refdir string = fullfile(fileparts(mfilename("fullpath")), "reference_pathogen");
        options.Verbose (1,1) logical = true;
        options.Atol (1,1) double = 1e-7;
        options.Rtol (1,1) double = 1e-6;
        options.RelFloor (1,1) double = 1e-12;
        % Variant knobs, mirroring PGOLD_LIGHT_RUN in compare_pathogen.jl. The
        % defaults reproduce reference_pathogen (dark, uniform seed); the light
        % variant needs DarkRespiration=0.1, LightConst=0.8 and a
        % phototroph-heavy seed, which together exercise the dark-respiration
        % light floor max(fdark, I*e^(1-I)).
        options.DarkRespiration (1,1) double = 0.0;
        options.LightConst double = [];   % [] -> default diel forcing
        options.SeedMatrix double = [0.05 0.05 0.02 0.01];
        options.SeedEnclosedLiquids (1,1) double = 1e-3;
    end

    % --- shared run spec (keep in sync with compare_pathogen.jl) -------------
    N       = 20;
    TSIM    = 1e-5;
    DT      = 1e-6;
    NFRAMES = 5;
    %          HET   PHO   POM  PAT   O2    IC    NH4   HPO4 DOM
    INFLOW = [1e-3, 1e-3, 0.0, 1e-5, 1e-2, 1e-2, 1e-5, 0.0, 1e-4] / 10;
    SEED_MATRIX = options.SeedMatrix;         % HET PHO POM PAT
    SEED_ENCL_L = options.SeedEnclosedLiquids; % enclosed liquids
    SEED_WATER  = 0.1;                        % enclosed water volume

    results = runPathogenGolden(N, TSIM, DT, NFRAMES, INFLOW, ...
                                SEED_MATRIX, SEED_ENCL_L, SEED_WATER, ...
                                options.DarkRespiration, options.LightConst);

    % --- reference metadata --------------------------------------------------
    meta = readMeta(fullfile(refdir, "meta.txt"));
    names = split(string(meta.names), ",").';
    flagsAgree = strtrim(string(meta.flag)) == string(results.Flag);

    if options.Verbose
        fprintf("MATLAB(legacy) flag: %s | MPC-SSF flag: %s\n\n", ...
                string(meta.flag), string(results.Flag));
    end

    % --- field-by-field diff -------------------------------------------------
    labels = strings(0, 1);
    maxAbs = zeros(0, 1);
    maxRel = zeros(0, 1);
    okFlag = false(0, 1);

    function add(label, a, b)
        assert(isequal(size(a), size(b)), ...
               "compare_pathogen:shape", "shape mismatch for %s: [%s] vs [%s]", ...
               label, num2str(size(a)), num2str(size(b)));
        absdiff = abs(a - b);
        ma = max([absdiff(:); 0]);
        sig = abs(b) > options.RelFloor;
        if any(sig(:))
            mr = max([absdiff(sig)./abs(b(sig)); 0]);
        else
            mr = 0;
        end
        labels(end+1, 1) = label;
        maxAbs(end+1, 1) = ma;
        maxRel(end+1, 1) = mr;
        okFlag(end+1, 1) = (ma <= options.Atol) || (mr <= options.Rtol);
    end

    depths = results.SandFilter.GridPoints.Centers(:);
    add("depths", depths, readmatrix(fullfile(refdir, "depths.csv")));
    add("times", results.Frames.Time(:), readmatrix(fullfile(refdir, "times.csv")));

    C = results.Frames.Concentrations;
    regions = ["Matrix", "Enclosed", "Flowing"];
    for i = 1:numel(names)
        for r = regions
            ref = readmatrix(fullfile(refdir, sprintf("conc_%s_%s.csv", names(i), lower(r))));
            cij = C{names(i), r};
            val = cij{1};
            if isscalar(val)   % liquids have no matrix region (stored as scalar 0)
                val = zeros(numel(depths), NFRAMES);
            end
            add(names(i) + "/" + lower(r), val, ref);
        end
    end
    cij = C{"Water", "Enclosed"};
    add("Water/enclosed", cij{1}, readmatrix(fullfile(refdir, "conc_Water_enclosed.csv")));

    diffTable = table(labels, maxAbs, maxRel, okFlag, ...
                      VariableNames=["Field", "MaxAbs", "MaxRel", "OK"]);

    if options.Verbose
        fprintf("field                         differences (MPC-SSF MATLAB vs legacy reference)\n");
        verdicts = repmat("FAIL", height(diffTable), 1);
        verdicts(diffTable.OK) = "ok";
        for i = 1:height(diffTable)
            fprintf("  %-28s maxabs=%.3e  maxrel=%.3e  %s\n", ...
                    diffTable.Field(i), diffTable.MaxAbs(i), diffTable.MaxRel(i), verdicts(i));
        end
    end

    out = struct("Match", all(diffTable.OK), ...
                 "WorstAbs", max(diffTable.MaxAbs), ...
                 "WorstRel", max(diffTable.MaxRel), ...
                 "FlagsAgree", flagsAgree, ...
                 "Table", diffTable);

    if options.Verbose
        fprintf("\nworst over all fields: maxabs=%.3e  maxrel=%.3e\n", out.WorstAbs, out.WorstRel);
        if out.Match
            fprintf("RESULT: MATCH -- MPC-SSF MATLAB reproduces the pathogen reference.\n");
        else
            fprintf("RESULT: MISMATCH beyond tolerance.\n");
            disp(diffTable(~diffTable.OK, :));
        end
        if ~flagsAgree
            fprintf("WARNING: flags disagree.\n");
        end
    end
end

function results = runPathogenGolden(N, tsim, dt, nframes, inflow, seedMatrix, seedEnclL, seedWater, darkResp, lightConst)
    % A constant light irradiation replaces the diel forcing when requested, so
    % the light factor is time-invariant and the dark-respiration floor is what
    % the comparison actually probes.
    if isempty(lightConst)
        filter = SandFilter(Temperature=20);   % degC -> theta^0 = 1 (mu = nominal)
    else
        filter = SandFilter(Temperature=20, LightIrradiation=@(t) 0*t + lightConst);
    end
    filter = filter.addGridPoints(N);
    model = modelPathogen(WaterFactor=1e-3, SandPathogen=0.1, DarkRespiration=darkResp);
    state = State(filter, model);

    % Seed the same uniform mature biofilm as the reference: matrix particles +
    % enclosed liquids + enclosed water (enclosed particles stay 0).
    Nc = numel(filter.GridPoints.Centers);
    state.GlobalConcentration.Matrix = repmat(seedMatrix, Nc, 1);
    state.GlobalConcentration.EnclosedLiquids(:) = seedEnclL;
    state.EnclosedWaterVolume(:) = seedWater;

    results = state.simulate(InflowConcentrations=inflow, SimulationTime=tsim, ...
                            TimeStep=dt, FrameNumber=nframes, ...
                            CloggingFraction=0.99, Quiet=true);
end

function meta = readMeta(path)
    meta = struct();
    lines = readlines(path);
    for line = lines.'
        if strlength(strtrim(line)) == 0
            continue
        end
        kv = split(line, "=");
        meta.(kv(1)) = strjoin(kv(2:end), "=");
    end
end
