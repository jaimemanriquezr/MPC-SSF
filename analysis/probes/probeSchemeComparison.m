function res = probeSchemeComparison(opts)
% PROBESCHEMECOMPARISON  Three-way comparison of the solver configurations.
%
%   A  Shin  + fully explicit Solver B      (the published/current scheme)
%   B  Shin  + partially implicit Solver B  (implicit dispersion)
%   C  Bailo + partially implicit Solver B
%
% All three run the identical model, mesh, inflow and duration, so cost,
% accuracy, what binds the timestep, and the free-running Solver A bounds are
% directly comparable.
%
% Accuracy is reported PAIRWISE (A-B, B-C, A-C) rather than against a reference:
% a converged reference would cost more than the comparison itself, and the
% pairwise differences already separate the two changes -- A-B isolates the
% dispersion splitting, B-C isolates the Cahn-Hilliard scheme.
%
% The free-running Solver A bounds are the reason C exists. The Shin mobility
% zeta_0*u(1-u) goes negative once u does, turning diffusion into anti-diffusion;
% Bailo's (x)^+ cannot. Whether that is enough to hold u in [0,1] once a reaction
% source is present is exactly what this measures -- Bailo's proof assumes a pure
% conservation law, so it does not transfer for free.

arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
    opts.MaxDt (1,1) double = 1e-2      % loose, so the CFL decides
    opts.Save (1,1) logical = true
end

here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

if opts.Days < 3
    error("probeSchemeComparison:inactiveRegime", ...
        "Days = %g is below the 3 d minimum: no biofilm has formed, so this " + ...
        "measures the startup transient. Phase 0 measured phib_sup = 0 at 0.3 d " + ...
        "against 0.1377 at 3 d.", opts.Days);
end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
inflow = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

cfg = struct( ...
    "tag",    {"A shin+explicit",  "B shin+implicitDisp", "C bailo+implicitDisp"}, ...
    "scheme", {"shin",             "shin",                "bailo"}, ...
    "impdis", {false,              true,                  true});

fprintf("\n=== scheme comparison: N = %d, kappa = %g, %.3g d ===\n\n", ...
    opts.NCells, opts.Kappa, opts.Days);

R = cell(1,3); wall = nan(1,3);
for k = 1:3
    t0 = tic;
    R{k} = simulate(State(f, m), InflowConcentrations=inflow, ...
        SimulationTime=opts.Days, TimeStep="adaptive", AdaptiveInitialDt=1e-8, ...
        AdaptiveMaxDt=opts.MaxDt, FrameNumber=12, ImplicitOsmosis=true, ...
        RecordCflBudget=true, ImplicitDispersion=cfg(k).impdis, ...
        CohesionScheme=cfg(k).scheme, Quiet=true);
    wall(k) = toc(t0);
    fprintf("  %-22s done: %s, %.0f s\n", cfg(k).tag, R{k}.Flag, wall(k));
end

regions = ["matrix","enclosed P","enclosed L","flowing P","flowing L"];
fprintf("\n%-22s %9s %9s %8s %7s  %s\n", "config", "steps", "wall s", "vs A", "flag", "binds");
steps = nan(1,3);
for k = 1:3
    b = R{k}.SimulationData.CflBudget; steps(k) = b.Steps;
    [~, kr] = max(b.RegionWins);
    fprintf("%-22s %9d %9.0f %7.2fx %7s  %s (%.0f%%)\n", cfg(k).tag, steps(k), wall(k), ...
        wall(1)/wall(k), R{k}.Flag, regions(kr), 100*b.RegionWins(kr)/steps(k));
end

fprintf("\n%-22s  %s\n", "binding-region terms", "w_v/dz   w_a/dz2   w_b      w_s");
for k = 1:3
    tf = R{k}.SimulationData.CflBudget.TermSums/steps(k);
    fprintf("%-22s  %-8.3f %-9.3f %-8.3f %-8.3f\n", cfg(k).tag, tf(1), tf(2), tf(3), tf(4));
end

fprintf("\nFREE-RUNNING Solver A (the component-elimination criterion)\n");
fprintf("%-22s %12s %12s %10s %10s\n", "config", "min u", "max u", "1st u<0", "1st NaN");
for k = 1:3
    pd = R{k}.SimulationData.PhibParDiag;
    fprintf("%-22s %12.3e %12.3e %10.4g %10.4g\n", cfg(k).tag, ...
        pd.minSeen, pd.maxSeen, pd.tFirstNeg, pd.tFirstNaN);
end

fprintf("\nPAIRWISE solution differences (max relative, all recorded fields)\n");
pairs = [1 2; 2 3; 1 3]; lbl = ["A-B (dispersion splitting)", "B-C (CH scheme)", "A-C (both)"];
dif = nan(1,3);
for k = 1:3
    dif(k) = maxrel(R{pairs(k,1)}, R{pairs(k,2)});
    fprintf("  %-28s %.4e\n", lbl(k), dif(k));
end

bn = [];
if isfield(R{3}.SimulationData, "BailoNewton")
    bn = R{3}.SimulationData.BailoNewton;
    fprintf("\nBailo semismooth Newton: mean %.2f it/step, max %d, " + ...
        "non-converged %.3f%%, zero-iter %.1f%%\n", ...
        bn.Mean, bn.Max, 100*bn.NonConvergedFrac, 100*bn.ZeroIterFrac);
    if bn.NonConvergedFrac > 0
        fprintf("  *** steps hit the 50-iteration cap: those did NOT solve the " + ...
            "implicit system; treat C as suspect ***\n");
    end
end

res = struct("tags", {[cfg.tag]}, "steps", steps, "wall", wall, ...
    "flags", {[R{1}.Flag, R{2}.Flag, R{3}.Flag]}, "pairwise", dif, ...
    "bailoNewton", bn, "NCells", opts.NCells, "Days", opts.Days);
if opts.Save
    save(fullfile(S, sprintf("schemecmp_n%d_%gd.mat", opts.NCells, opts.Days)), "res", "-v7");
end
end

function d = maxrel(a, b)
Ca = a.Frames.Concentrations; Cb = b.Frames.Concentrations; d = 0;
for nm = string(Ca.Properties.RowNames).'
    for ph = string(Ca.Properties.VariableNames)
        xa = Ca{nm, ph}{1}; xb = Cb{nm, ph}{1};
        if isnumeric(xa) && isequal(size(xa), size(xb)) && ~isscalar(xa)
            den = max(abs(xa(:)));
            if den > 0, d = max(d, max(abs(xa(:) - xb(:)))/den); end
        end
    end
end
end
