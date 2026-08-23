function res = probeImplicitDispersion(opts)
% PROBEIMPLICITDISPERSION  Does implicit dispersion deliver the predicted dt, and
% at what cost in accuracy?
%
% Phase 0 measured dispersion at 94% of the binding region's sum at N=500, with
% the matrix region's margin at 0.0688, predicting AT MOST ~14.5x on dt before
% advection takes over at dz/q. Two things can eat that:
%   - the dispersion coefficient and phiFlowing are LAGGED, so the splitting has
%     its own accuracy limit independent of stability;
%   - once dispersion leaves the CFL, whatever binds next may bind sooner than
%     the margin suggests.
% Both are measured here rather than assumed.
%
% Check A (accuracy): explicit vs implicit at the SAME small dt. Both are first
%   order in the split, so they should agree to O(dt) and converge together.
% Check B (timestep): raise MaxDt until the CFL binds, and compare the dt each
%   scheme sustains, plus steps and wall time for the same simulated span.

arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 0.5
    opts.MaxDt (1,1) double = 1e-2      % deliberately loose: let the CFL decide
    opts.Save (1,1) logical = true
end

here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
inflow = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

run = @(imp, maxdt) runOne(f, m, inflow, opts.Days, maxdt, imp);

fprintf("\n=== implicit dispersion, N = %d, kappa = %g, %.3g d ===\n", ...
    opts.NCells, opts.Kappa, opts.Days);

% ---- Check A: same small dt, do they agree? --------------------------------
small = 3e-6;
fprintf("\nA. accuracy at a common MaxDt = %g (both CFL-free)\n", small);
[re, te] = run(false, small);
[ri, ti] = run(true,  small);
dA = compareRuns(re, ri);
fprintf("   explicit flag %s (%.0f s), implicit flag %s (%.0f s)\n", re.Flag, te, ri.Flag, ti);
fprintf("   max rel difference over all recorded fields: %.4e\n", dA);

% ---- Check B: let the CFL decide -------------------------------------------
fprintf("\nB. timestep, MaxDt raised to %g so the CFL binds\n", opts.MaxDt);
[rE, tE] = run(false, opts.MaxDt);
[rI, tI] = run(true,  opts.MaxDt);
sE = rE.SimulationData.CflBudget.Steps;  sI = rI.SimulationData.CflBudget.Steps;
capE = rE.SimulationData.CflBudget.CapBound/sE;
capI = rI.SimulationData.CflBudget.CapBound/sI;
fprintf("   explicit: %8d steps, %6.0f s, cap-bound %5.1f%%, flag %s\n", sE, tE, 100*capE, rE.Flag);
fprintf("   implicit: %8d steps, %6.0f s, cap-bound %5.1f%%, flag %s\n", sI, tI, 100*capI, rI.Flag);
fprintf("   STEP REDUCTION %.2fx   (Phase 0 predicted at most 14.5x)\n", sE/sI);
fprintf("   wall-time change %.2fx  (a tridiagonal solve per component per step)\n", tE/tI);
regions = ["matrix","enclosed P","enclosed L","flowing P","flowing L"];
bI = rI.SimulationData.CflBudget;
fprintf("   what binds now: ");
[~, kk] = max(bI.RegionWins); fprintf("%s (%.0f%% of steps)\n", regions(kk), 100*bI.RegionWins(kk)/sI);
tf = bI.TermSums/sI;
fprintf("   binding-region shares: w_v/dz %.2f  w_a/dz2 %.2f  w_b %.2f  w_s %.2f\n", tf(1), tf(2), tf(3), tf(4));

res = struct("accuracyRelDiff", dA, "stepsExplicit", sE, "stepsImplicit", sI, ...
    "stepReduction", sE/sI, "wallExplicit", tE, "wallImplicit", tI, ...
    "capExplicit", capE, "capImplicit", capI, "termFracImplicit", tf, ...
    "NCells", opts.NCells, "Days", opts.Days);
if opts.Save
    save(fullfile(S, sprintf("impdisp_n%d_%gd.mat", opts.NCells, opts.Days)), "res", "-v7");
end
end

function [r, wall] = runOne(f, m, inflow, days, maxdt, implicitDisp)
t0 = tic;
r = simulate(State(f, m), InflowConcentrations=inflow, SimulationTime=days, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=maxdt, ...
    FrameNumber=6, ImplicitOsmosis=true, RecordCflBudget=true, ...
    ImplicitDispersion=implicitDisp, Quiet=true);
wall = toc(t0);
end

function d = compareRuns(a, b)
Ca = a.Frames.Concentrations; Cb = b.Frames.Concentrations;
d = 0;
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
