function out = probeSolverADrift(scheme, opts)
% PROBESOLVERADRIFT  One scheme's free-running Solver A drift, saved to data/.
%
% Runner half of plotSolverADrift: it runs ONE cohesion scheme and writes the
% arrays the plot needs. Split out so the three schemes can go to Slurm as an
% array (one task each) instead of running in sequence inside the plotter --
% they are completely independent, so the wall clock is the slowest scheme
% rather than their sum.
%
% Configuration is IDENTICAL to plotSolverADrift's (same filter, model, inflow,
% MaxDt, FrameNumber), because the point is to compare schemes against each
% other and against the existing 3 d / 20 d drift numbers.
%
% Writes data/solverA_drift_<scheme>_n<NCells>.mat holding one struct `run`.

arguments
    scheme (1,1) string {mustBeMember(scheme, ["shin", "matched", "bailo"])}
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 20
    opts.MaxDt (1,1) double = 3e-6
end
if opts.Days < 3
    error("probeSolverADrift:inactiveRegime", ...
        "Days = %g is below the 3 d minimum: no biofilm has formed, so this " + ...
        "measures the startup transient rather than the drift.", opts.Days);
end

here = fileparts(mfilename("fullpath"));            % .../analysis/probes
W = fileparts(fileparts(here));                      % repo root
S = fullfile(here, "data");
if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

fprintf("scheme %s: N = %d, kappa = %g, %.0f d, MaxDt = %g\n", ...
    scheme, opts.NCells, opts.Kappa, opts.Days, opts.MaxDt);
t0 = tic;
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=72, ImplicitOsmosis=true, TrackCFL=true, ...
    CohesionScheme=scheme, Quiet=true);
wall = toc(t0);
checkRunFlag(r, sprintf("scheme %s, N = %d, %g d", scheme, opts.NCells, opts.Days));
fprintf("  %-7s %s  %.0f s\n", scheme, r.Flag, wall);

n0 = f.GridZero;
rec.scheme = scheme;
rec.z      = f.GridPoints.Centers(1:n0);
rec.phiA   = r.SimulationData.PhibPar;                  % free-running Solver A
phiS       = biofilmFrac(r);
rec.phiS   = phiS(1:n0, :);                             % authoritative sum
rec.phiC   = r.SimulationData.PhibCH;                   % re-seeded, one-step
% PhibCH is written from frame 1; PhibPar only once the free state is seeded, so
% mask both by the free state's (narrower) validity window.
ok         = ~isnan(rec.phiA(1,:));
rec.tt     = r.Frames.Time(ok).';
rec.phiA   = rec.phiA(:,ok);
rec.phiS   = rec.phiS(:,ok);
rec.phiC   = rec.phiC(:,ok);
rec.D      = max(abs(rec.phiA - rec.phiS), [], 1);      % drift vs time
rec.flag   = r.Flag;
rec.wall   = wall;
rec.steps  = r.SimulationData.CflBudget.Steps;
rec.parDiag = r.SimulationData.PhibParDiag;
if isfield(r.SimulationData, "BailoNewton")
    rec.bailoNewton = r.SimulationData.BailoNewton;
end
rec.NCells = opts.NCells;  rec.Kappa = opts.Kappa;  rec.Days = opts.Days;

fn = fullfile(S, sprintf("solverA_drift_%s_n%d.mat", scheme, opts.NCells));
save(fn, "rec", "-v7");
fprintf("  wrote %s\n", fn);
out = rec;
end

function p = biofilmFrac(r)
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
p = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], p = p + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   p = p + C{nm,"Enclosed"}{1}/dL; end
end
