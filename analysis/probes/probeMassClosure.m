function out = probeMassClosure(opts)
% PROBEMASSCLOSURE  Zero-biology particle-mass closure with eps-weighted stocks.
%
% Closes the open question of 2026-08-21c ("transport creates particle mass"),
% which was explained as an eps-unweighted probe integral (JOURNAL 2026-08-22)
% but never re-run. simulate.m updates every stored concentration as
%   c += (dt/dz)*(epsFace_in*F_in - epsFace_out*F_out)/epsCenter + dt*rhs,
% so the conserved stock is  M = sum_z epsCenter * c * dz  summed over phases, and
% with all reaction rates zero  dM/dt = q*(c_in - c_out)  must hold for every
% particle species to solver tolerance. Attachment/detachment/transfer only move
% mass between phases of the same cell and cancel in the sum.
arguments
    opts.NCells (1,1) double = 100
    opts.Days (1,1) double = 10
    opts.MaxDt (1,1) double = 5e-5
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
rx = mp.Reactions;
for k = 1:numel(rx), rx(k).NominalRate = 0.0; end
m = Model(mp.Components, rx, Kappa=1e-6, Zeta0=5, Zeta1=0.27, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), ...
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, ...
    OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveMaxDt=opts.MaxDt, ImplicitOsmosis=true, ...
    FrameNumber=241, Quiet=true);
checkRunFlag(r, "probeMassClosure");

C = r.Frames.Concentrations; ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
eps = computePorosity(f, z);
q = f.InflowVelocity;
names = [m.Particles.Name];
out = struct("tag", "massclosure", "days", opts.Days, "NCells", opts.NCells, "species", names);
fprintf("MASS CLOSURE (zero biology, N=%d, %g d, q=%g):\n", opts.NCells, opts.Days, q);
fprintf("%6s %12s %12s %12s %12s %10s\n", "sp", "dM(eps)", "supply", "export", "residual", "res/supply");
for k = 1:numel(names)
    nm = names(k);
    cin = infl(k);
    total = C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1};   % N x frames
    M = sum(eps.*total, 1)*dz;                                                % eps-weighted stock
    cout = C{nm,"Flowing"}{1}(end, :);
    supply = q*cin*(ts(end) - ts(1));
    export = q*trapz(ts, cout(:));
    dM = M(end) - M(1);
    resid = dM - (supply - export);
    out.(nm) = struct("dM", dM, "supply", supply, "export", export, "residual", resid, ...
        "residualOverSupply", resid/max(supply, 1e-300), "stock", M, "cout", cout);
    fprintf("%6s %12.5g %12.5g %12.5g %12.3e %10.2e\n", nm, dM, supply, export, resid, resid/max(supply,1e-300));
end
out.ts = ts; out.eps = eps; out.z = z; out.SolverOptions = r.SolverOptions;
save(fullfile(S, "massclosure_n" + opts.NCells + "_" + opts.Days + "d.mat"), "out", "-v7.3");
end
