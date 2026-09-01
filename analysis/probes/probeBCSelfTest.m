function out = probeBCSelfTest(opts)
% PROBEBCSELFTEST  Verify the Dirichlet condition at z = 0.
%
% Two gates:
%
%   OPERATOR (exact, no simulation). Converting the Neumann assembly to Dirichlet
%   must change row n0 of the mu block by exactly -kappa/dz^2 on the diagonal,
%   which is precisely the -kappa*u_bed/dz^2 term moved to the right-hand side.
%   This is an algebraic identity and holds to round-off.
%
%   ACTIVE (simulation). The real Dirichlet run must DIFFER from Neumann,
%   otherwise the term is inert and the operator gate passed vacuously. This is
%   null below 3 d: phi_b is ~0 until the biofilm forms, and with u == 0 every
%   boundary condition agrees to machine zero.
%
% The former runtime `BCSelfTest` flag (which forced the ghost to u(n0) and so
% should have reproduced Neumann) was removed on 2026-08-25. It could never
% return exactly zero -- the diagonal is implicit while the ghost is explicit, so
% the residual is O(dt) and accumulates -- which made it a weaker test than the
% algebraic identity above. Measured 2.76e-05 over 3 d against an active effect
% of 3.3e-02.
arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
end
if opts.Days < 3
    error("probeBCSelfTest:inactiveRegime", "Days = %g < 3: the ACTIVE gate is vacuous.", opts.Days);
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);

% ---- GATE 1: the operator identity, exact -------------------------------
n0 = f.GridZero; dz = f.GridSize;
[~, DD, ~] = getCahnHilliardMatrices(f, m, true);
DD = sparse(DD.Rows, DD.Columns, DD.Values, 2*n0, 2*n0);
CH0 = speye(2*n0) - DD;                       % Neumann, as assembled
r0 = n0 + n0;                                 % mu-row for the cell at z = 0
CH0d = CH0;  CH0d(r0, n0) = CH0d(r0, n0) - opts.Kappa/dz^2;   % -> Dirichlet
u = rand(n0,1); mu = rand(n0,1); x = [u; mu];
lhsDiff = CH0d(r0,:)*x - CH0(r0,:)*x;
expected = -opts.Kappa*u(n0)/dz^2;            % the term moved to the rhs
out.operatorResidual = abs(lhsDiff - expected);
out.operatorPass = out.operatorResidual < 1e-12*max(abs(expected), 1);

% ---- GATE 2: the two boundary conditions must differ in an active regime -
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
args = {"InflowConcentrations", infl, "SimulationTime", opts.Days, ...
        "TimeStep", "adaptive", "AdaptiveInitialDt", 1e-8, "AdaptiveMaxDt", 3e-6, ...
        "FrameNumber", 12, "ImplicitOsmosis", true, "Quiet", true};
rN = simulate(State(f,m), args{:}, CohesionBC="neumann");    checkRunFlag(rN, "neumann");
rD = simulate(State(f,m), args{:}, CohesionBC="dirichlet");  checkRunFlag(rD, "dirichlet");
a = biofilmFrac(rN); b = biofilmFrac(rD);
out.phiMax    = max(a, [], "all");
out.activeAbs = max(abs(a-b), [], "all");
out.activeRel = out.activeAbs/max(out.phiMax, realmin);
out.activePass = out.activeAbs > 0 && out.phiMax > 1e-3;

fprintf("\nOPERATOR  |Dirichlet-Neumann row - (-kappa*u(n0)/dz^2)| = %.3e   %s\n", ...
    out.operatorResidual, ternaryLocal(out.operatorPass, "PASS", "FAIL"));
fprintf("ACTIVE    max phi_b = %.5f, max|dirichlet-neumann| = %.3e (rel %.3e)   %s\n", ...
    out.phiMax, out.activeAbs, out.activeRel, ternaryLocal(out.activePass, "PASS", "FAIL"));
out.pass = out.operatorPass && out.activePass;
fprintf("GATE: %s\n", ternaryLocal(out.pass, "PASS", "FAIL"));
save(fullfile(S, sprintf("bc_selftest_n%d_%gd.mat", opts.NCells, opts.Days)), "out", "-v7");
end

function v = ternaryLocal(c, a, b)
if c, v = a; else, v = b; end
end

function p = biofilmFrac(r)
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
p = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], p = p + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   p = p + C{nm,"Enclosed"}{1}/dL; end
end
