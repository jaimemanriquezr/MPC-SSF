function out = probeRefinement(ncells, bc, opts)
% PROBEREFINEMENT  One (N, boundary-condition) arm of the grid-refinement study.
%
% The 0-cell sits exactly on z = 0, where v_b is undefined (model.tex:148-152), so
% phi_b(0) is not obviously a convergent quantity: if the flux jump acts as a
% surface source, its cell average scales as 1/dz and has no continuum limit. This
% probe measures WHICH outputs converge, and reports every metric without ever
% prescribing a length -- the schmutzdecke thickness must be an output.
%
% METRICS
%   Sigma          eps-weighted areal density over the supernatant [m3/m2].
%                  Convergent by construction (it is the integral of the delta).
%   phiPeak        max phi_b. NOT expected to converge if a surface source forms.
%   L              Sigma/phiPeak, an equivalent thickness with no threshold in it.
%   frontZ         where phi_b crosses zeta_1 -- self-consistent, since that is
%                  where Psi'' changes sign and the layer physics changes character.
%   sharpness      width over which phi_b falls 0.9*peak -> 0.1*peak. This is what
%                  should collapse as zeta_1 rises and a real front forms.
%   zeroCellRatio  phi_b(n0)/phi_b(n0-1). THE DELTA DETECTOR: flat in N means no
%                  surface source; growing like 1/dz means phi_b(0) has no limit.
%   bedDepth       depth at which bed phi_b falls to 10% of its top value.
arguments
    ncells (1,1) double
    bc (1,1) string {mustBeMember(bc, ["neumann","dirichlet"])} = "neumann"
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
    opts.MaxDt (1,1) double = 3e-6
end
if opts.Days < 3
    error("probeRefinement:inactiveRegime", ...
        "Days = %g is below the 3 d minimum: phi_b is ~0 before then, so every " + ...
        "boundary condition and every mesh agree trivially and nothing is measured.", opts.Days);
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis"));

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(ncells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

fprintf("refinement: N = %d, BC = %s, kappa = %g, %.0f d\n", ncells, bc, opts.Kappa, opts.Days);
t0 = tic;
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=12, ImplicitOsmosis=true, CohesionBC=bc, Quiet=true);
wall = toc(t0);
checkRunFlag(r, sprintf("N = %d, BC = %s, %g d", ncells, bc, opts.Days));
fprintf("  %s N=%d  %s  %.0f s\n", bc, ncells, r.Flag, wall);

n0 = f.GridZero;  dz = f.GridSize;
z  = f.GridPoints.Centers;
eps = computePorosity(f, z);
P  = biofilmFrac(r);                            % FULL history, all cells x frames
p  = P(:, end);                                 % final profile
sup = 1:n0;                                      % supernatant + the 0-cell
zeta1 = m.CohesionSubModel.Zeta1;

rec.ncells = ncells;  rec.bc = bc;  rec.dz = dz;  rec.n0 = n0;
rec.z = z;  rec.phi = p;  rec.eps = eps;
% Frame history, so every quantity here can also be plotted against time rather
% than only reported at the final frame.
rec.ts     = r.Frames.Time(:).';
rec.phiT   = P;
rec.bedInt = sum(eps(z > 0).*P(z > 0, :), 1)*dz;
rec.supInt = sum(eps(sup).*P(sup, :), 1)*dz;
rec.bedPeakT = max(P(z > 0, :), [], 1);
rec.supPeakT = max(P(sup, :), [], 1);
rec.flag = r.Flag;  rec.wall = wall;  rec.Kappa = opts.Kappa;  rec.Days = opts.Days;
rec.Sigma   = sum(eps(sup).*p(sup))*dz;
rec.phiPeak = max(p(sup));
rec.L       = rec.Sigma/max(rec.phiPeak, realmin);
rec.zeroCellRatio = p(n0)/max(p(n0-1), realmin);
% front: highest z (deepest into the supernatant) where phi_b crosses zeta_1
ic = find(p(sup) >= zeta1, 1);
rec.frontZ = ternaryLocal(isempty(ic), NaN, z(max(ic,1)));
% sharpness: 0.9*peak -> 0.1*peak, measured on the supernatant profile
i9 = find(p(sup) >= 0.9*rec.phiPeak, 1);
i1 = find(p(sup) >= 0.1*rec.phiPeak, 1);
rec.sharpness = ternaryLocal(isempty(i9) || isempty(i1), NaN, abs(z(i9) - z(i1)));
% bed penetration: below z = 0, where phi_b falls to 10% of the first bed cell
bed = (n0+1):numel(z);
if ~isempty(bed) && p(bed(1)) > 0
    ib = find(p(bed) <= 0.1*p(bed(1)), 1);
    rec.bedDepth = ternaryLocal(isempty(ib), z(bed(end)), z(bed(ib)));
else
    rec.bedDepth = NaN;
end

fprintf("  Sigma %.5e  peak %.5f  L %.5f  ratio %.5f  front %+.4f  sharp %.4f  bed %.4f\n", ...
    rec.Sigma, rec.phiPeak, rec.L, rec.zeroCellRatio, rec.frontZ, rec.sharpness, rec.bedDepth);
fn = fullfile(S, sprintf("refine_n%d_%s_%gd.mat", ncells, bc, opts.Days));
save(fn, "rec", "-v7");
fprintf("  wrote %s\n", fn);
out = rec;
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
