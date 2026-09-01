function out = probeRegressionBisect(arm, opts)
% PROBEREGRESSIONBISECT  Which recent change made the biofilm leave the sand?
%
% Published Manriquez2026 Fig 9(a) at t = 20 d: biofilm IN the bed, phi_b ~ 0.55 at
% z = 0 decaying over ~10-20 cm, with a thin (~2 cm) supernatant layer. The current
% code inverts this -- bed 0.278 / supernatant 0.325 at 20 d, reaching 96.9% of all
% biomass in the supernatant by 90 d. Campos2002 and Demir2017 both put the biomass
% in the sand (top 1.5-2 cm), so the published behaviour is the correct one.
%
% zeta_1 cannot explain it: the published value was 0.005-0.01, the same as now.
% Three defaults changed in the last five days, each feeding mu and hence v_b:
%   kappa0    kappa was INERT before 2026-08-23 (block (2,1) was empty, so the solve
%             reduced to mu = Psi'(u) exactly). kappa = 0 reproduces that.
%   noupwind  IsUpwinded flipped false -> true on 2026-08-24.
%   rawlight  the normalised-light convention arrived 2026-08-20.
%
% SCORE against the published figure: bedPeak should be ~0.55 and supFrac small.
arguments
    arm (1,1) string {mustBeMember(arm, ["baseline","kappa0","noupwind","rawlight"])}
    opts.NCells (1,1) double = 100
    opts.Days (1,1) double = 20
    opts.MaxDt (1,1) double = 3e-6
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis"));

kappa = 1e-6;  upwind = true;  normLight = true;
switch arm
    case "kappa0",   kappa = 0;
    case "noupwind", upwind = false;
    case "rawlight", normLight = false;
end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=normLight);
m = Model(mp.Components, mp.Reactions, Kappa=kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

fprintf("bisect arm %s: kappa=%g upwind=%d normLight=%d, N=%d, %.0f d\n", ...
    arm, kappa, upwind, normLight, opts.NCells, opts.Days);
t0 = tic;
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=20, ImplicitOsmosis=true, IsUpwinded=upwind, Quiet=true);
wall = toc(t0);
checkRunFlag(r, sprintf("arm %s, N = %d, %g d", arm, opts.NCells, opts.Days));
fprintf("  %-9s %s  %.0f s\n", arm, r.Flag, wall);

z = f.GridPoints.Centers;
eps = computePorosity(f, z);
P = biofilmFrac(r); ts = r.Frames.Time(:).';
rec.arm = arm; rec.kappa = kappa; rec.upwind = upwind; rec.normLight = normLight;
rec.z = z; rec.ts = ts; rec.phi = P; rec.eps = eps; rec.flag = r.Flag; rec.wall = wall;
rec.NCells = opts.NCells; rec.Days = opts.Days;
bed = z > 0; sup = z <= 0;
w = eps(:).*P;
rec.bedPeak = max(P(bed, end));
rec.supPeak = max(P(sup, end));
rec.supFrac = sum(w(sup, end))/max(sum(w(:, end)), realmin);
zb = z(bed); b = w(bed, end);
rec.top2Frac = sum(b(zb <= 0.02))/max(sum(b), realmin);
fprintf("  t=%.0f d: bedPeak %.4f  supPeak %.4f  supFrac %.1f%%  top2/bed %.1f%%\n", ...
    ts(end), rec.bedPeak, rec.supPeak, 100*rec.supFrac, 100*rec.top2Frac);
fprintf("  PUBLISHED Fig 9(a) target: bedPeak ~0.55, supPeak ~0.15, supFrac small\n");
save(fullfile(S, sprintf("bisect_%s_n%d_%gd.mat", arm, opts.NCells, opts.Days)), "rec", "-v7");
out = rec;
end

function p = biofilmFrac(r)
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
p = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], p = p + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   p = p + C{nm,"Enclosed"}{1}/dL; end
end
