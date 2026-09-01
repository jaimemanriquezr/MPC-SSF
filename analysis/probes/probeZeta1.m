function out = probeZeta1(zeta1, opts)
% PROBEZETA1  One arm of the zeta_1 sweep.
%
% zeta_1 sets the equilibrium biofilm fraction (the single well of Psi sits at
% 3*zeta_1/2) and the spinodal (Psi'' = 3u(u - zeta_1) changes sign at u = zeta_1).
% At zeta_1 = 0.01 the well is at 0.015 while phi_b reaches ~0.32, so the potential
% is in its SPREADING branch everywhere and the layer smears (29 cm half-max at 20 d,
% N=100). Raising zeta_1 puts the dilute tail inside the spinodal band, where the
% potential is cohesive, which should produce a real front at an emergent height.
%
% zeta_1 MUST be set by reconstructing the Model: CahnHilliardModel rebuilds
% PotentialGradient from Zeta1 in its constructor, so assigning the field alone
% changes the number without changing the physics.
%
% SCORING TARGET (Campos2002, Demir2017): biomass in the SAND, with the top 2 cm
% holding ~50% of bed biomass. Today's baseline gives 10-16%.
arguments
    zeta1 (1,1) double
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
    opts.Frames (1,1) double = 12
    opts.MaxDt (1,1) double = 3e-6
    opts.Scheme (1,1) string = "shin"
    opts.ImplicitDispersion (1,1) logical = true
    % Detachment scaling. Baseline handle is 0.14*sqrt(|v|/18); DetachScale
    % multiplies the coefficient. 0 switches detachment off entirely, which is the
    % discriminating test for whether the bed peak-and-decline is driven by
    % surface erosion (Leg A) or by the supernatant layer starving the bed (Leg B).
    opts.DetachScale (1,1) double = 1.0
end
if opts.Days < 3
    error("probeZeta1:inactiveRegime", "Days = %g < 3: no biofilm has formed.", opts.Days);
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis"));

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=zeta1, ...
    DetachmentFunction=@(v) opts.DetachScale*0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
assert(abs(m.CohesionSubModel.Zeta1 - zeta1) < 1e-15, "Zeta1 did not take");
% Verify the potential regenerated, rather than trusting the constructor.
uchk = linspace(1e-3, 0.999, 51).';
assert(max(abs(m.CohesionSubModel.PotentialGradient(uchk) - uchk.^2.*(uchk - 1.5*zeta1))) < 1e-12, ...
    "PotentialGradient did not regenerate for Zeta1 = %g", zeta1);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

fprintf("zeta_1 = %g | well %.4f | spinodal %.4f | %s impDisp=%d detach=%.2g | N=%d, %.0f d, MaxDt=%g\n", ...
    zeta1, 1.5*zeta1, zeta1, opts.Scheme, opts.ImplicitDispersion, opts.DetachScale, opts.NCells, opts.Days, opts.MaxDt);
t0 = tic;
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=opts.Frames, ImplicitOsmosis=true, ImplicitDispersion=opts.ImplicitDispersion, ...
    CohesionScheme=opts.Scheme, Quiet=true);
wall = toc(t0);
checkRunFlag(r, sprintf("zeta_1 = %g, N = %d, %g d, MaxDt = %g", zeta1, opts.NCells, opts.Days, opts.MaxDt));

z = f.GridPoints.Centers; dz = f.GridSize; n0 = f.GridZero;
eps = computePorosity(f, z);
P = biofilmFrac(r); p = P(:, end); w = eps(:).*p;
% The sand bed starts at z = -delta (SandRoughness), not z = 0: that is where
% epsilon first drops below 1, i.e. where sand grains are first present. So bed
% mass integrates from -delta, and Campos2002's 0-2 cm sand core is the 2 cm
% BELOW -delta. The previous "bed = z > 0" was wrong twice over -- it started at
% the wrong depth AND excluded the cell centred exactly at z = 0, which carries
% the largest phi_b in the profile.
%
% At N = 100 this changes almost nothing, because dz = 1.00 cm > delta = 0.50 cm
% and NO cell falls in [-delta, 0). At N = 500 it halves the supernatant fraction
% (9.3% -> 4.5%). The bed/supernatant split is not resolvable at N = 100.
delta = f.SandRoughness;
bed = z >= -delta; sup = z < -delta; zb = z(bed); b = w(bed);
rec.zeta1 = zeta1; rec.z = z; rec.phi = p; rec.flag = r.Flag; rec.wall = wall;
rec.scheme = opts.Scheme; rec.impDisp = opts.ImplicitDispersion;
rec.NCells = opts.NCells; rec.Days = opts.Days; rec.Kappa = opts.Kappa;
rec.Sigma    = sum(w(sup))*dz;
rec.phiPeak  = max(p);
rec.bedPeak  = max(p(bed));
rec.supPeak  = max(p(sup));
rec.supFrac  = sum(w(sup))/max(sum(w), realmin);
rec.top2Frac = sum(b(zb <= -delta + 0.02))/max(sum(b), realmin);
rec.L        = rec.Sigma/max(rec.supPeak, realmin);
rec.zeroCellRatio = p(n0)/max(p(n0-1), realmin);
% TRAJECTORY. The final profile cannot say whether the bed peaks and then
% declines -- the failure mode that breaks Manriquez2026 Fig 9(a) and contradicts
% Campos2002's monotone logistic rise. Record the eps-weighted split every frame.
rec.ts     = r.Frames.Time(:).';
% Full profile history, not just the final frame. Needed to draw the
% Manriquez2026 Fig 9(a) style plot (one curve per day) -- the integrals alone
% cannot reconstruct the shape, which is why the first 20 d runs had to be redone.
rec.phiT   = P;
rec.bedInt = sum(eps(bed).*P(bed, :), 1)*dz;
rec.supInt = sum(eps(sup).*P(sup, :), 1)*dz;
rec.bedPeakT = max(P(bed, :), [], 1);
% EFFLUENT time series -- the quantity a plant operator and a referee both care
% about, and the only model output directly comparable with Bellamy/Bae/Elemo
% removal data. Taken from the deepest cell of the FLOWING phase.
C = r.Frames.Concentrations;
lNames = [r.Model.Liquids.Name];  pNames = [r.Model.Particles.Name];
rec.effNames = [lNames, pNames];
rec.effluent = zeros(numel(rec.effNames), numel(rec.ts));
for q = 1:numel(lNames)
    a = C{lNames(q), "Flowing"}{1};  rec.effluent(q, :) = a(end, :);
end
for q = 1:numel(pNames)
    a = C{pNames(q), "Flowing"}{1};  rec.effluent(numel(lNames)+q, :) = a(end, :);
end
rec.DetachScale = opts.DetachScale;
rec.influent = infl;
[~, kpk] = max(rec.bedInt);
rec.tBedPeak = rec.ts(kpk);
rec.bedDecline = rec.bedInt(kpk)/max(rec.bedInt(end), realmin);
fprintf("  bed peaks at t = %.1f d, then falls %.2fx by t = %.0f d\n", ...
    rec.tBedPeak, rec.bedDecline, rec.ts(end));
i9 = find(p(sup) >= 0.9*rec.supPeak, 1); i1 = find(p(sup) >= 0.1*rec.supPeak, 1);
if isempty(i9) || isempty(i1), rec.sharpness = NaN; else, rec.sharpness = abs(z(i9)-z(i1)); end
fprintf("  %s %.0fs | bedPeak %.4f supPeak %.4f | supFrac %5.1f%% | top2/bed %5.1f%% | sharp %.4f | ratio %.4f\n", ...
    r.Flag, wall, rec.bedPeak, rec.supPeak, 100*rec.supFrac, 100*rec.top2Frac, rec.sharpness, rec.zeroCellRatio);
% Filename carries MaxDt: two runs differing only in the cap are different
% results, and writing both to one name has already destroyed good data twice
% today. (Note MATLAB's %g prints 0.40 as "0.4" -- match that when globbing.)
fn = fullfile(S, sprintf("zeta1_%g_n%d_%gd_dt%g_det%g.mat", zeta1, opts.NCells, opts.Days, opts.MaxDt, opts.DetachScale));
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
