function probeEndog(season, withHet)
% Campos-style endogenous respiration probe (biomass self-oxidation, sign-changed):
%   PHO endog: k_ra = 0.276/d (Campos2006 avg, Brown & Barnwell), theta 1.08
%   HET endog: k_rb = 1.72/d  (Campos2006 avg, Bowie et al.),      theta 1.047
% Stoich (both, per unit biomass; reverse of phototroph growth row):
%   X -1, O2 -0.9301, IC +0.36, NH4 +0.06, HPO4 +0.01
% Constant day/night (Campos k_ra bundles excretion/photorespiration); O2-Monod protected.
% Baseline: corrected growth rates (mu_HET=1.008, mu_PHO=3.0), old half-sats,
% PG-excess respiration OFF, normalized light. 30 d, 100 cells (matches probeSeason).
% season: "summer" | "winter"; withHet: false = PHO only, true = PHO + HET.
W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(NormalizedLight=true);
rx = m.Reactions;
names = [rx.Name];
rx(names == "Heterotroph growth").NominalRate = 0.042*24;
rx(names == "Phototroph growth").NominalRate  = 0.125*24;
endogStoich = @(X) dictionary( ...
    [X, "O2", "IC", "NH4", "HPO4"], [-1.0, -0.9301, 0.36, 0.06, 0.01]);
phoEndog = Reaction(Name="Phototroph endogenous respiration", ...
    NominalRate=0.276, TemperatureCorrectionFactor=1.08, ...
    Order=dictionary("PHO", 1), ...
    HalfSaturationConstants=dictionary("O2", 3.00E-03), ...
    StoichiometricCoefficients=endogStoich("PHO"));
rx = [rx; phoEndog];
tag = sprintf("%s_pho", season);
if withHet
    hetEndog = Reaction(Name="Heterotroph endogenous respiration", ...
        NominalRate=1.72, TemperatureCorrectionFactor=1.047, ...
        Order=dictionary("HET", 1), ...
        HalfSaturationConstants=dictionary("O2", 3.00E-03), ...
        StoichiometricCoefficients=endogStoich("HET"));
    rx = [rx; hetEndog];
    tag = sprintf("%s_phohet", season);
end
m = Model(m.Components, rx, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
if season == "summer"
    f = SandFilter(Temperature=19);
    light = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
else
    f = SandFilter(Temperature=3);
    light = @(t) max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
end
f = f.addGridPoints(100);
f.LightIrradiation = light;
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=30.0, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=240, ImplicitOsmosis=true, Quiet=true);
C = r.Frames.Concentrations;
ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
o2 = C{"O2","Flowing"}{1}*1000;
nh4 = C{"NH4","Flowing"}{1}*1000;
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name]
    phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP;
end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end
mass = @(nm) sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}, 1)*dz;
lastDay = ts >= max(ts) - 1;
res.tag = tag; res.flag = r.Flag; res.tFinal = r.TimeFinal;
res.o2_out_mean = mean(o2(end, lastDay));
res.o2_out_range = [min(o2(end, lastDay)), max(o2(end, lastDay))];
res.o2_bed_min = min(o2(:, lastDay), [], "all");
res.nh4_out_mean = mean(nh4(end, lastDay));
res.phib_max = max(phiB(:, end)); res.phib_sup = max(phiB(z < 0, end));
res.pho_mass_end = mass("PHO"); res.pho_mass_end = res.pho_mass_end(end);
res.het_mass_end = mass("HET"); res.het_mass_end = res.het_mass_end(end);
save(fullfile(S, sprintf("endog_%s.mat", tag)), "res", "z", "ts", "o2", "nh4", "phiB", "-v7");
disp(res);
fprintf("ENDOG PROBE %s: flag=%s t=%.1f O2out=%.2f (%.2f-%.2f) bedmin=%.2f NH4=%.4f phib=%.3f sup=%.3f PHO=%.3f HET=%.3f\n", ...
    tag, r.Flag, r.TimeFinal, res.o2_out_mean, res.o2_out_range(1), res.o2_out_range(2), ...
    res.o2_bed_min, res.nh4_out_mean, res.phib_max, res.phib_sup, res.pho_mass_end, res.het_mass_end);
end
