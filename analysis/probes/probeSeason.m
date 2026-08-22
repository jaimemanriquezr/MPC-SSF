function probeSeason(season, respRate)
% season: "summer" (19 C) | "winter" (3 C); respRate: PG-excess respiration /d.
% Corrected growth rates (Campos2006 upper range): mu_HET=1.008/d, mu_PHO=3.0/d.
W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(PhototrophRespiration=respRate, PGExcess=true, NormalizedLight=true);
rx = m.Reactions;
names = [rx.Name];
rx(names == "Heterotroph growth").NominalRate = 0.042*24;
rx(names == "Phototroph growth").NominalRate  = 0.125*24;
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
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name]
    phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP;
end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end
mass = @(nm) sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}, 1)*dz;
lastDay = ts >= max(ts) - 1;
res.season = season; res.respRate = respRate; res.flag = r.Flag; res.tFinal = r.TimeFinal;
res.o2_out_mean = mean(o2(end, lastDay));
res.o2_bed_min = min(o2(:, lastDay), [], "all");
res.phib_max = max(phiB(:, end)); res.phib_sup = max(phiB(z < 0, end));
res.pho_mass_end = mass("PHO"); res.pho_mass_end = res.pho_mass_end(end);
res.het_mass_end = mass("HET"); res.het_mass_end = res.het_mass_end(end);
save(fullfile(S, sprintf("season_%s_r%g.mat", season, respRate)), ...
    "res", "z", "ts", "o2", "phiB", "-v7");
disp(res);
fprintf("SEASON PROBE %s r=%g: flag=%s t=%.1f O2out=%.2f phib=%.3f sup=%.3f PHO=%.3f HET=%.3f\n", ...
    season, respRate, r.Flag, r.TimeFinal, res.o2_out_mean, res.phib_max, ...
    res.phib_sup, res.pho_mass_end, res.het_mass_end);
end
