W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(m.Components, m.Reactions, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
f = SandFilter(Temperature=19);
f = f.addGridPoints(100);
f.LightIrradiation = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=10.0, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=240, ImplicitOsmosis=true, Quiet=true);
C = r.Frames.Concentrations;
ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
o2 = C{"O2","Flowing"}{1};                      % z x frames
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
nf = numel(ts);
globPart = zeros(numel(z), nf);
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name]
    gp = C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1};
    globPart = globPart + gp;
    phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP;
end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end
etaW = f.LightAttenuationEtaWater(:); etaS = f.LightAttenuationEtaSand(:);
save(fullfile(S, "probe_profiles.mat"), "z", "dz", "ts", "o2", "globPart", ...
    "phiB", "etaW", "etaS", "-v7");
fprintf("profiles probe: flag=%s t=%.2f O2out=%.2f mg/L, saved\n", r.Flag, r.TimeFinal, o2(end,end)*1000);
