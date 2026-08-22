W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
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
    FrameNumber=10, ImplicitOsmosis=true, Quiet=true);
C = r.Frames.Concentrations;
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
phiB = C{"Water","Enclosed"}{1}(:,end)/densityL;
for nm = [m.Particles.Name], phiB = phiB + (C{nm,"Matrix"}{1}(:,end) + C{nm,"Enclosed"}{1}(:,end))/densityP; end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}(:,end)/densityL; end
z = f.GridPoints.Centers(:);
[mx, i] = max(phiB);
o2out = C{"O2","Flowing"}{1}(end,end);
fprintf("norm probe: flag=%s t=%.2f phib_max=%.3f at z=%.3f  sup_phib_max=%.3f  O2out=%.2f mg/L\n", ...
    r.Flag, r.TimeFinal, mx, z(i), max(phiB(z<0)), o2out*1000);
