function probeReserve(variant, lightMode, pgInfluent)
% variant: "tracked" | "excess"; lightMode: "summer" | "covered"
% pgInfluent (tracked only): influent PG concentration (algae arrive charged)
if nargin < 3, pgInfluent = 0.0; end
W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(PhototrophRespiration=0.55, PGExcess=(variant == "excess"), ...
    NormalizedLight=true);
m = Model(m.Components, m.Reactions, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
f = SandFilter(Temperature=19);
f = f.addGridPoints(100);
summer = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
if lightMode == "covered", f.LightIrradiation = @(t) 0.01*summer(t);
else, f.LightIrradiation = summer; end
if variant == "tracked"
    infl = [2.68e-3, 1.00e-2, 0.0, 0.0, pgInfluent, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
else
    infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
end
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=10.0, ...
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
phoTot = sum(C{"PHO","Matrix"}{1} + C{"PHO","Enclosed"}{1} + C{"PHO","Flowing"}{1}, 1)*dz;
day10 = ts >= 9;
res.variant = variant; res.light = lightMode; res.flag = r.Flag; res.tFinal = r.TimeFinal;
res.o2_out_mean = mean(o2(end, day10)); res.o2_out_range = [min(o2(end, day10)), max(o2(end, day10))];
res.o2_bed_min = min(o2(:, day10), [], "all");
res.phib_max = max(phiB(:, end)); res.phib_sup = max(phiB(z < 0, end));
res.pho_mass_end = phoTot(end); res.pho_mass_growth = phoTot(end)/phoTot(1);
if variant == "tracked"
    pgTot = sum(C{"PG","Matrix"}{1} + C{"PG","Enclosed"}{1} + C{"PG","Flowing"}{1}, 1)*dz;
    quot = (C{"PG","Matrix"}{1}(:,end) + C{"PG","Enclosed"}{1}(:,end)) ./ ...
           max(C{"PHO","Matrix"}{1}(:,end) + C{"PHO","Enclosed"}{1}(:,end), 1e-30);
    res.pg_mass_end = pgTot(end);
    res.quot_range = [min(quot), max(quot)];
    res.quot_surface = quot(find(z >= 0, 1));
end
tag = ""; if pgInfluent > 0, tag = "_charged"; end
save(fullfile(S, sprintf("reserve_%s_%s%s.mat", variant, lightMode, tag)), ...
    "res", "z", "ts", "o2", "phiB", "-v7");
disp(res);
fprintf("RESERVE PROBE %s/%s: flag=%s O2out=%.2f phib=%.3f phoX=%.2f\n", ...
    variant, lightMode, r.Flag, res.o2_out_mean, res.phib_max, res.pho_mass_growth);
end
