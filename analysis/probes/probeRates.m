function probeRates(withRespiration)
% Corrected growth rates per Jaime 2026-08-20 (Campos2006 Table 3 upper range):
% mu_HET = 0.042/h * 24 = 1.008/d, mu_PHO = 0.125/h * 24 = 3.0/d.
% withRespiration: false = respiration off; true = PG-excess respiration 0.55/d.
W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
if withRespiration
    m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
    tag = "resp";
else
    m = pathogenModel(NormalizedLight=true);
    tag = "noresp";
end
rx = m.Reactions;
names = [rx.Name];
rx(names == "Heterotroph growth").NominalRate = 0.042*24;   % 1.008 /d
rx(names == "Phototroph growth").NominalRate  = 0.125*24;   % 3.0 /d
m = Model(m.Components, rx, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
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
o2 = C{"O2","Flowing"}{1}*1000;
dom = C{"DOM","Flowing"}{1}*1000;
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name]
    phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP;
end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end
mass = @(nm) sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}, 1)*dz;
hetTot = mass("HET"); phoTot = mass("PHO");
day10 = ts >= 9;
res.tag = tag; res.flag = r.Flag; res.tFinal = r.TimeFinal;
res.o2_out_mean = mean(o2(end, day10)); res.o2_out_range = [min(o2(end, day10)), max(o2(end, day10))];
res.o2_bed_min = min(o2(:, day10), [], "all");
res.dom_out_mean = mean(dom(end, day10));
res.phib_max = max(phiB(:, end)); res.phib_sup = max(phiB(z < 0, end));
res.het_mass_end = hetTot(end); res.pho_mass_end = phoTot(end);
save(fullfile(S, sprintf("rates_%s.mat", tag)), "res", "z", "ts", "o2", "dom", "phiB", "-v7");
disp(res);
fprintf("RATES PROBE %s: flag=%s O2out=%.2f (%.2f-%.2f) bedmin=%.2f DOMout=%.4f phib=%.3f HET=%.3f PHO=%.3f\n", ...
    tag, r.Flag, res.o2_out_mean, res.o2_out_range(1), res.o2_out_range(2), ...
    res.o2_bed_min, res.dom_out_mean, res.phib_max, res.het_mass_end, res.pho_mass_end);
end
