W = fileparts(fileparts(mfilename("fullpath")));
addpath(genpath(fullfile(W, "src")));

% Structure: off = original 5 reactions; on = 6, reversed stoichiometry, floor 0.
m0 = modelLund();
assert(length(m0.Reactions) == 5 && m0.Reactions(2).MinimumLightFactor == 0.01);
m = modelLund(PhototrophRespiration=0.276);
assert(length(m.Reactions) == 6);
resp = m.Reactions(6);  gro = m.Reactions(2);
assert(resp.Name == "Phototroph respiration" && ~resp.IsLightDependent);
assert(gro.MinimumLightFactor == 0.0);
for k = ["O2", "IC", "NH4", "HPO4"]
    assert(resp.StoichiometricCoefficients(k) == -gro.StoichiometricCoefficients(k));
end

% Dark-column behavior + anchor export for the Julia comparison.
infl = [0.0, 1.0e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.0e-5, 0.0, 0.0];
function r = darkrun(m, infl)
    f = SandFilter();
    f = f.addGridPoints(15);
    f.LightIrradiation = @(t) 0.0;
    r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=2e-3, ...
        TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
        FrameNumber=3, ImplicitOsmosis=true, Quiet=true);
end
rOn = darkrun(m, infl);  rOff = darkrun(m0, infl);
assert(rOn.Flag == "OK" && rOff.Flag == "OK");
o2tot = @(r) sum(r.Frames.Concentrations{"O2","Flowing"}{1}(:,end)) ...
           + sum(r.Frames.Concentrations{"O2","Enclosed"}{1}(:,end));
photot = @(r) sum(r.Frames.Concentrations{"PHO","Matrix"}{1}(:,end)) ...
            + sum(r.Frames.Concentrations{"PHO","Enclosed"}{1}(:,end)) ...
            + sum(r.Frames.Concentrations{"PHO","Flowing"}{1}(:,end));
assert(o2tot(rOn) < o2tot(rOff), "respiration must consume O2 in the dark");
assert(photot(rOn) < photot(rOff), "respiration must remove biomass");

% Anchor export: final-frame O2 (flowing, enclosed) and PHO (matrix) profiles.
out = [rOn.Frames.Concentrations{"O2","Flowing"}{1}(:,end), ...
       rOn.Frames.Concentrations{"O2","Enclosed"}{1}(:,end), ...
       rOn.Frames.Concentrations{"PHO","Matrix"}{1}(:,end)];
writematrix(out, "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad/anchor_matlab.csv");
fprintf("MATLAB RESPIRATION TESTS PASS  o2on=%.10e o2off=%.10e\n", o2tot(rOn), o2tot(rOff));
