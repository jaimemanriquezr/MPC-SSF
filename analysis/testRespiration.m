% Verification of the phototroph metabolic split, Wolf2007 (PHOBIA) r6 form
% with the internal polyglucose pool. Mirrors the Julia testset
% "phototroph respiration (Wolf2007 r6, polyglucose pool)" in runtests.jl.
W = fileparts(fileparts(mfilename("fullpath")));

% Off by default: original 5-reaction, 9-component model.
m0 = modelLund();
assert(length(m0.Reactions) == 5 && length(m0.Components) == 9);
assert(m0.Reactions(2).MinimumLightFactor == 0.01);

% On: PG appended as 5th particle; growth stores f; r6 grows PHO on PG.
f = 0.2;  Y = 0.63;
m = modelLund(PhototrophRespiration=0.55);
pNames = [m.Particles.Name];
assert(isequal(pNames, ["HET", "PHO", "POM", "PAT", "PG"]));
assert(length(m.Reactions) == 6);
gro = m.Reactions(2);
assert(gro.MinimumLightFactor == 0.0);
assert(abs(gro.StoichiometricCoefficients("PG") - f) < 1e-12);
assert(abs(gro.StoichiometricCoefficients("O2") - (0.9301 + 1.0667*f)) < 1e-12);
assert(abs(gro.StoichiometricCoefficients("IC") + (0.36 + 0.4*f)) < 1e-12);
resp = m.Reactions(6);
assert(resp.Name == "Phototroph respiration");
assert(abs(resp.LightInhibition - 1.0) < 1e-12);   % Monod light term, K in I_hat units
assert(abs(resp.HalfSaturationConstants("PG/PHO") - 0.005) < 1e-12);
assert(resp.StoichiometricCoefficients("PG") == -1.0);
assert(abs(resp.StoichiometricCoefficients("PHO") - Y) < 1e-12);          % biomass PRODUCED
assert(abs(resp.StoichiometricCoefficients("NH4") + 0.06*Y) < 1e-12);     % ammonia CONSUMED
assert(abs(resp.StoichiometricCoefficients("O2") + (1.0667 - 0.9301*Y)) < 1e-12);
assert(abs(resp.StoichiometricCoefficients("IC") - (0.4 - 0.36*Y)) < 1e-12);
% COD and carbon close.
assert(abs(1.0667 - (Y*0.9301 + abs(resp.StoichiometricCoefficients("O2")))) < 1e-12);
assert(abs(0.4 - (Y*0.36 + resp.StoichiometricCoefficients("IC"))) < 1e-12);

% Charge/discharge: light builds the PG pool, darkness drains it (r6 is
% dark-only via K/(K+I)).
infl = [0.0, 1.0e-2, 0.0, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.0e-5, 0.0, 0.0];
tswitch = 1.5e-3;
f2 = SandFilter();
f2 = f2.addGridPoints(15);
f2.LightIrradiation = @(t) double(t < tswitch);
r = simulate(State(f2, m), InflowConcentrations=infl, SimulationTime=3e-3, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=7, ImplicitOsmosis=true, Quiet=true);
assert(r.Flag == "OK");
ts = r.Frames.Time;
C = r.Frames.Concentrations;
pgm = sum(C{"PG","Matrix"}{1} + C{"PG","Enclosed"}{1} + C{"PG","Flowing"}{1}, 1);
[~, imid] = min(abs(ts - tswitch));
assert(pgm(imid) > pgm(1), "light must charge the PG pool");
assert(pgm(end) < pgm(imid), "darkness must drain the PG pool");
fprintf("PG pool: charge %.4g -> %.4g, discharge -> %.4g\n", pgm(1), pgm(imid), pgm(end));

% PG-in-excess variant: pool untracked (9 components), respiration = r6 row
% minus the PG column, light factor = complement 1 - Steele.
mx = modelLund(PhototrophRespiration=0.55, PGExcess=true);
assert(length(mx.Components) == 9 && length(mx.Reactions) == 6);
rx = mx.Reactions(6);
assert(abs(rx.LightInhibition - 1.0) < 1e-12);   % Monod light term (complement retired)
assert(~isKey(rx.StoichiometricCoefficients, "PG"));
assert(abs(rx.StoichiometricCoefficients("PHO") - 0.63) < 1e-12);
assert(abs(rx.StoichiometricCoefficients("NH4") + 0.06*0.63) < 1e-12);
assert(abs(rx.StoichiometricCoefficients("O2") + (1.0667 - 0.9301*0.63)) < 1e-12);

% Dark switch on the inhibition factor itself (PG-free micro-model): bright
% light suppresses an inhibited reaction; darkness leaves it on.
phoC = Particle(Name="PHO", Density=1.117e3, Dispersivity=1.2e-2, Transport=5.47);
o2C  = Liquid(Name="O2", Density=998.0, Dispersivity=1.2e-2, Transport=600.0);
respOnly = Reaction(Name="resp", NominalRate=0.55, Order=dictionary("PHO", 1.0), ...
    LightInhibition=8e-5/1.814e-2, ...
    HalfSaturationConstants=dictionary("O2", 3.0e-3), ...
    StoichiometricCoefficients=dictionary(["PHO", "O2"], [-1.0, -0.9301]));
mR = Model([phoC; o2C], respOnly, Kappa=1e-6, Zeta0=1e2, Zeta1=1e-2, DetachmentFunction=@(v) 0*v);
function r = lightrun(mR, lum)
    f = SandFilter();
    f = f.addGridPoints(15);
    f.LightIrradiation = @(t) lum;
    r = simulate(State(f, mR), InflowConcentrations=[1e-2, 9.1e-3], ...
        SimulationTime=2e-3, TimeStep=3e-6, FrameNumber=3, ...
        ImplicitOsmosis=true, Quiet=true);
end
rD = lightrun(mR, 0.0);  rB = lightrun(mR, 1.0);
o2t = @(r) sum(r.Frames.Concentrations{"O2","Flowing"}{1}(:,end)) ...
         + sum(r.Frames.Concentrations{"O2","Enclosed"}{1}(:,end));
assert(rD.Flag == "OK" && rB.Flag == "OK");
assert(o2t(rD) < o2t(rB), "dark must consume more O2 than bright");
fprintf("dark-switch: o2dark=%.10e o2bright=%.10e\n", o2t(rD), o2t(rB));
% Cross-implementation anchor (2026-08-19): the equivalent Julia runs match
% these outputs at golden tolerance; see the anchor check in the session log.
fprintf("MATLAB RESPIRATION (r6+PG) TESTS PASS\n");
