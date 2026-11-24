% TESTASYMMETRYOPTIONS  Unit checks for the five additive modelLund/pathogenModel
% temperature-asymmetry options: ThetaDHET, ThetaDPHO,
% HETCardinal, Grazing, GrazingTheta. All defaults are inert (NaN/NaN/
% [NaN NaN NaN]/0.0/1.08), so (i) checks that the no-argument call is unchanged.
%
% Runs locally; the smoke simulate() calls in (vi) are ~0.01 d at
% N = 50.

W = fileparts(fileparts(fileparts(mfilename("fullpath"))));   % tests/<suite>/<file> -> repo root
addpath(fullfile(W, "src"));

%% ------------------------------------------------------------------ (i)
fprintf("=== (i) modelLund() with no options is unchanged ===\n");
lund = mpcssf.presets.modelLund();
assert(numel(lund.Reactions) == 5, "expected 5 reactions, got %d", numel(lund.Reactions));
expectedNames = ["Heterotroph growth", "Phototroph growth", "Heterotroph death", ...
                  "Phototroph death", "Hydrolysis"];
assert(isequal([lund.Reactions.Name], expectedNames), "reaction names/order changed");
expectedRate  = [4.80, 2.00, 0.40, 0.276, 3.00];
expectedTheta = [1.0725, 1.047, 1.0725, 1.047, 1.0725];
assert(all(abs([lund.Reactions.NominalRate] - expectedRate) < 1e-12), "NominalRate changed");
assert(all(abs([lund.Reactions.TemperatureCorrectionFactor] - expectedTheta) < 1e-12), ...
    "TemperatureCorrectionFactor changed");
assert(lund.Reactions(1).TemperatureResponse == "exponential", "HET growth response changed");
fprintf("  PASS\n");

%% ----------------------------------------------------------------- (ii)
fprintf("=== (ii) ThetaDHET / ThetaDPHO override only their own reaction ===\n");
lDH = mpcssf.presets.modelLund(ThetaDHET=1.0);
thetaDH = [lDH.Reactions.TemperatureCorrectionFactor];
assert(abs(thetaDH(3) - 1.0) < 1e-12, "ThetaDHET did not reach heterotrophDeath");
assert(isequal(thetaDH([1 2 4 5]), expectedTheta([1 2 4 5])), "ThetaDHET touched another reaction");

lDP = mpcssf.presets.modelLund(ThetaDPHO=1.0);
thetaDP = [lDP.Reactions.TemperatureCorrectionFactor];
assert(abs(thetaDP(4) - 1.0) < 1e-12, "ThetaDPHO did not reach phototrophDeath");
assert(isequal(thetaDP([1 2 3 5]), expectedTheta([1 2 3 5])), "ThetaDPHO touched another reaction");
fprintf("  PASS\n");

%% ---------------------------------------------------------------- (iii)
fprintf("=== (iii) HETCardinal switches heterotroph GROWTH only ===\n");
lC = mpcssf.presets.modelLund(HETCardinal=[-5 37 45]);
rxG = lC.Reactions(1); rxD = lC.Reactions(3);
assert(rxG.TemperatureResponse == "cardinal", "heterotrophGrowth not switched to cardinal");
assert(isequal(rxG.CardinalTemperatures, [-5 37 45]), "CardinalTemperatures not set");
assert(abs(rxG.computeRate(20) - 4.8) < 1e-12, "cardinal anchor mu(20) failed");
assert(rxG.computeRate(2) < 4.8*1.0725^(2 - 20), "cardinal did not suppress mu(2) below exponential");
assert(rxD.TemperatureResponse == "exponential", "heterotrophDeath response changed by HETCardinal");
fprintf("  PASS\n");

%% ----------------------------------------------------------------- (iv)
fprintf("=== (iv) Grazing appends a sixth reaction ===\n");
lG = mpcssf.presets.modelLund(Grazing=0.8);
assert(numel(lG.Reactions) == 6, "expected 6 reactions with Grazing on, got %d", numel(lG.Reactions));
gz = lG.Reactions(6);
assert(gz.Name == "Heterotroph grazing", "sixth reaction misnamed: %s", gz.Name);
assert(abs(gz.NominalRate - 0.8) < 1e-12, "Grazing rate not forwarded");
assert(abs(gz.TemperatureCorrectionFactor - 1.08) < 1e-12, "GrazingTheta default not 1.08");
sc = gz.StoichiometricCoefficients;
assert(abs(sc("HET") - (-1.0)) < 1e-12, "HET coefficient wrong");
assert(abs(sc("POM") - 0.4*0.9123) < 1e-12, "POM coefficient wrong");
assert(abs(sc("IC") - 0.6*0.3600) < 1e-12, "IC coefficient wrong");
assert(abs(sc("O2") - (-0.6*1.0667)) < 1e-12, "O2 coefficient wrong");
assert(abs(sc("NH4") - 0.0653) < 1e-12, "NH4 coefficient wrong");
assert(abs(sc("HPO4") - 0.0209) < 1e-12, "HPO4 coefficient wrong");
hs = gz.HalfSaturationConstants;
assert(abs(hs("HET") - 2.00e-3) < 1e-15, "HET half-saturation wrong");

lG0 = mpcssf.presets.modelLund(Grazing=0.8, GrazingTheta=1.0);
assert(abs(lG0.Reactions(6).TemperatureCorrectionFactor - 1.0) < 1e-12, "GrazingTheta override not forwarded");
fprintf("  PASS\n");

%% ------------------------------------------------------------------ (v)
fprintf("=== (v) pathogenModel forwards all five options ===\n");
pm = mpcssf.presets.pathogenModel(Grazing=0.8, ThetaDHET=1.0, HETCardinal=[-5 37 45]);
rx = pm.Reactions;
assert(any([rx.Name] == "Heterotroph grazing"), "grazing row missing from pathogenModel");
hd = rx([rx.Name] == "Heterotroph death");
assert(abs(hd.TemperatureCorrectionFactor - 1.0) < 1e-12, "ThetaDHET not forwarded through pathogenModel");
hg = rx([rx.Name] == "Heterotroph growth");
assert(hg.TemperatureResponse == "cardinal", "HETCardinal not forwarded through pathogenModel");
% Defaults stay inert: no grazing row, no cardinal, unchanged thetas.
pmDefault = mpcssf.presets.pathogenModel();
rxDef = pmDefault.Reactions;
assert(~any([rxDef.Name] == "Heterotroph grazing"), "default pathogenModel gained a grazing row");
assert(rxDef([rxDef.Name] == "Heterotroph growth").TemperatureResponse == "exponential", ...
    "default pathogenModel heterotroph growth response changed");
fprintf("  PASS\n");

%% ------------------------------------------------------------------ (vi)
fprintf("=== (vi) smoke simulate with Grazing and with HETCardinal ===\n");
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
smokeArgsCommon = {"InflowConcentrations", infl, "SimulationTime", 0.01, ...
    "TimeStep", "adaptive", "AdaptiveInitialTimeStep", 1e-8, "AdaptiveMaxTimeStep", 3e-6, ...
    "ImplicitOsmosis", true, "CohesionBoundaryConditions", "neumann", "CohesionScheme", "shin", ...
    "FrameNumber", 3, "Quiet", true};

smokeLabels = ["grazing", "cardinal"];
smokeModels = {mpcssf.presets.pathogenModel(Grazing=0.8), mpcssf.presets.pathogenModel(HETCardinal=[-5 37 45])};
for si = 1:numel(smokeLabels)
    label = smokeLabels(si); m = smokeModels{si};
    f = mpcssf.SandFilter(Temperature=19, ...
        LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
    f = f.addGridPoints(50);
    r = simulate(mpcssf.State(f, m), smokeArgsCommon{:});
    assert(r.Flag == "OK", "%s smoke run did not finish OK: %s", label, r.Flag);
    C = r.Frames.Concentrations;
    for nm = [m.Particles.Name]
        for ph = ["Matrix", "Enclosed", "Flowing"]
            assert(all(isfinite(C{nm, ph}{1}(:))), "%s: non-finite %s/%s", label, nm, ph);
        end
    end
    for nm = [m.Liquids.Name]
        for ph = ["Enclosed", "Flowing"]
            assert(all(isfinite(C{nm, ph}{1}(:))), "%s: non-finite %s/%s", label, nm, ph);
        end
    end
    fprintf("  %s: PASS (Flag = %s)\n", label, r.Flag);
end

fprintf("MATLAB ASYMMETRY-OPTIONS TESTS PASS\n");
