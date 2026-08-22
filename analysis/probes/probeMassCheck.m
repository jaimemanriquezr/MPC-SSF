function probeMassCheck(season, ~)
% ZERO-BIOLOGY mass-bookkeeping check: all reaction rates set to 0; transport,
% attachment, detachment, CH untouched. dM must equal supply - export exactly;
% any deviation measures the bias of the budget-residual method used across
% the probe chain (probeBudget/probeLight/probeTheta growth residuals).
% DIAGNOSTIC theta_growth,PHO sweep (uncited values, threshold-finding only —
% not a candidate fix; the citable implementation of steep cold-suppression is
% the cardinal form, probeCardinal). Measures the flip threshold theta* that an
% honest CTMI triplet must beat. Baseline = probeLight at lightScale 1.
lightScale = 1.0;
% Winter light-curve sweep (10 d) on the DECIDED structure: field influent
% 5e-4, corrected kinetics, PHO-only endogenous respiration (k_rb dropped,
% Jaime 2026-08-20), theta_growth,PHO = 1.066 (Campos theta_kga, audit fix).
% lightScale multiplies the seasonal light curve (winter peaks 0.6/0.3/0.15
% at scale 1/0.5/0.25). Budget instrumentation kept.
phoIn = 5.0e-4;
% = probeDeath (field influent, corrected kinetics, endogenous respiration,
% source-faithful death; always PHO+HET) with two additions:
%   1. parameterized influent PHO (winter seasonality sweep: 5e-4 / 1e-4 / 2e-5)
%   2. PHO budget instrumentation: death and endogenous losses integrate exactly
%      (first-order, constant T), supply and effluent export come from the frames,
%      and GROWTH is the closure residual of
%      dM/dt = q(c_in - c_out) + growth - death - endog.
withHet = false;
% = probeFieldInfl (field influent 5e-4, corrected kinetics, endogenous respiration)
% with the SOURCE-FAITHFUL death structure (audit addendum 2026-08-20):
%   d_PHO 0.4 -> 0.09 /d (Wolf b_ina,PH; our 0.4 was Wolf's HETEROTROPH b_ina,H)
%   d_HET stays 2.0 /d (matches Campos avg k_db = 2.06)
%   death stoichiometry made O2-NEUTRAL (both sources; Lund rows released O2)
% 10 days (day-10 W/S phiB ratio is 90-95% of the day-30 value in this family).
% = probeCombined (corrected kinetics + endogenous respiration) with the influent
% PHO reduced to the Campos2006b Fig.1(a) field value: 5e-4 kg/m3 (~3-7 ug/L Chla,
% Kempton Park summer Run 1; Table B.1's 1e-2 is ~15x field strength).
% HET influent kept at 2.68e-3 (unchanged - only the algal load is field-anchored).
% COMBINED probe: corrected kinetics + Campos-style endogenous respiration.
%   - corrected growth rates: mu_HET = 1.008/d, mu_PHO = 3.0/d (Campos2006)
%   - corrected half-sats (audit 2026-08-20, source-faithful):
%       HET: O2 2.0e-4, DOM 4.0e-3, NH4 1.0e-6 (protective), HPO4 2.0e-5
%       PHO: IC 1.2e-3, NH4 2.0e-5, HPO4 2.0e-5;  hydrolysis POM/HET 0.1
%   - endogenous respiration (X -1, O2 -0.9301, IC +0.36, NH4 +0.06, HPO4 +0.01):
%       PHO k_ra = 0.276/d theta 1.08 (always); HET k_rb = 1.72/d theta 1.047 (withHet)
%   - PG-excess OFF, normalized light, 30 d, 100 cells (matches probeSeason/probeEndog)
here = fileparts(mfilename("fullpath"));            % .../analysis/probes
W = fileparts(fileparts(here));                      % repo root
S = fullfile(here, "data");
if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(NormalizedLight=true);
rx = m.Reactions;
names = [rx.Name];
iH = names == "Heterotroph growth";
iP = names == "Phototroph growth";
iY = names == "Hydrolysis";
rx(iH).NominalRate = 0.042*24;
rx(iP).NominalRate = 0.125*24;
for k = 1:numel(rx), rx(k).NominalRate = 0.0; end
rx(iH).HalfSaturationConstants = dictionary( ...
    ["O2", "DOM", "NH4", "HPO4"], [2.0E-04, 4.0E-03, 1.0E-06, 2.0E-05]);
rx(iP).HalfSaturationConstants = dictionary( ...
    ["IC", "NH4", "HPO4"], [1.2E-03, 2.0E-05, 2.0E-05]);
rx(iY).HalfSaturationConstants = dictionary("POM/HET", 0.1);
iDH = names == "Heterotroph death";
iDP = names == "Phototroph death";
rx(iDP).NominalRate = 0.09;
rx(iDH).StoichiometricCoefficients = dictionary( ...
    ["HET", "POM", "NH4", "HPO4"], [-1.0, 0.9123, 0.0653, 0.0209]);
rx(iDP).StoichiometricCoefficients = dictionary( ...
    ["PHO", "POM", "NH4", "HPO4"], [-1.0, 0.6316, 0.0221, 0.0037]);
endogStoich = @(X) dictionary( ...
    [X, "O2", "IC", "NH4", "HPO4"], [-1.0, -0.9301, 0.36, 0.06, 0.01]);
phoEndog = Reaction(Name="Phototroph endogenous respiration", ...
    NominalRate=0.276, TemperatureCorrectionFactor=1.08, ...
    Order=dictionary("PHO", 1), ...
    HalfSaturationConstants=dictionary("O2", 2.0E-04), ...
    StoichiometricCoefficients=endogStoich("PHO"));
rx = [rx; phoEndog];
tag = sprintf("%s_masscheck", season);
if withHet
    hetEndog = Reaction(Name="Heterotroph endogenous respiration", ...
        NominalRate=1.72, TemperatureCorrectionFactor=1.047, ...
        Order=dictionary("HET", 1), ...
        HalfSaturationConstants=dictionary("O2", 2.0E-04), ...
        StoichiometricCoefficients=endogStoich("HET"));
    rx = [rx; hetEndog];
end
m = Model(m.Components, rx, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
if season == "summer"
    f = SandFilter(Temperature=19);
    light = @(t) lightScale*max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
else
    f = SandFilter(Temperature=3);
    light = @(t) lightScale*max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
end
f = f.addGridPoints(100);
f.LightIrradiation = light;
infl = [2.68e-3, phoIn, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=10.0, ...
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
phoM = sum(C{"PHO","Matrix"}{1} + C{"PHO","Enclosed"}{1} + C{"PHO","Flowing"}{1}, 1)'*dz;
phoOut = C{"PHO","Flowing"}{1}(end, :)';
q = 7.2;  T = f.Temperature;
kDeath = 0.09*1.08^(T - 20);
kEndog = 0.276*1.08^(T - 20);
B.supply  = q*phoIn*(ts(end) - ts(1));
B.export  = q*trapz(ts, phoOut);
B.death   = kDeath*trapz(ts, phoM);
B.endog   = kEndog*trapz(ts, phoM);
B.dM      = phoM(end) - phoM(1);
B.growth  = B.dM - B.supply + B.export + B.death + B.endog;
lastDay = ts >= max(ts) - 1;
res.tag = tag; res.flag = r.Flag; res.tFinal = r.TimeFinal;
res.o2_out_mean = mean(o2(end, lastDay));
res.o2_out_range = [min(o2(end, lastDay)), max(o2(end, lastDay))];
res.o2_bed_min = min(o2(:, lastDay), [], "all");
res.nh4_out_mean = mean(nh4(end, lastDay));
res.phib_max = max(phiB(:, end)); res.phib_sup = max(phiB(z < 0, end));
res.pho_mass_end = mass("PHO"); res.pho_mass_end = res.pho_mass_end(end);
res.het_mass_end = mass("HET"); res.het_mass_end = res.het_mass_end(end);
res.budget = B;
save(fullfile(S, sprintf("mck_%s.mat", tag)), "res", "z", "ts", "o2", "nh4", "phiB", "phoM", "phoOut", "-v7");
disp(res);
fprintf("MASSCHECK PROBE %s: flag=%s O2out=%.2f PHO=%.4f HET=%.4f\n", ...
    tag, r.Flag, res.o2_out_mean, res.pho_mass_end, res.het_mass_end);
fprintf("MASSCHECK-BUDGET %s [kg/m2 over %g d]: supply=%.4f export=%.4f death=%.4f endog=%.4f dM=%.4f growth(residual)=%.4f (%.0f%% of supply)\n", ...
    tag, ts(end)-ts(1), B.supply, B.export, B.death, B.endog, B.dM, B.growth, 100*B.growth/max(B.supply,1e-12));
end
