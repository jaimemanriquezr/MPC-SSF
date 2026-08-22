function probePInventory(season)
% G-probe (plan: fixed-influent flip evaluation, step 1). Phosphorus inventory,
% boundary fluxes, recycle fluxes, and closure for summer/winter on the decided
% structure (field influent 5e-4, corrected kinetics, PHO-only endogenous
% respiration k_ra=0.276 theta 1.08, death 0.09/2.0 O2-neutral, theta_PHO 1.066).
%
% P contents per unit biomass (growth stoichiometry): P_PHO = 0.01, P_HET = 0.0141.
% Known non-closure candidates (quantified here, not hidden): HET death releases
% 0.0209 HPO4 per unit HET (> its 0.0141 content); POM carries no P state;
% PHO death releases 0.0037 (< content 0.01, remainder implicitly in POM).
% PHO endogenous respiration releases exactly its content (0.01) - closes.
W = "/Users/jaime/Research/SSF/code/MPC-SSF/.claude/worktrees/agent-a5f45f8cf63f54f77";
S = "/private/tmp/claude-501/-Users-jaime-Research-SSF/18cb4123-7947-42f0-8409-e994c38c4600/scratchpad";
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));
m = pathogenModel(NormalizedLight=true);
rx = m.Reactions;
names = [rx.Name];
iH = names == "Heterotroph growth";  iP = names == "Phototroph growth";
iY = names == "Hydrolysis";
iDH = names == "Heterotroph death";  iDP = names == "Phototroph death";
rx(iH).NominalRate = 0.042*24;
rx(iP).NominalRate = 0.125*24;
rx(iP).TemperatureCorrectionFactor = 1.066;
rx(iH).HalfSaturationConstants = dictionary( ...
    ["O2", "DOM", "NH4", "HPO4"], [2.0E-04, 4.0E-03, 1.0E-06, 2.0E-05]);
rx(iP).HalfSaturationConstants = dictionary( ...
    ["IC", "NH4", "HPO4"], [1.2E-03, 2.0E-05, 2.0E-05]);
rx(iY).HalfSaturationConstants = dictionary("POM/HET", 0.1);
rx(iDP).NominalRate = 0.09;
rx(iDH).StoichiometricCoefficients = dictionary( ...
    ["HET", "POM", "NH4", "HPO4"], [-1.0, 0.9123, 0.0653, 0.0209]);
rx(iDP).StoichiometricCoefficients = dictionary( ...
    ["PHO", "POM", "NH4", "HPO4"], [-1.0, 0.6316, 0.0221, 0.0037]);
phoEndog = Reaction(Name="Phototroph endogenous respiration", ...
    NominalRate=0.276, TemperatureCorrectionFactor=1.08, ...
    Order=dictionary("PHO", 1), ...
    HalfSaturationConstants=dictionary("O2", 2.0E-04), ...
    StoichiometricCoefficients=dictionary( ...
        ["PHO", "O2", "IC", "NH4", "HPO4"], [-1.0, -0.9301, 0.36, 0.06, 0.01]));
rx = [rx; phoEndog];
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
infl = [2.68e-3, 5.00e-4, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=30.0, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=240, ImplicitOsmosis=true, Quiet=true);
C = r.Frames.Concentrations;
ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
q = 7.2;  T = f.Temperature;
mass = @(nm) (sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}, 1)*dz)';
liqM = @(nm) (sum(C{nm,"Flowing"}{1} + C{nm,"Enclosed"}{1}, 1)*dz)';
outc = @(nm) C{nm,"Flowing"}{1}(end, :)';
phoM = mass("PHO"); hetM = mass("HET"); pomM = mass("POM"); hpo4M = liqM("HPO4");
P_PHO = 0.01; P_HET = 0.0141;
% temperature-corrected first-order rate constants
th = @(theta) theta^(T - 20);
kDP = 0.09*th(1.08); kDH = 2.0*th(1.08); kRA = 0.276*th(1.08);
% inventories
invP = @(i) hpo4M(i) + P_PHO*phoM(i) + P_HET*hetM(i);
% recycle fluxes (exact first-order integrals; Monod factors ~1, noted)
I = @(v, tEnd) trapz(ts(ts <= tEnd), v(ts <= tEnd));
res.season = season; res.flag = r.Flag;
for tEnd = [10, 30]
    tag = sprintf("d%d", tEnd);
    i = find(ts <= tEnd, 1, "last");
    B.release_deathPHO = 0.0037*kDP*I(phoM, tEnd);
    B.release_deathHET = 0.0209*kDH*I(hetM, tEnd);
    B.release_endog    = 0.01*kRA*I(phoM, tEnd);
    B.release = B.release_deathPHO + B.release_deathHET + B.release_endog;
    B.in_particulate = q*(P_PHO*infl(2) + P_HET*infl(1))*ts(i);
    B.out_dissolved  = q*I(outc("HPO4"), tEnd);
    B.out_particulate = q*(P_PHO*I(outc("PHO"), tEnd) + P_HET*I(outc("HET"), tEnd));
    B.out_POM_massflux = q*I(outc("POM"), tEnd);   % P-unquantified channel
    B.hpo4_inv = hpo4M(i); B.biomassP_inv = P_PHO*phoM(i) + P_HET*hetM(i);
    B.totalP_inv = invP(i); B.dInv = invP(i) - invP(1);
    B.closure_defect = B.dInv - (B.in_particulate - B.out_dissolved - B.out_particulate);
    B.tau_P = mean(hpo4M(ts <= tEnd))/(B.release/ts(i));   % days
    % gross PHO production (budget residual, as probeBudget)
    B.grossPHO = (phoM(i) - phoM(1)) - q*infl(2)*ts(i) + q*I(outc("PHO"), tEnd) ...
                 + kDP*I(phoM, tEnd) + kRA*I(phoM, tEnd);
    res.(tag) = B;
    fprintf("PINV %s %s: HPO4=%.3e biomP=%.4f totP=%.4f dInv=%.4f | in=%.4f outD=%.3e outP=%.3e | rel=%.4f (dPHO %.4f dHET %.4f end %.4f) | tau_P=%.2f d | grossPHO=%.4f | defect=%.4f\n", ...
        season, tag, B.hpo4_inv, B.biomassP_inv, B.totalP_inv, B.dInv, ...
        B.in_particulate, B.out_dissolved, B.out_particulate, ...
        B.release, B.release_deathPHO, B.release_deathHET, B.release_endog, ...
        B.tau_P, B.grossPHO, B.closure_defect);
end
save(fullfile(S, sprintf("pinv_%s.mat", season)), "res", "ts", "phoM", "hetM", "pomM", "hpo4M", "-v7");
fprintf("PINV PROBE %s: flag=%s PHO(30d)=%.4f\n", season, r.Flag, phoM(end));
end
