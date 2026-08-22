function probeSeason90(season, variant)
% 90-DAY seasonal pair at the manuscript's own resolution and metric, to test
% whether the winter/summer inversion survives in the quantity results.tex:165
% actually claims ("growth in summer being substantially greater",
% fig:seasons-results = biofilm VOLUME FRACTION after 90 d).
%
% The existing probe chain measured PHO MASS over 10-30 d at field influent and
% found W/S = 1.9-3.4; recomputing those same runs in phi_b gives only 1.05-1.18
% at the peak and ~1.0 at the sand surface. Three confounds separate the probe
% chain from the manuscript figure: duration (10-30 d vs 90 d), influent
% (5e-4 vs Table B.1's 1e-2), and kinetics (corrected vs published). This probe
% varies the last two at the manuscript's 90 d and reports BOTH metrics, so each
% confound is isolated:
%
%   variant "manuscript" : published kinetics  + Table B.1 influent (PHO 1e-2)
%                          -> must reproduce fig:seasons-results (summer > winter)
%   variant "corrected"  : corrected kinetics  + Table B.1 influent (PHO 1e-2)
%                          -> isolates the effect of the kinetics audit
%   variant "field"      : corrected kinetics  + field influent   (PHO 5e-4)
%                          -> isolates the effect of influent strength
%
% Fixed influent in ALL arms (no seasonal boundary): this is the E1/E2 protocol,
% and it is also Bellamy1985's protocol -- their Filter 6 was chilled to 5/2 C
% while Filter 1 ran at 17 C on the SAME Horsetooth raw water (Table 1), so a
% fixed-influent contrast is what the literature temperature data is comparable
% to. Bellamy Table 4 (coliform removal 97 -> 87%) and Bae2023 Fig 4 (DOC
% removal 43 -> 28%) constrain seasonal ACTIVITY, not standing stock; Bae2023
% Table 2 is the only standing-stock measurement and reads winter 1.44 vs summer
% 1.25 x1e9 cells/g (confounded by run time: 200 d vs 330 d since scraping).
%
% Manuscript configuration mirrored from analysis/manuscriptExperiments.m:
% E1/E3 = 90 d, summer 19 C / winter 3 C, clean start, NCells 100, MaxDt 3e-6,
% Kappa 1e-7, detachment 0.14*sqrt(|v|/18), NormalizedLight, Table B.1 influent
% (baseInfluent, PHO slot 1.00e-2).
%
% Every depth integral is reported RAW and POROSITY-WEIGHTED. The solver
% conserves sum(eps_i * c_i * dz) -- fluxes are porosity-weighted and divided by
% porosityCenters (simulate.m:498-513), so the conserved stock carries eps(z),
% which is 1 in the supernatant and 0.4 in the bed. The probe family's unweighted
% sum(c)*dz is what made the zero-biology check (probeMassCheck) read as
% transport creating mass; the eps-weighted integrals here are the correct ones.
arguments
    season (1,1) string {mustBeMember(season, ["summer", "winter"])}
    variant (1,1) string {mustBeMember(variant, ["manuscript", "corrected", "field"])}
end
here = fileparts(mfilename("fullpath"));            % .../analysis/probes
W = fileparts(fileparts(here));                      % repo root
S = fullfile(here, "data");
if ~isfolder(S), mkdir(S); end
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis"));

if variant == "field"
    phoIn = 5.0e-4;     % Campos2006b Fig.1(a) field value
else
    phoIn = 1.0e-2;     % Table B.1 (manuscriptExperiments/baseInfluent)
end

if variant == "manuscript"
    % Published stack exactly as E1/E3 runs it (manuscriptExperiments defaults).
    m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
    rx = m.Reactions;
else
    % Corrected stack (probeTheta at thetaG = 1.066, the decided structure):
    % corrected growth rates and half-sats, source-faithful O2-neutral death,
    % PHO-only endogenous respiration, PG-excess off.
    m = pathogenModel(NormalizedLight=true);
    rx = m.Reactions;
    names = [rx.Name];
    iH = names == "Heterotroph growth";
    iP = names == "Phototroph growth";
    iY = names == "Hydrolysis";
    rx(iH).NominalRate = 0.042*24;
    rx(iP).NominalRate = 0.125*24;
    rx(iP).TemperatureCorrectionFactor = 1.066;   % Campos theta_kga (audit fix)
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
    phoEndog = Reaction(Name="Phototroph endogenous respiration", ...
        NominalRate=0.276, TemperatureCorrectionFactor=1.08, ...
        Order=dictionary("PHO", 1), ...
        HalfSaturationConstants=dictionary("O2", 2.0E-04), ...
        StoichiometricCoefficients=dictionary( ...
            ["PHO", "O2", "IC", "NH4", "HPO4"], [-1.0, -0.9301, 0.36, 0.06, 0.01]));
    rx = [rx; phoEndog];
end

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

infl = [2.68e-3, phoIn, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
tag = sprintf("%s_%s", season, variant);
fprintf("SEASON90 start %s: phoIn=%g T=%g\n", tag, phoIn, f.Temperature);

r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=90.0, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=3e-6, ...
    FrameNumber=720, ImplicitOsmosis=true, Quiet=true);

C = r.Frames.Concentrations;
ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
poros = computePorosity(f, z);        % 1 in the supernatant, SandPorosity in the bed

o2 = C{"O2","Flowing"}{1}*1000;
nh4 = C{"NH4","Flowing"}{1}*1000;
densityL = mean([m.Liquids.Density]); densityP = mean([m.Particles.Density]);

% Biofilm volume fraction: the manuscript's plotted quantity (fig:seasons-results).
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name]
    phiB = phiB + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/densityP;
end
for nm = [m.Liquids.Name], phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end

% Depth integrals, raw and porosity-weighted (see header).
massRaw = @(nm) (sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}, 1)*dz)';
massEps = @(nm) (sum(poros.*(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}), 1)*dz)';

% Manuscript totalBiomass convention (manuscriptExperiments.m:414): attached
% phases only -- particles Matrix+Enclosed, liquids Enclosed.
bioRaw = zeros(size(ts)); bioEps = zeros(size(ts));
for nm = [m.Particles.Name]
    X = C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1};
    bioRaw = bioRaw + (sum(X, 1)*dz)';  bioEps = bioEps + (sum(poros.*X, 1)*dz)';
end
for nm = [m.Liquids.Name]
    X = C{nm,"Enclosed"}{1};
    bioRaw = bioRaw + (sum(X, 1)*dz)';  bioEps = bioEps + (sum(poros.*X, 1)*dz)';
end

phoM = massRaw("PHO");  phoMe = massEps("PHO");
hetM = massRaw("HET");  hetMe = massEps("HET");
phiInt = (sum(phiB, 1)*dz)';  phiIntE = (sum(poros.*phiB, 1)*dz)';

lastDay = ts >= max(ts) - 1;
res.tag = tag; res.season = season; res.variant = variant;
res.flag = r.Flag; res.tFinal = r.TimeFinal; res.phoIn = phoIn;
% Manuscript metric
res.phib_max = max(phiB(:, end));
res.phib_sup = max(phiB(z < 0, end));
res.phib_int_raw = phiInt(end);   res.phib_int_eps = phiIntE(end);
% Mass metrics
res.pho_mass_raw = phoM(end);     res.pho_mass_eps = phoMe(end);
res.het_mass_raw = hetM(end);     res.het_mass_eps = hetMe(end);
res.biomass_raw  = bioRaw(end);   res.biomass_eps  = bioEps(end);
% Chemistry (comparable to the O2 outputs that already stand)
res.o2_out_mean = mean(o2(end, lastDay));
res.o2_out_range = [min(o2(end, lastDay)), max(o2(end, lastDay))];
res.o2_bed_min = min(o2(:, lastDay), [], "all");
res.nh4_out_mean = mean(nh4(end, lastDay));

save(fullfile(S, sprintf("s90_%s.mat", tag)), "res", "z", "ts", "poros", "o2", "nh4", ...
    "phiB", "phoM", "phoMe", "hetM", "hetMe", "bioRaw", "bioEps", "phiInt", "phiIntE", "-v7");
disp(res);
% NB: "a" + "b", not ["a" "b"] -- the bracket form builds a STRING ARRAY and
% fprintf rejects it as a format ("Invalid file identifier"). That is what made
% every task of cosmos array 3530219 exit 1 AFTER completing its 90-day run and
% saving valid results.
fprintf("SEASON90 %s: flag=%s phib_max=%.4f phib_sup=%.4f phib_int(raw/eps)=%.4f/%.4f " + ...
        "PHO(raw/eps)=%.4f/%.4f HET(raw/eps)=%.4f/%.4f biomass(raw/eps)=%.4f/%.4f O2out=%.2f\n", ...
    tag, r.Flag, res.phib_max, res.phib_sup, res.phib_int_raw, res.phib_int_eps, ...
    res.pho_mass_raw, res.pho_mass_eps, res.het_mass_raw, res.het_mass_eps, ...
    res.biomass_raw, res.biomass_eps, res.o2_out_mean);
end
