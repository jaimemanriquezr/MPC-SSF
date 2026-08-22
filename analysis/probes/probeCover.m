function res = probeCover(idx, opts)
% PROBECOVER  Run one cell of the covered/uncovered driver sweep.
%
%   res = probeCover(idx)                     % idx into coverCases("all")
%   res = probeCover(idx, Stage="A")          % idx into coverCases("A")
%   res = probeCover(idx, Stage="smoke", Force=true)
%
% Writes  analysis/probes/data/cover_<tag>.mat        (cfg, res, profiles, limitation)
%         analysis/results/cover_sweep/profile_<tag>.csv
% and prints a one-line summary. Skips work if the .mat already exists unless
% Force=true.
%
% WHY THIS EXISTS -----------------------------------------------------------
% Campos2002 measured, on two full-scale Thames beds ripened in parallel from
% bare sand (Bed 9 uncovered / Bed 10 under a light-excluding membrane, 104 d):
%   0-2 cm biomass   60 ugC/g (126.8 peak) vs 15   -> ratio 4-8x
%   0-10 cm biomass  22.5 vs 5.6                   -> ratio 4x
%   schmutzdecke     present vs ABSENT
%   TOC/DOC removal  25/23% vs 23/23%              -> IDENTICAL
% Our E2 gives a 4% biomass difference; 0.97 on the corrected stack. This probe
% sweeps the drivers that could each independently gate that response.
%
% The limitation diagnostic (.claude/decisions/2026-08-22-liebig-limitation-
% diagnostic.md) already established that the null has TWO regionally separate
% causes on the corrected stack:
%   supernatant  lit (Ihat 0.06-0.79) but nutrient-throttled (min-Monod ~0.005,
%                binding on IC in the biofilm and HPO4 in the flowing phase)
%   top 0-2 cm   nutrient-replete (min-Monod 0.703) but Ihat = 1.4e-27
% so covering removes light where growth was already nutrient-capped, and
% removes nothing where nutrients would allow growth.
%
% CONVENTIONS ---------------------------------------------------------------
% * Clean start, no maturation. Matches Campos (both beds resanded and ripened
%   in parallel) AND sidesteps the getMature cache-key trap entirely: that cache
%   keys on NCells + smoke flag only (manuscriptExperiments.m:190-201), so it
%   ignores kinetics, temperature, influent and etaSand.
% * Every output path derives from coverTag(cfg). No hand-typed filenames.
% * Campos normalises per gram of dry sand, so the commensurable depth average
%   is the plain ARITHMETIC mean of the bulk-volumetric attached concentration,
%   not the eps-weighted one -- eps is constant inside the bed and the quantity
%   is a concentration, not a conserved stock. Conserved-mass integrals ARE
%   eps-weighted (see massEps below); both are reported and labelled.
%
% See also COVERCASES, COVERTAG, COLLECTCOVER, PROBELIMITDIAG.

arguments
    idx (1,1) double {mustBeInteger, mustBePositive}
    opts.Stage (1,1) string = "all"
    opts.Force (1,1) logical = false
end

here = fileparts(mfilename("fullpath"));       % .../analysis/probes
W    = fileparts(fileparts(here));             % repo root
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);

dataDir = fullfile(here, "data");
csvDir  = fullfile(W, "analysis", "results", "cover_sweep");
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(csvDir),  mkdir(csvDir);  end

cases = coverCases(opts.Stage);
assert(idx <= numel(cases), "idx %d exceeds %d cases in stage %s", idx, numel(cases), opts.Stage);
cfg = cases(idx);

outMat = fullfile(dataDir, "cover_" + cfg.tag + ".mat");
if isfile(outMat) && ~opts.Force
    fprintf("COVER skip %s (exists; Force=false)\n", cfg.tag);
    res = load(outMat).res; return
end

%% ---- model -------------------------------------------------------------
m = buildCoverModel(cfg.kinetics);

%% ---- filter ------------------------------------------------------------
% 19 C and the summer light curve in every cell: this is a LIGHT experiment,
% so temperature is held fixed and only the cover multiplier moves.
f = SandFilter(Temperature=19);
f.LightAttenuationCoeffSand = cfg.etaSand;
f = f.addGridPoints(cfg.ncells);
Ls = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);   % lightSummer, manuscriptExperiments.m:84
s  = cfg.cover;
f.LightIrradiation = @(t) s*Ls(t);

%% ---- influent ----------------------------------------------------------
% Dictionary rather than a bare 9-vector: index-order-proof if a preset ever
% inserts a component (modelLund inserts PG at slot 5 when PGExcess=false).
infl = dictionary( ...
    ["HET","PHO","POM","PAT","O2","IC","NH4","HPO4","DOM"], ...
    [2.68e-3, cfg.phoIn, 0.0, 0.0, 9.10e-3, cfg.icIn, cfg.nh4In, cfg.hpo4In, 1.75e-4]);

fprintf("COVER start %-42s cover=%.4g pho=%.3g hpo4=%.3g ic=%.3g eta=%g n=%d t=%gd kin=%s\n", ...
    cfg.tag, cfg.cover, cfg.phoIn, cfg.hpo4In, cfg.icIn, cfg.etaSand, cfg.ncells, cfg.tsim, cfg.kinetics);
tRun = tic;

r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=cfg.tsim, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=cfg.maxdt, ...
    FrameNumber=max(round(24*cfg.tsim), 12), ImplicitOsmosis=true, Quiet=true, ...
    RecordLimitation=true);

wall = toc(tRun)/60;

%% ---- geometry and fields ----------------------------------------------
C  = r.Frames.Concentrations;
ts = r.Frames.Time(:);
z  = f.GridPoints.Centers(:);
dz = f.GridSize;
poros = computePorosity(f, z);
densityL = mean([m.Liquids.Density]);  densityP = mean([m.Particles.Density]);

sup  = z < 0;                       % supernatant
bed  = z >= 0;
top2 = z >= 0 & z <= 0.02;          % Campos's 0-2 cm sampling layer
top10= z >= 0 & z <= 0.10;          % Campos's 0-10 cm core

% Attached biomass per component: matrix + enclosed, bulk-volumetric [kg/m3].
att = @(nm) C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1};

% Biofilm volume fraction, the E2 / fig:seasons-results quantity.
phiB = C{"Water","Enclosed"}{1}/densityL;
for nm = [m.Particles.Name], phiB = phiB + att(nm)/densityP; end
for nm = [m.Liquids.Name],   phiB = phiB + C{nm,"Enclosed"}{1}/densityL; end

%% ---- Campos-commensurable biomass -------------------------------------
% kappa converts bulk-volumetric COD to ugC per g dry sand:
%   rho_sand = (1 - eps0)*rho_quartz = 0.6*2650 = 1590 kg/m3
%   COD -> C  for C5H7O2N            = 0.531/1.42 = 0.374 kgC/kgCOD
%   kappa = 0.374/1590 * 1e6         = 235.2 (ugC/g) per (kg COD/m3)
% Sanity: Campos's 60 ugC/g <-> 0.255 kg COD/m3, right order for a ripened top.
KAPPA = 0.374/((1 - f.SandPorosity)*2650) * 1e6;

attHET = att("HET");  attPHO = att("PHO");  attPOM = att("POM");
live  = attHET + attPHO;                    % Campos assayed total microbial biomass
withP = live + attPOM;                      % bracketing variant incl. detrital POM
% Plain mean over the layer at the FINAL frame: per-gram-sand normalisation, so
% the arithmetic (not eps-weighted) mean is the commensurable one.
layerMean = @(M, mask) mean(M(mask, end), 1);

res = struct();
res.tag = cfg.tag;  res.stage = cfg.stage;  res.cfg = cfg;
res.flag = string(r.Flag);  res.tFinal = r.TimeFinal;  res.wallMin = wall;
res.kappa = KAPPA;

res.B_0_2       = KAPPA*layerMean(live,   top2);
res.B_0_10      = KAPPA*layerMean(live,   top10);
res.B_0_2_pom   = KAPPA*layerMean(withP,  top2);
res.B_0_10_pom  = KAPPA*layerMean(withP,  top10);
res.Bpho_0_2    = KAPPA*layerMean(attPHO, top2);
res.Bpho_0_10   = KAPPA*layerMean(attPHO, top10);
res.Bhet_0_2    = KAPPA*layerMean(attHET, top2);

[bpk, ipk]      = max(KAPPA*live(bed, end));
zbed            = z(bed);
res.Bpeak       = bpk;
res.zPeak       = zbed(ipk);

%% ---- E2 continuity -----------------------------------------------------
res.phib_0_2  = mean(phiB(top2, end));
res.phib_0_10 = mean(phiB(top10, end));
res.phib_max  = max(phiB(:, end));
res.phib_sup  = max(phiB(sup, end));

%% ---- schmutzdecke presence/absence ------------------------------------
% "No schmutzdecke at all" on the covered bed is a PRESENCE measurement and
% needs its own metric -- a ratio of biomass cannot express it.
supMat = sum(live(sup, end))*dz;                       % kg/m2 above the sand
res.mat_mass  = supMat;
thick = phiB(sup, end) > 0.05;
res.mat_depth = sum(thick)*dz;                         % m of biofilm-bearing supernatant
res.mat_present = double(res.mat_depth > 0);

%% ---- Campos's second result: removal ----------------------------------
% Last-10%-of-run means, so a ripening transient does not contaminate them.
tail = ts >= ts(end) - max(0.1*cfg.tsim, 1.0);
domIn = 1.75e-4;
domOut = C{"DOM","Flowing"}{1}(end, :);
res.doc_removal = 1 - mean(domOut(tail))/domIn;
tocIn  = domIn + 2.68e-3 + cfg.phoIn;
tocOut = domOut + C{"HET","Flowing"}{1}(end,:) + C{"PHO","Flowing"}{1}(end,:) + C{"POM","Flowing"}{1}(end,:);
res.toc_removal = 1 - mean(tocOut(tail))/tocIn;
res.o2_out = mean(C{"O2","Flowing"}{1}(end, tail))*1000;      % mg/L
res.o2_bed_min = min(C{"O2","Flowing"}{1}(:, tail), [], "all")*1000;

%% ---- limitation diagnostic reductions ---------------------------------
res.lim = reduceLimitation(r.Frames.Limitation, z, sup, top2, tail);

%% ---- photosynthesis partition -----------------------------------------
% Proxy: lightFactor * monod * X_PHO, normalised to a partition, so the omitted
% constants (mu(T), phi) cancel. Answers "where in the column does the model
% actually photosynthesise" head-on.
L = r.Frames.Limitation;
jP = find(L.ReactionNames == "Phototroph growth", 1);
if ~isempty(jP)
    prod = double(L.LightFactor(:, :, jP)) .* double(L.MonodFlowing(:, :, jP)) .* C{"PHO","Flowing"}{1} ...
         + double(L.LightFactor(:, :, jP)) .* double(L.MonodBiofilm(:, :, jP)) .* attPHO;
    pv = sum(prod(:, tail), 2, "omitnan");
    tot = sum(pv);
    res.psFrac_sup   = sum(pv(sup))/max(tot, realmin);
    res.psFrac_0_2   = sum(pv(top2))/max(tot, realmin);
    res.psFrac_deep  = sum(pv(z > 0.02))/max(tot, realmin);
else
    [res.psFrac_sup, res.psFrac_0_2, res.psFrac_deep] = deal(NaN);
end

%% ---- mass audit --------------------------------------------------------
% Raw and eps-weighted stocks. NOTE: this residual includes reaction-network
% non-closure (the death rows are known not to close, ~45% -- see the kinetic
% audit), so it is NOT an absolute conservation check. It exists so the two arms
% of a pair can be compared: collectCover rejects a pair whose residuals differ
% by more than 5%, which is what licenses trusting the ratio.
allNames = [m.Particles.Name, m.Liquids.Name];
stockRaw = 0; stockEps = 0;
for nm = allNames
    tot = att(nm) + C{nm,"Flowing"}{1};
    stockRaw = stockRaw + sum(tot(:, end))*dz;
    stockEps = stockEps + sum(poros.*tot(:, end))*dz;
end
res.stock_raw = stockRaw;
res.stock_eps = stockEps;
res.supplied  = f.InflowVelocity * sum(lookup(infl, allNames, FallbackValue=0)) * cfg.tsim;
res.mass_ratio = stockEps / max(res.supplied, realmin);

% CLOSURE RESIDUAL -- the quantity the pair guard actually needs.
%
% stock/supplied was the first attempt and it is NOT a numerical diagnostic: two
% arms that grow differently hold different stock for entirely physical reasons,
% so it flagged every Stage A pair (0.181 uncovered vs 0.082 covered) while
% nothing was wrong. Replaced by in - out - dStock, normalised by what was
% supplied, which is zero for a conservative scheme regardless of how much
% biomass each arm happens to carry.
%
% Caveat kept explicit: the reaction network is NOT mass-closed (the death rows
% are known to violate closure ~45%, see the kinetic audit), so this residual is
% not expected to vanish. Both arms run the identical network, so the DIFFERENCE
% between arms is still the right common-mode test -- that is all the pair guard
% claims.
outflux = 0;
for nm = allNames
    cOut = C{nm,"Flowing"}{1}(end, :);
    outflux = outflux + trapz(ts, cOut) * f.InflowVelocity;
end
res.outflux = outflux;
res.closure_residual = (res.supplied - outflux - stockEps) / max(res.supplied, realmin);

%% ---- save --------------------------------------------------------------
profiles = struct("z", z, "poros", poros, "phiB", phiB(:,end), ...
    "HET", attHET(:,end), "PHO", attPHO(:,end), "POM", attPOM(:,end), ...
    "HPO4", C{"HPO4","Flowing"}{1}(:,end), "NH4", C{"NH4","Flowing"}{1}(:,end), ...
    "IC", C{"IC","Flowing"}{1}(:,end), "O2", C{"O2","Flowing"}{1}(:,end), ...
    "Ihat", double(L.LightAttenuated(:,end)));
save(outMat, "cfg", "res", "profiles", "ts", "-v7");

writematrix([profiles.z, profiles.poros, profiles.phiB, profiles.HET, profiles.PHO, ...
             profiles.POM, profiles.HPO4, profiles.NH4, profiles.IC, profiles.O2, profiles.Ihat], ...
            fullfile(csvDir, "profile_" + cfg.tag + ".csv"));

fprintf("COVER done  %-42s flag=%s t=%.1f wall=%.1fmin | B02=%.3g B010=%.3g ugC/g | " + ...
        "mat=%.4g kg/m2 (%.3g m) | DOCrem=%.1f%% | O2out=%.2f | psSup=%.3f\n", ...
    cfg.tag, res.flag, res.tFinal, wall, res.B_0_2, res.B_0_10, ...
    res.mat_mass, res.mat_depth, 100*res.doc_removal, res.o2_out, res.psFrac_sup);
end


% =========================================================================
function m = buildCoverModel(kinetics)
% Corrected stack per .claude/decisions/2026-08-20-kinetic-parameter-audit.md,
% mirroring probeSeason90.m so the two are comparable. "manuscript" is the
% published stack exactly as manuscriptExperiments runs it.
if kinetics == "manuscript"
    m = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
    rx = m.Reactions;
else
    m = pathogenModel(NormalizedLight=true);
    rx = m.Reactions;
    nm = [rx.Name];
    iH = nm == "Heterotroph growth";  iP = nm == "Phototroph growth";
    iY = nm == "Hydrolysis";
    iDH = nm == "Heterotroph death";  iDP = nm == "Phototroph death";
    rx(iH).NominalRate = 0.042*24;                 % 1.008 /d  (Campos2006 upper)
    rx(iP).NominalRate = 0.125*24;                 % 3.0   /d
    rx(iP).TemperatureCorrectionFactor = 1.066;    % Campos theta_kga
    % Retire the dark floor EXPLICITLY. modelLund only zeroes
    % MinimumLightFactor inside its own PhototrophRespiration > 0 branches
    % (modelLund.m:119,146). The corrected stack builds respiration by appending
    % phoEndog by hand, so those branches are NOT taken and the 0.01 floor
    % survives -- meaning a fully covered filter (cover = 0) would still
    % photosynthesise at 1% of optimum, capping the very contrast this sweep
    % measures. The manuscript stack, which does pass PhototrophRespiration,
    % has no floor. Backwards, and caught by the smoke test.
    % NOTE: probeSeason90.m builds its "corrected"/"field" variants the same way
    % and therefore carries the floor; its committed data is affected.
    rx(iP).MinimumLightFactor = 0.0;
    rx(iH).HalfSaturationConstants = dictionary( ...
        ["O2","DOM","NH4","HPO4"], [2.0e-4, 4.0e-3, 1.0e-6, 2.0e-5]);
    rx(iP).HalfSaturationConstants = dictionary( ...
        ["IC","NH4","HPO4"], [1.2e-3, 2.0e-5, 2.0e-5]);
    rx(iY).HalfSaturationConstants = dictionary("POM/HET", 0.1);
    rx(iDP).NominalRate = 0.09;                    % Wolf b_ina,PH
    % O2-neutral death (RWQM1 sets Y_ALG,death to avoid consuming O2/nutrients).
    rx(iDH).StoichiometricCoefficients = dictionary( ...
        ["HET","POM","NH4","HPO4"], [-1.0, 0.9123, 0.0653, 0.0209]);
    rx(iDP).StoichiometricCoefficients = dictionary( ...
        ["PHO","POM","NH4","HPO4"], [-1.0, 0.6316, 0.0221, 0.0037]);
    phoEndog = Reaction(Name="Phototroph endogenous respiration", ...
        NominalRate=0.276, TemperatureCorrectionFactor=1.08, ...
        Order=dictionary("PHO", 1), ...
        HalfSaturationConstants=dictionary("O2", 2.0e-4), ...
        StoichiometricCoefficients=dictionary( ...
            ["PHO","O2","IC","NH4","HPO4"], [-1.0, -0.9301, 0.36, 0.06, 0.01]));
    rx = [rx; phoEndog];
end
m = Model(m.Components, rx, Kappa=1e-7, Zeta0=1e2, Zeta1=1e-2, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
    BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
end


% =========================================================================
function out = reduceLimitation(L, z, sup, top2, tail)
% Collapse the per-frame limitation arrays into scalars: for each of the two
% growth reactions x {supernatant, top 0-2 cm} x {biofilm, flowing}, the modal
% binding substrate, its occupancy fraction, and the mean value of the min.
out = struct();
for reac = ["Heterotroph growth", "Phototroph growth"]
    j = find(L.ReactionNames == reac, 1);
    if isempty(j), continue, end
    key0 = matlab.lang.makeValidName(extractBefore(reac + " ", " "));
    for regn = ["sup", "top2"]
        mask = sup; if regn == "top2", mask = top2; end
        for phase = ["Biofilm", "Flowing"]
            v  = L.(phase)(mask, tail, j);   v = v(:);   v = v(v > 0);
            mo = double(L.("Monod"+phase)(mask, tail, j));
            key = key0 + "_" + regn + "_" + lower(phase);
            if isempty(v)
                out.(key + "_species") = "none";
                out.(key + "_frac") = NaN;
            else
                u = unique(v);  cnt = arrayfun(@(x) sum(v==x), u);
                [mx, o] = max(cnt);
                out.(key + "_species") = L.Names(u(o));
                out.(key + "_frac") = mx/numel(v);
            end
            out.(key + "_monod") = mean(mo(:), "omitnan");
        end
    end
end
la = double(L.LightAttenuated);
out.Ihat_surface = max(la(1, tail));
out.Ihat_z0      = max(la(find(z >= 0, 1), tail));
iBed = find(z > 0, 1) + 1;
if iBed <= size(la,1), out.Ihat_bed1 = max(la(iBed, tail)); else, out.Ihat_bed1 = NaN; end
end
