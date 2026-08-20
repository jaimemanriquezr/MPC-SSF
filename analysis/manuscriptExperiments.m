function manuscriptExperiments(names, options)
% MANUSCRIPTEXPERIMENTS  Redo the Manriquez2026 manuscript experiments on the
% corrected model (PG-excess respiration + published transfer rates).
%
%   manuscriptExperiments("all")
%   manuscriptExperiments(["E4", "X1"], NCells=100)
%
% Protocols follow the HEAD scripts of the legacy slow-sand-filtration repo
% (see .claude/plans/2026-08-19 manuscript-driver plan and PARAMETERS.md):
%   E1 seasons (90 d summer 19C / winter 3C, clean start)
%   E2 covered vs uncovered (30 d, covered = 0.01x summer light)
%   E3 long-term outputs (from the E1 summer run: snapshots, effluent series,
%      (z,t) liquid fields, total mass)
%   E4 scraping (from the 30 d summer state; depths 0/5/15/25/50 cm, refill)
%   E5 marker pulses (day 30-32; Pulse = 1x nominal, BigPulse = 100x)
%   E6 controls (Clean = no marker; ConstantFeed = 5.3616e-3 steady)
%   E7/E8 (z,t) marker fields from E5   E9 HET/PHO log removals from E5/E6
%   E10 scraping x constant marker feed (depths 0/15/25/50 cm, 10 d)
%   X1 O2 day/night depth profiles + effluent diel series from the summer run
%
% Dependencies are resolved through a state cache in results/manuscript/cache:
% the 90 d summer run produces mature30_summer, consumed by E4/E5/E6/E10.
% Smoke=true shrinks every duration to minutes for end-to-end checks.
arguments
    names string = "all"
    options.Respiration (1,1) {mustBeNumeric} = 0.55;
    options.PGExcess (1,1) logical = true;
    options.NCells (1,1) {mustBeNumeric} = 100;
    options.MaxDt (1,1) {mustBeNumeric} = 3e-6;
    % Article value (Table 2) — the presets carry 1e-6, which is CH-unstable at
    % fine grids (dt 3e-6 blew up at 100 cells, t=5.9 d); the manuscript's 1e-7
    % is what allowed the legacy article runs at N=500, dt 1e-5.
    options.Kappa (1,1) {mustBeNumeric} = 1e-7;
    options.Smoke (1,1) logical = false;
    options.OutRoot (1,1) string = "";
end

here = fileparts(mfilename('fullpath'));
if options.OutRoot == ""
    options.OutRoot = fullfile(here, "results", "manuscript");
end
if ~isfolder(fullfile(options.OutRoot, "cache")), mkdir(fullfile(options.OutRoot, "cache")); end

% Durations (days). Smoke mode compresses everything to minutes of wall time.
if options.Smoke
    D = struct("long", 0.2, "matureAt", 0.1, "covered", 0.1, "postScrape", 0.05, ...
               "patStart", 0.1, "pulseOn", 0.12, "pulseOff", 0.14, "patLen", 0.08);
else
    D = struct("long", 90, "matureAt", 30, "covered", 30, "postScrape", 30, ...
               "patStart", 30, "pulseOn", 30, "pulseOff", 32, "patLen", 10);
end

order = ["E3", "E1", "E2", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "X1"];
if isscalar(names) && names == "all"
    names = order;
else
    [~, idx] = ismember(names, order);
    [~, srt] = sort(idx);
    names = names(srt);
end

for nm = names
    fprintf("\n----- %s (%s) -----\n", nm, string(datetime));
    tStage = tic;
    switch nm
        case "E1",  runSummer(options, D);  e1_winter(options, D);
        case "E2",  e2_covered(options, D);
        case "E3",  runSummer(options, D);
        case "E4",  e4_scraping(options, D);
        case "E5",  e5_pulses(options, D);
        case "E6",  e6_controls(options, D);
        case {"E7", "E8"}, e5_pulses(options, D);   % fields are written by E5
        case "E9",  e5_pulses(options, D); e6_controls(options, D);
        case "E10", e10_scrape_pat(options, D);
        case "X1",  runSummer(options, D);          % X1 files written with the summer run
        otherwise, error("unknown experiment id: %s", nm);
    end
    fprintf("----- %s done in %.1f min -----\n", nm, toc(tStage)/60);
end
end

% ============================ shared machinery ================================
function L = lightSummer()
L = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
end

function L = lightWinter()
L = @(t) max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
end

function infl = baseInfluent(pat)
% Table B.1 (9 components); PAT slot per experiment (biofilm runs use nominal 0).
infl = [2.68e-3, 1.00e-2, 0.0, pat, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
end

function m = theModel(opts, family)
% family = "biofilm" or "pathogen". DEVIATION from legacy (2026-08-20): the
% legacy repo ran pathogen experiments with detachment 1.4e-5*sqrt(v/18)
% (srun_pathogen), 1e-4 of the biofilm value -- under the corrected model the
% unshed biofilm grows from the mature phib ~0.5 to the 0.99 clog inside the
% 10-day pathogen window. Both families therefore use the biofilm detachment
% 0.14*sqrt(v/18); flagged for co-author review (changes PAT matrix-release
% tails vs the published figures).
m = pathogenModel(PhototrophRespiration=opts.Respiration, PGExcess=opts.PGExcess, NormalizedLight=true);
if family == "biofilm" || family == "pathogen"
    m = Model(m.Components, m.Reactions, Kappa=opts.Kappa, ...
        Zeta0=m.CohesionSubModel.Zeta0, Zeta1=m.CohesionSubModel.Zeta1, ...
        DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=m.WaterDensity, ...
        BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
    return
end
if opts.Kappa ~= m.CohesionSubModel.Kappa
    m = Model(m.Components, m.Reactions, Kappa=opts.Kappa, ...
        Zeta0=m.CohesionSubModel.Zeta0, Zeta1=m.CohesionSubModel.Zeta1, ...
        DetachmentFunction=m.DetachmentFunction, WaterDensity=m.WaterDensity, ...
        BiofilmPorosity=m.BiofilmPorosity, OsmosisRate=m.OsmosisRate);
end
end

function f = theFilter(opts, tempC, light)
f = SandFilter(Temperature=tempC);
f = f.addGridPoints(opts.NCells);
f.LightIrradiation = light;
end

function r = runSim(s, tsim, infl, nframes, opts)
r = simulate(s, InflowConcentrations=infl, SimulationTime=tsim, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=nframes, ImplicitOsmosis=true, Quiet=true);
assert(r.Flag == "OK", "run failed: flag=" + r.Flag);
end

function st = frameState(r, k)
% State snapshot of frame k as a re-homeable struct (finalState generalized).
C = r.Frames.Concentrations;
pNames = [r.Model.Particles.Name];  lNames = [r.Model.Liquids.Name];
st.pNames = pNames;  st.lNames = lNames;
st.Matrix = grab(C, pNames, "Matrix", k);
st.EnclosedParticles = grab(C, pNames, "Enclosed", k);
st.FlowingParticles = grab(C, pNames, "Flowing", k);
st.EnclosedLiquids = grab(C, lNames, "Enclosed", k);
st.FlowingLiquids = grab(C, lNames, "Flowing", k);
densityL = mean([r.Model.Liquids.Density]);
st.EnclosedWaterVolume = C{"Water", "Enclosed"}{1}(:, k)/densityL;
st.VelocityBiofilm = r.Frames.Velocity.Biofilm(:, k);
st.Time = r.Frames.Time(k);
end

function M = grab(C, names, vol, k)
M = [];
for nm = names
    v = C{nm, vol}{1};
    if isscalar(v), v = zeros(size(C{names(1), "Enclosed"}{1}, 1), size(C{names(1), "Enclosed"}{1}, 2)); end
    M = [M, v(:, k)]; %#ok<AGROW>
end
end

function s = rehome(f, m, st, t0)
s = State(f, m);
s.GlobalConcentration.Matrix = st.Matrix;
s.GlobalConcentration.EnclosedParticles = st.EnclosedParticles;
s.GlobalConcentration.FlowingParticles = st.FlowingParticles;
s.GlobalConcentration.EnclosedLiquids = st.EnclosedLiquids;
s.GlobalConcentration.FlowingLiquids = st.FlowingLiquids;
s.EnclosedWaterVolume = st.EnclosedWaterVolume;
s.Velocity.Biofilm = st.VelocityBiofilm;
s.Time = t0;
end

function st = scrapeState(st, zs, dz, centers)
% Port of legacy @SDresults/scrape.m (refill-sand family) / Julia scrape!:
% cells 0 <= z <= ceil(zs/dz + 1/2)*dz lose biofilm and enclosed water; the
% flowing suspension is emptied everywhere; velocities dropped; clock reset.
scraped = ceil(zs/dz + 1/2)*sign(zs)*dz;
cells = (centers >= 0) & (centers <= scraped) & (zs > 0);
if zs == 0, cells = false(size(centers)); end   % "0 cm" removes only supernatant biofilm
sup = centers < 0;
st.Matrix(cells | sup, :) = 0;
st.EnclosedParticles(cells | sup, :) = 0;
st.EnclosedLiquids(cells | sup, :) = 0;
st.EnclosedWaterVolume(cells | sup) = 0;
st.FlowingParticles(:) = 0;
st.FlowingLiquids(:) = 0;
st.VelocityBiofilm(:) = 0;
end

function [st, cachePath] = getMature(opts, D)
% mature30_summer: the E1/E3 summer run's state at t = matureAt, cached.
cachePath = fullfile(opts.OutRoot, "cache", sprintf("mature30_summer_N%d%s.mat", ...
    opts.NCells, ternary(opts.Smoke, "_smoke", "")));
if isfile(cachePath)
    st = load(cachePath).st;
    return
end
runSummer(opts, D);
st = load(cachePath).st;
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function writeProfileCsv(outdir, name, centers, M)
writematrix([centers(:), M], fullfile(outdir, name));
end

% ============================ E1/E3/X1: summer + winter ======================
function runSummer(opts, D)
% One pass produces: E1 summer phi_b outputs, all E3 outputs, X1 O2 diel
% outputs, and the mature30_summer cache. Skipped when the sentinel exists.
outdir = fullfile(opts.OutRoot, "E3_longterm");
sentinel = fullfile(outdir, "done.txt");
if isfile(sentinel), return, end
if ~isfolder(outdir), mkdir(outdir); end

nframes = max(round(24*D.long), 12);            % hourly frames
f = theFilter(opts, 19, lightSummer());
m = theModel(opts, "biofilm");
r = runSim(State(f, m), D.long, baseInfluent(0.0), nframes, opts);
centers = f.GridPoints.Centers(:);
dz = f.GridSize;
ts = r.Frames.Time(:);
C = r.Frames.Concentrations;

% --- E3: effluent series for the five liquids (hourly) -----------------------
liq = ["O2", "IC", "NH4", "HPO4", "DOM"];
eff = ts;
for nm = liq, eff = [eff, C{nm, "Flowing"}{1}(end, :)']; end %#ok<AGROW>
writematrix(eff, fullfile(outdir, "effluent_liquids.csv"));

% --- E3: (z,t) fields, daily resolution --------------------------------------
keep = 1:max(1, round(24)):length(ts);          % ~daily columns
for nm = liq
    writematrix([[0; centers], [ts(keep)'; C{nm, "Flowing"}{1}(:, keep)]], ...
        fullfile(outdir, "field_" + nm + ".csv"));
end

% --- E1/E3: phi_b snapshots + total biomass ----------------------------------
phiB = biofilmFraction(r);
snaps = unique(min([1 2 5 10 20, D.matureAt + [0 1 2 5 10 20]], D.long));
ks = arrayfun(@(d) nearestIdx(ts, d), snaps);
writeProfileCsv(outdir, "phib_snapshots_summer.csv", centers, phiB(:, ks));
writematrix(snaps(:), fullfile(outdir, "phib_snapshot_days.csv"));
mass = totalBiomass(r, dz);
writematrix([ts, mass(:)], fullfile(outdir, "total_biomass_summer.csv"));

% --- X1: O2 day/night profiles + diel effluent -------------------------------
xdir = fullfile(opts.OutRoot, "X1_o2_diel");
if ~isfolder(xdir), mkdir(xdir); end
lastDay = ts >= (D.long - 1);
o2 = C{"O2", "Flowing"}{1};
L = lightSummer();
lightVals = arrayfun(L, ts);
tsl = ts(lastDay); o2l = o2(:, lastDay); lv = lightVals(lastDay);
[~, kNoon] = max(lv); [~, kMid] = min(lv);
writeProfileCsv(xdir, "o2_profiles_noon_midnight.csv", centers, [o2l(:, kNoon), o2l(:, kMid)]);
writematrix([ts, o2(end, :)', lightVals], fullfile(xdir, "effluent_o2_diel.csv"));

% --- cache mature30_summer ----------------------------------------------------
st = frameState(r, nearestIdx(ts, D.matureAt)); %#ok<NASGU>
save(fullfile(opts.OutRoot, "cache", sprintf("mature30_summer_N%d%s.mat", ...
    opts.NCells, ternary(opts.Smoke, "_smoke", ""))), "st");
fid = fopen(sentinel, "w"); fprintf(fid, "%s\n", string(datetime)); fclose(fid);
end

function e1_winter(opts, D)
outdir = fullfile(opts.OutRoot, "E1_seasons");
if isfile(fullfile(outdir, "phib_final_winter.csv")), return, end
if ~isfolder(outdir), mkdir(outdir); end
f = theFilter(opts, 3, lightWinter());
r = runSim(State(f, theModel(opts, "biofilm")), D.long, baseInfluent(0.0), max(round(2*D.long), 6), opts);
centers = f.GridPoints.Centers(:);
phiB = biofilmFraction(r);
writeProfileCsv(outdir, "phib_final_winter.csv", centers, phiB(:, end));
writematrix([r.Frames.Time(:), totalBiomass(r, f.GridSize)'], ...
    fullfile(outdir, "total_biomass_winter.csv"));
% Summer final profile for the side-by-side figure (from the E3 run's data).
e3 = fullfile(opts.OutRoot, "E3_longterm", "phib_snapshots_summer.csv");
copyfile(e3, fullfile(outdir, "phib_snapshots_summer.csv"));
end

% ============================ E2: covered vs uncovered ========================
function e2_covered(opts, D)
outdir = fullfile(opts.OutRoot, "E2_covered");
if isfile(fullfile(outdir, "phib_final_covered.csv")), return, end
if ~isfolder(outdir), mkdir(outdir); end
Ls = lightSummer();
for v = ["uncovered", "covered"]
    scale = ternary(v == "covered", 0.01, 1.0);
    f = theFilter(opts, 19, @(t) scale*Ls(t));
    r = runSim(State(f, theModel(opts, "biofilm")), D.covered, baseInfluent(0.0), max(round(2*D.covered), 6), opts);
    phiB = biofilmFraction(r);
    writeProfileCsv(outdir, "phib_final_" + v + ".csv", f.GridPoints.Centers(:), phiB(:, end));
end
end

% ============================ E4: scraping (biofilm) ==========================
function e4_scraping(opts, D)
outdir = fullfile(opts.OutRoot, "E4_scraping");
if isfile(fullfile(outdir, "total_biomass.csv")), return, end
if ~isfolder(outdir), mkdir(outdir); end
st0 = getMature(opts, D);
f = theFilter(opts, 19, lightSummer());
centers = f.GridPoints.Centers(:);
dz = f.GridSize;
depths = [0 5 15 25 50]/100;                    % m (HEAD batch_scraping.m:20)
if opts.Smoke, depths = depths(1:2); end
massT = [];
for zs = depths
    st = scrapeState(st0, zs, dz, centers);
    s = rehome(f, theModel(opts, "biofilm"), st, 0.0);
    r = runSim(s, D.postScrape, baseInfluent(0.0), max(round(2*D.postScrape), 6), opts);
    tag = sprintf("GP%d", round(zs*100));
    phiB = biofilmFraction(r);
    ts = r.Frames.Time(:);
    ks = arrayfun(@(d) nearestIdx(ts, d), unique(min([0 1 2 5 10 20], D.postScrape)));
    writeProfileCsv(outdir, "phib_snapshots_" + tag + ".csv", centers, phiB(:, ks));
    massT = pad2(massT, [ts, totalBiomass(r, dz)']);
end
writematrix(massT, fullfile(outdir, "total_biomass.csv"));
end

% ============================ E5/E6: marker feeds =============================
function e5_pulses(opts, D)
runMarker(opts, D, "E5_pulses", ["Pulse", "BigPulse"]);
end

function e6_controls(opts, D)
runMarker(opts, D, "E6_controls", ["Clean", "ConstantFeed"]);
end

function runMarker(opts, D, dirname, variants)
outdir = fullfile(opts.OutRoot, dirname);
if isfile(fullfile(outdir, "done.txt")), return, end
if ~isfolder(outdir), mkdir(outdir); end
st0 = getMature(opts, D);
f = theFilter(opts, 19, lightSummer());
patNom = 5.3616e-3;
for v = variants
    base = baseInfluent(0.0);
    switch v
        case "Pulse",        infl = pulseHandle(base, patNom, D.pulseOn, D.pulseOff);  cref = patNom;
        case "BigPulse",     infl = pulseHandle(base, 100*patNom, D.pulseOn, D.pulseOff); cref = 100*patNom;
        case "Clean",        infl = base;                                              cref = NaN;
        case "ConstantFeed", infl = baseInfluent(patNom);                              cref = patNom;
    end
    s = rehome(f, theModel(opts, "pathogen"), st0, D.patStart);
    r = runSim(s, D.patLen, infl, max(round(144*D.patLen), 12), opts);
    C = r.Frames.Concentrations;
    ts = r.Frames.Time(:);
    patOut = C{"PAT", "Flowing"}{1}(end, :)';
    hetF = log10(2.68e-3./max(C{"HET", "Flowing"}{1}(end, :)', 1e-30));
    phoF = log10(1.00e-2./max(C{"PHO", "Flowing"}{1}(end, :)', 1e-30));
    patF = log10(max(cref, eps)./max(patOut, 1e-30));
    writematrix([ts, patOut, patF, hetF, phoF], fullfile(outdir, "effluent_" + v + ".csv"));
    % (z,t) marker fields (E7/E8): daily-or-denser columns
    keep = 1:max(1, round(length(ts)/48)):length(ts);
    centers = f.GridPoints.Centers(:);
    for vol = ["Flowing", "Matrix", "Enclosed"]
        writematrix([[0; centers], [ts(keep)'; C{"PAT", vol}{1}(:, keep)]], ...
            fullfile(outdir, "patfield_" + v + "_" + vol + ".csv"));
    end
end
fid = fopen(fullfile(outdir, "done.txt"), "w"); fprintf(fid, "%s\n", string(datetime)); fclose(fid);
end

function h = pulseHandle(base, amp, t0, t1)
pulsed = base; pulsed(4) = amp;
h = @(t) selectRow(t >= t0 && t < t1, pulsed, base);
end

function out = selectRow(c, a, b)
if c, out = a; else, out = b; end
end

% ============================ E10: scraping x marker ==========================
function e10_scrape_pat(opts, D)
outdir = fullfile(opts.OutRoot, "E10_scrape_pat");
if isfile(fullfile(outdir, "done.txt")), return, end
if ~isfolder(outdir), mkdir(outdir); end
st0 = getMature(opts, D);
f = theFilter(opts, 19, lightSummer());
centers = f.GridPoints.Centers(:);
dz = f.GridSize;
patNom = 5.3616e-3;
depths = [0 15 25 50]/100;                      % fig_pathogens.m:86-92
if opts.Smoke, depths = depths(1:2); end
for zs = depths
    st = scrapeState(st0, zs, dz, centers);
    s = rehome(f, theModel(opts, "pathogen"), st, 0.0);
    r = runSim(s, D.patLen, baseInfluent(patNom), max(round(24*D.patLen), 12), opts);
    ts = r.Frames.Time(:);
    patOut = r.Frames.Concentrations{"PAT", "Flowing"}{1}(end, :)';
    writematrix([ts, patOut, log10(max(patNom, eps)./max(patOut, 1e-30))], ...
        fullfile(outdir, sprintf("logremoval_GP%d.csv", round(zs*100))));
end
fid = fopen(fullfile(outdir, "done.txt"), "w"); fprintf(fid, "%s\n", string(datetime)); fclose(fid);
end

% ============================ small numerics helpers ==========================
function phiB = biofilmFraction(r)
C = r.Frames.Concentrations;
densityL = mean([r.Model.Liquids.Density]);
densityP = mean([r.Model.Particles.Density]);
phiB = C{"Water", "Enclosed"}{1}/densityL;
for nm = [r.Model.Particles.Name]
    phiB = phiB + (C{nm, "Matrix"}{1} + C{nm, "Enclosed"}{1})/densityP;
end
for nm = [r.Model.Liquids.Name]
    phiB = phiB + C{nm, "Enclosed"}{1}/densityL;
end
end

function mass = totalBiomass(r, dz)
C = r.Frames.Concentrations;
mass = 0;
for nm = [r.Model.Particles.Name]
    mass = mass + sum(C{nm, "Matrix"}{1} + C{nm, "Enclosed"}{1}, 1)*dz;
end
for nm = [r.Model.Liquids.Name]
    mass = mass + sum(C{nm, "Enclosed"}{1}, 1)*dz;
end
end

function k = nearestIdx(ts, d)
[~, k] = min(abs(ts - d));
end

function M = pad2(M, N)
% Column-concatenate matrices with possibly different row counts (pad NaN).
if isempty(M), M = N; return, end
h = max(size(M, 1), size(N, 1));
M(end+1:h, :) = NaN;  N(end+1:h, :) = NaN;
M = [M, N];
end
