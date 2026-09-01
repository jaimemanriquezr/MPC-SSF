function reproduceManriquez2026(what, opts)
% REPRODUCEMANRIQUEZ2026  Regenerate the Manriquez2026 results figures from the
% AUDITED WORKING SET (EXPERIMENTS.md E4/E7), not from the published parameters.
%
%   reproduceManriquez2026()                 % all figures from data on disk
%   reproduceManriquez2026("figures")
%   reproduceManriquez2026("scrape", Depth=0.04, Days=10)   % one scraping run
%
% The published figures came from the pre-audit model (pathogenModel presets,
% zeta0 = 1e2, sqrt detachment, Table B.1 influent HET 2.68e-3 / PHO 1.0e-2,
% 90 d summer / winter, N = 100). This script draws the SAME panels from the
% working set instead: 2x field influent (HET 3.0e-4 / PHO 1.0e-3), zeta0 = 1,
% zeta1 = 0.27, linear detachment, transfer x10, K_DOM 3e-4, K_HPO4 1e-6,
% N = 500, 19 C, and the E7 chain that now reaches t = 104 d. Every deviation
% from the published figure is listed in
% analysis/results/figures/repro/README.md.
%
% Data sources (all already on disk; nothing here starts a long run):
%   analysis/probes/data/chain/chain_fld2x_lit_leg{1..11}.mat   uncovered arm
%   analysis/probes/data/chain/chain_fld2x_dark_leg{1..11}.mat  covered arm
%   analysis/results/figures/repro/data/scrape_GP*.mat          scraping runs
%                                                (written by the "scrape" mode)
% The winter arm has no counterpart in the working set and is NOT reproduced.
arguments
    what (1,1) string = "all"
    opts.LitTag  (1,1) string = "fld2x_lit"
    opts.DarkTag (1,1) string = "fld2x_dark"
    % scrape mode
    opts.Depth (1,1) double = 0.0      % m, scraping depth
    opts.Days  (1,1) double = 10       % post-scrape regrowth
    opts.FromLeg (1,1) double = 3      % leg whose end state is scraped (leg 3 = 30 d)
    opts.NCells (1,1) double = 500
end
here = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(fileparts(here), "src")));
addpath(here); addpath(fullfile(here, "probes"));
outDir = fullfile(here, "results", "figures", "repro");
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(fullfile(outDir, "data")), mkdir(fullfile(outDir, "data")); end

switch what
    case "scrape", runScrape(opts, here, outDir); return
    case {"all", "figures"}
    otherwise, error("reproduceManriquez2026:what", "what must be all|figures|scrape");
end

lit  = loadChain(here, opts.LitTag);
dark = loadChain(here, opts.DarkTag);
fprintf("lit  %s: %d legs, t = %g d\n", opts.LitTag, lit.legs, lit.ts(end));
fprintf("dark %s: %d legs, t = %g d\n", opts.DarkTag, dark.legs, dark.ts(end));

fig1_light(outDir);
fig2_seasons(lit, outDir);
fig3_roofed(lit, dark, outDir);
fig4_evolution(lit, outDir);
fig5_outflow(lit, outDir);
fig6_fields(lit, outDir);
figX1_o2diel(lit, dark, outDir);
fig7_scrape(outDir);
fprintf("figures written to %s\n", outDir);
end

% ===================== figure 1: light forcing ===========================
function fig1_light(outDir)
% fig:seasons-light -- I(t)/I_opt for summer and winter.
% The published curves are manuscriptExperiments' lightSummer/lightWinter; the
% working set uses the SandFilter curve wired in probeChain, which is the
% manuscript's Table-B.1 curve. All three are analytic, so this figure is
% reproduced exactly, with the working-set curve added.
t = linspace(0, 3, 3000);
Ls = @(t) max(.5*(sin(2*pi*(t - 0.3)) + 1) - 0.2, 0);
Lw = @(t) max(.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
Lp = @(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50);
f = figure("Visible", "off", "Position", [0 0 720 320]);
plot(t, Ls(t), "-", LineWidth=1.6); hold on
plot(t, Lw(t), "--", LineWidth=1.6);
plot(t, Lp(t), ":", LineWidth=2.0);
xlabel("time (d)"); ylabel("I(t)/I_{opt}");
legend(["summer (published)", "winter (published)", "working set (E4/E7, 19 \circC)"], ...
    Location="northeast"); grid on; ylim([0 1.05]);
title("fig:seasons-light — effective light irradiation");
saveBoth(f, outDir, "fig1_seasons_light");
end

% ===================== figure 2: seasons =================================
function fig2_seasons(lit, outDir)
% fig:seasons-results -- phi_b after 90 d, full filter + zoom, summer vs winter.
% Summer = the E7 lit arm at t = 90 d. Winter is NOT reproduced (it needs a
% fresh 90 d run at 3 C, deliberately not started).
[phi90, t90] = profileAt(lit, 90);
[phi104, t104] = profileAt(lit, lit.ts(end));
f = figure("Visible", "off", "Position", [0 0 980 380]);
tl = tiledlayout(f, 1, 2, TileSpacing="compact");
nexttile; plot(lit.z, phi90, "-", LineWidth=1.4); hold on
plot(lit.z, phi104, "--", LineWidth=1.0);
xlabel("depth z (m)"); ylabel("\phi_b"); grid on; xlim([lit.z(1) lit.z(end)]);
legend([sprintf("summer, working set, t = %.0f d", t90), sprintf("t = %.0f d", t104)], Location="north");
title("entire filter");
nexttile; plot(lit.z, phi90, "-", LineWidth=1.4); hold on
plot(lit.z, phi104, "--", LineWidth=1.0);
xlabel("depth z (m)"); ylabel("\phi_b"); grid on; xlim([-0.02 0.10]);
title("upper layers (zoom)");
text(0.03, 0.6*max(phi90), "winter: not reproduced" + newline + "(needs a fresh 90 d run at 3 \circC)", ...
    FontAngle="italic");
title(tl, "fig:seasons-results — \phi_b after 90 d (summer only)");
saveBoth(f, outDir, "fig2_seasons_results");
end

% ===================== figure 3: covered vs uncovered ====================
function fig3_roofed(lit, dark, outDir)
% fig:roofed-results -- phi_b after 30 d, covered vs uncovered, three zooms.
% Uncovered = lit arm (LightScale 1), covered = dark arm (LightScale 0; the
% published "covered" was 0.1-1 % of full light, not exactly zero).
[pl, t1] = profileAt(lit, 30); pd = profileAt(dark, 30);
f = figure("Visible", "off", "Position", [0 0 1200 360]);
tl = tiledlayout(f, 1, 3, TileSpacing="compact");
xl = {[lit.z(1) lit.z(end)], [-0.05 0.20], [-0.008 0.02]};
nm = ["entire filter", "upper layers", "sand surface (ultra zoom)"];
for k = 1:3
    nexttile; plot(lit.z, pl, "-", LineWidth=1.4); hold on
    plot(dark.z, pd, "--", LineWidth=1.4);
    xlabel("depth z (m)"); ylabel("\phi_b"); grid on; xlim(xl{k}); title(nm(k));
    if k == 1, legend(["uncovered (LightScale 1)", "covered (LightScale 0)"], Location="north"); end
end
title(tl, sprintf("fig:roofed-results — \\phi_b at t = %.0f d, covered vs uncovered", t1));
saveBoth(f, outDir, "fig3_roofed_results");
end

% ===================== figure 4: evolution + mass ========================
function fig4_evolution(lit, outDir)
% fig:month-change-full / -zoom / -mass. The published panels show 20 d of
% evolution and a mass curve that flattens at ~30 d; the working set runs to
% 104 d, so the snapshot list is extended and the mass curve covers the lot.
days = [1 2 5 10 20 30 60 90 104];
days = days(days <= lit.ts(end));
f = figure("Visible", "off", "Position", [0 0 980 380]);
tl = tiledlayout(f, 1, 2, TileSpacing="compact");
cols = parula(numel(days));
for k = 1:2
    nexttile; hold on
    for j = 1:numel(days)
        p = profileAt(lit, days(j));
        plot(lit.z, p, LineWidth=1.2, Color=cols(j,:));
    end
    xlabel("depth z (m)"); ylabel("\phi_b"); grid on
    if k == 1, xlim([lit.z(1) lit.z(end)]); title("entire filter");
        legend(compose("%g d", days), Location="northwest", NumColumns=2);
    else, xlim([-0.02 0.20]); title("upper layers (zoom)"); end
end
title(tl, "fig:month-change-full / -zoom — \phi_b evolution, working set");
saveBoth(f, outDir, "fig4_month_change_profiles");

f = figure("Visible", "off", "Position", [0 0 720 360]);
plot(lit.ts, lit.mTot, "-", LineWidth=1.5); hold on
plot(lit.ts, lit.mBed, "--", LineWidth=1.2);
plot(lit.ts, lit.mRough, ":", LineWidth=1.6);
plot(lit.ts, lit.mSup, "-.", LineWidth=1.0);
xlabel("time (d)"); ylabel("biomass (kg/m^2)"); grid on
legend(["total column", "sand bed", "roughness layer", "supernatant"], Location="southeast");
title("fig:month-change-mass — total biofilm mass, working set");
saveBoth(f, outDir, "fig5_month_change_mass");
end

% ===================== figure 5: effluent liquids ========================
function fig5_outflow(lit, outDir)
% fig:1d-outflow-liquids -- the five liquids at the outflow.
nm = ["O2" "IC" "NH4" "HPO4" "DOM"];
f = figure("Visible", "off", "Position", [0 0 1100 560]);
tl = tiledlayout(f, 2, 3, TileSpacing="compact");
for k = 1:numel(nm)
    nexttile; j = find(lit.effNames == nm(k));
    plot(lit.ts, lit.eff(j, :)*1e3, LineWidth=1.0); hold on
    yline(lit.influent(liquidIndex(nm(k)))*1e3, "--", Color=[.5 .5 .5]);
    xlabel("time (d)"); ylabel(nm(k) + " (g/m^3)"); grid on; title(nm(k));
end
nexttile; axis off
text(0, 0.5, ["dashed: influent" newline "hourly frames" newline ...
    "diel oscillation is the light cycle"], FontAngle="italic");
title(tl, "fig:1d-outflow-liquids — effluent concentrations, working set");
saveBoth(f, outDir, "fig6_outflow_liquids");
end

function j = liquidIndex(nm)
% position in the 9-element Influent vector [HET PHO POM PAT O2 IC NH4 HPO4 DOM]
j = find(["HET" "PHO" "POM" "PAT" "O2" "IC" "NH4" "HPO4" "DOM"] == nm);
end

% ===================== figure 6: (z,t) liquid fields =====================
function fig6_fields(lit, outDir)
% fig:2d-plots -- (z,t) fields of the five liquids in the flowing suspension.
nm = ["O2" "IC" "NH4" "HPO4" "DOM"];
f = figure("Visible", "off", "Position", [0 0 1200 620]);
tl = tiledlayout(f, 2, 3, TileSpacing="compact");
for k = 1:numel(nm)
    nexttile;
    F = lit.fields.(nm(k));                       % nz x nt, daily columns
    imagesc(lit.fieldT, lit.z, F*1e3); set(gca, YDir="reverse");   % z grows downwards
    xlabel("time (d)"); ylabel("depth z (m)"); c = colorbar; c.Label.String = "g/m^3";
    title(nm(k));
end
nexttile; axis off
text(0, 0.5, ["flowing suspension" newline "daily columns" newline ...
    sprintf("t = 0 .. %.0f d", lit.ts(end))], FontAngle="italic");
title(tl, "fig:2d-plots — liquid components in the flowing suspension");
saveBoth(f, outDir, "fig7_2d_liquid_fields");
end

% ===================== figure X1: O2 diel ================================
function figX1_o2diel(lit, dark, outDir)
% X1 in analysis/manuscriptExperiments.m: O2 depth profiles at noon and midnight
% of the last simulated day, plus the diel effluent series. The covered arm is
% added because it is the control that shows how much of the diel signal is
% photosynthetic.
L = @(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50);
t = lit.o2LastT; lv = arrayfun(L, t);
last = t >= t(end) - 1;
tl_ = t(last); o2 = lit.o2Last(:, last); lvl = lv(last);
[~, kN] = max(lvl); [~, kM] = min(lvl);
f = figure("Visible", "off", "Position", [0 0 1000 380]);
tl = tiledlayout(f, 1, 2, TileSpacing="compact");
nexttile; plot(o2(:, kN)*1e3, lit.z, "-", LineWidth=1.5); hold on
plot(o2(:, kM)*1e3, lit.z, "--", LineWidth=1.5);
set(gca, YDir="reverse"); xlabel("O_2 (g/m^3)"); ylabel("depth z (m)"); grid on
legend([sprintf("noon (t = %.2f d)", tl_(kN)), sprintf("midnight (t = %.2f d)", tl_(kM))], Location="southwest");
title("O_2 depth profile, last day");
nexttile;
plot(tl_, o2(end, :)*1e3, "-", LineWidth=1.2); hold on
kd = dark.o2LastT >= dark.o2LastT(end) - 1;
plot(dark.o2LastT(kd), dark.o2Last(end, kd)*1e3, "--", LineWidth=1.2);
yyaxis right; plot(tl_, lvl, ":", LineWidth=1.4); ylabel("I(t)/I_{opt}");
yyaxis left; xlabel("time (d)"); ylabel("effluent O_2 (g/m^3)"); grid on
legend(["uncovered", "covered", "light"], Location="best");
title("diel effluent O_2");
title(tl, sprintf("X1 — O_2 diel structure at t = %.0f d (working set)", lit.ts(end)));
saveBoth(f, outDir, "fig10_X1_o2_diel");
end

% ===================== figure 7: scraping ================================
function fig7_scrape(outDir)
% fig:month-scrape-0/4/8/12 and -mass, from the scrape runs if they exist.
g = dir(fullfile(outDir, "data", "scrape_GP*.mat"));
if isempty(g)
    fprintf("no scrape_GP*.mat in %s -- scraping figures skipped\n", fullfile(outDir, "data"));
    return
end
S = struct([]); cm = [];
for k = 1:numel(g)
    D = load(fullfile(g(k).folder, g(k).name), "sc");
    S(k).sc = D.sc; cm(k) = D.sc.depth*100; %#ok<AGROW>
end
[cm, o] = sort(cm); S = S(o);
f = figure("Visible", "off", "Position", [0 0 1200 620]);
tl = tiledlayout(f, 2, 2, TileSpacing="compact");
for k = 1:numel(S)
    nexttile; sc = S(k).sc;
    ks = round(linspace(1, numel(sc.ts), min(6, numel(sc.ts))));
    hold on
    for j = ks, plot(sc.z, sc.phi(:, j), LineWidth=1.1); end
    xlabel("depth z (m)"); ylabel("\phi_b"); grid on; xlim([-0.02 0.25]);
    legend(compose("%.1f d", sc.ts(ks)), Location="northeast");
    title(sprintf("scraped %.0f cm", cm(k)));
end
title(tl, "fig:month-scrape — \phi_b after scraping a 30 d filter");
saveBoth(f, outDir, "fig8_month_scrape_profiles");

f = figure("Visible", "off", "Position", [0 0 720 360]);
hold on
for k = 1:numel(S)
    plot(S(k).sc.ts, S(k).sc.mass, LineWidth=1.4);
end
xlabel("time after scraping (d)"); ylabel("biomass (kg/m^2)"); grid on
legend(compose("scraped %.0f cm", cm), Location="southeast");
title("fig:month-scrape-mass — regrowth after scraping");
saveBoth(f, outDir, "fig9_month_scrape_mass");
end

% ===================== scraping runs =====================================
function runScrape(opts, here, outDir)
% One scraping experiment: take the end state of leg `FromLeg` of the lit chain
% (leg 3 = t = 30 d), remove the biofilm above `Depth`, run `Days` of regrowth.
%
% probeChain has no scraping option and must not be edited, so the filter and
% model are rebuilt here with the E4 settings, VERBATIM from probeChain.m
% (2026-08-27). If probeChain's construction changes, this copy must follow.
src = fullfile(here, "probes", "data", "chain", ...
    sprintf("chain_%s_leg%d.mat", opts.LitTag, opts.FromLeg));
D = load(src, "snap", "rec");
st = D.snap; t0 = D.rec.tEnd;
fprintf("scrape: from %s (t = %g d), depth %.3f m, %g d of regrowth\n", src, t0, opts.Depth, opts.Days);

[f, m, common] = workingSetSetup(opts.NCells);
z = f.GridPoints.Centers(:); dz = f.GridSize;
st = scrapeState(st, opts.Depth, dz, z);
s = rehomeState(f, m, st, 0.0);
tic; r = simulate(s, common{:}, SimulationTime=opts.Days, AdaptiveInitialDt=1e-8);
wall = toc;
checkRunFlag(r, sprintf("scrape %.0f cm", opts.Depth*100));

C = r.Frames.Concentrations;
dL = mean([m.Liquids.Density]); dP = mean([m.Particles.Density]);
phi = C{"Water","Enclosed"}{1}/dL;
for nm = [m.Particles.Name], phi = phi + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [m.Liquids.Name],   phi = phi + C{nm,"Enclosed"}{1}/dL; end
eps = computePorosity(f, z);
tot = 0;
for nm = ["HET" "PHO" "POM"], tot = tot + C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1}; end
sc = struct("depth", opts.Depth, "fromT", t0, "ts", r.Frames.Time(:).', "z", z, ...
    "eps", eps(:), "dz", dz, "delta", f.SandRoughness, "phi", phi, ...
    "mass", sum(eps(:).*tot, 1)*dz, "wall", wall, "flag", string(r.Flag), ...
    "effO2", C{"O2","Flowing"}{1}(end, :)); %#ok<NASGU>
fn = fullfile(outDir, "data", sprintf("scrape_GP%d.mat", round(opts.Depth*100)));
save(fn, "sc", "-v7.3");
fprintf("scrape %.0f cm: %s, %.0f s -> %s\n", opts.Depth*100, r.Flag, wall, fn);
end

function st = scrapeState(st, zs, dz, centers)
% Port of the private scrapeState in analysis/manuscriptExperiments.m, itself a
% port of the legacy @SDresults/scrape.m (refill-sand family): cells
% 0 <= z <= ceil(zs/dz + 1/2)*dz lose biofilm and enclosed water, the flowing
% suspension is emptied everywhere, velocities dropped, clock reset.
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

function [f, m, common] = workingSetSetup(N)
% E4/E7 working set, copied from analysis/probes/probeChain.m (do not diverge).
f = SandFilter(Temperature=19, LightAttenuationCoeffSand=1500, SandRoughness=5e-3, ...
    LightIrradiation=@(t) 1*0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(N);
mp = pathogenModel(PhototrophRespiration=0.0, RespirationForm="reichert", NormalizedLight=true);
rx = mp.Reactions; rx([rx.Name] == "Phototroph growth").MinimumLightFactor = 0.0;
iH = find([rx.Name] == "Heterotroph growth"); H = rx(iH).HalfSaturationConstants;
H("DOM") = 3e-4; rx(iH).HalfSaturationConstants = H;
iP = find([rx.Name] == "Phototroph growth"); H = rx(iP).HalfSaturationConstants;
H("HPO4") = 1e-6; rx(iP).HalfSaturationConstants = H;
mp.Reactions = rx;
cs = mp.Components;
for q = 1:numel(cs)
    if isa(cs(q), "Liquid"), cs(q).TransportRate = 10*cs(q).TransportRate; end
end
mp.Components = cs;
m = Model(mp.Components, mp.Reactions, Kappa=1e-6, Zeta0=1, Zeta1=0.27, ...
    DetachmentFunction=@(v) 0.14*abs(v)/18, WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
common = {"InflowConcentrations", [3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3], ...
    "TimeStep", "adaptive", "AdaptiveMaxDt", 5e-5, "ImplicitOsmosis", true, ...
    "CohesionBC", "neumann", "CohesionScheme", "shin", "TransferForm", "constant", ...
    "FrameNumber", 241, "Quiet", true};
end

% ===================== chain reading =====================================
function A = loadChain(here, tag)
d = fullfile(here, "probes", "data", "chain");
g = dir(fullfile(d, sprintf("chain_%s_leg*.mat", tag)));
g = g(~contains(string({g.name}), "ABORTED"));
if isempty(g), error("reproduceManriquez2026:noChain", "no legs for %s", tag); end
n = arrayfun(@(x) str2double(regexp(string(x.name), "leg(\d+)\.mat", "tokens", "once")), g);
[~, o] = sort(n); g = g(o);
A.tag = tag; A.legs = numel(g);
ts = []; phi = []; eff = []; mS = []; mR = []; mB = [];
F = struct(); FT = [];
liq = ["O2" "IC" "NH4" "HPO4" "DOM"];
for nm = liq, F.(nm) = []; end
prevEnd = -inf;
for k = 1:numel(g)
    D = load(fullfile(g(k).folder, g(k).name), "rec", "results_py");
    rec = D.rec; py = D.results_py;
    t = rec.ts(:).';
    if k > 1 && abs(t(1) - 2*prevEnd) < abs(t(1) - prevEnd), t = t - prevEnd; end
    keep = 1:numel(t); if k > 1, keep = 2:numel(t); end
    z = py.z(:); eps = py.porosity(:); dz = py.dz; delta = py.delta;
    pN = string(strsplit(py.particleNames, '|')); lN = string(strsplit(py.liquidNames, '|'));
    tot = 0;
    for nm = ["HET" "PHO" "POM"], j = pN == nm; tot = tot + py.particles(:,:,j,1) + py.particles(:,:,j,2); end
    w = eps .* tot;
    sup = z < -delta; rough = z >= -delta & z < 0; bed = z >= 0;
    ts  = [ts,  t(keep)];                   %#ok<AGROW>
    phi = [phi, rec.phiT(:, keep)];         %#ok<AGROW>
    eff = [eff, rec.effluent(:, keep)];     %#ok<AGROW>
    mS  = [mS,  sum(w(sup, keep), 1)*dz];   %#ok<AGROW>
    mR  = [mR,  sum(w(rough, keep), 1)*dz]; %#ok<AGROW>
    mB  = [mB,  sum(w(bed, keep), 1)*dz];   %#ok<AGROW>
    % daily columns of the flowing liquid fields, for the (z,t) figure
    kd = keep(1:24:numel(keep));
    for nm = liq
        j = lN == nm;
        F.(nm) = [F.(nm), py.liquids(:, kd, j, 2)];
    end
    FT = [FT, t(kd)];                       %#ok<AGROW>
    if k == numel(g)                        % hourly O2 field of the last leg, for X1
        A.o2Last = py.liquids(:, :, lN == "O2", 2); A.o2LastT = t;
    end
    prevEnd = t(end);
    A.z = z; A.eps = eps; A.dz = dz; A.delta = delta;
    A.effNames = rec.effNames; A.influent = rec.influent;
end
A.ts = ts; A.phi = phi; A.eff = eff;
A.mSup = mS; A.mRough = mR; A.mBed = mB; A.mTot = mS + mR + mB;
A.fields = F; A.fieldT = FT;
end

function [p, tGot] = profileAt(A, tWant)
[~, k] = min(abs(A.ts - tWant));
p = A.phi(:, k); tGot = A.ts(k);
end

function saveBoth(f, outDir, name)
exportgraphics(f, fullfile(outDir, name + ".png"), Resolution=150);
exportgraphics(f, fullfile(outDir, name + ".pdf"), ContentType="vector");
close(f);
fprintf("  %s.{png,pdf}\n", name);
end
