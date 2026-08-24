function compareSchemeValues(opts)
% COMPARESCHEMEVALUES  Compare the SOLUTION FIELDS of the solver configurations,
% not just their cost. Runs A (shin + explicit Solver B) and B (shin + implicit
% dispersion), optionally C (bailo + implicit dispersion), keeps the full frames,
% and plots the values against each other.
%
% The cost probes save only summary statistics, so this exists to answer "do the
% schemes actually give the same answer" in the terms a modeller cares about --
% phi_b profiles, biomass, effluent -- rather than a single max-relative number.

arguments
    opts.NCells (1,1) double = 100
    opts.Kappa (1,1) double = 1e-6
    opts.Days (1,1) double = 3
    opts.MaxDt (1,1) double = 1e-2
    opts.WithBailo (1,1) logical = false
    opts.OutDir (1,1) string = ""
end
if opts.Days < 3
    error("compareSchemeValues:inactiveRegime", ...
        "Days = %g is below the 3 d minimum: no biofilm forms.", opts.Days);
end
here = fileparts(mfilename("fullpath"));
W = fileparts(here);
addpath(genpath(fullfile(W,"src")));
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

names  = ["A shin+explicit", "B shin+implicitDisp"];
scheme = ["shin", "shin"];  impdis = [false, true];
if opts.WithBailo
    names(3) = "C bailo+implicitDisp"; scheme(3) = "bailo"; impdis(3) = true;
end

R = cell(1, numel(names));
for k = 1:numel(names)
    t0 = tic;
    R{k} = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
        TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
        FrameNumber=72, ImplicitOsmosis=true, ImplicitDispersion=impdis(k), ...
        CohesionScheme=scheme(k), Quiet=true);
    fprintf("  %-22s %s  %.0f s\n", names(k), R{k}.Flag, toc(t0));
end

z  = f.GridPoints.Centers(:);
dz = f.GridSize;
t  = R{1}.Frames.Time(:);
phi = cellfun(@(r) biofilmFrac(r), R, UniformOutput=false);
mass = cellfun(@(r) totalBio(r, dz), R, UniformOutput=false);

% Figures are created invisible: on a headless cluster node figure() takes the
% GUI path and MATLAB dies with a fatal error in PFGuiApplication/GuiThread --
% job 3531774 completed all three 3 d simulations at N=500 and then crashed in
% plotting, losing the tables. Locally there is a display so it never showed up.
set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.5, ...
    "defaultFigureVisible", "off");
% A and B agree to ~4e-03, so drawn identically B hides A entirely and the figure
% looks like a single curve. Distinct widths/styles so the overlap is visibly an
% overlap rather than a missing series.
lw = [3.0 1.6 1.6]; ls = ["-" "--" ":"];
cols = lines(3);
fig = figure(Position=[80 80 1150 760], Visible="off");
tl = tiledlayout(fig, 2, 3, TileSpacing="compact", Padding="compact");

nexttile; hold on                                     % phi_b final profile
for k = 1:numel(R), plot(phi{k}(:,end), z, Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k), DisplayName=names(k)); end
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off");
xlabel("\phi_b [-]"); ylabel("z [m]"); title(sprintf("\\phi_b at t = %.1f d", opts.Days));
legend(Location="southeast", FontSize=8); grid on; box on

nexttile; hold on                                     % phi_b, supernatant zoom
sup = z < 0;
for k = 1:numel(R), plot(phi{k}(sup,end), z(sup), Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k)); end
set(gca, YDir="reverse"); xlabel("\phi_b [-]"); ylabel("z [m]");
title("supernatant (z < 0), where cohesion acts"); xlim([0 inf]); grid on; box on

nexttile; hold on                                     % difference profile
for k = 2:numel(R)
    plot(phi{k}(:,end) - phi{1}(:,end), z, Color=cols(k,:), DisplayName=names(k)+" - A");
end
set(gca, YDir="reverse"); xline(0, "k:", HandleVisibility="off");
xlabel("\Delta\phi_b vs A"); ylabel("z [m]"); title("pointwise difference from A");
legend(Location="southeast", FontSize=8); grid on; box on

nexttile; hold on                                     % biomass
for k = 1:numel(R), plot(t, mass{k}, Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k), DisplayName=names(k)); end
xlabel("t [d]"); ylabel("total biomass [kg/m^2]"); title("total biomass");
legend(Location="northwest", FontSize=8); grid on; box on

nexttile; hold on                                     % effluent O2
for k = 1:numel(R)
    o2 = R{k}.Frames.Concentrations{"O2","Flowing"}{1};
    plot(t, o2(end,:), Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k), DisplayName=names(k));
end
xlabel("t [d]"); ylabel("effluent O_2 [kg/m^3]"); title("effluent O_2");
legend(Location="best", FontSize=8); grid on; box on

nexttile                                              % per-component max rel diff
C1 = R{1}.Frames.Concentrations;
nm = string(C1.Properties.RowNames); ph = string(C1.Properties.VariableNames);
D = zeros(numel(nm), numel(R)-1);
for k = 2:numel(R)
    Ck = R{k}.Frames.Concentrations;
    for i = 1:numel(nm)
        d = 0;
        for j = 1:numel(ph)
            xa = C1{nm(i),ph(j)}{1}; xb = Ck{nm(i),ph(j)}{1};
            if isnumeric(xa) && isequal(size(xa),size(xb)) && ~isscalar(xa)
                den = max(abs(xa(:)));
                if den > 0, d = max(d, max(abs(xa(:)-xb(:)))/den); end
            end
        end
        D(i,k-1) = d;
    end
end
b = barh(categorical(nm, nm), D); set(gca, XScale="log");
for k = 1:size(D,2), b(k).DisplayName = names(k+1)+" vs A"; end
xlabel("max relative difference"); title("per component"); legend(Location="best", FontSize=8);
grid on; box on

title(tl, sprintf("Solver configurations, values compared (N = %d, \\kappa = %g, %.0f d)", ...
    opts.NCells, opts.Kappa, opts.Days));
out = fullfile(opts.OutDir, sprintf("fig6_scheme_values_n%d.png", opts.NCells));
exportgraphics(fig, out, Resolution=200);
fprintf("\nwrote %s\n", out);

% ---------- per-component profiles and v_b ----------------------------------
% phi_b is a sum, so agreement there can hide cancelling per-component errors.
% Plot each component separately, and v_b, which is the ONLY thing Solver A hands
% to Solver B and therefore the channel through which the CH scheme can act.
parts = [R{1}.Model.Particles.Name];  liqs = [R{1}.Model.Liquids.Name];
panels = [ compose("%s|Matrix", parts), compose("%s|Flowing", liqs), "Water|Enclosed" ];
np = numel(panels);
fig2 = figure(Position=[60 60 1250 780], Visible="off");
tl2 = tiledlayout(fig2, 3, ceil((np+1)/3), TileSpacing="compact", Padding="compact");
devC = nan(1, np);
for i = 1:np
    pp = split(panels(i), "|"); nm = pp(1); ph = pp(2);
    nexttile; hold on
    xa = R{1}.Frames.Concentrations{nm, ph}{1};
    for k = 1:numel(R)
        xk = R{k}.Frames.Concentrations{nm, ph}{1};
        plot(xk(:,end), z, Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k), ...
             DisplayName=names(k));
    end
    xc = R{end}.Frames.Concentrations{nm, ph}{1};
    den = max(abs(xa(:,end)));
    if den > 0, devC(i) = max(abs(xc(:,end) - xa(:,end)))/den; end
    set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off");
    title(sprintf("%s / %s   (C-A %.1e)", nm, ph, devC(i)), FontSize=9);
    xlabel("[kg/m^3]"); if mod(i-1, ceil((np+1)/3)) == 0, ylabel("z [m]"); end
    grid on; box on
end

nexttile; hold on                                   % v_b, Solver A's only output
zf = f.GridPoints.Boundaries(2:end-1);
for k = 1:numel(R)
    vb = R{k}.Frames.Velocity.Biofilm;
    plot(vb(:,end), zf, Color=cols(k,:), LineWidth=lw(k), LineStyle=ls(k), ...
         DisplayName=names(k));
end
set(gca, YDir="reverse"); yline(0, "k:", HandleVisibility="off"); xline(0, "k:", HandleVisibility="off");
xlabel("v_b [m/d]"); title("v_b -- Solver A's only output", FontSize=9);
legend(Location="best", FontSize=7); grid on; box on
title(tl2, sprintf("Per-component profiles at t = %.0f d (N = %d, \kappa = %g)", ...
    opts.Days, opts.NCells, opts.Kappa));
out2 = fullfile(opts.OutDir, sprintf("fig7_components_n%d.png", opts.NCells));
exportgraphics(fig2, out2, Resolution=200);
fprintf("wrote %s\n", out2);

fprintf("\nC vs A, per component (max relative on the final profile):\n");
[sv, si] = sort(devC, "descend", MissingPlacement="last");
for i = 1:np
    if ~isnan(sv(i)), fprintf("   %-16s %.3e\n", panels(si(i)), sv(i)); end
end

fprintf("\n%-22s %12s %12s %12s\n", "config", "max phi_b", "biomass(end)", "vs A (max rel)");
for k = 1:numel(R)
    if k == 1, s = "-"; else, s = sprintf("%.3e", max(D(:,k-1))); end
    fprintf("%-22s %12.5f %12.5f %12s\n", names(k), max(phi{k}(:,end)), mass{k}(end), s);
end
close all
end

function p = biofilmFrac(r)
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
p = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], p = p + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   p = p + C{nm,"Enclosed"}{1}/dL; end
end

function mss = totalBio(r, dz)
C = r.Frames.Concentrations; mss = 0;
for nm = [r.Model.Particles.Name]
    mss = mss + sum(C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1}, 1)*dz;
end
end
