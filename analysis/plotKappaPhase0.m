function plotKappaPhase0(opts)
% PLOTKAPPAPHASE0  Figures for the 2026-08-23 kappa fix and Bailo Phase 0.
%
% NOTE: none of these are Bailo results. The Bailo scheme is NOT implemented.
% Every figure measures the EXISTING solver -- Phase 0 exists precisely to decide
% whether building Bailo is justified, and on this evidence three of its four
% justifications are dead. Bailo results would be the Phase 1 gates (mass,
% bounds, monotone energy, second-order accuracy, and bounds holding at 100x the
% explicit dt limit); none have been run.
%
%   plotKappaPhase0
%
% Produces, into analysis/results/figures/:
%   fig1_cfl_margins   which CFL region binds, and by how much, vs mesh
%   fig2_cohesion_role cohesive share of |v_b|, speedup ceiling X, cap-bound frac
%   fig3_dt_scaling    0e: Solver A phi_b vs diagnostic sum, first order in dt
%   fig4_mms           MMS convergence of the restored kappa operator
%   fig5_kappa_sizing  interface width vs mesh: which (kappa, N) are resolved
%
% Data: analysis/probes/data/cflbudget_*.mat (probeCflBudget) and
%       analysis/probes/data/mms_cahnhilliard.mat (probeMMSCahnHilliard).

arguments
    opts.OutDir (1,1) string = ""
end

here = fileparts(mfilename("fullpath"));
D = fullfile(here, "probes", "data");
if opts.OutDir == "", opts.OutDir = fullfile(here, "results", "figures"); end
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

set(groot, "defaultAxesFontSize", 11, "defaultLineLineWidth", 1.6);
regions = ["matrix", "enclosed P", "enclosed L", "flowing P", "flowing L"];
cols = lines(5);

% ---- load the margin sweep (N = 100, 200, 500 at MaxDt = 3e-6) --------------
Ns = [100 200 500];
share = nan(numel(Ns), 5); shareMax = nan(numel(Ns), 5);
mNoCoh = nan(1, numel(Ns)); capF = mNoCoh; Xm = mNoCoh; cohV = mNoCoh; cohVmax = mNoCoh;
termF = nan(numel(Ns), 4);
for k = 1:numel(Ns)
    f = fullfile(D, sprintf("cflbudget_n%d_k1e-06_3d_dt3e-06.mat", Ns(k)));
    if ~isfile(f), warning("missing %s", f); continue, end
    r = load(f).res; b = r.budget; n = b.Steps;
    share(k,:)    = b.RegionShareSum/n;
    shareMax(k,:) = b.RegionShareMax;
    mNoCoh(k) = b.MatrixShareNoCoh/n;
    capF(k)   = r.capFrac;
    Xm(k)     = r.Xmean;
    cohV(k)   = b.VbCohSum/n;  cohVmax(k) = b.VbCohMax;
    termF(k,:) = r.termFrac;
end

% ============================ FIG 1: CFL margins =============================
f1 = figure(Position=[100 100 980 400]);
tl = tiledlayout(f1, 1, 2, TileSpacing="compact", Padding="compact");

nexttile; hold on
for j = 1:5
    plot(Ns, share(:,j), "-o", Color=cols(j,:), DisplayName=regions(j));
end
plot(Ns, mNoCoh, "--s", Color=cols(1,:)*0.5, DisplayName="matrix, cohesion removed");
yline(1, "k:", "binding", LabelHorizontalAlignment="left", HandleVisibility="off");
set(gca, XScale="log"); xticks(Ns); xticklabels(string(Ns));
xlabel("cells N"); ylabel("region total / binding total");
title("Which region sets the timestep"); ylim([0 1.08]);
legend(Location="east", FontSize=8); grid on; box on

nexttile; hold on
b = bar(categorical(string(Ns), string(Ns)), termF, "stacked");
names = ["w_v/dz (advection+cohesion)", "w_a/dz^2 (dispersion)", ...
         "w_b (exchange)", "w_s (ecology)"];
for j = 1:4, b(j).DisplayName = names(j); end
xlabel("cells N"); ylabel("share of the binding region's sum");
title("What the binding region is made of"); ylim([0 1]);
legend(Location="southoutside", FontSize=8, NumColumns=2); grid on; box on
title(tl, "CURRENT scheme: CFL budget (\kappa = 10^{-6}, \zeta_0 = 10^2, 3 d, active regime)");
subtitle(tl, "Bailo Phase 0 decision data -- the Bailo scheme is NOT implemented; nothing here is a Bailo result", ...
    FontSize=9, Interpreter="none");
exportgraphics(f1, fullfile(opts.OutDir, "fig1_cfl_margins.png"), Resolution=200);

% ======================= FIG 2: cohesion's role ==============================
f2 = figure(Position=[100 100 980 380]);
tl = tiledlayout(f2, 1, 3, TileSpacing="compact", Padding="compact");

nexttile; hold on
plot(Ns, cohV, "-o", DisplayName="mean");
plot(Ns, cohVmax, "--s", DisplayName="max");
set(gca, XScale="log"); xticks(Ns); xticklabels(string(Ns));
xlabel("cells N"); ylabel("cohesive share of max|v_b|");
title("Cohesion IS most of v_b"); ylim([0 1]); legend(Location="northwest"); grid on; box on

nexttile
plot(Ns, Xm, "-o", Color=[0.85 0.2 0.2]);
yline(1, "k:", HandleVisibility="off");
set(gca, XScale="log"); xticks(Ns); xticklabels(string(Ns));
xlabel("cells N"); ylabel("X = dt_{CFL}(no cohesion)/dt_{CFL}");
title("...but the speedup ceiling is 1"); ylim([0.98 1.10]); grid on; box on

nexttile
plot(Ns, 100*capF, "-o", Color=[0.2 0.4 0.8]);
set(gca, XScale="log"); xticks(Ns); xticklabels(string(Ns));
xlabel("cells N"); ylabel("% of steps pinned to MaxDt");
title("CFL only binds at fine meshes"); ylim([0 100]); grid on; box on
title(tl, "Current scheme: cohesion is dominant in v_b yet irrelevant to dt");
subtitle(tl, "measured on the EXISTING solver; Bailo is unimplemented", FontSize=9, Interpreter="none");
exportgraphics(f2, fullfile(opts.OutDir, "fig2_cohesion_role.png"), Resolution=200);

% ===================== FIG 3: 0e first-order dt scaling ======================
dts = [3e-6, 7.5e-7]; errs = nan(1,2);
for k = 1:2
    f = fullfile(D, sprintf("cflbudget_n100_k1e-06_3d_dt%g.mat", dts(k)));
    if isfile(f), errs(k) = load(f).res.phibConsistencyAbs; end
end
f3 = figure(Position=[100 100 500 430]);
plot(dts, errs, "o-", MarkerFaceColor="auto", DisplayName="measured"); hold on
plot(dts, errs(1)*(dts/dts(1)), "k--", DisplayName="slope 1 (first order)");
% loglog() alone left these linear -- only two points spanning well under a
% decade, so MATLAB's automatic ruler picked a linear multiplier style. Force it.
set(gca, XScale="log", YScale="log");
xlim([5e-7 4e-6]); ylim([4e-6 4e-5]);
xlabel("MaxDt [d]"); ylabel("max |\phi_b^{SolverA} - \phi_b^{sum}|");
title(sprintf("0e: splitting offset is first order in dt (ratio %.3f)", errs(1)/errs(2)));
subtitle("two points; dt ratio 4.00, error ratio 3.997", FontSize=9, Interpreter="none");
legend(Location="southeast"); grid on; box on
exportgraphics(f3, fullfile(opts.OutDir, "fig3_dt_scaling.png"), Resolution=200);

% =========================== FIG 4: MMS ======================================
fm = fullfile(D, "mms_cahnhilliard.mat");
if isfile(fm)
    m = load(fm).res;
    f4 = figure(Position=[100 100 520 420]);
    loglog(m.dz, m.errInterior, "o-", DisplayName="interior"); hold on
    loglog(m.dz, m.errAll, "s--", DisplayName="incl. Neumann rows");
    ref = m.errInterior(1)*(m.dz/m.dz(1)).^2;
    loglog(m.dz, ref, "k:", DisplayName="slope 2");
    set(gca, XDir="reverse");
    xlabel("\Deltaz [m]"); ylabel("max error in \kappa-operator");
    title(sprintf("MMS: restored \\kappa operator, order %.3f", m.orderInterior(end)));
    legend(Location="southeast"); grid on; box on
    exportgraphics(f4, fullfile(opts.OutDir, "fig4_mms.png"), Resolution=200);
end

% ====================== FIG 5: kappa / mesh sizing ===========================
z1 = 1e-2; psi2 = 3*0.3*(0.3 - z1);          % Psi'' at phi = 0.3
kappas = [1e-7 1e-6 1e-5]; Nsweep = round(logspace(1.5, 3.7, 60));
dzs = 1./(Nsweep + 0.5);
f5 = figure(Position=[100 100 560 420]); hold on
for k = 1:numel(kappas)
    ell = sqrt(kappas(k)/psi2);
    plot(Nsweep, dzs/ell, "-", DisplayName=sprintf("\\kappa = 10^{%d}", log10(kappas(k))));
end
yline(1, "k--", "resolved", HandleVisibility="off");
for N = [100 200 500], xline(N, ":", "N="+string(N), HandleVisibility="off", FontSize=8); end
set(gca, XScale="log", YScale="log");
xlabel("cells N"); ylabel("\Deltaz / L_{int}");
subtitle("L_int = sqrt(kappa/Psi_zz),  Psi_zz = 3*phi*(phi-zeta_1) = 0.26 at phi = 0.3", ...
    FontSize=9, Interpreter="none");
title("Which (\kappa, N) actually resolve the interface");
legend(Location="southwest"); grid on; box on
exportgraphics(f5, fullfile(opts.OutDir, "fig5_kappa_sizing.png"), Resolution=200);

fprintf("wrote figures to %s\n", opts.OutDir);
d = dir(fullfile(opts.OutDir, "*.png"));
for k = 1:numel(d), fprintf("   %s (%.0f kB)\n", d(k).name, d(k).bytes/1024); end
close all
end
