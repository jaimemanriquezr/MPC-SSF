function res = probeCflBudget(opts)
% PROBECFLBUDGET  Phase 0 of the Bailo investigation: is the scheme worth building?
%
%   probeCflBudget
%   probeCflBudget(NCells=500, Kappa=1e-6, MatureDays=3)
%
% See .claude/plans/2026-08-23-bailo-scheme.md, "Phase 0". Answers:
%
%   0a  Is dt actually CFL-bound, or pinned to AdaptiveMaxDt? If the cap binds,
%       every other number here is meaningless. (This exact trap produced a
%       worthless kappa-fix dt measurement on 2026-08-23.)
%   0b  Which region attains the max in
%           dt_CFL = cflFactor / max(w_v/dz + w_a/dz^2 + w_b + w_s),
%       and what share of that region's SUM does each term carry? The four
%       weights add within a region, so this is a share-of-sum, not a max.
%   0c  THE DECISION NUMBER. v_b = q - zeta_0(1-phi)grad(mu). Bailo would treat
%       the cohesive term implicitly but cannot touch q, so the best it can ever
%       do is delete the cohesive part of vbmax. X = dt_CFL(advective v_b only)
%       / dt_CFL(actual) is the resulting speedup ceiling. X ~ 1 => no speed
%       case, whatever the scheme.
%
% REGIME. Must run where cohesion is active, i.e. with biofilm in the
% supernatant. A clean start has phi_b ~ 0, the cohesive term is inert, and 0c
% returns X ~ 1 for reasons that have nothing to do with the scheme. MatureDays
% defaults to 3 (phib_sup ~ 0.14, verified active).
%
% See also PROBEMMSCAHNHILLIARD, MANUSCRIPTEXPERIMENTS.

arguments
    opts.NCells (1,1) double = 500
    opts.Kappa (1,1) double = 1e-6
    opts.MatureDays (1,1) double = 3
    opts.MaxDt (1,1) double = 3e-6
    opts.Respiration (1,1) double = 0.55
    opts.CohesionScheme (1,1) string = "shin"
    opts.PositiveMobility (1,1) logical = false
    opts.Save (1,1) logical = true
end

here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end
% Below 3 d there is no biofilm formation, so every number here describes the
% startup transient rather than the model. This produced four wrong readings in
% one session (kappa at 0.05 d, CFL budget at 0.3 d, implicit dispersion at
% 0.02 d, w_s spread at 0.3 d), each a plausible-looking number rather than an
% obvious error. Refuse rather than warn: a warning scrolls past.
if opts.MatureDays < 3
    error("probeCflBudget:inactiveRegime", ...
        "MatureDays = %g d is below the 3 d minimum: no biofilm has formed, so " + ...
        "this measures the startup transient. Phase 0 measured phib_sup = 0 at " + ...
        "0.3 d against 0.1377 at 3 d. Pass MatureDays >= 3, or use a " + ...
        "state-independent operator-level check instead.", opts.MatureDays);
end


f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=opts.Respiration, PGExcess=true, NormalizedLight=true);
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=mp.CohesionSubModel.Zeta0, Zeta1=mp.CohesionSubModel.Zeta1, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), WaterDensity=mp.WaterDensity, ...
    BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);

% Table B.1, exactly as manuscriptExperiments/baseInfluent(0) builds it, so this
% probe measures E3's regime and not a regime of its own invention.
%          HET       PHO      POM  PAT  O2       IC       NH4      HPO4 DOM
inflow = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];

fprintf("\n=== Phase 0: CFL budget ===\n");
fprintf("N = %d, kappa = %g, zeta_0 = %g, %.2f d, MaxDt = %g, PositiveMobility = %d\n\n", ...
    opts.NCells, opts.Kappa, m.CohesionSubModel.Zeta0, opts.MatureDays, opts.MaxDt, opts.PositiveMobility);

t0 = tic;
% Same call as manuscriptExperiments/runSim, plus the budget recorder.
r = simulate(State(f, m), InflowConcentrations=inflow, ...
    SimulationTime=opts.MatureDays, TimeStep="adaptive", ...
    AdaptiveInitialDt=1e-8, AdaptiveMaxDt=opts.MaxDt, ...
    FrameNumber=24, ImplicitOsmosis=true, RecordCflBudget=true, ...
    PositiveMobility=opts.PositiveMobility, CohesionScheme=opts.CohesionScheme, Quiet=true);
wall = toc(t0);

b = r.SimulationData.CflBudget;
n = b.Steps;
regions = ["matrix", "enclosed P", "enclosed L", "flowing P", "flowing L"];
terms   = ["w_v/dz (advection+cohesion)", "w_a/dz^2 (dispersion)", ...
           "w_b (exchange)", "w_s (ecology)"];

% --- regime check: was cohesion active at all? -------------------------------
phiB = biofilmFractionLocal(r);
sup  = f.GridPoints.Centers < 0;
phibSup = max(phiB(sup, end));
fprintf("REGIME  phib_sup(final) = %.4f   flag = %s\n", phibSup, r.Flag);
if phibSup < 1e-3
    fprintf("  *** phib_sup ~ 0: cohesion INACTIVE. Numbers below are a null test. ***\n");
end

fprintf("\n0a  steps = %d, wall = %.1f s\n", n, wall);
fprintf("    dt pinned to AdaptiveMaxDt: %.1f%% of steps  (CFL-bound: %.1f%%)\n", ...
    100*b.CapBound/n, 100*(1 - b.CapBound/n));

fprintf("\n0b  region attaining the max:\n");
for k = 1:5
    fprintf("      %-12s %6.2f%%\n", regions(k), 100*b.RegionWins(k)/n);
end
fprintf("    MARGIN -- each region's total as a fraction of the binding one:\n");
for k = 1:5
    fprintf("      %-12s mean %7.4f   worst-case %7.4f\n", regions(k), ...
        b.RegionShareSum(k)/n, b.RegionShareMax(k));
end
fprintf("      matrix region with cohesion removed: mean %7.4f\n", b.MatrixShareNoCoh/n);
fprintf("      (matrix/binding scales ~dz^-2, so x25 from N=100 to N=500)\n");
fprintf("    mean share of that region's sum:\n");
for k = 1:4
    fprintf("      %-28s %6.2f%%\n", terms(k), 100*b.TermSums(k)/n);
end

if opts.CohesionScheme == "bailo"
    bn = r.SimulationData.BailoNewton;
    fprintf("\n0h  Bailo semismooth Newton cost\n");
    fprintf("      mean %.2f it/step   max %d   non-converged %.3f%%   zero-iter %.1f%%\n", ...
        bn.Mean, bn.Max, 100*bn.NonConvergedFrac, 100*bn.ZeroIterFrac);
    fprintf("      (zero-iter steps are those where the residual already met tolerance,\n");
    fprintf("       i.e. nothing was happening anywhere; a high fraction means the mean\n");
    fprintf("       understates the cost where the term is actually active)\n");
    if bn.NonConvergedFrac > 0
        fprintf("      *** %d steps hit the 50-iteration cap: those did NOT solve the\n", ...
            round(bn.NonConvergedFrac*bn.Steps));
        fprintf("          implicit system. Treat all results from this run as suspect. ***\n");
    end
end

fprintf("\n0g  w_s spread ACROSS CELLS (is a few depleted cells setting dt?)\n");
fprintf("      per-cell ecology weight in the enclosed-liquid region, mean over steps:\n");
fprintf("        p50 %.4e   p90 %.4e   p99 %.4e   max %.4e\n", ...
    b.WsP50Sum/n, b.WsP90Sum/n, b.WsP99Sum/n, b.WsCellMaxSum/n);
% Concentration as a RATIO OF TIME-MEANS, not a mean of per-step ratios. The
% latter is dominated by startup steps where the median is ~1e-12 but nonzero,
% giving meaningless 1e17 values while the time-averaged profile is flat.
concRatio = (b.WsCellMaxSum/n)/max(b.WsP50Sum/n, realmin);
conc90    = (b.WsP90Sum/n)/max(b.WsP50Sum/n, realmin);
fprintf("      max/median %.2f   p90/median %.2f   (ratios of time-means)\n", concRatio, conc90);
fprintf("      bound actually used / per-cell max: %.3f", (b.WsBoundSum/n)/max(b.WsCellMaxSum/n, eps));
fprintf("   (>1 = the two maxima land in different cells)\n");
if concRatio > 10
    fprintf("      => HIGHLY concentrated: dt is set by a small minority of cells.\n");
    fprintf("         A tighter bound or local time-stepping would be far cheaper\n");
    fprintf("         than making reactions implicit.\n");
else
    fprintf("      => broadly distributed: the reaction stiffness is a genuine global\n");
    fprintf("         constraint, so only implicit reactions would remove it.\n");
end

fprintf("\n0c  DECISION NUMBER -- speedup ceiling if cohesion were fully implicit\n");
fprintf("      X = dt_CFL(advective v_b) / dt_CFL(actual)\n");
fprintf("      mean %.3f   min %.3f   max %.3f\n", b.XSum/n, b.XMin, b.XMax);
fprintf("      cohesive share of max|v_b|: mean %.3f   max %.3f\n", ...
    b.VbCohSum/n, b.VbCohMax);
if b.XSum/n < 1.05
    fprintf("      => X ~ 1: NO speed case for Bailo. Justify on bound preservation only.\n");
else
    fprintf("      => X = %.2fx: a real speed case exists.\n", b.XSum/n);
end

% --- 0e: do Solver A's phi_b and the diagnostic component sum agree? ---------
% Gates the component-elimination design (plan, "The payoff that makes
% boundedness structural"): the elimination is only legitimate if these two
% routes to phi_b are the same equation. Any mismatch would be dumped entirely
% into the reconstructed component.
n0 = f.GridZero;
phibCH  = r.SimulationData.PhibCH;          % Solver A, cells 1..n0
phibSum = phiB(1:n0, :);                    % diagnostic sum, same cells
ok = ~isnan(phibCH(1,:));                   % frames actually written
dAbs = abs(phibCH(:,ok) - phibSum(:,ok));
% Scale by the DOMAIN maximum, not pointwise: phi_b is ~0 over most of the
% supernatant early on, and a pointwise denominator manufactures huge relative
% errors out of nothing.
den  = max(max(abs(phibSum(:,ok)), [], "all"), 1e-12);
fprintf("\n0e  Solver A phi_b vs diagnostic component sum (cells 1..%d)\n", n0);
fprintf("      max abs diff  %.4e\n", max(dAbs, [], "all"));
fprintf("      max rel diff  %.4e   (scaled by max|phi_b| = %.4e)\n", max(dAbs, [], "all")/den, den);
fprintf("      final-frame max abs %.4e\n", max(dAbs(:,end)));
% A single-dt threshold is the WRONG test here: the two routes are coupled by a
% first-order operator splitting, so an O(dt) mismatch is expected at any finite
% dt and says nothing about whether they are the same equation. Measured
% 2026-08-23 at N=100, 3 d: MaxDt 3e-6 -> 2.7281e-05, MaxDt 7.5e-7 -> 6.8254e-06,
% a ratio of 3.997 against a dt ratio of 4. First order, so it converges away and
% the two routes ARE the same equation. Report the number; judge it by scaling.
fprintf("      judge by dt-scaling, not by this number alone: O(dt) splitting\n");
fprintf("      offset is expected. See slurm/cfl_0e_dtscale.sbatch.\n");

% --- 0f: ACCUMULATED drift of a free-running Solver A ------------------------
% 0e re-seeds u from Solver B every step, so it can only ever show a one-step
% difference. The elimination design makes Solver A's phi_b authoritative for
% the WHOLE run, so what matters is whether an un-reseeded CH state tracks the
% component sum over thousands of steps. That is this.
phibPar = r.SimulationData.PhibPar;
dPar = abs(phibPar(:,ok) - phibSum(:,ok));
tt = r.Frames.Time(ok); tt = tt(:).';
fprintf("\n0f  FREE-RUNNING Solver A phi_b vs component sum (never re-seeded)\n");
fprintf("      max abs drift  %.4e   (%.3f%% of max|phi_b|)\n", ...
    max(dPar, [], "all"), 100*max(dPar, [], "all")/den);
fprintf("      final-frame    %.4e   (%.3f%% of max|phi_b|)\n", ...
    max(dPar(:,end)), 100*max(dPar(:,end))/den);
fprintf("      one-step (0e) for comparison: %.4e -- ratio %.1fx\n", ...
    max(dAbs, [], "all"), max(dPar, [], "all")/max(max(dAbs, [], "all"), eps));
pd = r.SimulationData.PhibParDiag;
fprintf("      free state range over ALL steps: min %.4e  max %.4e\n", pd.minSeen, pd.maxSeen);
fprintf("      first phi_b < 0 at t = %.4g d\n", pd.tFirstNeg);
fprintf("      first phi_b > 1 at t = %.4g d\n", pd.tFirstAbove1);
fprintf("      first non-finite at t = %.4g d  (dead = %d)\n", pd.tFirstNaN, pd.dead);
if ~isnan(pd.tFirstNeg)
    fprintf("      => UNSTABLE. mobility zeta_0*u(1-u) < 0 once u < 0, so diffusion\n");
    fprintf("         becomes anti-diffusion and the blow-up self-reinforces. Solver A\n");
    fprintf("         CANNOT carry phi_b on its own with the present scheme.\n");
end
nshow = min(8, numel(tt));
idx = round(linspace(1, numel(tt), nshow));
fprintf("      drift growth:  %s\n", strjoin(compose("t=%.2f:%.1e", tt(idx).', max(dPar(:,idx)).'), "  "));

res = struct("budget", b, "phibSup", phibSup, ...
    "phibConsistencyAbs", max(dAbs, [], "all"), ...
    "phibConsistencyRel", max(dAbs, [], "all")/den, ...
    "phibParDriftAbs", max(dPar, [], "all"), ...
    "phibParDriftFinal", max(dPar(:,end)), ...
    "wsP50", b.WsP50Sum/n, "wsP90", b.WsP90Sum/n, "wsP99", b.WsP99Sum/n, ...
    "wsCellMax", b.WsCellMaxSum/n, "wsConc", concRatio, "wsBound", b.WsBoundSum/n, "phibParDiag", r.SimulationData.PhibParDiag, "flag", r.Flag, "wall", wall, ...
    "NCells", opts.NCells, "Kappa", opts.Kappa, "MatureDays", opts.MatureDays, ...
    "Xmean", b.XSum/n, "capFrac", b.CapBound/n, ...
    "regionFrac", b.RegionWins/n, "termFrac", b.TermSums/n, ...
    "regionShare", b.RegionShareSum/n, "regionShareMax", b.RegionShareMax, ...
    "matrixShareNoCoh", b.MatrixShareNoCoh/n);

if opts.Save
    tag = sprintf("cflbudget_n%d_k%g_%gd_dt%g%s", opts.NCells, opts.Kappa, opts.MatureDays, opts.MaxDt, ...
        ternaryStr(opts.PositiveMobility, "_posmob", ""));
    save(fullfile(S, tag + ".mat"), "res", "-v7");
    fprintf("\nwrote %s\n", fullfile(S, tag + ".mat"));
end
end

function out = ternaryStr(c, a, b)
if c, out = a; else, out = b; end
end

function phiB = biofilmFractionLocal(r)
% Mirror of manuscriptExperiments/biofilmFraction, kept local so this probe does
% not depend on that file's private functions.
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
