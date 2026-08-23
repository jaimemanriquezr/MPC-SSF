function res = validateBailo()
% VALIDATEBAILO  Phase 1 gates for the Bailo prototype (bailoCH1D).
%
% Gate 1  mass conserved to machine precision
% Gate 2  |phi| <= 1 at every cell and step
% Gate 3  discrete energy F_Delta monotone non-increasing, Eq. (2.2)/(2.3)
% Gate 4  second-order convergence to the EXACT steady state, Eq. (4.1)
% Gate 5  UNCONDITIONALITY -- gates 1-3 still hold as dt is raised far past the
%         explicit limit. This is the whole point; a prototype that only works
%         at small dt has proved nothing.
%
% See .claude/plans/2026-08-23-bailo-scheme.md, Phase 1.

here = fileparts(mfilename("fullpath")); addpath(here);
SV = "newton";  PM = 60;
% Newton, not Picard. Picard diverged at M >= 400 (2000-iteration cap) and made
% the Gate 4 and Gate 5 tables meaningless; Newton converges in 2-4 iterations at
% every mesh tried. PM = 60 is a generous cap that also detects failure.

fprintf("\n=== Gate 4: convergence to the exact steady state (eps = 0.1) ===\n");
Ms = [50 100 200 400 800];
errs = nan(size(Ms)); dxs = errs;
for k = 1:numel(Ms)
    o = bailoCH1D(Eps=0.1, M=Ms(k), Solver=SV, PicardMax=PM, Quiet=true);
    errs(k) = o.errSteady; dxs(k) = o.dx;
    if k == 1
        fprintf("  M=%4d dx=%.3e err=%.4e    -\n", Ms(k), dxs(k), errs(k));
    else
        p = log(errs(k-1)/errs(k))/log(dxs(k-1)/dxs(k));
        fprintf("  M=%4d dx=%.3e err=%.4e  order %.3f\n", Ms(k), dxs(k), errs(k), p);
    end
end
ordSteady = log(errs(1:end-1)./errs(2:end))./log(dxs(1:end-1)./dxs(2:end));

fprintf("\n=== Gates 1-3 across the paper's range of eps ===\n");
epss = [1 0.1 0.01];
tab = nan(numel(epss), 4);
for k = 1:numel(epss)
    o = bailoCH1D(Eps=epss(k), M=200, Solver=SV, PicardMax=PM, Quiet=true);
    tab(k,:) = [o.massDrift, o.boundViolation, o.maxEnergyIncrease, o.picardMean];
    fprintf("  eps=%-6.3g mass %.2e | bounds %+.2e | max dE %+.2e | Picard %.1f\n", ...
        epss(k), tab(k,1), tab(k,2), tab(k,3), tab(k,4));
end

fprintf("\n=== Gate 5: UNCONDITIONALITY (eps = 0.1, M = 200) ===\n");
base = 0.1*0.1^2;                 % the dt used above
mults = [1 10 100 1000];
unc = nan(numel(mults), 5);
for k = 1:numel(mults)
    o = bailoCH1D(Eps=0.1, M=200, Dt=base*mults(k), Solver=SV, PicardMax=PM, Quiet=true);
    unc(k,:) = [o.massDrift, o.boundViolation, o.maxEnergyIncrease, o.picardMean, o.errSteady];
    fprintf("  dt=%9.3g (%5gx) mass %.2e | bounds %+.2e | max dE %+.2e | Picard %5.1f | err %.3e\n", ...
        base*mults(k), mults(k), unc(k,1), unc(k,2), unc(k,3), unc(k,4), unc(k,5));
end

fprintf("\n--- verdicts ---\n");
fprintf("  Gate 1 mass        : %d  (all drifts < 1e-9)\n", all(tab(:,1) < 1e-9) && all(unc(:,1) < 1e-9));
fprintf("  Gate 2 bounds      : %d  (all violations < 1e-8, i.e. at solve roundoff)\n", ...
    all(tab(:,2) < 1e-8) && all(unc(:,2) < 1e-8));
fprintf("  Gate 3 energy      : %d  (no increase beyond 1e-10)\n", ...
    all(tab(:,3) < 1e-10) && all(unc(:,3) < 1e-10));
fprintf("  Gate 4 second order: %d  (asymptotic %.3f)\n", abs(ordSteady(end) - 2) < 0.2, ordSteady(end));
fprintf("  Gate 5 uncondition : %d  (bounds+energy hold at 1000x dt, Newton CONVERGED)\n", ...
    unc(end,2) < 1e-8 && unc(end,3) < 1e-10 && unc(end,4) < PM);
fprintf("  (Newton iteration counts above are the cost: compare against the ONE\n");
fprintf("   linear solve the present Solver A does per step.)\n");

res = struct("Ms", Ms, "dx", dxs, "errSteady", errs, "orderSteady", ordSteady, ...
    "eps", epss, "epsTable", tab, "dtMults", mults, "uncond", unc);
save(fullfile(here, "validateBailo.mat"), "res", "-v7");
end
