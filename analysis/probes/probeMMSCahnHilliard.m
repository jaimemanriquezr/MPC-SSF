function res = probeMMSCahnHilliard(opts)
% PROBEMMSCAHNHILLIARD  Manufactured-solution verification of the kappa operator.
%
%   probeMMSCahnHilliard
%   probeMMSCahnHilliard(NList=[20 40 80 160 320 640], Kappa=1e-7)
%
% WHAT IS BEING VERIFIED ----------------------------------------------------
% The chemical-potential row of the Cahn-Hilliard mixed system (Eq. 29):
%
%     mu_j = -kappa * d2u/dz2 |_j  +  Psi'(u_j)
%
% assembled as block (2,1) of the 2*n0 system. Until 2026-08-23 that block was
% EMPTY -- the kappa operator sat in block (1,1), whose solution is discarded, so
% mu collapsed to Psi'(u) and kappa had no effect on any output
% (.claude/decisions/2026-08-23-cahn-hilliard-kappa-inert.md). This probe checks
% that the restored operator converges at the design order.
%
% METHOD --------------------------------------------------------------------
% Take a smooth manufactured u(z) whose exact second derivative is known, apply
% the ASSEMBLED block-(2,1) operator to it, and compare against -kappa*u''.
% Refining the mesh must give second-order convergence.
%
% The boundary condition matters. getCahnHilliardMatrices clamps the diffusion
% indices to [1, n0], which imposes homogeneous NEUMANN (zero-flux) at both faces
% of the CH domain -- the condition manuscript Appendix A.1 adopted in place of
% the Diehl2025 Dirichlet one. The manufactured solution must therefore satisfy
% u'= 0 at both faces, else the boundary rows carry an O(1) consistency error
% that has nothing to do with the operator being tested. Hence
%
%     u(z) = a + b*cos(m*pi*zeta),   zeta = (z - zLo)/(zHi - zLo)
%
% on the CH domain [zLo, zHi] = [boundaries(1), boundaries(n0+1)], which has zero
% derivative at both ends for any integer m.
%
% Two error norms are reported: over the interior only, and over the whole
% domain including the two boundary rows. If the clamping implements Neumann
% correctly both are second order; if it does not, the boundary norm degrades to
% first order or worse. That contrast is the point.
%
% See also GETCAHNHILLIARDMATRICES, PROBELIMITDIAG.

arguments
    opts.NList (1,:) double = [20 40 80 160 320 640]
    opts.Kappa (1,1) double = 1e-7
    opts.Amplitude (1,1) double = 0.30     % b; mean level a below
    opts.Mean (1,1) double = 0.45          % a: keeps u well inside (0,1)
    opts.Modes (1,1) double = 3            % m
    opts.Save (1,1) logical = true
end

here = fileparts(mfilename("fullpath"));
W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

kappa = opts.Kappa;  a = opts.Mean;  b = opts.Amplitude;  m = opts.Modes;

fprintf("\n=== MMS: Cahn-Hilliard kappa operator ===\n");
fprintf("kappa = %g, u(z) = %g + %g*cos(%d*pi*zeta), Neumann at both faces\n\n", ...
    kappa, a, b, m);
fprintf("%-7s %-12s %-13s %-8s %-13s %-8s\n", ...
    "N", "dz", "err interior", "order", "err all", "order");

nN = numel(opts.NList);
[dzs, eInt, eAll] = deal(nan(1, nN));

for k = 1:nN
    N = opts.NList(k);
    f = SandFilter();
    f = f.addGridPoints(N);
    n0 = f.GridZero;
    mdl = Model(Component.empty, Reaction.empty, ...
        Kappa=kappa, Zeta0=1e2, Zeta1=1e-2);

    [~, DD, ~] = getCahnHilliardMatrices(f, mdl, false);
    A = sparse(DD.Rows, DD.Columns, DD.Values, 2*n0, 2*n0);
    B21 = A(n0+1:end, 1:n0);          % the operator under test: kappa*L = -kappa*Delta_h

    zc  = f.GridPoints.Centers(1:n0);
    zLo = f.GridPoints.Boundaries(1);
    zHi = f.GridPoints.Boundaries(n0+1);
    Ldom = zHi - zLo;
    zeta = (zc - zLo)/Ldom;

    u      = a + b*cos(m*pi*zeta);
    d2u    = -b*(m*pi/Ldom)^2 * cos(m*pi*zeta);
    exact  = -kappa*d2u;               % the kappa contribution to mu
    approx = B21*u;                    % what the assembled operator produces

    err  = abs(approx - exact);
    dzs(k)  = f.GridSize;
    eInt(k) = max(err(2:end-1));
    eAll(k) = max(err);

    if k == 1
        fprintf("%-7d %-12.4e %-13.4e %-8s %-13.4e %-8s\n", N, dzs(k), eInt(k), "-", eAll(k), "-");
    else
        oI = log(eInt(k-1)/eInt(k))/log(dzs(k-1)/dzs(k));
        oA = log(eAll(k-1)/eAll(k))/log(dzs(k-1)/dzs(k));
        fprintf("%-7d %-12.4e %-13.4e %-8.3f %-13.4e %-8.3f\n", N, dzs(k), eInt(k), oI, eAll(k), oA);
    end
end

ordInt = log(eInt(1:end-1)./eInt(2:end))./log(dzs(1:end-1)./dzs(2:end));
ordAll = log(eAll(1:end-1)./eAll(2:end))./log(dzs(1:end-1)./dzs(2:end));

fprintf("\nasymptotic order (last pair): interior %.3f, all %.3f\n", ordInt(end), ordAll(end));
fprintf("VERDICT interior second order : %d\n", abs(ordInt(end) - 2) < 0.15);
fprintf("VERDICT Neumann rows preserved: %d\n", abs(ordAll(end) - 2) < 0.15);

% Linearity in kappa: the operator is exactly linear, so doubling kappa must
% double the error. A cheap independent check that nothing else crept in.
f = SandFilter(); f = f.addGridPoints(160); n0 = f.GridZero;
e2 = nan(1,2); ks = [kappa, 2*kappa];
for i = 1:2
    mdl = Model(Component.empty, Reaction.empty, Kappa=ks(i), Zeta0=1e2, Zeta1=1e-2);
    [~, DD, ~] = getCahnHilliardMatrices(f, mdl, false);
    A = sparse(DD.Rows, DD.Columns, DD.Values, 2*n0, 2*n0);
    zc = f.GridPoints.Centers(1:n0);
    zLo = f.GridPoints.Boundaries(1); zHi = f.GridPoints.Boundaries(n0+1);
    Ldom = zHi - zLo; zeta = (zc - zLo)/Ldom;
    u = a + b*cos(m*pi*zeta);
    d2u = -b*(m*pi/Ldom)^2 * cos(m*pi*zeta);
    e2(i) = max(abs(A(n0+1:end,1:n0)*u - (-ks(i)*d2u)));
end
fprintf("kappa doubling -> error ratio %.6f (exact 2.0)\n", e2(2)/e2(1));

res = struct("NList", opts.NList, "dz", dzs, "errInterior", eInt, "errAll", eAll, ...
    "orderInterior", ordInt, "orderAll", ordAll, "kappa", kappa, ...
    "kappaLinearity", e2(2)/e2(1));

if opts.Save
    save(fullfile(S, "mms_cahnhilliard.mat"), "res", "-v7");
    writematrix([opts.NList(:), dzs(:), eInt(:), eAll(:)], ...
        fullfile(S, "mms_cahnhilliard.csv"));
    fprintf("\nwrote %s/mms_cahnhilliard.{mat,csv}\n", S);
end
end
