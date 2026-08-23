function out = bailoCH1D(opts)
% BAILOCH1D  Bailo et al. (2023) bound-preserving finite-volume Cahn-Hilliard, 1-D.
%
% Standalone prototype -- NOT wired into the solver. Implements Eq. (2.1a)-(2.1i)
% of documents/Bailo2023.pdf (arXiv:2105.05351v2).
%
%   phi_t + (F_{i+1/2} - F_{i-1/2})/dx = 0                                (2.1a)
%   F_{i+1/2} = M(phi_i, phi_{i+1})(u_{i+1/2})^+ + M(phi_{i+1}, phi_i)(u_{i+1/2})^-
%   u_{i+1/2} = -(xi_{i+1} - xi_i)/dx                                     (2.1c)
%   M(x,y)    = M0 (1+x)^+ (1-y)^+                                        (2.1d)
%   xi_i      = Hc'(phi_i^{n+1}) - He'(phi_i^n) - eps2*(Lap phi)_i^{n+1}  (2.1e)
%   F_{1/2} = F_{M+1/2} = 0,  ghosts phi_0 = phi_1, phi_{M+1} = phi_M     (2.1h,i)
%
% Everything at n+1 except He', which is the explicit (expansive) half of the
% convex splitting. The wetting term W is zero here: we have no contact-angle
% physics, and our z = 0 boundary is homogeneous Neumann.
%
% WHY BOUND PRESERVATION HOLDS, and what it rests on: the contradiction proof
% (paper p. 6) sums (2.1a) over a contiguous run of cells with phi > 1 and kills
% the two dangerous flux terms using only (1 - y)^+ = 0. It uses NO property of
% the potential. That is what makes it applicable to our very non-standard Psi.
%
% NONLINEAR SOLVE. The paper does not specify one. Used here: Picard with the
% mobility and the upwind direction lagged, so each iterate is a linear
% pentadiagonal solve (xi couples 3 cells, the flux difference couples 5).
% Lagging is safe at convergence -- the converged iterate satisfies the fully
% implicit equations, which is what the proofs require. Iteration counts are
% reported because they are the cost that must be weighed against the present
% single linear solve.
%
% VARIABLE CONVENTION. Bailo works on phi in [-1,1]. This file stays in that
% convention deliberately: the validation case of section 4.1 and its exact
% steady state are stated there, so reproducing them literally is an INDEPENDENT
% check of the implementation. The SSF model lives on (0,1); the map is
% phi_B = 2u - 1 with M0 = zeta_0, eps2 = kappa/4, H(phi_B) = Psi((1+phi_B)/2),
% and natively M(x,y) = zeta_0 (x)^+ (1-y)^+. See
% .claude/plans/2026-08-23-bailo-scheme.md.
%
% Usage:
%   out = bailoCH1D(Eps=0.1, M=200)                 % deep-quench validation case
%   out = bailoCH1D(..., Potential="ssf", Zeta1=1e-2)

arguments
    opts.Eps (1,1) double = 0.1          % interface parameter (eps, not eps^2)
    opts.M (1,1) double = 200            % cells
    opts.M0 (1,1) double = 1.0           % mobility scale
    opts.Dt (1,1) double = NaN           % NaN -> 0.1*eps^2 ((1,1) rejects [])
    opts.TFinal (1,1) double = NaN       % NaN -> 20*eps^2 (paper's choice)
    opts.Potential (1,1) string = "deepquench"   % or "ssf"
    opts.Zeta1 (1,1) double = 1e-2       % only for Potential="ssf"
    opts.PicardTol (1,1) double = 1e-12
    opts.PicardMax (1,1) double = 100
    opts.Quiet (1,1) logical = false
end

eps2 = opts.Eps^2;
if isnan(opts.TFinal), TFinal = 20*eps2; else, TFinal = opts.TFinal; end
if isnan(opts.Dt),     dt = 0.1*eps2;      else, dt = opts.Dt;         end

% --- domain and datum, section 4.1 -------------------------------------------
L = 1.5*pi*opts.Eps;                       % Omega = [-3*pi*eps/2, 3*pi*eps/2]
M = opts.M;  dx = 2*L/M;
x = -L + ((1:M)' - 0.5)*dx;

switch opts.Potential
    case "deepquench"
        % H(phi) = (1-phi^2)/2 is CONCAVE, so the splitting is Hc = 0,
        % He = phi^2/2 - 1/2, both convex. Hence Hc' = 0 and He' = phi.
        Hc_p = @(p) zeros(size(p));
        He_p = @(p) p;
        Hfun = @(p) (1 - p.^2)/2;
        phi0 = -ones(M,1);
        inner = abs(x) <= pi*opts.Eps/2;
        phi0(inner) = cos(x(inner)/opts.Eps) - 1;
    case "ssf"
        % Our potential mapped into Bailo's variable: H(phi_B) = Psi((1+phi_B)/2)
        % with Psi(u) = u^4/4 - 0.5*zeta_1*u^3. Convexity deficit of Psi is only
        % -3*zeta_1^2/4 = -7.5e-05, so the splitting constant is tiny; see the
        % plan. Split in u: Psi_c = Psi + a*u^2/2, Psi_e = a*u^2/2.
        z1 = opts.Zeta1;  a = 3*z1^2/4;
        u_of = @(p) (1 + p)/2;
        % dH/dphi_B = Psi'(u)/2
        Hc_p = @(p) ( u_of(p).^2.*(u_of(p) - 1.5*z1) + a*u_of(p) )/2;
        He_p = @(p) ( a*u_of(p) )/2;
        Hfun = @(p) u_of(p).^4/4 - 0.5*z1*u_of(p).^3;
        phi0 = -ones(M,1);
        inner = abs(x) <= pi*opts.Eps/2;
        phi0(inner) = cos(x(inner)/opts.Eps) - 1;
    otherwise
        error("bailoCH1D:potential", "unknown Potential '%s'", opts.Potential);
end

% --- discrete Laplacian with the ghost values of (2.1i) ----------------------
% phi_0 := phi_1 and phi_{M+1} := phi_M, i.e. homogeneous Neumann.
e = ones(M,1);
Lap = spdiags([e, -2*e, e], -1:1, M, M);
Lap(1,1) = -1; Lap(M,M) = -1;              % ghost rows
Lap = Lap/dx^2;

Mfun = @(xx, yy) opts.M0 * max(1 + xx, 0) .* max(1 - yy, 0);

phi = phi0;
nSteps = max(1, round(TFinal/dt));
energy = nan(nSteps+1, 1);  mass = energy;  minphi = energy;  maxphi = energy;
picard = nan(nSteps, 1);
energy(1) = discreteEnergy(phi, Hfun, eps2, dx);
mass(1) = sum(phi)*dx;  minphi(1) = min(phi);  maxphi(1) = max(phi);

for n = 1:nSteps
    phin = phi;
    Hep = He_p(phin);                       % explicit half, frozen this step
    w = phin;                                % Picard iterate
    converged = false;
    for k = 1:opts.PicardMax
        % lagged mobility + lagged upwind direction -> linear system in wNew
        xi_w = Hc_p(w) - Hep - eps2*(Lap*w);
        uface = -(xi_w(2:end) - xi_w(1:end-1))/dx;         % M-1 interior faces
        Mp = Mfun(w(1:end-1), w(2:end));                    % for u >= 0
        Mm = Mfun(w(2:end),   w(1:end-1));                  % for u <  0
        Aface = Mp.*(uface >= 0) + Mm.*(uface < 0);         % chosen coefficient

        % F = Aface * uface, uface = -(xi_{i+1} - xi_i)/dx, xi linear in wNew
        % (Hc' lagged for the "ssf" potential; identically zero for deepquench)
        % build the operator: div(F) acting on the unknown v, with the
        % eps2*Lap part of xi kept implicit and the Hc'/He' part moved to rhs
        Dface = spdiags([-e(1:M-1), e(1:M-1)], [0 1], M-1, M)/dx;   % face gradient
        Uop = -Dface*(-eps2*Lap);                 % uface as a linear map on v
        Fop = spdiags(Aface, 0, M-1, M-1)*Uop;
        Div = spdiags([-e, e], [-1 0], M, M-1)/dx;  % no-flux ends, (2.1h)
        Aop = speye(M)/dt + Div*Fop;
        % constant part of xi (Hc'(w) - He'(phi^n)) also drives the flux
        xi_const = Hc_p(w) - Hep;
        u_const = -(xi_const(2:end) - xi_const(1:end-1))/dx;
        rhs = phin/dt - Div*(Aface.*u_const);
        wNew = Aop \ rhs;

        dw = norm(wNew - w, inf);
        w = wNew;
        if dw < opts.PicardTol, converged = true; picard(n) = k; break, end
    end
    if ~converged, picard(n) = opts.PicardMax; end
    phi = w;

    energy(n+1) = discreteEnergy(phi, Hfun, eps2, dx);
    mass(n+1) = sum(phi)*dx;  minphi(n+1) = min(phi);  maxphi(n+1) = max(phi);
end

% --- exact steady state, Eq. (4.1) -- deep quench only -----------------------
phiInf = [];
if opts.Potential == "deepquench"
    phiInf = -ones(M,1);
    ii = abs(x) <= pi*opts.Eps;
    phiInf(ii) = (1/pi)*(1 + cos(x(ii)/opts.Eps)) - 1;
end

out = struct("x", x, "dx", dx, "phi", phi, "phi0", phi0, "phiInf", phiInf, ...
    "energy", energy, "mass", mass, "minphi", minphi, "maxphi", maxphi, ...
    "picard", picard, "dt", dt, "nSteps", nSteps, "eps", opts.Eps, ...
    "massDrift", abs(mass(end) - mass(1))/max(abs(mass(1)), eps), ...
    "boundViolation", max(max(maxphi) - 1, -min(minphi) - 1), ...
    "energyIncreases", sum(diff(energy) > 0), ...
    "maxEnergyIncrease", max([0; diff(energy)]), ...
    "picardMean", mean(picard), "picardMax", max(picard));
if ~isempty(phiInf), out.errSteady = max(abs(phi - phiInf)); end

if ~opts.Quiet
    fprintf("eps=%.3g M=%4d dt=%.3g steps=%d | mass drift %.2e | bounds viol %.2e | " + ...
        "energy up %d/%d (max %.2e) | Picard %.1f (max %d)", ...
        opts.Eps, M, dt, nSteps, out.massDrift, out.boundViolation, ...
        out.energyIncreases, nSteps, out.maxEnergyIncrease, out.picardMean, out.picardMax);
    if ~isempty(phiInf), fprintf(" | err vs phi_inf %.3e", out.errSteady); end
    fprintf("\n");
end
end

function E = discreteEnergy(phi, Hfun, eps2, dx)
% F_Delta, Eq. (2.2), without the wetting terms.
grad = (phi(2:end) - phi(1:end-1))/dx;
E = sum(Hfun(phi))*dx + sum((eps2/2)*grad.^2)*dx;
end
