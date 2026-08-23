function [convection, diffusion, diffusionMobility] = getCahnHilliardMatrices(filter, model, upwindConvection)
    arguments
        filter  SandFilter
        model  Model
        upwindConvection = true
    end
    depthCenters = filter.GridPoints.Centers;
    depthBoundaries = filter.GridPoints.Boundaries;

    porosityCenters = computePorosity(filter, depthCenters);
    porosityBoundaries = computePorosity(filter, depthBoundaries);

    % Use the filter's own accessors rather than re-deriving them: these were
    % duplicated here, so a fix to either had to be made twice.
    dz = filter.GridSize;
    n0 = filter.GridZero;

    [convection, diffusion, diffusionMobility] = deal(struct('Rows', [], 'Columns', [], 'Values', []));
    if upwindConvection
        convectionStencil = [-1; 1; 0; 0];    %% UPWIND!!!
    else
        convectionStencil = [-1; 1; -1; 1]/2; %% CENTERED!!!
    end
    diffusionStencil = [1; -1; -1; 1];

    [localColumns, localRows] = meshgrid(1:2);
    % CONVECTION MATRIX
    q = filter.InflowVelocity;
    convection.Rows = localRows(:) + (0:n0-2);
    convection.Columns = localColumns(:) + (0:n0-2);
    convectionValues = q * (convectionStencil .* ones(1, n0-1)) ./ porosityCenters(convection.Rows);
    convection.Values = convectionValues / dz;

    % DIFFUSION MATRIX
    % This is the kappa*L operator of the CHEMICAL POTENTIAL equation
    %   mu_j = -kappa*(Du_{j+1/2} - Du_{j-1/2})/dz + Psi'(u_j)          (Eq. 29)
    % i.e. mu = -kappa*Delta_h u + Psi'(u) = kappa*L*u + Psi'(u), L = -Delta_h.
    % So it belongs in the mu row -- block (2,1) of the mixed system -- giving
    % mu - kappa*L*u = Psi'(u). Hence the "+ n0" on Rows, mirroring the "+ n0"
    % that diffusionMobility applies to its Columns to reach block (1,2).
    %
    % It was previously assembled with Rows in 1..n0, i.e. block (1,1). Because
    % every operator here has Rows <= n0, that left block (2,1) EMPTY and block
    % (2,2) a bare identity, so the solve reduced to mu = Psi'(u) exactly and
    % KAPPA HAD NO EFFECT ON ANY OUTPUT. Confirmed at operator level: with random
    % u and random mobility, max|muCH - Psi'(u)| was 0.000e+00 for any state, and
    % kappa ranked 27/27 with exactly zero sensitivity in the log-OAT campaign.
    % See .claude/decisions/2026-08-23-cahn-hilliard-kappa-inert.md.
    %
    % The index clamping to [1, n0] stays: it imposes the homogeneous Neumann
    % condition d(phi_b)/dz = 0 at z = 0 that manuscript Appendix A.1 adopted in
    % place of the Diehl2025 Dirichlet condition. numerical-method.md:250-262
    % transcribes the superseded Dirichlet variant, including a -kappa*u_1/dz^2
    % term in the mu right-hand side, which must NOT be reinstated here.
    kappa = model.CohesionSubModel.Kappa;
    diffusion.Rows = min(max(1, localRows(:) + (0:n0) - 1), n0) + n0;
    diffusion.Columns = min(max(1, localColumns(:) + (0:n0) - 1), n0);
    diffusionValues = kappa * (diffusionStencil.*ones(1, n0+1));
    diffusion.Values = diffusionValues / dz^2;

    % WEIGHTED DIFFUSION MATRIX
    lambda = porosityBoundaries(2:n0).';
    diffusionMobility.Rows = localRows(:) + (0:n0-2);
    diffusionMobility.Columns = localColumns(:) + (0:n0-2) + n0;
    diffusionMobilityValues = -(diffusionStencil .* lambda) ./ porosityCenters(diffusionMobility.Rows);
    diffusionMobility.Values = diffusionMobilityValues / dz^2;
end