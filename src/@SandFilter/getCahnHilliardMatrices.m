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

    dz = mean(diff(depthBoundaries));
    n0 = find(abs(depthCenters) < dz/2);

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
    kappa = model.CohesionSubModel.Kappa;
    diffusion.Rows = min(max(1, localRows(:) + (0:n0) - 1), n0);
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