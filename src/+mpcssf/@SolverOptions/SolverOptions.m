classdef SolverOptions
    properties
        CFLFactor double = .99
        AdaptiveVelocityFactor double = 0.0
        AdaptiveTimeTolerance double = 1e-3
        AdaptiveInitialTimeStep double = 5e-5;
        AdaptiveMaxTimeStep double = 5e-5;

        UpwindedVelocity logical = true
        ImplicitOsmosis logical = true
        ImplicitDispersion logical = true
        ConvexSplitting logical = false

    end

    properties ( Hidden)
        % Flowing<->enclosed liquid exchange.
        % "constant": b = Liquid.TransportRate
        % "film": b = D/L_f^2 with L_f the coating thickness eps*phi_b / (6(1-eps)/d)
        TransferForm string = "constant"
        TransferDiffusivity double = 1.0e-4     % m^2/d (Wolf2007: O2 1.73e-4, CO2 1.65e-4, cations 1.15e-4)
        GrainDiameter double = 3.0e-4           % m (Campos2002 Table 1, d10)
        MatThickness double = 1.0e-3            % m (Wolf2007 Table VI, L_f,max)
        CohesionScheme string = "shin"
        CohesionBoundaryConditions string = "neumann"
    end

    methods
        function obj = SolverOptions(parameters)
            arguments
                parameters.CFLFactor (1,1) {mustBeNumeric} = .99;
                parameters.AdaptiveVelocityFactor (1,1) {mustBeNumeric} = 0;
                parameters.AdaptiveTimeTolerance (1,1) {mustBeNumeric} = 1e-3;
                parameters.AdaptiveInitialTimeStep (1,1) {mustBeNumeric} = 5e-5;
                parameters.AdaptiveMaxTimeStep (1,1) {mustBeNumeric} = 5e-5;

                parameters.UpwindedVelocity (1,1) = true;
                parameters.ImplicitOsmosis (1,1) = true;
                parameters.ImplicitDispersion (1,1) logical = true;
                parameters.ConvexSplitting (1,1) logical = false;

                parameters.CohesionScheme string {mustBeMember(parameters.CohesionScheme, ["shin","matched","bailo"])} = "shin";
                parameters.CohesionBoundaryConditions string {mustBeMember(parameters.CohesionBoundaryConditions, ["neumann","dirichlet"])} = "neumann";
            end

            for field = string(fieldnames(parameters)).'
                obj.(field) = parameters.(field);
            end
        end

        function s = summary(obj)
            % One-line provenance string for logs and figure captions.
            s = sprintf("scheme=%s bc=%s upwind=%d impOsm=%d impDisp=%d cvx=%d " + ...
                "maxDt=%g initDt=%g cfl=%g", obj.CohesionScheme, ...
                obj.CohesionBoundaryConditions, obj.UpwindedVelocity, ...
                obj.ImplicitOsmosis, obj.ImplicitDispersion, obj.ConvexSplitting, ...
                obj.AdaptiveMaxTimeStep, obj.AdaptiveInitialTimeStep, obj.CFLFactor);
        end

        function obj = override(obj, options)
            arguments
                obj mpcssf.SolverOptions
                options struct
            end
            % Copies every non-empty field of `options` that names a property.
            for field = string(fieldnames(options)).'
                if isprop(obj, field) && ~isempty(options.(field))
                    obj.(field) = options.(field);
                end
            end
        end
    end
end
