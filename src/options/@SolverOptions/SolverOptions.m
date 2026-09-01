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

        CohesionScheme string = "shin"
        CohesionBoundaryConditions string = "neumann"
        % Flowing<->enclosed liquid exchange. "constant": b = Liquid.TransportRate
        % (preset 600/d). "film": b = D/L_f^2 with L_f the coating thickness
        % eps*phi_b/a_s on grains of diameter GrainDiameter (a_s = 6(1-eps)/d), and
        % the mat thickness (capped at MatThickness) where eps = 1; advanced with
        % the exact relaxation over the step. See decisions/2026-08-26-film-transfer.md.
        TransferForm string = "constant"
        TransferDiffusivity double = 1.0e-4     % m^2/d (Wolf2007: O2 1.73e-4, CO2 1.65e-4, cations 1.15e-4)
        GrainDiameter double = 3.0e-4           % m (Campos2002 Table 1, d10)
        MatThickness double = 1.0e-3            % m (Wolf2007 Table VI, L_f,max)
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
                obj SolverOptions
                options struct
            end
            % `aliases` maps the LEGACY simulate() names onto the canonical
            % property names. 39 files under analysis/ still pass the legacy
            % names; once they are migrated, delete this map and the lookup.
            aliases = dictionary( ...
                ["AdaptiveInitialDt",       "AdaptiveMaxDt", ...
                 "IsUpwinded",              "CohesionBC"], ...
                ["AdaptiveInitialTimeStep", "AdaptiveMaxTimeStep", ...
                 "UpwindedVelocity",        "CohesionBoundaryConditions"]);
            for field = string(fieldnames(options)).'
                name = field;
                if isKey(aliases, field), name = aliases(field); end
                % NOTE: options.(field), not options.field -- the latter looks up
                % a literal field called "field" and throws.
                if isprop(obj, name) && ~isempty(options.(field))
                    obj.(name) = options.(field);
                end
            end
        end
    end
end
