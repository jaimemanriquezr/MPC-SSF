classdef Model
    properties
        Components (:, 1)

        WaterDensity double = 998
        BiofilmPorosity double
        OsmosisRate double
        DetachmentFunction function_handle

        Reactions =   mpcssf.Reaction.empty;
        CohesionSubModel =  mpcssf.CahnHilliardModel()
    end

    properties (Dependent, Hidden)
        Particles
        Liquids

        StoichiometricCoefficients
        HalfSaturationConstants
        Order

        StoichiometricMatrixParticles
        StoichiometricMatrixLiquids
    end

    methods
        function obj = Model(components, reactions, input)
            arguments
                components = mpcssf.Component.empty;
                reactions = mpcssf.Reaction.empty;
                input.Components =  mpcssf.Component.empty;
                input.Kappa double = [];
                input.Zeta0 double = [];
                input.Zeta1 double = [];
                input.MobilityFunction function_handle = @(u) u .* (1 - u);
                % No default: its absence signals that CahnHilliardModel should
                % build the published Zeta1-dependent handle. See below.
                input.PotentialGradient function_handle

                input.WaterDensity double = 998;
                input.DetachmentFunction function_handle = @(v, qnom) sqrt(abs(v)/qnom);
                input.OsmosisRate double = 1e-5;
                input.BiofilmPorosity double = .99;

                input.Reactions = mpcssf.Reaction.empty;
            end
            % Positional args override the name-value defaults only when
            % actually supplied; otherwise a name-value call
            % (Model(Components=..., Reactions=...), as modelLund uses) would be
            % clobbered by the empty positional defaults.
            if ~isempty(components)
                input.Components = components;
            end
            if ~isempty(reactions)
                input.Reactions = reactions;
            end
            % Preset was removed 2026-09-22: its only non-empty branches called
            % presets.modelLund(), and no +presets package exists, so they could
            % only ever have errored. Zero callers repo-wide. What follows is the
            % former `case string.empty` body, which is the path every caller took.
                propNames = string(properties(obj));
                for field = propNames(~contains(propNames, "SubModel")).'
                    obj.(field) = input.(field);
                end
                for superField = propNames(contains(propNames, "SubModel")).'
                    for field = string(fieldnames(obj.(superField)).')
                        if ~isfield(input, field)
                            continue    % optional arg genuinely not supplied
                        end
                        inputValue = input.(field);
                        if ~isempty(inputValue)
                            obj.(superField).(field) = inputValue;
                        end
                    end
                end
                % PotentialGradient closes over Zeta1, and the loop above
                % assigns Zeta1 field-by-field AFTER construction -- which
                % leaves the handle built from the old Zeta1. Rebuild through
                % the constructor so the published form is the single source
                % of truth. Skipped when the caller supplied a handle.
                if ~isfield(input, "PotentialGradient")
                    c = obj.CohesionSubModel;
                    obj.CohesionSubModel = mpcssf.CahnHilliardModel(Kappa=c.Kappa, ...
                        Zeta0=c.Zeta0, Zeta1=c.Zeta1, ...
                        MobilityFunction=c.MobilityFunction);
                end
        end

        function p = get.Particles(obj)
            p = obj.Components(arrayfun(@(C)isa(C, 'mpcssf.Particle'), obj.Components));
        end

        function l = get.Liquids(obj)
            l = obj.Components(arrayfun(@(C)isa(C, 'mpcssf.Liquid'), obj.Components));
        end

        function sigma = get.StoichiometricCoefficients(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Components);
        end

        function sigma = get.StoichiometricMatrixParticles(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Particles);
        end

        function sigma = get.StoichiometricMatrixLiquids(obj)
            sigma = obj.Reactions.lookupStoichiometricCoefficients(obj.Liquids);
        end

        function K = get.HalfSaturationConstants(obj)
            K = obj.Reactions.lookupHalfSaturationConstants(obj.Components);
        end

        function p = get.Order(obj)
            p = obj.Reactions.lookupOrder(obj.Components);
        end

        mu = computeReactionRates(obj, temperature);
    end
end
