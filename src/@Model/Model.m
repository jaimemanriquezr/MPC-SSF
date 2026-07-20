classdef Model
    properties
        Components (:, 1)

        WaterDensity double = 998
        BiofilmPorosity double
        OsmosisRate double
        DetachmentFunction function_handle

        Reactions =   Reaction.empty;
        CohesionSubModel =  CahnHilliardModel()
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
                components = Component.empty;
                reactions = Reaction.empty;
                input.Components =  Component.empty;
                input.Kappa double = [];
                input.Zeta0 double = [];
                input.Zeta1 double = [];
                input.MobilityFunction function_handle = @(u) u .* (1 - u);
                input.PotentialGradient function_handle = @(u) 0.25*((u.^2).*((1 - u).^2));

                input.WaterDensity double = 998;
                input.DetachmentFunction function_handle = @(v, qnom) sqrt(abs(v)/qnom);
                input.OsmosisRate double = 1e-5;
                input.BiofilmPorosity double = .99;

                input.Reactions = Reaction.empty;
                input.Preset = string.empty;
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
            switch input.Preset
                case string.empty
                    propNames = string(properties(obj));
                    for field = propNames(~contains(propNames, "SubModel")).'
                        obj.(field) = input.(field);
                    end
                    for superField = propNames(contains(propNames, "SubModel")).'
                        for field = string(fieldnames(obj.(superField)).')
                            inputValue = input.(field);
                            if ~isempty(inputValue)
                                obj.(superField).(field) = inputValue;
                            end
                        end
                    end
                case {"Lund", "Rosenqvist"}
                    obj =  presets.modelLund();
                otherwise
                    error("Unknown model preset.");
            end
        end

        function p = get.Particles(obj)
            p = obj.Components(arrayfun(@(C)isa(C, 'Particle'), obj.Components));
        end

        function l = get.Liquids(obj)
            l = obj.Components(arrayfun(@(C)isa(C, 'Liquid'), obj.Components));
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
