classdef CahnHilliardModel
    properties
        Kappa double = [];
        Zeta0 double = [];
        Zeta1 double = [];
        MobilityFunction function_handle = @(u) u .* (1 - u);
        PotentialGradient function_handle = @(u) 0.25*((u.^2).*((1 - u).^2));
    end
    methods
        function obj = CahnHilliardModel(input)
            arguments
                input.Kappa double = [];
                input.Zeta0 double = [];
                input.Zeta1 double = [];
                input.MobilityFunction function_handle = @(u) u .* (1 - u);
                input.PotentialGradient function_handle = @(u) 0.25*((u.^2).*((1 - u).^2));
            end
            for field = string(fieldnames(input).')
                obj.(field) = input.(field);
            end
        end
    end
end