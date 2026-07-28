classdef CahnHilliardModel
    properties
        Kappa double = [];
        Zeta0 double = [];
        Zeta1 double = [];
        MobilityFunction function_handle = @(u) u .* (1 - u);
        % Published chemical-potential gradient: psi(u) = u^3(u/2 - Zeta1)/2, so
        % dpsi/du = u^2(u - 3*Zeta1/2). The legacy .mat models store exactly this
        % as a baked handle (SDparameters*.mat with Zeta1 = 1/200, thesis_model.mat
        % with 1/100). Rebuilt in the constructor once Zeta1 is known; this
        % declaration is the Zeta1 = 0 degenerate case.
        PotentialGradient function_handle = @(u) u.^2 .* u;
    end
    methods
        function obj = CahnHilliardModel(input)
            arguments
                input.Kappa double = [];
                input.Zeta0 double = [];
                input.Zeta1 double = [];
                input.MobilityFunction function_handle = @(u) u .* (1 - u);
                % No default: its absence is how we detect that the caller did not
                % supply a handle and the published form should be built below.
                input.PotentialGradient function_handle
            end
            for field = string(fieldnames(input).')
                obj.(field) = input.(field);
            end
            % MATLAB cannot reference a sibling argument in an `arguments` default,
            % so the Zeta1-dependent handle is built here instead.
            if ~isfield(input, "PotentialGradient")
                zeta1 = obj.Zeta1;
                if isempty(zeta1)
                    zeta1 = 0;
                end
                obj.PotentialGradient = @(u) u.^2 .* (u - 1.5*zeta1);
            end
        end
    end
end