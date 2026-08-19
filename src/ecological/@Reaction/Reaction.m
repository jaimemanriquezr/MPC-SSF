classdef Reaction
    properties
        Name string
        NominalRate double {mustBeScalarOrEmpty}
        TemperatureCorrectionFactor double {mustBeScalarOrEmpty}
        Order dictionary
        HalfSaturationConstants dictionary
        StoichiometricCoefficients dictionary
        EfficiencyBiofilm double = 1.0
        EfficiencyFlowing double = 1.0
        IsLightDependent logical
        MinimumLightFactor double
        OptimalLightFactor double
        % Dark-switch (Wolf2007 PHOBIA r6): when > 0, the reaction's light
        % factor is K/(K + I_local), I_local in optimal-intensity units --
        % active in darkness, suppressed in light. Mutually exclusive with
        % IsLightDependent. 0.0 disables (default; goldens unchanged).
        LightInhibition double = 0.0
    end

    methods
        function obj = Reaction(input)
            arguments
                input.Name string = ""
                input.NominalRate double = 1.0
                input.TemperatureCorrectionFactor double = 1.0;
                input.Order dictionary = dictionary("", 0);
                input.HalfSaturationConstants dictionary = dictionary("", 0)
                input.StoichiometricCoefficients dictionary = dictionary("", 0)
                input.EfficiencyBiofilm double = 1.0;
                input.EfficiencyFlowing double = 1.0;
                input.IsLightDependent logical = false
                input.MinimumLightFactor double = 0.0;
                input.OptimalLightFactor double = 0.0;
                input.LightInhibition double = 0.0;
            end
            for field = string(fieldnames(input).')
                if isa(input.(field), "dictionary") && ~isa(input.(field).keys, "string")
                    input_field = dictionary([input.(field).keys.Name].', input.(field).values);
                else
                    input_field = input.(field);
                end
                obj.(field) = input_field;
            end
        end
        
        mu = computeRate(obj, temperature);
    end
end