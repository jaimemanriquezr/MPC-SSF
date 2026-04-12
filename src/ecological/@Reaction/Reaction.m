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
    end

    methods
        function obj = Reaction(input)
            arguments
                input.Name string = ""
                input.NominalRate double = 0.0
                input.TemperatureCorrectionFactor double = 1.0;
                input.Order dictionary = dictionary("", 0);
                input.HalfSaturationConstants dictionary = dictionary("", 0)
                input.StoichiometricCoefficients dictionary = dictionary("", 0)
                input.EfficiencyBiofilm double = 1.0;
                input.EfficiencyFlowing double = 1.0;
                input.IsLightDependent logical = false
                input.MinimumLightFactor double = 0.0;
                input.OptimalLightFactor double = 0.0;
            end
            for field = string(fieldnames(input).')
                obj.(field) = input.(field);
            end
        end
        
        mu = computeRate(obj, temperature);
    end
end