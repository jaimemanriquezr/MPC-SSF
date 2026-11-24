classdef Component < matlab.mixin.Heterogeneous
    properties
        Name string = ""
        Density double = 0.0 % [kg/m^3]
        Dispersivity double = 0.0 % [m]
        TransportRate double = 0.0 % [1/day]
    end
    methods
        function obj = Component(input)
            arguments
                input.Name string = string.empty
                input.Density double = []
                input.Dispersivity double = []
                input.TransportRate double = []
            end
            for field = string(fieldnames(input).')
                input_value = input.(field);
                if ~isempty(input_value)
                    obj.(field) = input_value;
                end
            end
        end
    end

    methods (Static)
        function displayHomogeneousNonScalarObject(objArray, dimStr, cName)
            headerStr = [dimStr,' ',cName,' array:'];
            fprintf('%s\n',headerStr);

            propNames = string(properties(objArray));
            propValues = cellfun(@(s) [objArray.(s)], propNames, "UniformOutput", false);
            objNames = propValues{propNames == "Name"};
            propValues(propNames == "Name") = [];
            propNames(propNames == "Name") = [];
            propValues = cell2mat(propValues);
            
            nameLength = sprintf("%%%is", max(max(strlength(objNames)), 4) + 2);
            propNameLengthString = sprintf(" %%%is", max((strlength(propNames)), 10) + 2);
            propNameLengthDouble = sprintf(" %%%i.4e", max((strlength(propNames)), 10) + 2);
            fprintf(nameLength + propNameLengthString + "\n", ["Name"; propNames]);
            for i = 1:length(objNames)
                fprintf(nameLength + propNameLengthDouble + "\n", objNames(i), propValues(:, i));
            end
        end
    end
end