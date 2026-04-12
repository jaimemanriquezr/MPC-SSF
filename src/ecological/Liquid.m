classdef Liquid <  Component & matlab.mixin.CustomDisplay
    methods
        function obj = Liquid(input)
            arguments
                input.Name string = string.empty
                input.Density double = []
                input.Dispersivity double = []
                input.TransportRate double = []
            end
            input_args = namedargs2cell(input);
            obj = obj@ Component(input_args{:});
        end
    end

    methods (Access = protected)
        function displayNonScalarObject(objArray)
            dimStr = matlab.mixin.CustomDisplay.convertDimensionsToString(objArray);
            cName = matlab.mixin.CustomDisplay.getClassNameForHeader(objArray);
             Component.displayHomogeneousNonScalarObject(objArray, dimStr, cName);
        end
    end
end