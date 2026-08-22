classdef Particle < Component & matlab.mixin.CustomDisplay
    properties
        AttachmentMatrix double = 0.0 % [1/day]
        AttachmentSand double = 0.0 % [1/day]
        Attenuation double = 0.0 % [m^2/kg]
        % Scales attachment to BARE SAND relative to AttachmentSand; the biofilm
        % share of the flowing-attachment factor is unaffected. 1.0 (the default)
        % reproduces the undifferentiated behaviour of every existing preset.
        SandAttachmentFactor double = 1.0 % [-]
    end
    methods
        function obj = Particle(input)
            arguments
                input.Name string = string.empty
                input.Density double = []
                input.Dispersivity double = []
                input.AttachmentSand double = []
                input.AttachmentMatrix double = []
                input.Transport double = []
                input.Attenuation double = []
                input.SandAttachmentFactor double = []
            end
            extra_fields = ["AttachmentSand", "AttachmentMatrix", "Attenuation", ...
                            "SandAttachmentFactor"];
            input_args = namedargs2cell(rmfield(input, extra_fields));
            obj = obj@ Component(input_args{:});
            for field = extra_fields
                input_value = input.(field);
                if ~isempty(input_value)
                    obj.(field) = input_value;
                end
            end
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