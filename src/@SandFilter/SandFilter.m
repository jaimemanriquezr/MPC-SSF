classdef SandFilter
    properties
        Height double
        Depth double

        SandPorosity double
        SandRoughness double

        InflowVelocity double

        Temperature double
        LightIrradiation function_handle
        LightAttenuationCoeffWater double
        LightAttenuationCoeffSand double

        GridPoints
    end
    
    properties (Dependent, Hidden)
        GridSize
        GridZero
        
        LightAttenuationEtaWater
        LightAttenuationEtaSand
    end       

    methods
        function obj = SandFilter(input)
            arguments
                input.Height {mustBeNumeric, mustBeScalarOrEmpty} = 1;
                input.Depth {mustBeNumeric, mustBeScalarOrEmpty} = 1;
                input.SandPorosity {mustBeNumeric, mustBeScalarOrEmpty} = .4;
                input.SandRoughness {mustBeNumeric, mustBeScalarOrEmpty} = 5e-3;
                input.InflowVelocity (1,1) {mustBeNumeric} = 0.3 * 24;
                input.Temperature (1,1) {mustBeNumeric} = 15 + 273;
                input.LightIrradiation = @(t) 0.8 * max(sin(2*pi*(t - 13/48)) + 31/50, 0) / (1 + 31/50);
                input.LightAttenuationCoeffWater = 0.32;
                input.LightAttenuationCoeffSand = 1500;
                input.GridPoints = [];
            end

            for field = string(fieldnames(input).')
                obj.(field) = input.(field);
            end
        end
        
        function dz = get.GridSize(obj)
            dz = mean(diff(obj.GridPoints.Boundaries));
        end

        function n0 = get.GridZero(obj)
            n0 = find(abs(obj.GridPoints.Centers) < obj.GridSize/2);
        end

        function eta = get.LightAttenuationEtaWater(obj)
            eta = (obj.GridPoints.Centers + obj.Height);
            eta = obj.LightAttenuationCoeffWater*eta;
        end

        function eta = get.LightAttenuationEtaSand(obj)
            eta = zeros(size(obj.GridPoints.Centers));
            roughIndex = obj.GridPoints.Centers >= -obj.SandRoughness;
            packedIndex = obj.GridPoints.Centers >= 0;
            roughIndex(packedIndex) = false;

            etaSandDelta1 = (1 - obj.SandPorosity)*(obj.GridPoints.Centers(roughIndex) + obj.SandRoughness);
            etaSandDelta2 = (1 + (obj.GridPoints.Centers(roughIndex) - obj.SandRoughness)/(2*obj.SandRoughness));
            eta(roughIndex) = etaSandDelta1.*etaSandDelta2;
            eta(packedIndex) = (1 - obj.SandPorosity)*(obj.SandRoughness/2 + obj.GridPoints.Centers(packedIndex));

            eta = obj.LightAttenuationCoeffSand*eta;
        end

        function obj = set.GridPoints(obj, points)
            if isempty(points)
                obj.GridPoints = struct.empty;
            else
                points = points(:);
                centers = 0.5*(points(1:end-1) + points(2:end));
                obj.GridPoints = struct("Boundaries", points, "Centers", centers);
            end
        end
        
        function obj = addGridPoints(obj, intervalNumber)
            dz = obj.Depth / (intervalNumber + 1/2);
            lowerGrid = [obj.Depth+dz, obj.Depth:(-dz):(-obj.Height)];
            if all(abs(lowerGrid) > 1e-16)
                obj.GridPoints = fliplr(lowerGrid);
            else % a computational boundary is approximately at z = 0
                obj = addGridPoints(obj, intervalNumber+1);
            end
        end
        
        function plot(obj, axisHandle)
            arguments
                obj SandFilter
                axisHandle matlab.graphics.axis.Axes = axes(figure());
            end
            axisHandle.NextPlot = "add";

            H = obj.Height;
            B = obj.Depth;
            % plot(axisHandle, [0, 0], [-H, B], LineStyle="-", Color="k", LineWidth=3.0);
            z = linspace(-H, B, 1e3);
            plot(axisHandle, obj.computePorosity(z), z, Color="k", LineWidth=3.0);
            xlim(axisHandle, [0 1.1]);
            ylabel(axisHandle, "Depth $z$ [m]", Interpreter="latex", FontSize=20);
            xlabel(axisHandle, "Porosity $\tilde{\varepsilon}(z)$ [-]", Interpreter="latex", FontSize=20);
            axisHandle.YDir = "reverse";
            axis(axisHandle, "square");
        end

        % OVERRIDING FUNCTIONS
        function is_eq = isequal(self, other)
            is_eq = true;
            for field = string(fieldnames(self).')
                self_field = self.(field);
                other_field = other.(field);

                if isequal(class(self_field),'function_handle')
                    self_field = sym(self_field);
                    other_field = sym(other_field);
                end

                if ~isequal(self_field, other_field)
                    is_eq = false;
                    break;
                end
            end
        end

    end
end
