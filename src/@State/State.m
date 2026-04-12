classdef State
    properties
        SandFilter  SandFilter
        Model  Model
        
        Time double = 0.0
        GlobalConcentration
        EnclosedWaterVolume
        Velocity
    end
    
    properties (Hidden, Dependent)
        GlobalConcentrationBiofilm
        GlobalConcentrationFlowing
        
        VolumeFractions
    end

    methods
        function obj = State(filter, model)
            arguments
                filter  SandFilter =  SandFilter()
                model  Model =  Model()
            end
            obj.SandFilter = filter;
            obj.Model = model;

            if ~isempty(obj.SandFilter.GridPoints)
                N = length(obj.SandFilter.GridPoints.Centers);
            else
                N = [];
            end

            obj.Velocity.Biofilm = zeros(N-1, 1);
            obj.Velocity.Flowing = zeros(N+1, 1);
            
            kP = length(obj.Model.Particles);
            obj.GlobalConcentration.Matrix = zeros(N, kP);
            obj.GlobalConcentration.EnclosedParticles = zeros(N, kP);
            obj.GlobalConcentration.FlowingParticles = zeros(N, kP);

            kL = length(obj.Model.Liquids);
            obj.GlobalConcentration.EnclosedLiquids = zeros(N, kL);
            obj.GlobalConcentration.FlowingLiquids = zeros(N, kL);
            
            obj.EnclosedWaterVolume = zeros(N, 1);
        end
        
        function C = get.GlobalConcentrationBiofilm(obj)
            globalC = obj.GlobalConcentration;
            C = [globalC.Matrix, globalC.EnclosedParticles, globalC.EnclosedLiquids];
        end

        function C = get.GlobalConcentrationFlowing(obj)
            globalC = obj.GlobalConcentration;
            C = [globalC.FlowingParticles, globalC.FlowingLiquids];
        end

        function phi = get.VolumeFractions(obj)
            density = struct("Matrix", [obj.Model.Particles.Density], ...
                             "EnclosedParticles", [obj.Model.Particles.Density], ...
                             "EnclosedLiquids", [obj.Model.Liquids.Density], ...
                             "EnclosedWater", obj.Model.DensityWater, ...
                             "FlowingParticles", [obj.Model.Particles.Density], ...
                             "FlowingLiquids", [obj.Model.Liquids.Density]);
            phi = struct("Matrix", [], "EnclosedParticles", [], "EnclosedLiquids", [], ...
                                         "FlowingParticles", [], "FlowingLiquids", []);
            
            fieldnames = ["Enclosed", "Flowing"] + ["Particles"; "Liquids"];
            for field = ["Matrix", fieldnames(:)']
                phi.(field) = obj.GlobalConcentration.(field) ./ density.(field);
            end
        end
    end
end