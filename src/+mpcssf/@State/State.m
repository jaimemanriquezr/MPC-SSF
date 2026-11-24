classdef State
    properties
        SandFilter  mpcssf.SandFilter
        Model  mpcssf.Model
        
        Time double = 0.0
        GlobalConcentration
        EnclosedWaterVolume
        Velocity
    end
    
    properties (Hidden, Dependent)
        GlobalConcentrationBiofilm
        GlobalConcentrationFlowing
    end

    methods
        function obj = State(filter, model)
            arguments
                filter  mpcssf.SandFilter =  mpcssf.SandFilter()
                model  mpcssf.Model =  mpcssf.Model()
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

        % get.VolumeFractions removed 2026-09-22: it read obj.Model.DensityWater,
        % and Model has no such property (it is WaterDensity), so it threw on any
        % access. Zero callers repo-wide; Results.getVolumeFractions is the live one.
    end
end