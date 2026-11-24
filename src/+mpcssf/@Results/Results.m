classdef Results
    properties
        SandFilter
        Model
        Frames
        TimeStart
        TimeFinal
        Flag
        SimulationData
        % The exact solver configuration this run used. Carried WITH the data so
        % a result is self-describing: scheme, boundary condition, timestep caps
        % and the implicit-treatment switches travel with the frames rather than
        % living in a filename or a lab notebook. Set by simulate().
        SolverOptions mpcssf.SolverOptions = mpcssf.SolverOptions.empty
    end

    methods
        function obj = Results(filter,model)
            arguments
                filter  mpcssf.SandFilter =  mpcssf.SandFilter();
                model  mpcssf.Model = mpcssf.Model
            end

            obj.SandFilter = filter;
            obj.Model = model;
            obj.Frames = struct( ...
                'Time', [],...
                'Concentrations',[], ...
                'Velocity',[]);
            obj.TimeStart = 0;
            obj.TimeFinal = 0;
            obj.Flag = "UNINITIATED";
            obj.SimulationData = [];
            obj.SolverOptions = mpcssf.SolverOptions.empty;
        end


        [phi, z, t] = getVolumeFractions(obj);
        plotConcentrations(obj);
        plotConcentration(obj);

        function plot(obj, variable, options)
            arguments
                obj mpcssf.Results
                variable = "biofilm"
                options.Parent = []   % used by "biofilm"; see resolveAxes
            end
            switch variable
                case {"biofilm", "phi_b"}
                    plotBiofilm(obj, obj.TimeFinal, struct.empty, Parent=options.Parent);
                case {"concentrations"}
                    plotConcentrations(obj)
            end
        end

        % concatenate/plus removed 2026-09-22: both read lower-case fields
        % (.flag, .model, .time_start) that this class does not define -- it uses
        % Flag, Model, TimeFinal -- so any call threw. Zero callers repo-wide.

        bf_line = plotBiofilm(obj, plot_options, t, options)



    end
end
