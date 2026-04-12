classdef Results
    properties
        SandFilter
        Model
        Frames
        TimeStart
        TimeFinal
        Flag
        SimulationData
    end

    methods
        function obj = Results(filter,model)
            arguments
                filter  SandFilter =  SandFilter();
                model  Model = Model
            end

            obj.SandFilter = filter;
            obj.Model = model;
            obj.Frames = struct( ...
                'time', [],...
                'concentrations',[], ...
                'velocity',[]);
            obj.TimeStart = 0;
            obj.TimeFinal = 0;
            obj.Flag = "UNINITIATED";
            obj.SimulationData = [];
        end

        function obj = plus(obj1, obj2)
            obj = concatenate(obj1, obj2);
        end

        function obj = concatenate(obj,new_obj)
            if isequal(obj.Flag,"UNINITIATED")
                obj = new_obj;
                return
            end

            if isequal(new_obj.flag,"UNINITIATED")
                return
            end


            if isequal(obj.SandFilter, new_obj.SandFilter) && isequal(obj.Model, new_obj.model)
                if isequal(obj.TimeFinal, new_obj.time_start)
                    obj.Frames.time = [obj.Frames.time; new_obj.frames.time(2:end)];
                    for field = "concentrations"
                        for i = 1:size(obj.Frames.(field),1)
                            for j = 1:size(obj.Frames.(field),2)
                                if isequal(obj.Frames.(field){i,j}, {0})
                                    continue;
                                end
                                obj.Frames.(field){i,j} = ...
                                    {[obj.Frames.(field){i,j}{:}, new_obj.frames.(field){i,j}{:}(:,2:end,:)]};
                            end
                        end
                    end

                    obj.Frames.velocity.biofilm = [obj.Frames.velocity.biofilm, new_obj.frames.velocity.biofilm(:, 2:end)];
                    obj.Frames.velocity.flowing = [obj.Frames.velocity.flowing, new_obj.frames.velocity.flowing(:, 2:end)];

                    obj.TimeFinal = new_obj.time_final;
                    obj.Flag = new_obj.flag;

                    field_names = string(fieldnames(obj.SimulationData));
                    for fname = field_names.'
                        if isa(obj.SimulationData.(fname), "function_handle")
                            continue;
                        end
                        if fname == "biofilm_velocity"
                            obj.SimulationData.(fname) = [obj.SimulationData.(fname); new_obj.simulation_data.(fname)];
                        else
                            obj.SimulationData.(fname) = [obj.SimulationData.(fname), new_obj.simulation_data.(fname)];
                        end
                    end

                else
                    error("Initial time of concatenated object does not coincide with final time of the original.")
                end
            else
                error("Filter or model not compatible between results objects.")
            end
        end

        bf_line = plotBiofilm(obj, plot_options, t, options)



    end
end
