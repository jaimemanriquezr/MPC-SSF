function plotConcentrations(obj, options)
% PLOTCONCENTRATIONS  One figure per component and phase (space-time surfaces).
%
%   plotConcentrations(r)                                   % figures only
%   plotConcentrations(r, Export=true, filename_append="run1_")   % also write PNGs
%
% For a single panel inside an app or a layout, use plotConcentration(..., Parent=ax).
arguments
    obj mpcssf.Results
    options.sub_phases = [ {["Particles" ,"Matrix"];
                            ["Particles", "Enclosed"];
                            ["Liquids", "Enclosed"];
                            ["Particles", "Flowing"];
                            ["Liquids", "Flowing"]}, {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'}];
    options.Export (1,1) logical = false   % write <filename_append><phase>_<name>.png
    options.filename_append = "";
end
sub_phases = options.sub_phases;

[T, Z] = meshgrid(obj.Frames.Time, obj.SandFilter.GridPoints.Centers);

for s = 1:length(sub_phases(:,2))
    component_type = sub_phases{s,1}(1);
    volume_type = sub_phases{s,1}(2);
    phase_name = sub_phases{s,2};
    names = [obj.Model.(component_type).Name];
    for i = 1:length(names)
        component_name = names(i);
        Ci_gamma = obj.Frames.Concentrations{component_name, volume_type}{:};
        ax = axes(figure(), 'Units', 'centimeters', 'Position', [3 3 7 7]);
        surf(ax, T, Z, Ci_gamma);
        shading(ax, 'interp');
        title(ax, sprintf("$\\rm %s$", component_name), 'Interpreter', 'latex')
        c = colorbar(ax);
        xlabel(ax, "Time $t$ [days]", 'Interpreter', 'latex');
        ylabel(ax, "Depth $z$ [m]", 'Interpreter', 'latex');

        axis(ax, 'square')
        set(ax, 'XLimitMethod', 'tight')
        set(ax, 'YLim', [-obj.SandFilter.Height, obj.SandFilter.Depth])
        set(ax, 'YDir', 'reverse')
        view(ax, 2)
        colormap(ax, 'turbo')
        c.Label.String = "Concentration [kg/m$^3$]";
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 11;

        if options.Export
            filename = options.filename_append + sprintf("%s_%s.png", phase_name, component_name);
            exportgraphics(ax, filename, 'Resolution', 300)
        end
    end
end
end
