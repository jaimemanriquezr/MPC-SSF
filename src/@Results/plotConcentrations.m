function plotConcentrations(obj, options)
arguments
    obj Results
    options.sub_phases = [ {["Particles" ,"Matrix"]; 
                            ["Particles", "Enclosed"]; 
                            ["Liquids", "Enclosed"]; 
                            ["Particles", "Flowing"]; 
                            ["Liquids", "Flowing"]}, {'M'; 'Pe'; 'Le'; 'Pf'; 'Lf'}];
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
        sp_axes(i) = axes(figure(), 'Units', 'centimeters', 'Position', [3 3 7 7]);
        surf(sp_axes(i),T,Z,Ci_gamma);
        shading interp;
        title(sp_axes(i), sprintf("$\\rm %s$",component_name), ...
            'Interpreter','latex')
        c(i) = colorbar;
        ylim(sp_axes(i), [-1 1]);
        xlabel("Time $t$ [days]", 'Interpreter', 'latex');
        ylabel("Depth $z$ [m]", 'Interpreter', 'latex');

        axis(sp_axes(i), 'square')
        set(sp_axes(i),'XLimitMethod','tight')
        set(sp_axes(i),'YLim', [-obj.SandFilter.Height, obj.SandFilter.Depth])
        set(sp_axes(i), 'YDir', 'reverse')
        view(sp_axes(i), 2)
        colormap turbo
        c(i).Label.String = "Concentration [kg/m$^3$]";
        c(i).Label.Interpreter = 'latex';
        c(i).Label.FontSize = 11;
    
        filename = options.filename_append + sprintf("%s_%s.png", phase_name, component_name);
        exportgraphics(sp_axes(i), filename, 'Resolution', 300)
    end
end
end