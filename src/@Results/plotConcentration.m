function plotConcentration(obj, component_name, volume_type, options)
arguments
    obj Results
    component_name
    volume_type
    options.AxesHandle = axes(figure());
    options.Colormap = "turbo";
end
ax = options.AxesHandle;

[T, Z] = meshgrid(obj.Frames.Time, obj.SandFilter.GridPoints.Centers);
C = obj.Frames.Concentrations{component_name, volume_type}{:};
if isscalar(C) || isempty(C)
    C = 0*T;
end
surf(ax, T, Z, C, LineStyle='none');
%title(ax, sprintf("$\\rm %s$",component_name), Interpreter='latex');
colorbar(ax);
xlabel("Time $t$ [days]", Interpreter='latex', FontSize=20);
ylabel("Depth $z$ [m]", Interpreter='latex', FontSize=20);
colormap(ax, options.Colormap);
set(ax, YDir='reverse', XLimitMethod='tight', YLimitMethod='tight');
view(ax, 2);
end