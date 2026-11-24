function plotConcentration(obj, component_name, volume_type, options)
arguments
    obj mpcssf.Results
    component_name
    volume_type
    options.Parent = []   % axes/uiaxes to draw into, or a container; see resolveAxes
    options.Colormap = "turbo";
end
ax = mpcssf.internal.resolveAxes(options.Parent);

[T, Z] = meshgrid(obj.Frames.Time, obj.SandFilter.GridPoints.Centers);
C = obj.Frames.Concentrations{component_name, volume_type}{:};
if isscalar(C) || isempty(C)
    C = 0*T;
end
surf(ax, T, Z, C, LineStyle='none');
colorbar(ax);
xlabel(ax, "Time $t$ [days]", Interpreter='latex', FontSize=20);
ylabel(ax, "Depth $z$ [m]", Interpreter='latex', FontSize=20);
colormap(ax, options.Colormap);
set(ax, YDir='reverse', XLimitMethod='tight', YLimitMethod='tight');
view(ax, 2);
