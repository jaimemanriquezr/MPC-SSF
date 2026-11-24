function bf_line = plotBiofilm(obj, t, plot_options, options)
    arguments
        obj (1,1)  mpcssf.Results
        t (1,1) = obj.TimeFinal;
        plot_options struct = struct.empty;

        options.Parent = []   % axes/uiaxes to draw into, or a container; see resolveAxes
        options.FontSize = 12;
        options.Grid = true;

        options.Phase = "Biofilm";
        options.DepthUnits = "m";
    end
    switch options.Phase
        case {"Biofilm", "biofilm", "b"}
            phi = getVolumeFractions(obj).Biofilm;
            phi_label = "Biofilm volume fraction $\phi_{\rm b}$";
        case {"Matrix", "matrix", "M"}
            phi = getVolumeFractions(obj).Matrix;
            phi_label = "Matrix volume fraction $\phi_{\rm M}$";
        case {"Enclosed", "enclosed", "e"}
            phi = getVolumeFractions(obj).Enclosed;
            phi_label = "Enclosed volume fraction $\phi_{\rm e}$";
    end

    created = ~isa(options.Parent, "matlab.graphics.axis.Axes");
    ax = mpcssf.internal.resolveAxes(options.Parent, true, ...
        {"Size", Inf, "YLabel", phi_label, "YReverse", true});
    if created
        ax.YLim = [-obj.SandFilter.Height, obj.SandFilter.Depth];
    end

    switch options.DepthUnits
        case 'm'
            z_f = 1;
        case 'cm'
            z_f = 100;
        case 'mm'
            z_f = 1000;
    end
    z = z_f * obj.SandFilter.GridPoints.Centers;
    t = min(t, obj.TimeFinal);
    frame = find((obj.Frames.Time - t) >= 0, 1);
    if isempty(frame)
        frame = find(obj.Frames.Time > 0, 1, 'last');
    end
    phi = phi(:, frame);

    if isempty(plot_options)
        bf_line = plot(ax, phi, z, DisplayName=sprintf("T = %g", t));
    else
        plot_args = namedargs2cell(plot_options);
        bf_line = plot(ax, phi, z, plot_args{:}, DisplayName=sprintf("T = %g", t));
    end
end
