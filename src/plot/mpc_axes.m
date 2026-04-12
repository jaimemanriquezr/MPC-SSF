function ax = mpc_axes(options)
    arguments
        options.AxisHandle = matlab.graphics.axis.Axes.empty;
        options.AxisParent = matlab.ui.Figure.empty;
        options.Grid = true;

        options.Size = 700
        options.SizeUnits = "points"

        options.ParentPosition = [0 0 1.5 1];
        options.AxisPosition = [.15 .1 .8 .8];

        options.FontSize = 14
        options.LegendFontSize = 10;

        options.XLabel = "Depth $z$ [m]";
        options.YLabel = "";
        options.Legend = "show";
        options.LegendLocation = "Southeast";
        options.YReverse = false;
        options.YScale = "linear";

        options.AxisRatio = "square";
    end

    ax = options.AxisHandle;
    if isempty(ax)
        if isempty(options.AxisParent)
            ax = axes(figure(), "NextPlot", "add");
        else
            ax = axes(options.AxisParent, "NextPlot", "add");
        end
    end

    [ax.Units, ax.Parent.Units] = deal(options.SizeUnits);
    if ~isinf(options.Size)
        if ~isa(ax.Parent, "matlab.graphics.layout.TiledChartLayout")
            ax.Parent.Position = options.Size * options.ParentPosition;
            ax.Position = options.Size * options.AxisPosition;
        elseif ~isa(ax.Parent.Parent, "matlab.graphics.layout.TiledChartLayout")
            ax.Parent.Parent.Position = options.Size * options.ParentPosition;
            ax.Parent.Position = options.Size * options.AxisPosition;
        else
            ax.Parent.Parent.Parent.Position = options.Size * options.ParentPosition;
            % ax.Parent.Parent.Position = options.Size * options.AxesPosition;
        end
    end
    axis(ax, options.AxisRatio)

    if options.YReverse
        set(ax, "YDir", "reverse");
    end
    if options.Grid
        grid(ax, "on")
    end
    set(ax, "TickLabelInterpreter", "latex");
    set(ax, "YScale", options.YScale)

    y_label = ylabel(ax, options.YLabel);
    x_label = xlabel(ax, options.XLabel);
    ax_legend = legend(ax, options.Legend);

    set([y_label, x_label], ...
        "interpreter", "latex", ...
        "FontSize", options.FontSize);
    set(ax_legend, 'FontSize', options.LegendFontSize)
    set(ax_legend, 'Interpreter', 'latex')

    if options.Legend ~= "off"
        ax_legend.Location = options.LegendLocation;
    end
end
