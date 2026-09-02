function h_line = plotSupernatantHeight(obj, plot_options, options)
% PLOTSUPERNATANTHEIGHT  Equivalent supernatant biofilm height against time.
%
%   plotSupernatantHeight(obj)
%   plotSupernatantHeight(obj, struct(LineWidth=2), HeightUnits="mm")
%   plotSupernatantHeight(obj, [], AxisHandle=ax, DisplayName="Winter")
%
% Draws the time series returned by getSupernatantHeight; see that file for what the
% measure is and why it is an integral rather than a threshold.
arguments
    obj (1,1) Results
    plot_options struct = struct.empty;

    options.AxisParent = matlab.ui.Figure.empty;
    options.AxisHandle = matlab.graphics.axis.Axes.empty;
    options.FontSize (1,1) double = 12;
    options.Grid (1,1) logical = true;

    options.Reference (1,1) double {mustBePositive} = 0.3;
    options.SurfaceDepth (1,1) double = 0;
    options.HeightUnits (1,1) string ...
        {mustBeMember(options.HeightUnits, ["m", "cm", "mm"])} = "mm";
    options.DisplayName (1,1) string = "";
end
[h, t] = getSupernatantHeight(obj, ...
    Reference=options.Reference, SurfaceDepth=options.SurfaceDepth);

switch options.HeightUnits
    case "m",  h_f = 1;
    case "cm", h_f = 100;
    case "mm", h_f = 1000;
end

ax = options.AxisHandle;
if isempty(ax)
    ax = mpc_axes(AxisParent=options.AxisParent);
end
name = options.DisplayName;
if name == "", name = sprintf("$\\phi_{\\rm ref} = %g$", options.Reference); end

% Guard the caller's hold state: this is a time series meant to be overlaid one
% scenario at a time, so it must not clear an axes it was handed.
wasHeld = ishold(ax);
hold(ax, "on");
if isempty(plot_options)
    h_line = plot(ax, t, h_f * h, DisplayName=name);
else
    plot_args = namedargs2cell(plot_options);
    h_line = plot(ax, t, h_f * h, plot_args{:}, DisplayName=name);
end
if ~wasHeld, hold(ax, "off"); end

xlabel(ax, "Time [d]");
ylabel(ax, sprintf("Supernatant biofilm height [%s]", options.HeightUnits));
set(ax, FontSize=options.FontSize);
grid(ax, matlab.lang.OnOffSwitchState(options.Grid));
end
