function ax = resolveAxes(parent, styled, axesOptions)
% RESOLVEAXES  The axes a plot method draws into, given its Parent= option.
%
%   parent empty            -> new figure (styled with mpc_axes when `styled`)
%   parent axes or uiaxes   -> that axes, untouched (no resizing, no new figure)
%   parent container        -> new axes inside it (figure, uifigure, panel, tab,
%                              tiled layout); never resizes the container
%
% axesOptions are mpc_axes name-value pairs (labels, YReverse, ...), used only
% when the axes is created here and `styled` is true.
arguments
    parent = []
    styled (1,1) logical = false
    axesOptions cell = {}
end
if isa(parent, "matlab.graphics.axis.Axes")
    ax = parent;
elseif isempty(parent)
    if styled
        ax = mpcssf.internal.mpc_axes(axesOptions{:});
    else
        ax = axes(figure());
    end
elseif styled
    ax = mpcssf.internal.mpc_axes(axesOptions{:}, AxisParent=parent, Size=Inf);
else
    ax = axes(parent);
end
end
