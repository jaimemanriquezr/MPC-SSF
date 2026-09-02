function plotBiofilmComposition(obj, time, options)
arguments
    obj Results
    time = obj.TimeFinal

    options.DepthLimits = [0 obj.SandFilter.Depth];
    options.Phase (1,1) string {mustBeMember(options.Phase, ["Biofilm", "Matrix"])} = "Biofilm";
    options.AxisHandle = matlab.graphics.axis.Axes.empty;
    options.FontSize = 12;
    options.Grid = true;
end
z = obj.SandFilter.GridPoints.Centers;
t_frame = find((obj.Frames.Time - min(time, obj.TimeFinal)) >= 0, 1);
kp = length(obj.Model.Particles);

cell2arr = @(C) cell2mat(permute(cellfun(@(X) X(:, t_frame), C, UniformOutput=0), [3 2 1]));
particle_concentration = cell2arr(obj.Frames.Concentrations{[obj.Model.Particles.Name], :});

% 1: matrix, 2: enclosed, 3: flowing
switch options.Phase
    case "Biofilm"
        particle_biofilm_mass = squeeze(sum(particle_concentration(:, 1:2, :), 2));
    case "Matrix"
        particle_biofilm_mass = squeeze(sum(particle_concentration(:, 1, :), 2));
end
particle_biofilm_volfrac = particle_biofilm_mass ./ [obj.Model.Particles.Density];
total_biofilm_volfrac = sum(particle_biofilm_volfrac, 2);
particle_biofilm_relfrac = cumsum(particle_biofilm_volfrac, 2) ./ total_biofilm_volfrac;

if isempty(options.AxisHandle)
    ax = mpc_axes();
else
    ax = options.AxisHandle;
end
% The bands are CUMULATIVE fractions drawn largest-first, so each one overplots
% the previous. Without hold the axes are cleared on every area() call and only
% the last (smallest) band survives.
wasHeld = ishold(ax);
hold(ax, "on");
for i = 0:kp-1
    particle_name = obj.Model.Particles(end-i).Name;
    area(ax, z, particle_biofilm_relfrac(:, end-i), DisplayName=particle_name)
end
if ~wasHeld, hold(ax, "off"); end

xlim(ax, options.DepthLimits);
ylim(ax, [0 1]);
xlabel(ax, "Depth [m]");
ylabel(ax, "Relative volume fraction of " + lower(options.Phase));
set(ax, FontSize=options.FontSize);
grid(ax, matlab.lang.OnOffSwitchState(options.Grid));
% legend() with no handle targets gca, which is NOT ax when the caller passed
% AxisHandle (a subplot, say) -- the legend would land on someone else's axes.
legend(ax);
end
