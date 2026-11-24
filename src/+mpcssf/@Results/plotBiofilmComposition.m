function plotBiofilmComposition(obj, time, options)
arguments
    obj mpcssf.Results
    time = obj.TimeFinal

    options.DepthLimits = [0 obj.SandFilter.Depth];
    options.Phase (1,1) string {mustBeMember(options.Phase, ["Biofilm", "Matrix"])} = "Biofilm";
    options.Parent = []   % axes/uiaxes to draw into, or a container; see resolveAxes
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

ax = mpcssf.internal.resolveAxes(options.Parent, true);
% The bands are CUMULATIVE fractions drawn largest-first, so each one overplots
% the previous. Without hold the axes are cleared on every area() call and only
% the last (smallest) band survives.
% Colours are keyed to species name (not loop/draw order) so each band gets its
% intended colour regardless of the largest-first draw order.
speciesColors = containers.Map( ...
    ["PHO", "HET", "POM", "PAT"], ...
    {[0.20 0.55 0.25], [0.75 0.20 0.18], [0.55 0.55 0.55], [0.20 0.40 0.70]});
fallbackColor = [0.40 0.40 0.40];
% Species with no presence in the plotted depth window are dropped (band + legend)
% rather than drawn as a dead entry. Tested over the DepthLimits crop only, since a
% species can be zero in the crop but present elsewhere in the domain. Threshold is
% relative to a volume fraction (O(1) quantity): 1e-6 is far below anything physically
% resolved, so this only catches species that are genuinely absent here, not merely small.
zeroTol = 1e-6;
inWindow = z >= options.DepthLimits(1) & z <= options.DepthLimits(2);
present = max(abs(particle_biofilm_volfrac(inWindow, :)), [], 1) > zeroTol;
wasHeld = ishold(ax);
hold(ax, "on");
for i = 0:kp-1
    if ~present(end-i)
        continue
    end
    particle_name = obj.Model.Particles(end-i).Name;
    if isKey(speciesColors, particle_name)
        faceColor = speciesColors(particle_name);
    else
        faceColor = fallbackColor;
    end
    area(ax, z, particle_biofilm_relfrac(:, end-i), DisplayName=particle_name, FaceColor=faceColor)
end
if ~wasHeld, hold(ax, "off"); end

xlim(ax, options.DepthLimits);
ylim(ax, [0 1]);
xlabel(ax, "Depth [m]");
ylabel(ax, "Relative volume fraction of " + lower(options.Phase));
set(ax, FontSize=options.FontSize);
grid(ax, matlab.lang.OnOffSwitchState(options.Grid));
% legend() with no handle targets gca, which is NOT ax when the caller passed
% Parent (a subplot, say) -- the legend would land on someone else's axes.
legend(ax);
end
