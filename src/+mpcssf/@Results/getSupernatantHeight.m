function [h, t] = getSupernatantHeight(obj, options)
% GETSUPERNATANTHEIGHT  Equivalent height of the biofilm standing in the supernatant.
%
%   [h, t] = getSupernatantHeight(obj)
%   [h, t] = getSupernatantHeight(obj, Reference=0.3, SurfaceDepth=0)
%
% Returns, for every frame, the thickness the biofilm held ABOVE the sand surface
% would have if it were compacted into a free-standing layer of volume fraction
% `Reference`:
%
%       h(t) = ( \int_{z < z_s} eps(z) phi_b(z,t) dz ) / Reference          [m]
%
% with z_s = `SurfaceDepth` (z = 0, the top of the sand, by default), eps the
% porosity and phi_b the biofilm volume fraction from getVolumeFractions. The
% porosity weight makes eps*phi_b a fraction of TOTAL volume, so the integral is a
% biofilm volume per unit filter area -- the convention of the manuscript's
% bed integrals. Do not drop it: between -SandRoughness and 0 the
% porosity ramps from 1 down to SandPorosity, and that is exactly the band the
% supernatant biofilm occupies in the fld2x chains.
%
% WHY AN EQUIVALENT HEIGHT AND NOT A THRESHOLD HEIGHT. The obvious measure -- the
% most negative z at which phi_b exceeds some threshold -- is unusable on this grid.
% At N = 500 (dz ~ 2 mm) the supernatant biofilm of the working-set chains lives in
% one or two cells, so a threshold height can only ever return 0, 2 or 4 mm: it
% cannot separate the scenarios, and it flips between 4 mm and 2 mm on the same
% frame as the threshold moves from 1e-3 to 1e-2. The integral measure is continuous
% in time and threshold-free.
%
% CHOICE OF `Reference`. 0.3 is the biofilm volume fraction found at the sand
% surface in the working set (0.271 dark, 0.271 cov01, 0.306 winter, 0.324 summer at
% their final frames), so h reads as "the supernatant biofilm compacted to the
% density it has where it grew from". It is only a scale: h is linear in
% 1/Reference, so the ORDERING of scenarios does not depend on this number.
%
% CAVEAT worth carrying into any caption: on the fld2x chains h is of order 1 mm,
% below the cell size dz = 1.998 mm. The supernatant biofilm is not resolved at
% N = 500; h is a sub-cell quantity read off one or two cells.
% METHOD = "threshold" (added 2026-09-05). h(t) = -min{ z < z_s : phi_b > Threshold },
% i.e. the top of the biofilm measured up from the sand surface. The equivalent
% measure above was adopted because on the fld2x working set the supernatant film
% occupied one or two cells and a threshold height could only return 0, 2 or 4 mm.
% That objection is regime-dependent, NOT permanent: on the 2026-09-05 working set
% (AttachmentMatrix = 2.735e4) the mat spans ~8 cells and the threshold height is
% stable to one cell across phi* = 1e-3 ... 0.1 (9.99/13.99/15.98 mm at 10/20/30 d).
% Check that stability before quoting it; if the film is thin, use "equivalent".
% Threshold defaults to Zeta1 = 0.27, the free-energy phase boundary: Psi''(u) =
% 3u(u - Zeta1) changes sign there, so it separates the dense diffusive biofilm
% phase from the dilute spinodal tail rather than being an arbitrary cutoff.
arguments
    obj (1,1) mpcssf.Results
    options.Reference (1,1) double {mustBePositive} = 0.3;
    options.SurfaceDepth (1,1) double = 0;
    options.Method (1,1) string {mustBeMember(options.Method, ["equivalent","threshold"])} = "equivalent";
    options.Threshold (1,1) double {mustBeNonnegative} = 0.27;
end
[phi, z, t] = getVolumeFractions(obj);
z = z(:);
eps = computePorosity(obj.SandFilter, z);
dz = obj.SandFilter.GridSize;

sup = z < options.SurfaceDepth;
if options.Method == "threshold"
    % Sub-cell crossing: the topmost cell ABOVE the threshold is i, the cell above
    % it (smaller z) is below it, so phi crosses phi* between them. Interpolating
    % linearly removes the dz staircase -- without it h can only take multiples of
    % the cell size and jumps 2 mm at a time.
    h = zeros(1, size(phi.Biofilm, 2));
    for j = 1:numel(h)
        pbj = phi.Biofilm(:, j);
        i = find(sup & pbj > options.Threshold, 1, "first");
        if isempty(i), continue, end
        if i > 1
            pHi = pbj(i-1);                 % above the crossing, below phi*
            pLo = pbj(i);                   % at or below the crossing, above phi*
            f = (options.Threshold - pHi)/max(pLo - pHi, realmin);
            f = min(max(f, 0), 1);
            zc = z(i-1) + f*(z(i) - z(i-1));
        else
            zc = z(i);
        end
        h(j) = options.SurfaceDepth - zc;
    end
else
    h = sum(eps(sup) .* phi.Biofilm(sup, :), 1) * dz / options.Reference;
end
t = t(:).';
end
