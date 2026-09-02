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
% biofilm volume per unit filter area -- the same convention as `rec.bedInt` in
% analysis/probes/probeChain.m. Do not drop it: between -SandRoughness and 0 the
% porosity ramps from 1 down to SandPorosity, and that is exactly the band the
% supernatant biofilm occupies in the fld2x chains.
%
% WHY AN EQUIVALENT HEIGHT AND NOT A THRESHOLD HEIGHT. The obvious measure -- the
% most negative z at which phi_b exceeds some threshold -- is unusable on this grid.
% At N = 500 (dz ~ 2 mm) the supernatant biofilm of the working-set chains lives in
% one or two cells, so a threshold height can only ever return 0, 2 or 4 mm: it
% cannot separate the scenarios, and it flips between 4 mm and 2 mm on the same
% frame as the threshold moves from 1e-3 to 1e-2. The integral measure is continuous
% in time and threshold-free. See .claude/plans/2026-09-02-supernatant-height.md.
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
arguments
    obj (1,1) Results
    options.Reference (1,1) double {mustBePositive} = 0.3;
    options.SurfaceDepth (1,1) double = 0;
end
[phi, z, t] = getVolumeFractions(obj);
z = z(:);
eps = computePorosity(obj.SandFilter, z);
dz = obj.SandFilter.GridSize;

sup = z < options.SurfaceDepth;
h = sum(eps(sup) .* phi.Biofilm(sup, :), 1) * dz / options.Reference;
t = t(:).';
end
