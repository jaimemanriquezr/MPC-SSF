function out = probeMassClosureIsolation(opts)
% PROBEMASSCLOSUREISOLATION  Isolation battery for the PAT export non-closure.
%
% probeMassClosure (2026-08-27, audit Finding 5) leaves PAT with a residual of
% ~5e-5 of supply while HET/PHO close at 1e-11. Each Variant here switches off
% one mechanism; the one that collapses the PAT residual to the HET/PHO level
% names the term. Also records the per-frame residual curve
%   R(t) = [M(t) - M(t0)] - [q*cin*(t - t0) - q*int_t0^t cout dt]
% so the leak's timing (front passage vs steady accumulation) is visible.
% See .claude/plans/2026-09-01-pat-export-closure.md.
arguments
    opts.Variant (1,1) string {mustBeMember(opts.Variant, ...
        ["baseline", "explicitDispersion", "noDispersion", ...
         "noTransfer", "noAttach", "explicitOsmosis", "noOsmosis"])} = "baseline"
    opts.NCells (1,1) double = 100
    opts.Days (1,1) double = 10
    opts.MaxDt (1,1) double = 5e-5
    % probeMassClosure feeds no PAT (influent slot 4 is zero), so its outflow
    % term is never exercised -- the very gap Finding 5 flagged. Feed it here.
    opts.PatInflow (1,1) double = 1.75e-4
    % The budget integrates q*c_out over the saved frames; a residual that
    % falls ~4x per doubling here is trapz error on the breakthrough front,
    % not a solver leak.
    opts.NFrames (1,1) double = 241
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);
S = fullfile(here, "data"); if ~isfolder(S), mkdir(S); end

f = SandFilter(Temperature=19, ...
    LightIrradiation=@(t) 0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
mp = pathogenModel(PhototrophRespiration=0.55, PGExcess=true, NormalizedLight=true);
rx = mp.Reactions;
for k = 1:numel(rx), rx(k).NominalRate = 0.0; end

cs = mp.Components;
osmosisRate = mp.OsmosisRate;
implicitDispersion = true;
implicitOsmosis = true;
switch opts.Variant
    case "explicitDispersion"
        implicitDispersion = false;
    case "noDispersion"
        for i = 1:length(cs), cs(i).Dispersivity = 0; end
    case "noTransfer"
        for i = 1:length(cs)
            if isa(cs(i), "Particle") && cs(i).Name == "PAT", cs(i).TransportRate = 0; end
        end
    case "noAttach"
        % TransportRate=0 does not touch matrix attachment; PAT is the only
        % species holding full flowing concentration over the whole bed, so an
        % attachment-term phase imbalance would show for PAT alone.
        for i = 1:length(cs)
            if isa(cs(i), "Particle") && cs(i).Name == "PAT"
                cs(i).AttachmentMatrix = 0; cs(i).AttachmentSand = 0;
            end
        end
    case "explicitOsmosis"
        implicitOsmosis = false;
    case "noOsmosis"
        osmosisRate = 0;
end

m = Model(cs, rx, Kappa=1e-6, Zeta0=5, Zeta1=0.27, ...
    DetachmentFunction=@(v) 0.14*sqrt(abs(v)/18), ...
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, ...
    OsmosisRate=osmosisRate);
infl = [2.68e-3, 1.00e-2, 0.0, opts.PatInflow, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
r = simulate(State(f, m), InflowConcentrations=infl, SimulationTime=opts.Days, ...
    TimeStep="adaptive", AdaptiveMaxDt=opts.MaxDt, ...
    ImplicitOsmosis=implicitOsmosis, ImplicitDispersion=implicitDispersion, ...
    FrameNumber=opts.NFrames, Quiet=true);
checkRunFlag(r, "probeMassClosureIsolation:" + opts.Variant);

C = r.Frames.Concentrations; ts = r.Frames.Time(:);
z = f.GridPoints.Centers(:); dz = f.GridSize;
eps = computePorosity(f, z);
q = f.InflowVelocity;
% Audit liquids too: they exit at ~half the influent concentration, so they
% exercise the outflow numerics far harder than any particle — and they carry
% no attachment, which makes them a discriminator between a generic
% flowing-transport defect and a particle-only term. (probeMassClosure only
% ever audited particles; HET/PHO export ~2e-6 of supply, so their clean
% residuals say nothing about export-proportional errors.)
names = [m.Particles.Name, m.Liquids.Name];
kPart = numel(m.Particles);
out = struct("tag", "massclosure_" + opts.Variant, "days", opts.Days, ...
    "NCells", opts.NCells, "species", names);
fprintf("MASS CLOSURE ISOLATION variant=%s (N=%d, %g d, q=%g):\n", ...
    opts.Variant, opts.NCells, opts.Days, q);
fprintf("%6s %12s %12s %12s %12s %10s\n", "sp", "dM(eps)", "supply", "export", "residual", "res/supply");
for k = 1:numel(names)
    nm = names(k);
    cin = infl(k);
    if k <= kPart
        total = C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1}; % N x frames
    else
        total = C{nm,"Enclosed"}{1} + C{nm,"Flowing"}{1};   % liquids have no matrix phase
    end
    M = sum(eps.*total, 1)*dz;                                                % eps-weighted stock
    cout = C{nm,"Flowing"}{1}(end, :);
    supplyCurve = q*cin*(ts - ts(1));
    exportCurve = q*cumtrapz(ts, cout(:));
    residCurve = (M(:) - M(1)) - (supplyCurve - exportCurve);
    resid = residCurve(end);
    supply = supplyCurve(end);
    out.(nm) = struct("dM", M(end) - M(1), "supply", supply, "export", exportCurve(end), ...
        "residual", resid, "residualOverSupply", resid/max(supply, 1e-300), ...
        "stock", M, "cout", cout, "residCurve", residCurve);
    fprintf("%6s %12.5g %12.5g %12.5g %12.3e %10.2e\n", nm, M(end) - M(1), ...
        supply, exportCurve(end), resid, resid/max(supply,1e-300));
end
% PAT leak timing: quartile checkpoints of the residual curve
rc = out.PAT.residCurve;
idx = round(linspace(1, numel(ts), 5));
fprintf("PAT R(t) at t = "); fprintf("%.2f ", ts(idx)); fprintf("d:\n  ");
fprintf("%12.4e ", rc(idx)); fprintf("\n");
out.ts = ts; out.eps = eps; out.z = z; out.SolverOptions = r.SolverOptions;
save(fullfile(S, "massclosure_iso_" + opts.Variant + "_n" + opts.NCells + ...
    "_" + opts.Days + "d_f" + opts.NFrames + ".mat"), "out", "-v7.3");
end
