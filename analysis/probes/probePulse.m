function out = probePulse(outTag, opts)
% PROBEPULSE  PAT challenge pulse on a MATURE filter from a probeChain snapshot.
%
%   probePulse("e4lit_nom", Snapshot="chain_fld2x_lit_leg6.mat")
%
% Loads `snap` from analysis/probes/data/chain/<Snapshot>, rebuilds the SAME
% filter + model probeChain used to grow it (working set: zeta0 = 1, zeta1 =
% 0.27, LINEAR detachment x1, constant liquid transfer x10, K_DOM 3e-4,
% K_HPO4 1e-6, eta_sand 1500, delta 5 mm, audited pathogenModel), then runs
% TPost days with a PAT pulse on [PulseT0, PulseT1).
%
% THE PULSE. The snapshots were grown with influent PAT = 0, so the baseline
% influent PAT stays 0 here and the pulse is a clean spike from an empty
% column. Peak concentration
%       C_in,peak = PulseFactor * CRef,   CRef = 5.36e-3 kg/m3
% where CRef is the Manriquez2026 Table B.1 marker influent (results.tex:21) --
% the same nominal logOatCampaign.m:35 uses (`patIn = influent(4)`), and its
% pulse scenario multiplies it by 10 on [0.1, 0.3) d (logOatCampaign.m:37).
% PulseFactor = 10 reproduces that; PulseFactor = 100 is the BigPulse.
% Log removal is L(t) = log10(C_in,peak / c_PAT,out(t)), c_out floored at 1e-30.
%
% MODEL PROVENANCE. pathogenModel is used with LundFlowingInert at its default
% (true). That is NOT a choice made here: probeChain.m:77 builds its model from
% pathogenModel() as well, so the E1-E7 chain snapshots were themselves grown
% with the Lund reactions inert in the flowing suspension. Running the pulse
% any other way would change the host model mid-chain. See
% .claude/decisions/2026-08-27-pathogen-model-audit.md.
%
% TEMPERATURE. Temperature=3 gives theta-corrected kinetics only. The biofilm in
% the snapshot was GROWN AT 19 C; a 3 C arm is a cold-shock of a summer-grown
% filter, not a winter filter. Recorded in `rec.grownAtTemperature`.
%
% Saves analysis/probes/data/pulse/pulse_<outTag>.mat with `rec` (metrics),
% `results` (whole Results object) and `results_py`.
arguments
    outTag (1,1) string
    opts.Snapshot (1,1) string = "chain_fld2x_lit_leg6.mat"
    % --- pulse ---
    opts.CRef (1,1) double = 5.36e-3        % Manriquez Table B.1 PAT influent
    opts.PulseFactor (1,1) double = 10      % 10 = nominal, 100 = BigPulse
    opts.PulseT0 (1,1) double = 0.1
    opts.PulseT1 (1,1) double = 0.3
    opts.TPost (1,1) double = 3.0
    % Hydraulic surge applied FOR THE WHOLE RUN (logOatCampaign "flowstep"
    % multiplies InflowVelocity, which simulate cannot vary in time).
    opts.FlowSurge (1,1) double = 1.0
    opts.Temperature (1,1) double = 19
    % --- host model: must match the chain that produced the snapshot ---
    opts.NCells (1,1) double = 500
    opts.Zeta0 (1,1) double = 1
    opts.Zeta1 (1,1) double = 0.27
    opts.Kappa (1,1) double = 1e-6
    opts.DetachForm (1,1) string = "linear"
    opts.DetachScale (1,1) double = 1.0
    opts.LightScale (1,1) double = 1.0
    opts.KDOM (1,1) double = 3e-4
    opts.KHPO4 (1,1) double = 1e-6
    opts.TransferScale (1,1) double = 10
    opts.EtaSand (1,1) double = 1500
    opts.Delta (1,1) double = 5e-3
    opts.Influent (1,9) double = [3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3]
    % --- marker-rate overrides (fig:pulse-outflow reproduction) ---
    opts.InactivationRate (1,1) double = NaN   % NaN = preset 0.02
    opts.BacterivoryRate (1,1) double = NaN    % NaN = preset 8.0
    % --- solver ---
    opts.MaxDt (1,1) double = 5e-5
    % PAT SandAttachmentFactor; NaN = preset (0, i.e. the marker cannot attach to
    % bare sand). 0.06 is Schijven2013's T4 sticking efficiency (geometric mean).
    opts.SandPathogen (1,1) double = NaN
    % HET/PHO sand attachment factor; NaN = preset (1.0).
    opts.SandBiomass (1,1) double = NaN
    opts.NFrames (1,1) double = 145            % 3 d at 30 min
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
addpath(genpath(fullfile(W, "src"))); addpath(fullfile(W, "analysis")); addpath(here);
S = fullfile(here, "data", "pulse"); if ~isfolder(S), mkdir(S); end

snapFile = fullfile(here, "data", "chain", opts.Snapshot);
assert(isfile(snapFile), "probePulse:noSnapshot", "snapshot not found: %s", snapFile);
D = load(snapFile, "snap", "rec");
snap = D.snap;

% ---- filter and model: probeChain's construction, verbatim ------------------
f = SandFilter(Temperature=opts.Temperature, LightAttenuationCoeffSand=opts.EtaSand, ...
    SandRoughness=opts.Delta, ...
    LightIrradiation=@(t) opts.LightScale*0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(opts.NCells);
% The grid has more rows than NCells (supernatant + bed); compare against the
% actual grid so a snapshot from a different NCells cannot be rehomed silently.
assert(size(snap.Matrix, 1) == numel(f.GridPoints.Centers), "probePulse:gridMismatch", ...
    "snapshot has %d rows, this grid has %d (NCells = %d)", ...
    size(snap.Matrix, 1), numel(f.GridPoints.Centers), opts.NCells);
f.InflowVelocity = opts.FlowSurge*f.InflowVelocity;

mpArgs = {"NormalizedLight", true};
if ~isnan(opts.SandPathogen), mpArgs = [mpArgs, {"SandPathogen", opts.SandPathogen}]; end
if ~isnan(opts.SandBiomass), mpArgs = [mpArgs, {"SandBiomass", opts.SandBiomass}]; end
mp = pathogenModel(mpArgs{:});
rx = mp.Reactions;
rx([rx.Name] == "Phototroph growth").MinimumLightFactor = 0.0;   % as probeChain
iH = find([rx.Name] == "Heterotroph growth"); H = rx(iH).HalfSaturationConstants;
H("DOM") = opts.KDOM; rx(iH).HalfSaturationConstants = H;
iP = find([rx.Name] == "Phototroph growth"); Hp = rx(iP).HalfSaturationConstants;
Hp("HPO4") = opts.KHPO4; rx(iP).HalfSaturationConstants = Hp;
if ~isnan(opts.InactivationRate), rx([rx.Name] == "Inactivation").NominalRate = opts.InactivationRate; end
if ~isnan(opts.BacterivoryRate),  rx([rx.Name] == "Bacterivory").NominalRate  = opts.BacterivoryRate;  end
mp.Reactions = rx;
if opts.TransferScale ~= 1
    cs = mp.Components;
    for q = 1:numel(cs)
        if isa(cs(q), "Liquid"), cs(q).TransportRate = opts.TransferScale*cs(q).TransportRate; end
    end
    mp.Components = cs;
end
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, Zeta0=opts.Zeta0, Zeta1=opts.Zeta1, ...
    DetachmentFunction=detachmentFunction(opts.DetachForm, opts.DetachScale), ...
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, OsmosisRate=mp.OsmosisRate);

% ---- influent: baseline PAT 0, spiked on [PulseT0, PulseT1) -----------------
base = opts.Influent;  base(4) = 0.0;
cPeak = opts.PulseFactor*opts.CRef;
pulsed = base;  pulsed(4) = cPeak;
t0Snap = snap.Time;
inflow = @(t) tern(t >= t0Snap + opts.PulseT0 && t < t0Snap + opts.PulseT1, pulsed, base);

s = rehomeState(f, m, snap, t0Snap);
fprintf("pulse %s: snap=%s t0=%g d  N=%d  x%g pulse (peak %.4g kg/m3) on [%g,%g) d  T=%g C  surge=%g  light=%g  Tpost=%g d\n", ...
    outTag, opts.Snapshot, t0Snap, opts.NCells, opts.PulseFactor, cPeak, ...
    opts.PulseT0, opts.PulseT1, opts.Temperature, opts.FlowSurge, opts.LightScale, opts.TPost);

tic;
r = simulate(s, InflowConcentrations=inflow, SimulationTime=opts.TPost, ...
    TimeStep="adaptive", AdaptiveInitialDt=opts.MaxDt, AdaptiveMaxDt=opts.MaxDt, ...
    ImplicitOsmosis=true, CohesionBC="neumann", CohesionScheme="shin", ...
    TransferForm="constant", FrameNumber=opts.NFrames, Quiet=true);
wall = toc;
if r.Flag ~= "OK"
    % Keep the trajectory up to the abort under a name no normal load will pick
    % up, exactly as probeChain.m does. Without this a CLOGGED pulse produced
    % nothing at all and the run had to be repeated to see even where it died --
    % which is how the x100 arm of job 3549870 was lost. Frames past the abort
    % are zeros; trim on ts > 0. checkRunFlag still errors below, so no normal
    % pulse_<tag>.mat is written from a partial run.
    aborted = struct("flag", r.Flag, "tFinal", r.TimeFinal, "tStart", t0Snap, "wall", wall);
    results = r; results_py = flattenForPython(r); %#ok<NASGU>
    try
        rec = pulseMetrics(r, f, m, opts, outTag, cPeak, t0Snap, wall);   %#ok<NASGU>
    catch metricsErr
        % The zero-padded frames past the abort break interp1 in pulseMetrics
        % (job 3549875_1, 2026-08-27). Keep the trajectory; the metrics are
        % meaningless for a partial run anyway.
        rec = struct("tag", outTag, "flag", r.Flag, "metricsError", metricsErr.message); %#ok<NASGU>
    end
    fnAbort = fullfile(S, "pulse_" + outTag + "_ABORTED.mat");
    save(fnAbort, "results", "results_py", "rec", "aborted", "-v7.3");
    fprintf("  %s at t = %g d -- partial trajectory saved to %s\n", r.Flag, r.TimeFinal, fnAbort);
end
checkRunFlag(r, "probePulse " + outTag);

rec = pulseMetrics(r, f, m, opts, outTag, cPeak, t0Snap, wall);
results = r; results_py = flattenForPython(r); %#ok<NASGU>
fn = fullfile(S, "pulse_" + outTag + ".mat");
save(fn, "results", "results_py", "rec", "-v7.3");
fprintf("  %s  %.0f s  L(1d)=%.2f  L(3d)=%.2f  peak c_out=%.4g  Lmin=%.2f  effO2 %.4g -> %.4g  biomass %.5g -> %.5g\n", ...
    r.Flag, wall, rec.L1d, rec.L3d, rec.peakOut, rec.Lmin, rec.effO2(1), rec.effO2(end), ...
    rec.totalBiomass(1), rec.totalBiomass(end));
fprintf("  -> %s\n", fn);
out = rec;
end

% =============================================================================
function rec = pulseMetrics(r, f, m, opts, outTag, cPeak, t0Snap, wall)
z = f.GridPoints.Centers(:); dz = f.GridSize; delta = f.SandRoughness;
ep = computePorosity(f, z);
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);

% phi_b profile history (probeChain's legMetrics convention)
P = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], P = P + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   P = P + C{nm,"Enclosed"}{1}/dL; end

sup  = z < -delta;              % supernatant
rough = z >= -delta & z < 0;    % roughness layer (bare-sand attachment zone)
bed  = z >= 0;                  % sand bed

rec.tag = outTag; rec.wall = wall; rec.flag = string(r.Flag);
rec.snapshot = opts.Snapshot; rec.tSnapshot = t0Snap;
rec.grownAtTemperature = 19;    % every E1-E7 chain ran at 19 C
rec.Temperature = opts.Temperature; rec.FlowSurge = opts.FlowSurge;
rec.PulseFactor = opts.PulseFactor; rec.CRef = opts.CRef; rec.cPeak = cPeak;
rec.PulseT0 = opts.PulseT0; rec.PulseT1 = opts.PulseT1;
rec.LightScale = opts.LightScale; rec.NCells = opts.NCells;
rec.Zeta0 = opts.Zeta0; rec.Zeta1 = opts.Zeta1; rec.DetachForm = opts.DetachForm;
rec.KDOM = opts.KDOM; rec.KHPO4 = opts.KHPO4; rec.TransferScale = opts.TransferScale;
rec.EtaSand = opts.EtaSand; rec.Delta = opts.Delta; rec.Influent = opts.Influent;
rec.InactivationRate = opts.InactivationRate; rec.BacterivoryRate = opts.BacterivoryRate;
rec.MaxDt = opts.MaxDt;
rec.z = z; rec.eps = ep; rec.dz = dz; rec.delta = delta; rec.n0 = f.GridZero;

ts = r.Frames.Time(:).';
rec.ts = ts;  rec.tRel = ts - t0Snap;      % time since the snapshot
rec.phiT = P;
rec.phib = P(:, end);

% --- effluent ------------------------------------------------------------
lN = [r.Model.Liquids.Name]; pN = [r.Model.Particles.Name];
rec.effNames = [lN, pN];
rec.effluent = zeros(numel(rec.effNames), numel(ts));
for q = 1:numel(lN), a = C{lN(q), "Flowing"}{1}; rec.effluent(q,:) = a(end,:); end
for q = 1:numel(pN), a = C{pN(q), "Flowing"}{1}; rec.effluent(numel(lN)+q,:) = a(end,:); end
rec.effO2  = rec.effluent(lN == "O2", :);
rec.effPAT = rec.effluent(numel(lN) + find(pN == "PAT"), :);

% --- log removal ---------------------------------------------------------
cOut = max(rec.effPAT, 1e-30);
rec.L = log10(cPeak./cOut);
rec.peakOut = max(rec.effPAT);
rec.Lmin = min(rec.L);                       % worst-case removal = at breakthrough
rec.tPeakOut = rec.tRel(find(rec.effPAT == rec.peakOut, 1));
% A CLOGGED run leaves zero-padded frames past the abort (duplicate ts), which
% interp1 rejects; interpolate on the strictly increasing prefix only.
ok = [true, diff(rec.tRel) > 0];
rec.L1d = interp1(rec.tRel(ok), rec.L(ok), min(1.0, max(rec.tRel(ok))));
rec.L3d = rec.L(end);

% --- PAT depth profiles (matrix + enclosed + flowing, global) ------------
rec.patProfile = C{"PAT","Matrix"}{1} + C{"PAT","Enclosed"}{1} + C{"PAT","Flowing"}{1};

% --- biomass, by region and total ---------------------------------------
% eps-weighted integral of phi_b, the convention closed by probeMassClosure.
w = ep.*P;
rec.supInt   = sum(w(sup, :),   1)*dz;
rec.roughInt = sum(w(rough, :), 1)*dz;
rec.bedInt   = sum(w(bed, :),   1)*dz;
rec.totalBiomass = rec.supInt + rec.roughInt + rec.bedInt;
rec.supFrac = rec.supInt./max(rec.totalBiomass, realmin);
% top 2 cm of the sand bed (Campos2002 sampling depth), particulate carbon only
top2 = z >= 0 & z < 0.02;
Xc = zeros(numel(z), numel(ts));
for nm = [r.Model.Particles.Name], Xc = Xc + C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1}; end
rec.top2cmMass = sum(ep(top2).*Xc(top2, :), 1)*dz;      % kg/m2 over the top 2 cm
rec.top2cmDepth = 0.02;
% phototroph mass above the sand (lit/dark contrast)
Xpho = C{"PHO","Matrix"}{1} + C{"PHO","Enclosed"}{1};
rec.phoAboveSand = sum(ep(z < 0).*Xpho(z < 0, :), 1)*dz;
rec.phoBed = sum(ep(bed).*Xpho(bed, :), 1)*dz;
% profile maximum location, for the scoring rubric
[~, imax] = max(P(:, end));
rec.zPhibMax = z(imax);
rec.SolverOptions = r.SolverOptions.summary();
end

function f = detachmentFunction(form, scale)
switch form
    case "sqrt",   f = @(v) scale*0.14*sqrt(abs(v)/18);
    case "linear", f = @(v) scale*0.14*abs(v)/18;
    otherwise, error("probePulse:detachForm", "DetachForm must be sqrt or linear, got %s", form);
end
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end

function P = flattenForPython(r)
P.time = r.Frames.Time(:).';
P.flag = char(r.Flag);
C = r.Frames.Concentrations;
pN = [r.Model.Particles.Name];  lN = [r.Model.Liquids.Name];
P.particleNames = char(join(pN, "|"));
P.liquidNames   = char(join(lN, "|"));
P.phaseNames    = char("Matrix|Enclosed|Flowing");
nz = numel(r.SandFilter.GridPoints.Centers);  nt = numel(P.time);
P.particles = zeros(nz, nt, numel(pN), 3);
for q = 1:numel(pN)
    P.particles(:,:,q,1) = C{pN(q), "Matrix"}{1};
    P.particles(:,:,q,2) = C{pN(q), "Enclosed"}{1};
    P.particles(:,:,q,3) = C{pN(q), "Flowing"}{1};
end
P.liquids = zeros(nz, nt, numel(lN), 2);
for q = 1:numel(lN)
    P.liquids(:,:,q,1) = C{lN(q), "Enclosed"}{1};
    P.liquids(:,:,q,2) = C{lN(q), "Flowing"}{1};
end
P.water = C{"Water", "Enclosed"}{1};
P.velocityBiofilm = r.Frames.Velocity.Biofilm;
P.z = r.SandFilter.GridPoints.Centers(:).';
P.porosity = computePorosity(r.SandFilter, r.SandFilter.GridPoints.Centers(:)).';
P.dz = r.SandFilter.GridSize;
P.n0 = r.SandFilter.GridZero;
P.delta = r.SandFilter.SandRoughness;
P.opt_summary = char(r.SolverOptions.summary());
end
