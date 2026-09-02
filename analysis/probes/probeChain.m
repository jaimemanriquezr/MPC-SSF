function out = probeChain(tag, opts)
% PROBECHAIN  Run a long simulation in <=10 d legs, saving after every leg.
%
% Why chained: a 20 d run at N=500 is ~4 h and a 1e-6 arm is ~13 h. Run in one
% block, a wall-clock kill or a mid-run abort yields nothing. Run in legs, each
% leg writes its own result file, so a long job produces partial results and can
% be resumed by resubmitting -- completed legs are detected and skipped.
%
% Chaining is only legitimate if the snapshot carries the whole state; see
% probeRestart.m, which differences a chained run against a continuous one. It
% reports NO zeroed fields and ~1e-4 relative agreement, the residue being the
% adaptive stepper re-ramping at each leg boundary (seed it with ContinuationDt
% to suppress most of that).
%
% Saves per leg:  chain_<tag>_leg<k>.mat   with `rec` (metrics + trajectories +
% effluent + profile history for that leg) and `snap` (the restart snapshot).
arguments
    tag (1,1) string
    opts.Zeta1 (1,1) double = 0.27
    % Cohesion strength: lambda = Zeta0*phi_b*(1-phi_b), so Zeta0 scales the WHOLE
    % cohesive flux and hence the effective diffusivity D = lambda*Psi''. The three
    % presets disagree by four orders of magnitude -- modelPathogen 1e0,
    % pathogenModel 1e2 (the manuscript value, results.tex:79), modelLund 1e6 --
    % so it must be recorded with every result, not inferred from the preset.
    opts.Zeta0 (1,1) double = 100
    opts.DetachScale (1,1) double = 1.0
    opts.NCells (1,1) double = 500
    opts.Days (1,1) double = 20
    opts.LegDays (1,1) double = 10
    opts.MaxDt (1,1) double = 3e-6
    opts.Kappa (1,1) double = 1e-6
    opts.CohesionBC (1,1) string = "neumann"
    opts.Scheme (1,1) string = "shin"
    opts.FramesPerLeg (1,1) double = 41
    % Multiplies the diel light curve; 0 = covered/dark filter.
    opts.LightScale (1,1) double = 1.0
    % Influent [HET PHO POM PAT O2 IC NH4 HPO4 DOM] in kg/m3. Default is the
    % manuscript Table B.1; the "field" set (2026-08-25) is HET 1.5e-4, PHO 5e-4,
    % HPO4 5e-6 (Campos2006b Fig. 1(a); Chan2018 counts at 0.3 um^3/cell).
    opts.Influent (1,9) double = [2.68e-3, 1.00e-2, 0.0, 0.0, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]
    % Phototroph respiration rate (RWQM1 form), 0 = off. Since 2026-08-26 the
    % preset carries Campos k_ra (0.276/d) inside "Phototroph death" and no
    % respiration reaction; > 0 here re-adds RWQM1 (10) on top, for comparison only.
    opts.Respiration (1,1) double = 0.0
    % Heterotroph DOM half-saturation (kg/m3); NaN = preset value (4e-3, Wolf K_S,H,SS).
    opts.KDOM (1,1) double = NaN
    % Sand light-attenuation coefficient (1/m); preset 1500.
    opts.EtaSand (1,1) double = 1500
    % Multiplies every liquid's flowing->enclosed TransportRate (preset 600/d,
    % DOM 300/d). At 1 the enclosed phase runs 20-500x below the flowing phase
    % (kd_1e4_lit, 2026-08-26) and heterotroph growth is transfer-limited.
    opts.TransferScale (1,1) double = 1.0
    % "constant" (TransportRate x TransferScale) or "film" (D/L_f^2, see SolverOptions).
    opts.TransferForm (1,1) string = "constant"
    % Phototroph HPO4 half-saturation (kg P/m3); NaN = preset 2e-5 (Reichert2001).
    % Campos2006 Table 3 ksp: 1e-6 .. 5e-5, avg 2.55e-5.
    opts.KHPO4 (1,1) double = NaN
    % Detachment velocity dependence: "sqrt" (manuscript, 0.14*sqrt(|v|/18)) or
    % "linear" (0.14*|v|/18; same value at v = 18 m/d, grows with pore closure).
    opts.DetachForm (1,1) string = "sqrt"
    % Phototroph growth and loss rates (1/d at 20 C); NaN = preset (2.0 Reichert; 0.276 Campos k_ra).
    opts.MuPHO (1,1) double = NaN
    opts.DPHO (1,1) double = NaN
    % Heterotroph growth rate (1/d at 20 C); NaN = preset 2.0 (Reichert k_gro,H).
    opts.MuHET (1,1) double = NaN
    % Roughness-layer thickness delta (m): the lit zone above z = 0 with bare-sand
    % attachment (eps ramps 1 -> 0.4 over it). SandFilter default 5e-3.
    opts.Delta (1,1) double = 5e-3
    % Filter temperature (C). Every E1-E9 chain ran at 19; 3 = manuscript winter.
    opts.Temperature (1,1) double = 19
    % Diel light forcing: "chain" = the Table B.1 curve every chain run used;
    % "winter" = the manuscript's winter forcing (probeSeason90 / fig:seasons-light).
    % Biofilm self-shading nu_P [m^2/kg]; NaN = preset (52 since b1fcd03).
    % 0.094 reproduces the pre-b1fcd03 Gallegos2000 coefficient.
    opts.Attenuation (1,1) double = NaN
    opts.LightForm (1,1) string {mustBeMember(opts.LightForm, ["chain", "winter"])} = "chain"
end
here = fileparts(mfilename("fullpath")); W = fileparts(fileparts(here));
addpath(genpath(fullfile(W,"src"))); addpath(fullfile(W,"analysis")); addpath(here);
S = fullfile(here, "data", "chain"); if ~isfolder(S), mkdir(S); end

if opts.LightForm == "winter"
    lightFn = @(t) opts.LightScale*max(0.5*(sin(2*pi*(t - 0.2)) + 1) - 0.4, 0);
else
    lightFn = @(t) opts.LightScale*0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50);
end
f = SandFilter(Temperature=opts.Temperature, LightAttenuationCoeffSand=opts.EtaSand, SandRoughness=opts.Delta, ...
    LightIrradiation=lightFn);
f = f.addGridPoints(opts.NCells);
mpArgs = {"PhototrophRespiration", opts.Respiration, "RespirationForm", "reichert", ...
          "NormalizedLight", true};
if ~isnan(opts.Attenuation), mpArgs = [mpArgs, {"Attenuation", opts.Attenuation}]; end
mp = pathogenModel(mpArgs{:});
% modelLund only retires the 1 % dark-growth floor inside its respiration
% branches; with Respiration = 0 it would come back. The floor is an artefact term
% (decision 2026-08-18-phototroph-respiration-rate), so it is zeroed unconditionally.
rx = mp.Reactions; rx([rx.Name] == "Phototroph growth").MinimumLightFactor = 0.0;
if ~isnan(opts.KDOM)
    iH = find([rx.Name] == "Heterotroph growth"); H = rx(iH).HalfSaturationConstants;
    H("DOM") = opts.KDOM; rx(iH).HalfSaturationConstants = H;
end
if ~isnan(opts.KHPO4)
    iP = find([rx.Name] == "Phototroph growth"); H = rx(iP).HalfSaturationConstants;
    H("HPO4") = opts.KHPO4; rx(iP).HalfSaturationConstants = H;
end
if ~isnan(opts.MuPHO), rx([rx.Name] == "Phototroph growth").NominalRate = opts.MuPHO; end
if ~isnan(opts.MuHET), rx([rx.Name] == "Heterotroph growth").NominalRate = opts.MuHET; end
if ~isnan(opts.DPHO),  rx([rx.Name] == "Phototroph death").NominalRate  = opts.DPHO;  end
mp.Reactions = rx;
if opts.TransferScale ~= 1
    cs = mp.Components;
    for q = 1:numel(cs)
        if isa(cs(q), "Liquid"), cs(q).TransportRate = opts.TransferScale*cs(q).TransportRate; end
    end
    mp.Components = cs;
end
m = Model(mp.Components, mp.Reactions, Kappa=opts.Kappa, ...
    Zeta0=opts.Zeta0, Zeta1=opts.Zeta1, ...
    DetachmentFunction=detachmentFunction(opts.DetachForm, opts.DetachScale), ...
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, ...
    OsmosisRate=mp.OsmosisRate);
infl = opts.Influent;
common = {"InflowConcentrations", infl, "TimeStep", "adaptive", ...
    "AdaptiveMaxDt", opts.MaxDt, "ImplicitOsmosis", true, ...
    "CohesionBC", opts.CohesionBC, "CohesionScheme", opts.Scheme, "TransferForm", opts.TransferForm, ...
    "FrameNumber", opts.FramesPerLeg, "Quiet", true};

fprintf("chain %s: zeta0=%g zeta1=%g detach=%g N=%d BC=%s MaxDt=%g light=%g PHOin=%g HETin=%g resp=%g KDOM=%g etaSand=%g transfer=x%g form=%s KHPO4=%g detach=%s muPHO=%g muHET=%g dPHO=%g delta=%g | to t = %g d in legs of %g d\n", ...
    tag, opts.Zeta0, opts.Zeta1, opts.DetachScale, opts.NCells, opts.CohesionBC, opts.MaxDt, ...
    opts.LightScale, opts.Influent(2), opts.Influent(1), opts.Respiration, opts.KDOM, opts.EtaSand, opts.TransferScale, opts.TransferForm, opts.KHPO4, opts.DetachForm, opts.MuPHO, opts.MuHET, opts.DPHO, opts.Delta, opts.Days, opts.LegDays);

% Loop on TIME REACHED, not on a precomputed leg count. The earlier version used
% nLeg = ceil(Days/LegDays) and skipped any leg whose file existed, which is only
% correct when every leg is exactly LegDays long. Extending a 3 d run (one short
% leg) to 10 d computed nLeg = 1, found leg 1 on disk, skipped it and returned
% having run nothing -- silently, reporting tEnd = 3.
s = State(f, m); t = 0; L = 0;
while t < opts.Days - 1e-9
    L = L + 1;
    fn = fullfile(S, sprintf("chain_%s_leg%d.mat", tag, L));
    if isfile(fn)                              % resume: skip completed legs
        % Files written before 2026-08-25 hold only rec+snap; newer ones also
        % hold the full Results. Only rec and snap are needed to resume.
        D = load(fn, "snap", "rec");
        % Resuming into a run configured differently would silently mix two
        % physics settings under one tag. Refuse rather than guess.
        for chk = ["zeta0","zeta1","detach","NCells","Kappa","lightScale","influent","respiration","KDOM","etaSand","transferScale","transferForm","KHPO4","detachForm","muPHO","dPHO","deltaOpt","muHET"]
            % Legs written before LightScale existed (2026-08-26) ran at full light.
            if chk == "lightScale" && ~isfield(D.rec, "lightScale"), was = 1.0;
            elseif chk == "respiration" && ~isfield(D.rec, "respiration"), was = 0.1;   % pre-2026-08-25 legs ran RWQM1 at 0.1
            elseif chk == "KDOM" && ~isfield(D.rec, "KDOM"), was = NaN;
            elseif chk == "etaSand" && ~isfield(D.rec, "etaSand"), was = 1500;
            elseif chk == "transferScale" && ~isfield(D.rec, "transferScale"), was = 1.0;
            elseif chk == "transferForm" && ~isfield(D.rec, "transferForm"), was = "constant";
            elseif chk == "KHPO4" && ~isfield(D.rec, "KHPO4"), was = NaN;
            elseif chk == "detachForm" && ~isfield(D.rec, "detachForm"), was = "sqrt";
            elseif chk == "muPHO" && ~isfield(D.rec, "muPHO"), was = NaN;
            elseif chk == "dPHO" && ~isfield(D.rec, "dPHO"), was = NaN;
            elseif chk == "deltaOpt" && ~isfield(D.rec, "deltaOpt"), was = 5e-3;
            elseif chk == "muHET" && ~isfield(D.rec, "muHET"), was = NaN;
            else, was = D.rec.(chk); end
            now = struct("zeta0",opts.Zeta0, "zeta1",opts.Zeta1, "detach",opts.DetachScale, ...
                         "NCells",opts.NCells, "Kappa",opts.Kappa, "lightScale",opts.LightScale, ...
                         "influent",opts.Influent, "respiration",opts.Respiration, ...
                         "KDOM",opts.KDOM, "etaSand",opts.EtaSand, "transferScale",opts.TransferScale, "transferForm",opts.TransferForm, "KHPO4",opts.KHPO4, "detachForm",opts.DetachForm, "muPHO",opts.MuPHO, "dPHO",opts.DPHO, "deltaOpt",opts.Delta, "muHET",opts.MuHET).(chk);
            if ~isequaln(was, now)
                error("probeChain:parameterMismatch", ...
                    "%s leg %d on disk has %s = %s but this call asks for %s. " + ...
                    "Use a different tag rather than resuming into it.", tag, L, chk, mat2str(was), mat2str(now));
            end
        end
        snap = D.snap; t = D.rec.tEnd;
        s = rehomeState(f, m, snap, t);
        fprintf("  leg %d already done (t = %g d) -- skipped\n", L, t);
        continue
    end
    legLen = min(opts.LegDays, opts.Days - t);
    dt0 = 1e-8; if L > 1, dt0 = opts.MaxDt; end   % seed near the working dt
    t0 = tic;
    r = simulate(s, common{:}, SimulationTime=legLen, AdaptiveInitialDt=dt0);
    if r.Flag ~= "OK"
        % Keep the trajectory up to the abort under a name no resume will pick
        % up. Frames past the abort are zeros; Python callers trim on ts > 0.
        % checkRunFlag then errors as before, so no normal leg file is written.
        aborted = struct("flag", r.Flag, "tFinal", r.TimeFinal, "leg", L, ...
            "tStart", t, "wall", toc(t0));
        results = r; results_py = flattenForPython(r);
        fnAbort = fullfile(S, sprintf("chain_%s_leg%d_ABORTED.mat", tag, L));
        save(fnAbort, "results", "results_py", "aborted", "-v7.3");
        fprintf("  leg %d  %s at t = %g d -- partial trajectory saved to %s\n", L, r.Flag, r.TimeFinal, fnAbort);
    end
    checkRunFlag(r, sprintf("%s leg %d (t = %g..%g d)", tag, L, t, t+legLen));
    wall = toc(t0);
    results = r;              % saved whole -- see the note on the save() below
    rec = legMetrics(r, f, m, t, opts, tag, L, wall);
    snap = stateFromFrame(r, numel(r.Frames.Time));
    if ~isempty(snap.allZero)
        error("probeChain:zeroedField", "snapshot zeroed: %s", strjoin(snap.allZero, ", "));
    end
    % Save the WHOLE Results object, not just derived metrics.
    %
    % legMetrics keeps a hand-picked subset (phi_b, the bed/supernatant
    % integrals, effluent). Storing only that discards Frames.Concentrations --
    % every species profile except phi_b -- and results.SolverOptions, the
    % provenance object. It also stores phiT, which is itself derived from
    % Frames, so the file held a lossy projection of the thing that was thrown
    % away. Any new question (the O2/IC profiles needed to show the
    % substrate-starvation front directly, for instance) then costs a re-run.
    %
    % ~2-4 MB per leg compressed against 0.26 MB, which is nothing beside the
    % hours of cosmos time that produced it. `rec` is kept alongside for
    % convenience, but `results` is the source of truth and rec is regenerable
    % from it via legMetrics.
    %
    % -v7.3 (HDF5): -v7 caps individual arrays at 2 GB, which N=1000 with many
    % frames could approach, and -v7.3 also allows partial loading.
    % `results` is the MATLAB source of truth: the whole Results object, with
    % Frames.Concentrations (every species profile) and SolverOptions. It
    % deserialises only where @Results/@SolverOptions are on the path, and it
    % is NOT readable from Python -- MATLAB objects serialise as MCOS
    % references, which is why rec.flag came back as an MCOS blob in earlier
    % scipy analysis.
    %
    % `results_py` is the same run flattened to plain arrays and char, so the
    % Python plotting/analysis path can read it straight out of the HDF5 file.
    % It is DERIVED -- if the two ever disagree, `results` is correct.
    results_py = flattenForPython(r);
    save(fn, "results", "results_py", "rec", "snap", "-v7.3");
    fprintf("  leg %d  %s  %.0f s  t=%g..%g d  bedInt %.5f  supFrac %.1f%%  effO2 %.5f  -> %s\n", ...
        L, r.Flag, wall, t, t+legLen, rec.bedInt(end), 100*rec.supFrac, rec.effO2(end), fn);
    t = t + legLen;
    s = rehomeState(f, m, snap, t);
end
out = struct("tag", tag, "legs", L, "tEnd", t);
end

function rec = legMetrics(r, f, ~, ~, opts, tag, L, wall)   % model + tStart unused:
% r.Model is used, and simulate returns ABSOLUTE frame times (see rec.ts below)
z = f.GridPoints.Centers; dz = f.GridSize; n0 = f.GridZero;
delta = f.SandRoughness;
eps = computePorosity(f, z);
C = r.Frames.Concentrations;
dL = mean([r.Model.Liquids.Density]); dP = mean([r.Model.Particles.Density]);
P = C{"Water","Enclosed"}{1}/dL;
for nm = [r.Model.Particles.Name], P = P + (C{nm,"Matrix"}{1} + C{nm,"Enclosed"}{1})/dP; end
for nm = [r.Model.Liquids.Name],   P = P + C{nm,"Enclosed"}{1}/dL; end
bed = z >= -delta; sup = z < -delta;
p = P(:, end); w = eps(:).*p;
rec.tag = tag; rec.leg = L; rec.wall = wall; rec.flag = r.Flag;
rec.zeta1 = opts.Zeta1; rec.zeta0 = opts.Zeta0; rec.detach = opts.DetachScale; rec.NCells = opts.NCells;
rec.MaxDt = opts.MaxDt; rec.CohesionBC = opts.CohesionBC; rec.Kappa = opts.Kappa;
rec.lightScale = opts.LightScale;
rec.z = z; rec.eps = eps; rec.dz = dz; rec.delta = delta; rec.n0 = n0;
% simulate() returns ABSOLUTE frame times: rehomeState sets State.Time = t0 and
% the loop runs `while t < timeStart + simulationTime`. Adding tStart again
% double-counts, which put leg 2 of a 20 d chain at t = 20..30.
rec.ts = r.Frames.Time(:).';   rec.tEnd = rec.ts(end);
rec.phiT = P;                                   % full profile history
rec.bedInt = sum(eps(bed).*P(bed, :), 1)*dz;
rec.supInt = sum(eps(sup).*P(sup, :), 1)*dz;
rec.bedPeakT = max(P(bed, :), [], 1);
rec.supPeakT = max(P(sup, :), [], 1);
rec.supFrac = sum(w(sup))/max(sum(w), realmin);
rec.Sigma = sum(w(sup))*dz;
% EFFLUENT (deepest cell, flowing phase) for every species, plus O2 broken out.
lN = [r.Model.Liquids.Name]; pN = [r.Model.Particles.Name];
rec.effNames = [lN, pN];
rec.effluent = zeros(numel(rec.effNames), numel(rec.ts));
for q = 1:numel(lN), a = C{lN(q), "Flowing"}{1}; rec.effluent(q,:) = a(end,:); end
for q = 1:numel(pN), a = C{pN(q), "Flowing"}{1}; rec.effluent(numel(lN)+q,:) = a(end,:); end
rec.effO2 = rec.effluent(lN == "O2", :);
rec.influent = opts.Influent;
rec.respiration = opts.Respiration;
rec.KDOM = opts.KDOM; rec.etaSand = opts.EtaSand; rec.transferScale = opts.TransferScale; rec.transferForm = opts.TransferForm; rec.KHPO4 = opts.KHPO4; rec.detachForm = opts.DetachForm; rec.muPHO = opts.MuPHO; rec.muHET = opts.MuHET; rec.dPHO = opts.DPHO; rec.deltaOpt = opts.Delta;
end

function P = flattenForPython(r)
% FLATTENFORPYTHON  Plain-array view of a Results object, readable by h5py.
%
% MATLAB objects and tables serialise as MCOS references that neither h5py nor
% scipy can follow, so every field here is a numeric array or char. Species
% profiles are stacked into one 3-D array with a matching name list, since HDF5
% has no table equivalent.
P.time      = r.Frames.Time(:).';
P.flag      = char(r.Flag);
C           = r.Frames.Concentrations;
pN          = [r.Model.Particles.Name];
lN          = [r.Model.Liquids.Name];
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
P.water     = C{"Water", "Enclosed"}{1};
P.velocityBiofilm = r.Frames.Velocity.Biofilm;
P.z         = r.SandFilter.GridPoints.Centers(:).';
P.porosity  = computePorosity(r.SandFilter, r.SandFilter.GridPoints.Centers(:)).';
P.dz        = r.SandFilter.GridSize;
P.n0        = r.SandFilter.GridZero;
P.delta     = r.SandFilter.SandRoughness;
% SolverOptions, flattened -- the provenance must survive into Python too.
o = r.SolverOptions;
P.opt_CohesionScheme  = char(o.CohesionScheme);
P.opt_CohesionBC      = char(o.CohesionBoundaryConditions);
P.opt_UpwindedVelocity = double(o.UpwindedVelocity);
P.opt_ImplicitOsmosis  = double(o.ImplicitOsmosis);
P.opt_ImplicitDispersion = double(o.ImplicitDispersion);
P.opt_ConvexSplitting  = double(o.ConvexSplitting);
P.opt_AdaptiveMaxTimeStep = o.AdaptiveMaxTimeStep;
P.opt_AdaptiveInitialTimeStep = o.AdaptiveInitialTimeStep;
P.opt_CFLFactor        = o.CFLFactor;
P.opt_summary          = char(o.summary());
end

function f = detachmentFunction(form, scale)
% Detachment rate k_det(v_f) in 1/d. Both forms equal scale*0.14 at the nominal
% pore velocity 18 m/d (manuscript Table rhs-parameters); they differ in how the
% rate grows as pores close (v_f = q/(eps(1-phi_b))).
switch form
    case "sqrt",   f = @(v) scale*0.14*sqrt(abs(v)/18);
    case "linear", f = @(v) scale*0.14*abs(v)/18;
    otherwise, error("probeChain:detachForm", "DetachForm must be sqrt or linear, got %s", form);
end
end
