function outdir = logOatCampaign(disturbance, options)
% LOGOATCAMPAIGN  Logarithmic OAT sensitivity campaign, MATLAB anchor run.
%
%   outdir = logOatCampaign("pulse")
%
% Mirror of julia/analysis/log_oat_sensitivity.jl main() with the 2026-08-18
% design decisions applied:
%   * `velocity` EXCLUDED (fixed operating condition, 7.2 m/d; see
%     .claude/decisions/2026-08-18-velocity-fixed-operating-condition.md).
%   * `beta_porosity` perturbed WITHIN [0.95, 0.99]: runs at beta = 0.98 and
%     0.95 (gap 0.02 and 0.05 from nominal 0.99); the log-sensitivity uses the
%     secant denominator ln(0.05/0.02) = ln 2.5 instead of 2 ln 2.
% Everything else: parameter theta_i x2 and x1/2 around nominal, removal curve
% L(t) = log10(Cref/c_PAT_out(t)), measures I_rms/I_max/D_min/I_min/asymmetry
% over the post-disturbance window. Scenarios: "startup" (clean IC, ripening),
% "pulse" (mature + 10x PAT spike on [0.1, 0.3)), "flowstep" (mature + x1.4
% hydraulic surge). See analysis/PARAMETERS.md for the full parameter table.
arguments
    disturbance (1,1) string {mustBeMember(disturbance, ["startup", "pulse", "flowstep"])}
    options.TMature (1,1) {mustBeNumeric} = 3.0;
    options.TPost (1,1) {mustBeNumeric} = 1.5;
    options.NCells (1,1) {mustBeNumeric} = 30;
    options.NFrames (1,1) {mustBeNumeric} = 60;
    options.MaxDt (1,1) {mustBeNumeric} = 3e-6;
    options.ParamFilter string = string.empty;   % run a subset (smoke tests)
    options.SkipParams string = string.empty;    % exclude parameters by name
    % Baseline-model options, forwarded to pathogenModel (defaults reproduce
    % the 2026-08-18 anchor campaign's published-model baseline).
    options.Respiration (1,1) {mustBeNumeric} = 0.0;
    options.PGExcess (1,1) logical = false;
    options.OutTag (1,1) string = "";            % suffix for the results dir
    % ---- 2026-08-27: baseline on the AUDITED preset and the E1-E7 working set.
    % Preset = "campaign" reproduces the 2026-08-18 anchor exactly (nothing
    % below is applied). Preset = "workingset" switches the baseline to the E4
    % configuration -- audited pathogenModel, zeta0 = 1, LINEAR detachment,
    % constant liquid transfer x10, K_DOM 3e-4, K_HPO4,PHO 1e-6, eta_sand 1500,
    % delta 5 mm, T = 19 C, 2x field influent -- and moves every kinetic nominal
    % in campaignParams to the value that configuration actually runs at.
    % See .claude/decisions/2026-08-27-pathogen-model-audit.md.
    options.Preset (1,1) string ...
        {mustBeMember(options.Preset, ["campaign", "workingset"])} = "campaign";
    % Path to a probeChain leg file. "" = build the mature state here with a
    % TMature pre-run (the 2026-08-18 behaviour). A snapshot makes the anchor a
    % 60 d mature filter instead of a 3 d one, and removes the pre-run cost.
    % NCells must match the snapshot's grid.
    options.Snapshot (1,1) string = "";
    % Working-set knobs; only read when Preset = "workingset" (or when set
    % explicitly). Mirror probeChain.m's options of the same name.
    options.Zeta0 (1,1) {mustBeNumeric} = NaN;         % NaN -> preset default
    options.Zeta1 (1,1) {mustBeNumeric} = NaN;
    options.DetachForm (1,1) string = "";              % "" -> preset default
    options.KDOM (1,1) {mustBeNumeric} = NaN;
    options.KHPO4 (1,1) {mustBeNumeric} = NaN;
    options.TransferScale (1,1) {mustBeNumeric} = NaN;
    options.EtaSand (1,1) {mustBeNumeric} = NaN;
    options.Delta (1,1) {mustBeNumeric} = NaN;
    options.Temperature (1,1) {mustBeNumeric} = NaN;
    options.LightScale (1,1) {mustBeNumeric} = 1.0;
    options.Influent (1,9) {mustBeNumeric} = zeros(1,9);   % all-zero -> preset default
end

cfg = resolveConfig(options);
influent = cfg.influent;
patIn = cfg.patIn;
pulseFactor = 10.0;  pulseT0 = 0.1;  pulseT1 = 0.3;  flowSurge = 1.4;
% Absolute start time: a chain snapshot restarts at its own t (60 d for E4), and
% simulate returns ABSOLUTE frame times, so every window below is offset by it.
cfg.tStart = 0.0;
if cfg.snapshot ~= ""
    Dt = load(cfg.snapshot, "snap");  cfg.tStart = Dt.snap.Time;
end
pulseT0 = cfg.tStart + pulseT0;  pulseT1 = cfg.tStart + pulseT1;
td = cfg.tStart;
if disturbance == "pulse", td = pulseT0; end

params = campaignParams(cfg);
if ~isempty(options.ParamFilter)
    params = params(ismember([params.name], options.ParamFilter));
    % A typo in ParamFilter used to leave `params` empty, run only the baseline
    % and exit 0 -- an array task that produced a measures.csv with no rows.
    kept = string.empty; if ~isempty(params), kept = [params.name]; end
    missing = setdiff(options.ParamFilter, kept);
    assert(isempty(missing), "logOatCampaign:unknownParam", ...
        "ParamFilter names not in the design: %s", strjoin(missing, ", "));
end
if ~isempty(options.SkipParams)
    params = params(~ismember([params.name], options.SkipParams));
end

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, "results", "oat", "log_oat_" + disturbance + options.OutTag);
if ~isfolder(outdir), mkdir(outdir); end

fprintf("Log-OAT (%s): %d params, tmature=%g tpost=%g ncells=%d\n", ...
        disturbance, length(params), options.TMature, options.TPost, options.NCells);

% ---- shared mature state (nominal model); startup uses a clean IC ----------
mature = [];
if disturbance ~= "startup"
    if cfg.snapshot ~= ""
        D = load(cfg.snapshot, "snap");
        mature = D.snap;                         % rehomeState-compatible struct
        fprintf("mature state: snapshot %s at t = %g d\n", cfg.snapshot, mature.Time);
    else
        mature = buildMature(cfg, options);
    end
end

% ---- baseline ----------------------------------------------------------------
[t0, L0, flag0, tf0] = runChallenge(mature, [], NaN, disturbance, influent, ...
    patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options, cfg);
Lmin0 = min(L0(t0 >= td));
writematrix([t0(:), L0(:)], fullfile(outdir, "L0.csv"));
fprintf("baseline: flag=%s  Lmean=%.3f  Lmin=%.3f  (t_final=%.3f)\n", ...
        flag0, mean(L0(t0 >= td)), Lmin0, tf0);

% ---- parameter loop ------------------------------------------------------------
header = ["param", "block", "I_rms", "I_max", "D_min", "I_min", "asymmetry", ...
          "flag_plus", "flag_minus", "clog_driver", "source"];
rows = cell(length(params), length(header));
for i = 1:length(params)
    p = params(i);
    [tp, Lp, fp, tfp] = runChallenge(mature, p, p.plus, disturbance, influent, ...
        patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options, cfg);
    [~, Lm, fm, tfm] = runChallenge(mature, p, p.minus, disturbance, influent, ...
        patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options, cfg);

    tvalid = min([tf0, tfp, tfm]);
    win = (t0 >= td) & (t0 <= tvalid + 1e-9);
    isClog = fp ~= "OK" || fm ~= "OK";

    s = (Lp - Lm)/p.denom;
    sw = s(win);
    Irms = sqrt(mean(sw.^2));  Imax = max(abs(sw));
    Lminp = min(Lp(win));  Lminm = min(Lm(win));
    Dmin = max(Lminp - Lmin0, Lminm - Lmin0);
    Imin = (Lminp - Lminm)/p.denom;
    dp = Lp(win) - L0(win);  dm = Lm(win) - L0(win);
    asym = sqrt(sum((dp + dm).^2))/(sqrt(sum((dp - dm).^2)) + 1e-30);

    rows(i, :) = {p.name, p.block, Irms, Imax, Dmin, Imin, asym, fp, fm, isClog, p.source};
    writematrix([tp(:), Lp(:) - L0(:), Lm(:) - L0(:), s(:)], ...
                fullfile(outdir, "curves_" + p.name + ".csv"));
    fprintf("  %-16s [%-9s] I_rms=%.3e  I_max=%.3e  D_min=%+.3e  asym=%.2f  (%s/%s)%s\n", ...
            p.name, p.block, Irms, Imax, Dmin, asym, fp, fm, ternary(isClog, "  CLOG", ""));
end

measures = cell2table(rows, 'VariableNames', cellstr(header));
writetable(measures, fullfile(outdir, "measures.csv"));
fprintf("wrote %s\n", fullfile(outdir, "measures.csv"));
end

% =============================================================================
function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function cfg = resolveConfig(o)
% Collapse Preset + explicit overrides into one config struct. "campaign"
% reproduces the 2026-08-18 anchor bit for bit; "workingset" is E4.
here = fileparts(mfilename('fullpath'));
switch o.Preset
    case "campaign"
        d = struct("zeta0", 1e2, "zeta1", 1e-2, "detachForm", "campaign", ...
                   "kdom", NaN, "khpo4", NaN, "transferScale", 1, ...
                   "etaSand", 1500, "delta", 5e-3, "temperature", 15, ...
                   "influent", [2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4]);
    case "workingset"
        % E4 (EXPERIMENTS.md 2026-08-26): 2x field influent, PAT slot carries
        % the Table B.1 marker so the pulse has a baseline to be scaled from.
        d = struct("zeta0", 1, "zeta1", 0.27, "detachForm", "linear", ...
                   "kdom", 3e-4, "khpo4", 1e-6, "transferScale", 10, ...
                   "etaSand", 1500, "delta", 5e-3, "temperature", 19, ...
                   "influent", [3.0e-4, 1.0e-3, 0, 5.36e-3, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3]);
end
pick = @(v, dv) ternary(isnan(v), dv, v);
cfg.preset        = o.Preset;
cfg.zeta0         = pick(o.Zeta0, d.zeta0);
cfg.zeta1         = pick(o.Zeta1, d.zeta1);
cfg.detachForm    = ternary(o.DetachForm == "", d.detachForm, o.DetachForm);
cfg.kdom          = pick(o.KDOM, d.kdom);
cfg.khpo4         = pick(o.KHPO4, d.khpo4);
cfg.transferScale = pick(o.TransferScale, d.transferScale);
cfg.etaSand       = pick(o.EtaSand, d.etaSand);
cfg.delta         = pick(o.Delta, d.delta);
cfg.temperature   = pick(o.Temperature, d.temperature);
cfg.lightScale    = o.LightScale;
cfg.influent      = ternary(all(o.Influent == 0), d.influent, o.Influent);
cfg.patIn         = cfg.influent(4);
assert(cfg.patIn > 0, "logOatCampaign:noMarker", ...
    "influent PAT is 0; the removal curve L = log10(Cref/c_out) would be undefined");
cfg.respiration   = o.Respiration;
cfg.pgExcess      = o.PGExcess;
if o.Snapshot == ""
    cfg.snapshot = "";
elseif isfile(o.Snapshot)
    cfg.snapshot = o.Snapshot;
else
    cfg.snapshot = fullfile(here, "probes", "data", "chain", o.Snapshot);
    assert(isfile(cfg.snapshot), "logOatCampaign:noSnapshot", "not found: %s", cfg.snapshot);
end
end

function [f, m] = buildBase(cfg, ncells)
% The baseline filter and model. For Preset = "workingset" this is exactly
% probeChain.m's construction with the E4 options.
f = SandFilter(Temperature=cfg.temperature, LightAttenuationCoeffSand=cfg.etaSand, ...
    SandRoughness=cfg.delta, ...
    LightIrradiation=@(t) cfg.lightScale*0.8*max(sin(2*pi*(t - 13/48)) + 31/50, 0)/(1 + 31/50));
f = f.addGridPoints(ncells);
mp = pathogenModel(PhototrophRespiration=cfg.respiration, PGExcess=cfg.pgExcess, ...
                   NormalizedLight=(cfg.preset == "workingset"), ...
                   Zeta0=cfg.zeta0, DetachForm=cfg.detachForm);
if cfg.preset == "workingset"
    rx = mp.Reactions;
    rx([rx.Name] == "Phototroph growth").MinimumLightFactor = 0.0;
    if ~isnan(cfg.kdom)
        i = find([rx.Name] == "Heterotroph growth"); H = rx(i).HalfSaturationConstants;
        H("DOM") = cfg.kdom; rx(i).HalfSaturationConstants = H;
    end
    if ~isnan(cfg.khpo4)
        i = find([rx.Name] == "Phototroph growth"); H = rx(i).HalfSaturationConstants;
        H("HPO4") = cfg.khpo4; rx(i).HalfSaturationConstants = H;
    end
    mp.Reactions = rx;
    if cfg.transferScale ~= 1
        cs = mp.Components;
        for q = 1:numel(cs)
            if isa(cs(q), "Liquid"), cs(q).TransportRate = cfg.transferScale*cs(q).TransportRate; end
        end
        mp.Components = cs;
    end
end
m = Model(mp.Components, mp.Reactions, Kappa=mp.CohesionSubModel.Kappa, ...
    Zeta0=cfg.zeta0, Zeta1=cfg.zeta1, DetachmentFunction=mp.DetachmentFunction, ...
    WaterDensity=mp.WaterDensity, BiofilmPorosity=mp.BiofilmPorosity, ...
    OsmosisRate=mp.OsmosisRate);
end

function mature = buildMature(cfg, options)
[f, m] = buildBase(cfg, options.NCells);
r = simulate(State(f, m), InflowConcentrations=cfg.influent, ...
    SimulationTime=options.TMature, TimeStep="adaptive", AdaptiveInitialDt=1e-8, ...
    AdaptiveMaxDt=options.MaxDt, FrameNumber=5, ImplicitOsmosis=true, Quiet=true);
assert(r.Flag == "OK", "maturation failed: flag=" + r.Flag);
mature = finalState(r);
fprintf("mature state: t=%.3f  phib_max=%.4f\n", r.TimeFinal, max(mature.phiBiofilm));
end

function mature = finalState(r)
% Extract the last frame of a Results object as re-homeable arrays.
C = r.Frames.Concentrations;
pNames = [r.Model.Particles.Name];  lNames = [r.Model.Liquids.Name];
kP = length(pNames);  kL = length(lNames);
N = size(C{pNames(1), "Matrix"}{1}, 1);
mature.Matrix = zeros(N, kP);  mature.EnclosedParticles = zeros(N, kP);
mature.FlowingParticles = zeros(N, kP);
mature.EnclosedLiquids = zeros(N, kL);  mature.FlowingLiquids = zeros(N, kL);
for j = 1:kP
    mature.Matrix(:, j) = C{pNames(j), "Matrix"}{1}(:, end);
    mature.EnclosedParticles(:, j) = C{pNames(j), "Enclosed"}{1}(:, end);
    mature.FlowingParticles(:, j) = C{pNames(j), "Flowing"}{1}(:, end);
end
for j = 1:kL
    mature.EnclosedLiquids(:, j) = C{lNames(j), "Enclosed"}{1}(:, end);
    mature.FlowingLiquids(:, j) = C{lNames(j), "Flowing"}{1}(:, end);
end
densityL = mean([r.Model.Liquids.Density]);
mature.EnclosedWaterVolume = C{"Water", "Enclosed"}{1}(:, end)/densityL;
mature.VelocityBiofilm = r.Frames.Velocity.Biofilm(:, end);
densityP = mean([r.Model.Particles.Density]);
mature.phiBiofilm = mature.EnclosedWaterVolume ...
    + sum(mature.Matrix + mature.EnclosedParticles, 2)/densityP ...
    + sum(mature.EnclosedLiquids, 2)/densityL;
end

function [ts, L, flag, tFinal] = runChallenge(mature, p, value, disturbance, ...
    influent, patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options, cfg)
% One run: perturbed filter/model, scenario disturbance, removal curve.
[f, m] = buildBase(cfg, options.NCells);
patMult = 1.0;
if ~isempty(p)
    if p.name == "influent_PAT"
        patMult = value;
    else
        [f, m] = p.apply(f, m, value);
    end
end
if disturbance == "flowstep"
    f.InflowVelocity = flowSurge*f.InflowVelocity;
end
s = State(f, m);
if disturbance ~= "startup"
    s.GlobalConcentration.Matrix = mature.Matrix;
    s.GlobalConcentration.EnclosedParticles = mature.EnclosedParticles;
    s.GlobalConcentration.FlowingParticles = mature.FlowingParticles;
    s.GlobalConcentration.EnclosedLiquids = mature.EnclosedLiquids;
    s.GlobalConcentration.FlowingLiquids = mature.FlowingLiquids;
    s.EnclosedWaterVolume = mature.EnclosedWaterVolume;
    s.Velocity.Biofilm = mature.VelocityBiofilm;
    % A chain snapshot carries an absolute time; keeping it puts the diel light
    % curve at the right phase (finalState has no Time field, so guard).
    if isfield(mature, "Time"), s.Time = mature.Time; end
end

base = influent;  base(4) = patIn*patMult;
if disturbance == "pulse"
    cref = pulseFactor*patIn*patMult;
    pulsed = base;  pulsed(4) = cref;
    inflow = @(t) ternary(t >= pulseT0 && t < pulseT1, pulsed, base);
else
    cref = patIn*patMult;
    inflow = base;
end

r = simulate(s, InflowConcentrations=inflow, SimulationTime=options.TPost, ...
    TimeStep="adaptive", AdaptiveInitialDt=1e-8, AdaptiveMaxDt=options.MaxDt, ...
    FrameNumber=options.NFrames, ImplicitOsmosis=true, Quiet=true);

centers = f.GridPoints.Centers;
assert(centers(end) > 0, "grid orientation: expected effluent at the last cell");
cOut = max(r.Frames.Concentrations{"PAT", "Flowing"}{1}(end, :), 1e-30);
ts = r.Frames.Time(:);
L = log10(max(cref, eps)./cOut(:));
flag = string(r.Flag);
tFinal = r.TimeFinal;
end

% =============================================================================
function params = campaignParams(cfg)
% The 26-parameter design (velocity excluded). Mirrors PARAMS in
% log_oat_sensitivity.jl; nominals cited there and in PARAMETERS.md.
% Fields: name, block, plus/minus (the two perturbed VALUES), denom (the
% log-sensitivity denominator), source, apply(f, m, v) -> [f, m].
%
% NOMINALS ARE NOT DECORATIVE. Each perturbed run sets the parameter to the
% ABSOLUTE value 2*nom or nom/2, so a nominal that disagrees with what the
% baseline model actually runs at makes the campaign perturb around the wrong
% point AND silently change the baseline. The "campaign" column below is the
% 2026-08-18 published-model set; the "workingset" column is what the audited
% preset on the E4 configuration really carries
% (.claude/decisions/2026-08-25-kinetics-wolf-reichert.md,
%  .claude/decisions/2026-08-25-half-saturations-corrected.md,
%  .claude/decisions/2026-08-27-pathogen-model-audit.md).
ln2 = log(2);
P = @(name, block, nom, src, apply) struct("name", string(name), "block", string(block), ...
    "plus", 2*nom, "minus", 0.5*nom, "denom", 2*ln2, "source", string(src), "apply", apply);
ws = cfg.preset == "workingset";
pick = @(campaign, workingset) ternary(ws, workingset, campaign);

params = [ ...
    P("temperature",   "forcing",   cfg.temperature, "Campos2006/Manriquez seasonal",   {@(f,m,v) setFilter(f, m, "Temperature", v)}), ...
    P("influent_PAT",  "forcing",   1.0,    "Manriquez B.1 PAT (challenge multiplier)", {@(f,m,v) deal(f, m)}), ...
    P("dispersivity",  "transport", 0.012,  "Schijven2013 T1 dispersivity",             {@(f,m,v) setParticles(f, m, "Dispersivity", v)}), ...
    P("transport_P",   "transport", 5.47,   "Lund phase transfer",                      {@(f,m,v) setParticles(f, m, "TransportRate", v)}), ...
    P("attach_sand",   "transport", 547.0,  "Diehl2025/Lund b_sand",                    {@(f,m,v) setParticles(f, m, "AttachmentSand", v)}), ...
    P("sand_pathogen", "transport", 0.06,   "Schijven2013 T4 sticking alpha (geo-mean)", {@(f,m,v) setPAT(f, m, "SandAttachmentFactor", v)}), ...
    P("mu_HET",        "kinetics",  pick(0.0181, 2.0),  pick("Campos2006 kgmaxa", "Reichert2001 k_gro,H,aer"),   {@(f,m,v) setRate(f, m, "Heterotroph growth", v)}), ...
    P("mu_PHO",        "kinetics",  pick(5.5, 2.0),     pick("Campos2006 kgmaxb", "Reichert2001 k_gro,ALG"),      {@(f,m,v) setRate(f, m, "Phototroph growth", v)}), ...
    P("d_HET",         "kinetics",  pick(2.0, 0.40),    pick("Campos2006 T3 kdb", "Wolf2007 b_ina,H"),            {@(f,m,v) setRate(f, m, "Heterotroph death", v)}), ...
    P("d_PHO",         "kinetics",  pick(0.4, 0.276),   pick("Wolf2007 b_ina,H (slid)", "Campos2006 T3 k_ra"),    {@(f,m,v) setRate(f, m, "Phototroph death", v)}), ...
    P("hydrolysis",    "kinetics",  pick(0.09, 3.0),    pick("Campos2006 T3 kh", "Reichert2001 k_hyd / Wolf k_h"), {@(f,m,v) setRate(f, m, "Hydrolysis", v)}), ...
    P("theta_growth",  "kinetics",  pick(0.047, 0.0725), pick("Campos2006 T2 1.047 (dev)", "Reichert2001 beta_H 0.07 -> 1.0725 (dev)"), {@(f,m,v) setTheta(f, m, ["Heterotroph growth", "Phototroph growth"], v)}), ...
    P("theta_death",   "kinetics",  pick(0.066, 0.08),  pick("Campos2006 T2 1.066 (dev)", "Campos2006 theta_kra 1.08 (dev)"),           {@(f,m,v) setTheta(f, m, ["Heterotroph death", "Phototroph death"], v)}), ...
    P("K_O2_HET",      "kinetics",  pick(3.0e-3, 2.0e-4), "Reichert2001 K_O2,H = 0.2 g/m3 (audited)",             {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "O2", v)}), ...
    P("K_DOM_HET",     "kinetics",  pick(2.0e-4, cfg.kdom), pick("Campos2006 T3 ksCd", "oligotrophic K_DOM (UNCITED, decision pending)"), {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "DOM", v)}), ...
    P("K_HPO4_HET",    "kinetics",  pick(1.4e-8, 2.0e-5), "Reichert2001 K_HPO4,H = 0.02 gP/m3 (audited)",         {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "HPO4", v)}), ...
    P("K_O2_PAT",      "pathogen",  pick(3.0e-3, 2.0e-4), "audited HET row (2026-08-27); manuscript gives none",  {@(f,m,v) setHalfSat(f, m, "MarkerGrowth", "O2", v)}), ...
    P("K_pred",        "pathogen",  2.0e-3, "thesis_model; tab:eco-parameters gives none",                        {@(f,m,v) setHalfSat(f, m, "Bacterivory", "HET", v)}), ...
    P("marker_growth", "pathogen",  0.2,    "Manriquez B.4 mu_PAT",                     {@(f,m,v) setRate(f, m, "MarkerGrowth", v)}), ...
    P("inactivation",  "pathogen",  0.02,   "Schijven2013 mul,mus; Manriquez B.4",      {@(f,m,v) setRate(f, m, "Inactivation", v)}), ...
    P("bacterivory",   "pathogen",  8.0,    "Manriquez B.4 p_PAT",                      {@(f,m,v) setRate(f, m, "Bacterivory", v)}), ...
    P("zeta_0",        "biofilm",   cfg.zeta0, "Cahn-Hilliard cohesion",                {@(f,m,v) setModelArg(f, m, "Zeta0", v)}), ...
    P("kappa",         "biofilm",   1.0e-6, "Cahn-Hilliard interfacial width (sub-grid)", {@(f,m,v) setModelArg(f, m, "Kappa", v)}), ...
    P("zeta_1",        "biofilm",   cfg.zeta1, "Cahn-Hilliard stable-fraction parameter", {@(f,m,v) setModelArg(f, m, "Zeta1", v)}), ...
    P("detach_scale",  "biofilm",   1.0,    "detachment-law scale",                     {@(f,m,v) setDetachScale(f, m, v, cfg.detachForm)}), ...
    P("light_att_water", "light",   0.32,   "Lund supernatant optical depth",           {@(f,m,v) setFilter(f, m, "LightAttenuationCoeffWater", v)}), ...
    P("light_att_sand",  "light",   cfg.etaSand, "Lund sand-bed optical depth",         {@(f,m,v) setFilter(f, m, "LightAttenuationCoeffSand", v)}), ...
    P("attenuation_P",   "light",   0.094,  "Lund biofilm self-shading (all particles)", {@(f,m,v) setParticles(f, m, "Attenuation", v)}), ...
    ];

% beta_porosity: constrained scheme (2026-08-18) -- beta in [0.95, 0.99] only.
% Runs at gap = 0.02 (beta 0.98) and gap = 0.05 (beta 0.95); secant over ln(gap).
beta = P("beta_porosity", "biofilm", 0.0, "Lund beta=0.99; constrained to [0.95, 0.99]", ...
         {@(f,m,gap) setModelArg(f, m, "BiofilmPorosity", 1 - gap)});
beta.plus = 0.02;  beta.minus = 0.05;  beta.denom = log(0.05/0.02);
params = [params, beta];
end

% ---- appliers (all rebuild through the Model constructor so the cohesion
% potential regenerates when Zeta1 changes; see log_oat_sensitivity.jl notes) --
function [f, m] = setFilter(f, m, field, v)
f.(field) = v;
end

function [f, m] = setParticles(f, m, field, v)
c = m.Components;
for i = 1:length(c)
    if isa(c(i), "Particle"), c(i).(field) = v; end
end
m = rebuild(m, "Components", c);
end

function [f, m] = setPAT(f, m, field, v)
c = m.Components;
for i = 1:length(c)
    if isa(c(i), "Particle") && c(i).Name == "PAT", c(i).(field) = v; end
end
m = rebuild(m, "Components", c);
end

function [f, m] = setRate(f, m, name, v)
r = m.Reactions;
for i = 1:length(r)
    if r(i).Name == name, r(i).NominalRate = v; end
end
m = rebuild(m, "Reactions", r);
end

function [f, m] = setTheta(f, m, names, dev)
r = m.Reactions;
for i = 1:length(r)
    if ismember(r(i).Name, names), r(i).TemperatureCorrectionFactor = 1 + dev; end
end
m = rebuild(m, "Reactions", r);
end

function [f, m] = setHalfSat(f, m, rxName, key, v)
r = m.Reactions;
for i = 1:length(r)
    if r(i).Name == rxName
        K = r(i).HalfSaturationConstants;  K(key) = v;
        r(i).HalfSaturationConstants = K;
    end
end
m = rebuild(m, "Reactions", r);
end

function [f, m] = setModelArg(f, m, arg, v)
m = rebuild(m, arg, v);
end

function [f, m] = setDetachScale(f, m, v, form)
% Scales whichever law the BASELINE uses; hardcoding the campaign law here made
% detach_scale silently switch the model's detachment form on the working set.
switch form
    case "campaign", g = @(vel) v*1.4e-5*sqrt(abs(vel)/18);
    case "sqrt",     g = @(vel) v*0.14*sqrt(abs(vel)/18);
    case "linear",   g = @(vel) v*0.14*abs(vel)/18;
end
m = rebuild(m, "DetachmentFunction", g);
end

function m2 = rebuild(m, arg, v)
% Reconstruct the Model with one argument overridden. Constructing (rather than
% assigning fields) keeps PotentialGradient consistent with Zeta1.
args = struct("Components", {m.Components}, "Reactions", {m.Reactions}, ...
    "Kappa", m.CohesionSubModel.Kappa, "Zeta0", m.CohesionSubModel.Zeta0, ...
    "Zeta1", m.CohesionSubModel.Zeta1, "DetachmentFunction", m.DetachmentFunction, ...
    "WaterDensity", m.WaterDensity, "BiofilmPorosity", m.BiofilmPorosity, ...
    "OsmosisRate", m.OsmosisRate);
args.(arg) = v;
nv = namedargs2cell(rmfield(args, ["Components", "Reactions"]));
m2 = Model(args.Components, args.Reactions, nv{:});
end
