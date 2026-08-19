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
end

influent = [2.68e-3, 1.00e-2, 0.0, 5.36e-3, 9.10e-3, 6.23e-3, 2.00e-5, 0.0, 1.75e-4];
patIn = influent(4);
pulseFactor = 10.0;  pulseT0 = 0.1;  pulseT1 = 0.3;  flowSurge = 1.4;
td = 0.0;
if disturbance == "pulse", td = pulseT0; end

params = campaignParams();
if ~isempty(options.ParamFilter)
    params = params(ismember([params.name], options.ParamFilter));
end
if ~isempty(options.SkipParams)
    params = params(~ismember([params.name], options.SkipParams));
end

here = fileparts(mfilename('fullpath'));
outdir = fullfile(here, "results", "log_oat_" + disturbance + options.OutTag);
if ~isfolder(outdir), mkdir(outdir); end

fprintf("Log-OAT (%s): %d params, tmature=%g tpost=%g ncells=%d\n", ...
        disturbance, length(params), options.TMature, options.TPost, options.NCells);

% ---- shared mature state (nominal model); startup uses a clean IC ----------
mature = [];
if disturbance ~= "startup"
    mature = buildMature(options.NCells, options.TMature, options.MaxDt, influent, options.Respiration, options.PGExcess);
end

% ---- baseline ----------------------------------------------------------------
[t0, L0, flag0, tf0] = runChallenge(mature, [], NaN, disturbance, influent, ...
    patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options);
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
        patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options);
    [~, Lm, fm, tfm] = runChallenge(mature, p, p.minus, disturbance, influent, ...
        patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options);

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

function mature = buildMature(ncells, tmature, maxDt, influent, resp, pgx)
f = SandFilter();
f = f.addGridPoints(ncells);
r = simulate(State(f, pathogenModel(PhototrophRespiration=resp, PGExcess=pgx)), ...
    InflowConcentrations=influent, ...
    SimulationTime=tmature, TimeStep="adaptive", AdaptiveInitialDt=1e-8, ...
    AdaptiveMaxDt=maxDt, FrameNumber=5, ImplicitOsmosis=true, Quiet=true);
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
    influent, patIn, pulseFactor, pulseT0, pulseT1, flowSurge, options)
% One run: perturbed filter/model, scenario disturbance, removal curve.
f = SandFilter();
f = f.addGridPoints(options.NCells);
m = pathogenModel(PhototrophRespiration=options.Respiration, PGExcess=options.PGExcess);
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
function params = campaignParams()
% The 26-parameter design (velocity excluded). Mirrors PARAMS in
% log_oat_sensitivity.jl; nominals cited there and in PARAMETERS.md.
% Fields: name, block, plus/minus (the two perturbed VALUES), denom (the
% log-sensitivity denominator), source, apply(f, m, v) -> [f, m].
ln2 = log(2);
P = @(name, block, nom, src, apply) struct("name", string(name), "block", string(block), ...
    "plus", 2*nom, "minus", 0.5*nom, "denom", 2*ln2, "source", string(src), "apply", apply);

params = [ ...
    P("temperature",   "forcing",   15.0,   "Campos2006/Manriquez seasonal",            {@(f,m,v) setFilter(f, m, "Temperature", v)}), ...
    P("influent_PAT",  "forcing",   1.0,    "Manriquez B.1 PAT (challenge multiplier)", {@(f,m,v) deal(f, m)}), ...
    P("dispersivity",  "transport", 0.012,  "Schijven2013 T1 dispersivity",             {@(f,m,v) setParticles(f, m, "Dispersivity", v)}), ...
    P("transport_P",   "transport", 5.47,   "Lund phase transfer",                      {@(f,m,v) setParticles(f, m, "TransportRate", v)}), ...
    P("attach_sand",   "transport", 547.0,  "Diehl2025/Lund b_sand",                    {@(f,m,v) setParticles(f, m, "AttachmentSand", v)}), ...
    P("sand_pathogen", "transport", 0.06,   "Schijven2013 T4 sticking alpha (geo-mean)", {@(f,m,v) setPAT(f, m, "SandAttachmentFactor", v)}), ...
    P("mu_HET",        "kinetics",  0.0181, "Campos2006 kgmaxa",                        {@(f,m,v) setRate(f, m, "Heterotroph growth", v)}), ...
    P("mu_PHO",        "kinetics",  5.5,    "Campos2006 kgmaxb",                        {@(f,m,v) setRate(f, m, "Phototroph growth", v)}), ...
    P("d_HET",         "kinetics",  2.0,    "Campos2006 T3 kdb",                        {@(f,m,v) setRate(f, m, "Heterotroph death", v)}), ...
    P("hydrolysis",    "kinetics",  0.09,   "Campos2006 T3 kh",                         {@(f,m,v) setRate(f, m, "Hydrolysis", v)}), ...
    P("theta_growth",  "kinetics",  0.047,  "Campos2006 T2 theta_growth=1.047 (dev)",   {@(f,m,v) setTheta(f, m, ["Heterotroph growth", "Phototroph growth"], v)}), ...
    P("theta_death",   "kinetics",  0.066,  "Campos2006 T2 theta_death=1.066 (dev)",    {@(f,m,v) setTheta(f, m, ["Heterotroph death", "Phototroph death"], v)}), ...
    P("K_O2_HET",      "kinetics",  3.0e-3, "O2 half-sat (Reichert)",                   {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "O2", v)}), ...
    P("K_DOM_HET",     "kinetics",  2.0e-4, "Campos2006 T3 ksCd",                       {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "DOM", v)}), ...
    P("K_HPO4_HET",    "kinetics",  1.4e-8, "Campos2006 T3 ksp (P-limiting)",           {@(f,m,v) setHalfSat(f, m, "Heterotroph growth", "HPO4", v)}), ...
    P("marker_growth", "pathogen",  0.2,    "Manriquez B.4 mu_PAT",                     {@(f,m,v) setRate(f, m, "MarkerGrowth", v)}), ...
    P("inactivation",  "pathogen",  0.02,   "Schijven2013 mul,mus; Manriquez B.4",      {@(f,m,v) setRate(f, m, "Inactivation", v)}), ...
    P("bacterivory",   "pathogen",  8.0,    "Manriquez B.4 p_PAT",                      {@(f,m,v) setRate(f, m, "Bacterivory", v)}), ...
    P("zeta_0",        "biofilm",   1.0e2,  "Cahn-Hilliard cohesion",                   {@(f,m,v) setModelArg(f, m, "Zeta0", v)}), ...
    P("kappa",         "biofilm",   1.0e-6, "Cahn-Hilliard interfacial width (sub-grid)", {@(f,m,v) setModelArg(f, m, "Kappa", v)}), ...
    P("zeta_1",        "biofilm",   1.0e-2, "Cahn-Hilliard stable-fraction parameter",  {@(f,m,v) setModelArg(f, m, "Zeta1", v)}), ...
    P("detach_scale",  "biofilm",   1.0,    "detachment-law scale",                     {@(f,m,v) setDetachScale(f, m, v)}), ...
    P("light_att_water", "light",   0.32,   "Lund supernatant optical depth",           {@(f,m,v) setFilter(f, m, "LightAttenuationCoeffWater", v)}), ...
    P("light_att_sand",  "light",   1500.0, "Lund sand-bed optical depth",              {@(f,m,v) setFilter(f, m, "LightAttenuationCoeffSand", v)}), ...
    P("attenuation_P",   "light",   0.094,  "Lund biofilm self-shading (all particles)", {@(f,m,v) setParticles(f, m, "Attenuation", v)}), ...
    ];

% beta_porosity: constrained scheme (2026-08-18) -- beta in [0.95, 0.99] only.
% Runs at gap = 0.02 (beta 0.98) and gap = 0.05 (beta 0.95); secant over ln(gap).
beta = P("beta_porosity", "biofilm", NaN, "Lund beta=0.99; constrained to [0.95, 0.99]", ...
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

function [f, m] = setDetachScale(f, m, v)
m = rebuild(m, "DetachmentFunction", @(vel) v*1.4e-5*sqrt(abs(vel)/18));
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
