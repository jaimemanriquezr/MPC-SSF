function [score, metrics, notes] = scoreRun(files, opts)
% SCORERUN  Score one chained run against the literature rubric.
%
%   [score, metrics, notes] = scoreRun("fld2x_lit")
%   [score, metrics, notes] = scoreRun("fld2x_lit", Pair="fld2x_dark")
%   [score, metrics, notes] = scoreRun(["...leg1.mat" "...leg2.mat"], Out="tag")
%   [score, metrics, notes] = scoreRun(D)      % D = load(legfile) or a rec struct
%
% INPUT
%   files  - one of
%              * a tag, e.g. "fld2x_lit": expands to
%                analysis/probes/data/chain/chain_<tag>_leg*.mat, leg-sorted;
%              * a string array of .mat leg paths (sorted by leg number here);
%              * a struct with fields `rec` and (optionally) `results_py`, i.e.
%                the result of load() on one leg file, or an array of them.
%            A `results` (Results object) field is accepted and flattened, so a
%            fresh in-memory run can be scored without saving.
%
% OPTIONS
%   Pair        tag / files / struct for the opposite light arm (gives the
%               lit-dark contrast axis). Default "" = contrast not assessed.
%   Label       name used in the report (default: the tag).
%   Out         if non-empty, write the Markdown table to
%               analysis/results/scores/<Out>.md. "" = do not write.
%   WetFactor   wet:dry mass ratio of the biomass state variable (default 4).
%               f_dry is UNRESOLVED (EXPERIMENTS.md E1 caveat 1); both readings
%               are always reported, this only sets the ratio between them.
%   CarbonPerCOD  kg C per kg of state variable, default 0.374 (= 0.531/1.42 for
%               C5H7O2N, as in analysis/probes/probeCover.m).
%   TopDepth    depth of the Campos sampling layer, default 0.02 m.
%   SteadyDays  window for the steadiness tests, default 10 d.
%
% OUTPUT
%   score    struct: .profile (0-4) .biomass (0-3) .oxygen (0-3) .total (0-10)
%            plus .justification, one line per axis.
%   metrics  every measured number (see the field comments below).
%   notes    the Markdown report as a string (also written if Out is set).
%
% CONVENTIONS (state them with any number taken from here)
%   * Biomass per unit filter area, kg/m2:  int eps(z) * (Matrix + Enclosed) dz
%     over HET+PHO+POM. The porosity weight is required because the stored
%     concentrations are per unit PORE volume: phi_b = sum(att)/rho_P and the
%     bed integral that probeChain reports is int eps*phi_b dz.
%   * Campos-commensurable areal density, ug C per g dry sand:
%         kappa = CarbonPerCOD / ((1 - eps0) * rho_quartz) * 1e6
%               = 0.374 / (0.6 * 2650) * 1e6 = 235.2  (ug C/g) per (kg/m3)
%     applied to the eps-weighted layer mean at the final frame. This is the
%     "dry" reading (state variable read as dry organic matter); the "wet"
%     reading divides it by WetFactor.
%     NOTE: the E4 table in EXPERIMENTS.md is ~7 % higher than this because it
%     used 0.4 kg C/kg instead of 0.374; see analysis/results/scores/README-ish
%     note in the E4 validation block of the progress log.
%   * H/H0 is the Kozeny-Carman bed head loss of analysis/headlossKozenyCarman.py,
%     evaluated on z >= delta only (the roughness ramp has no clean-bed value).
%   * O2 in mg/L = 1000 x the stored kg/m3. "Consumed" = influent minus the mean
%     effluent over the final leg.
arguments
    files
    opts.Pair = ""
    opts.Label (1,1) string = ""
    opts.Out (1,1) string = ""
    % Which reading the Campos biomass point is scored on. "wet" (default) is
    % the reading EXPERIMENTS.md E4 compares with Campos; "dry" reads the state
    % variable as dry organic matter. The score is CONDITIONAL on this choice —
    % it is reported in every table, and both readings are always printed.
    opts.Convention (1,1) string {mustBeMember(opts.Convention, ["wet" "dry"])} = "wet"
    % Ratio between the two readings. NOT DERIVED: it is the factor quoted in
    % EXPERIMENTS.md E1 caveat 1 ("yields are ~4x too strong per kg of state
    % variable and the ug C/g figures move by the same factor") and used in the
    % E1/E4 tables (140 dry <-> 35 wet). f_dry is unresolved; treat 4 as a
    % placeholder with no independent provenance.
    opts.WetFactor (1,1) double = 4.0
    opts.CarbonPerCOD (1,1) double = 0.374
    opts.TopDepth (1,1) double = 0.02
    opts.SteadyDays (1,1) double = 10
end

A = readChain(files);
metrics = chainMetrics(A, opts);
metrics.label = opts.Label; if metrics.label == "", metrics.label = A.tag; end

% ---- lit/dark contrast --------------------------------------------------
metrics.hasPair = false;
if ~(isstring(opts.Pair) && isscalar(opts.Pair) && opts.Pair == "")
    B = readChain(opts.Pair);
    mb = chainMetrics(B, opts);
    metrics.pair = mb;
    metrics.hasPair = true;
    metrics.ratioPhoAbove = metrics.phoAbove / max(mb.phoAbove, realmin);
    metrics.ratioTop2     = metrics.top2Dry / max(mb.top2Dry, realmin);
    metrics.ratioBed      = metrics.massBed / max(mb.massBed, realmin);
    metrics.ratioBedInt   = metrics.bedInt  / max(mb.bedInt,  realmin);
end

metrics.convention = opts.Convention;
metrics.wetFactor = opts.WetFactor;
[score, metrics] = rubric(metrics, opts);
notes = report(metrics, score, opts);

if opts.Out ~= ""
    here = fileparts(mfilename("fullpath"));
    d = fullfile(here, "results", "scores");
    if ~isfolder(d), mkdir(d); end
    fid = fopen(fullfile(d, opts.Out + ".md"), "w");
    fprintf(fid, "%s", notes); fclose(fid);
    fprintf("wrote %s\n", fullfile(d, opts.Out + ".md"));
end
end

% ======================== loading ========================================
function A = readChain(files)
% Returns A: tag, ts, series (per frame), and the final-frame fields.
if isstruct(files)
    L = files(:).';
    for k = 1:numel(L), L(k) = normaliseLeg(L(k)); end
else
    fs = resolveFiles(files);
    L = struct("rec", cell(1, numel(fs)), "py", cell(1, numel(fs)));
    for k = 1:numel(fs)
        D = load(fs(k), "rec", "results_py");
        L(k).rec = normaliseRec(D.rec); L(k).py = D.results_py;
    end
end
A.tag = string(L(1).rec.tag);
A.legs = numel(L);

py = L(1).py;
A.z = py.z(:); A.dz = py.dz; A.eps = py.porosity(:); A.delta = py.delta;
A.influent = L(1).rec.influent;
A.lightScale = L(1).rec.lightScale;
A.pN = string(strsplit(py.particleNames, '|'));
A.lN = string(strsplit(py.liquidNames, '|'));

ts = []; mSup = []; mRough = []; mBed = []; effO2 = []; phiLast = [];
prevEnd = -inf;
for k = 1:numel(L)
    rec = L(k).rec; py = L(k).py;
    t = rec.ts(:).';
    % Guard inherited from analysis/headlossKozenyCarman.py: a few legacy legs
    % were written with tStart added twice.
    if k > 1 && abs(t(1) - 2*prevEnd) < abs(t(1) - prevEnd), t = t - prevEnd; end
    keep = 1:numel(t); if k > 1, keep = 2:numel(t); end
    z = py.z(:); eps = py.porosity(:); dz = py.dz; delta = py.delta;
    tot = attached(py, A.pN, ["HET" "PHO" "POM"]);
    sup = z < -delta; rough = z >= -delta & z < 0; bed = z >= 0;
    w = eps .* tot;
    ts    = [ts,    t(keep)];                        %#ok<AGROW>
    mSup  = [mSup,  sum(w(sup, keep), 1)*dz];        %#ok<AGROW>
    mRough= [mRough,sum(w(rough, keep), 1)*dz];      %#ok<AGROW>
    mBed  = [mBed,  sum(w(bed, keep), 1)*dz];        %#ok<AGROW>
    effO2 = [effO2, rec.effO2(keep)];                %#ok<AGROW>
    prevEnd = t(end);
    if k == numel(L)
        A.phiHist = rec.phiT;  A.phiHistT = t;
        A.lastPy = py; A.lastRec = rec;
    end
    phiLast = rec.phiT(:, end);
end
A.ts = ts; A.mSup = mSup; A.mRough = mRough; A.mBed = mBed;
A.mTot = mSup + mRough + mBed;
A.effO2 = effO2; A.phi = phiLast;
end

function rec = normaliseRec(rec)
% probeChain and probePulse write DIFFERENT rec schemas. probePulse capitalises
% Influent/LightScale, has no tEnd, and reports bed/roughness/supernatant
% integrals as time series. Map its names onto probeChain's so one scorer reads
% both; anything already in probeChain's schema passes through untouched.
if ~isfield(rec, "influent") && isfield(rec, "Influent"), rec.influent = rec.Influent; end
if ~isfield(rec, "lightScale") && isfield(rec, "LightScale"), rec.lightScale = rec.LightScale; end
if ~isfield(rec, "tEnd"), rec.tEnd = rec.ts(end); end
if ~isfield(rec, "supInt"), rec.supInt = 0; end
if ~isfield(rec, "tag"), rec.tag = "unnamed"; end
end

function s = normaliseLeg(s)
if isfield(s, "rec"), s.rec = normaliseRec(s.rec); end
if ~isfield(s, "py")
    if isfield(s, "results_py"), s.py = s.results_py;
    elseif isfield(s, "results"), s.py = flattenLite(s.results);
    else, error("scoreRun:noFields", "struct input needs results_py (or results)");
    end
end
if isfield(s, "rec") == 0, error("scoreRun:noRec", "struct input needs a rec field"); end
end

function P = flattenLite(r)
C = r.Frames.Concentrations;
pN = [r.Model.Particles.Name]; lN = [r.Model.Liquids.Name];
P.particleNames = char(join(pN, "|")); P.liquidNames = char(join(lN, "|"));
nz = numel(r.SandFilter.GridPoints.Centers); nt = numel(r.Frames.Time);
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
P.z = r.SandFilter.GridPoints.Centers(:).';
P.porosity = computePorosity(r.SandFilter, r.SandFilter.GridPoints.Centers(:)).';
P.dz = r.SandFilter.GridSize; P.delta = r.SandFilter.SandRoughness;
P.n0 = r.SandFilter.GridZero; P.time = r.Frames.Time(:).';
end

function M = attached(py, pN, names)
M = 0;
for nm = names
    j = find(pN == nm);
    if isempty(j), continue, end
    M = M + py.particles(:,:,j,1) + py.particles(:,:,j,2);
end
end

function fs = resolveFiles(files)
files = string(files);
if isscalar(files) && ~endsWith(files, ".mat")
    here = fileparts(mfilename("fullpath"));
    d = fullfile(here, "probes", "data", "chain");
    g = dir(fullfile(d, sprintf("chain_%s_leg*.mat", files)));
    g = g(~contains(string({g.name}), "ABORTED"));
    if isempty(g), error("scoreRun:noFiles", "no legs for tag %s in %s", files, d); end
    fs = string(fullfile(d, {g.name}));
else
    fs = files(:).';
end
% Sort by leg number where there is one. Files that are not chain legs
% (probePulse output, for instance) keep the order they were given in.
n = 1:numel(fs);
for k = 1:numel(fs)
    tok = regexp(fs(k), "leg(\d+)\.mat$", "tokens", "once");
    if ~isempty(tok), n(k) = str2double(tok(1)); end
end
[~, o] = sort(n); fs = fs(o);
end

% ======================== metrics ========================================
function M = chainMetrics(A, opts)
z = A.z; eps = A.eps; dz = A.dz; delta = A.delta;
M.tag = A.tag; M.legs = A.legs; M.tEnd = A.ts(end); M.ts = A.ts;
M.lightScale = A.lightScale; M.influent = A.influent;

% --- areal biomass, kg/m2, split by compartment ---
M.massSup = A.mSup(end); M.massRough = A.mRough(end); M.massBed = A.mBed(end);
M.massTot = M.massSup + M.massRough + M.massBed;
M.massSeries = A.mTot;
M.supFrac = M.massSup / max(M.massTot, realmin);
% probeChain's own bed/supernatant integrals of eps*phi_b (units m), kept so a
% score can be checked directly against the outcome tables in EXPERIMENTS.md.
M.bedInt = A.lastRec.bedInt(end); M.supInt = A.lastRec.supInt(end);

% --- Campos layer, ug C per g dry sand ---
kappa = opts.CarbonPerCOD / (0.6*2650) * 1e6;
M.kappa = kappa;
py = A.lastPy;
live = attached(py, A.pN, ["HET" "PHO"]);
top2 = z >= 0 & z <= opts.TopDepth;
M.top2Dry = kappa * mean(eps(top2) .* live(top2, end));
M.top2Wet = M.top2Dry / opts.WetFactor;
M.top2Cells = sum(top2);

% --- phototrophs above the sand (roughness layer), kg/m2 ---
pho = attached(py, A.pN, "PHO");
rough = z >= -delta & z < 0;
M.phoAbove = sum(eps(rough) .* pho(rough, end)) * dz;

% --- final-frame phi_b profile and its qualitative flags ---
phi = A.phi; M.phi = phi; M.z = z;
win = z >= -delta & z <= opts.TopDepth;
[~, iPk] = max(phi);
M.zPeak = z(iPk); M.phiPeak = phi(iPk);
M.peakAtSurface = win(iPk);
% no secondary maximum / upward migration below the peak: the largest rise of
% phi_b going downwards, relative to the peak
below = iPk:numel(phi);
d = diff(phi(below));
M.maxRiseBelow = max([0; d]) / max(M.phiPeak, realmin);
M.monotoneBelow = M.maxRiseBelow < 0.01;
% steadiness of the profile over the last SteadyDays
M.steadyWindow = opts.SteadyDays;
tt = A.phiHistT; P = A.phiHist;
M.profileChange = NaN; M.longEnough = (A.ts(end) - A.ts(1)) >= opts.SteadyDays;
if tt(end) - tt(1) >= opts.SteadyDays - 1e-9
    [~, j0] = min(abs(tt - (tt(end) - opts.SteadyDays)));
    M.profileChange = sum(abs(P(:, end) - P(:, j0))) / max(sum(abs(P(:, end))), realmin);
elseif M.longEnough
    % window spans more than one leg: fall back on the areal-mass change
    [~, j0] = min(abs(A.ts - (A.ts(end) - opts.SteadyDays)));
    M.profileChange = abs(A.mTot(end) - A.mTot(j0)) / max(A.mTot(end), realmin);
end
M.profileSteady = M.longEnough && M.profileChange < 0.01;

% --- monotone growth of the whole-column biomass ---
dm = diff(A.mTot);
M.maxDropFrac = max([0, -min([dm, 0])]) / max(max(A.mTot), realmin);
M.massMonotone = M.maxDropFrac < 0.01;

% --- head loss (Kozeny-Carman, z >= delta) ---
sel = z >= delta;
e0 = eps(sel); e = max(e0 .* (1 - phi(sel)), 1e-9);
M.HH0 = mean(((e0./e).^3) .* (((1 - e)./(1 - e0)).^2));

% --- effluent O2, mg/L, over the final leg ---
% Window: the LAST SIMULATED DAY, not the last leg. A leg that starts from a
% clean column (leg 1) begins with the deepest cell at O2 = 0 and stays near
% zero for the first residence time, so a leg-wide minimum would report the
% flush-through, not the filter. One day also contains exactly one diel cycle,
% which is what "the effluent minimum" means operationally.
win = A.ts >= A.ts(end) - 1;
o2 = A.effO2(win) * 1000;
M.o2Window = min(1, A.ts(end) - A.ts(1));
M.o2Mean = mean(o2); M.o2Min = min(o2); M.o2Max = max(o2);
M.o2In = A.influent(5) * 1000;
M.o2Consumed = M.o2In - M.o2Mean;
end

% ======================== rubric =========================================
function [S, M] = rubric(M, ~)
j = strings(0);

% ---- Profile, 0-4 ----
p = 0;
p1 = M.peakAtSurface;   p = p + p1;
p2 = M.monotoneBelow;   p = p + p2;
p3 = M.supFrac < 0.01;  p = p + p3;
p4 = M.profileSteady;   p = p + p4;
S.profileParts = [p1 p2 p3 p4];
j(1) = sprintf("Profile %d/4: peak at z = %.4f m (%s), largest rise below peak %.2f %% of peak (%s), supernatant %.3f %% of total (%s), %s.", ...
    p, M.zPeak, tf(p1, "in the roughness layer / top 2 cm", "BELOW the sand surface"), ...
    100*M.maxRiseBelow, tf(p2, "monotone", "secondary maximum"), 100*M.supFrac, ...
    tf(p3, "< 1 %", ">= 1 %"), steadyPhrase(M));

% ---- Total biomass, 0-3 ----
b = 0;
b1 = M.massMonotone; b = b + b1;
% Campos2002 60 ug C/g at d97, factor 2 => [30, 120]. Scored on ONE convention
% (opts.Convention), not on whichever passes: awarding the point to whichever
% reading happened to land in the window made the axis unfalsifiable.
inWin = @(x) x >= 30 && x <= 120;
M.top2Scored = M.top2Wet; if M.convention == "dry", M.top2Scored = M.top2Dry; end
b2 = inWin(M.top2Scored); b = b + b2;
if M.hasPair, b3 = M.ratioPhoAbove >= 1.5; else, b3 = false; end
b = b + b3;
S.biomassParts = [b1 b2 b3];
j(2) = sprintf("Biomass %d/3: whole-column %.4f kg/m2 (sup %.5f + roughness %.5f + bed %.5f), " + ...
    "largest drop %.2f %% (%s, and this axis cannot fail while a filter is simply ripening); " + ...
    "top-2 cm %.0f ug C/g dry / %.0f wet, SCORED ON THE %s READING vs Campos2002 60 at d97 (%s, " + ...
    "conditional on f_dry = %g); %s.", ...
    b, M.massTot, M.massSup, M.massRough, M.massBed, 100*M.maxDropFrac, ...
    tf(b1, "monotone", "not monotone"), M.top2Dry, M.top2Wet, upper(M.convention), ...
    tf(b2, "inside the factor-2 window", "outside the window"), M.wetFactor, ...
    contrastPhrase(M, b3));

% ---- Oxygen, 0-3 ----
c = M.o2Consumed;
if c >= 2 && c <= 5,  o1 = 2;
elseif (c >= 1 && c < 2) || (c > 5 && c <= 6), o1 = 1;
else, o1 = 0; end
o2p = M.o2Min >= 3;
o = o1 + o2p;
S.oxygenParts = [o1 o2p];
j(3) = sprintf("O2 %d/3: consumed %.2f mg/L (%s Elemo2024 2-5), effluent mean %.2f, minimum %.2f mg/L (%s 3 mg/L floor).", ...
    o, c, tf(o1 == 2, "inside", tf(o1 == 1, "in the 1-2 / 5-6 margin of", "outside")), ...
    M.o2Mean, M.o2Min, tf(o2p, "above the", "BELOW the"));

S.profile = p; S.biomass = b; S.oxygen = o; S.total = p + b + o;
S.justification = j;
end

function s = tf(c, a, b), if c, s = string(a); else, s = string(b); end, end

function s = contrastPhrase(M, awarded)
% Built outside the sprintf because MATLAB evaluates both branches of a
% ternary-style call, and M.ratioPhoAbove does not exist without a pair.
if ~M.hasPair
    s = "lit/dark contrast NOT ASSESSED (no pair given)";
else
    s = sprintf("phototrophs above the sand lit/dark %.2fx (%s)", ...
        M.ratioPhoAbove, tf(awarded, ">= 1.5x", "< 1.5x"));
end
end

function s = steadyPhrase(M)
if ~M.longEnough
    s = sprintf("run is %.1f d, shorter than the %g d steadiness window -- profile-steadiness point not available (cap 3/4)", ...
        M.tEnd - M.ts(1), M.steadyWindow);
elseif M.profileSteady
    s = sprintf("profile changed %.2f %% over the last %g d (< 1 %%)", 100*M.profileChange, M.steadyWindow);
else
    s = sprintf("profile changed %.2f %% over the last %g d (>= 1 %%)", 100*M.profileChange, M.steadyWindow);
end
end

% ======================== report =========================================
function s = report(M, S, opts)
L = strings(0);
L(end+1) = sprintf("# scoreRun — %s", M.label);
L(end+1) = "";
L(end+1) = sprintf("Tag `%s`, %d leg(s), t = %.4g d, LightScale %g, influent HET/PHO %.3g / %.3g kg/m3.", ...
    M.tag, M.legs, M.tEnd, M.lightScale, M.influent(1), M.influent(2));
if M.hasPair
    L(end+1) = sprintf("Paired against `%s` (t = %.4g d) for the light contrast.", M.pair.tag, M.pair.tEnd);
end
L(end+1) = "";
L(end+1) = sprintf("**Score %d/10** — profile %d/4, biomass %d/3, oxygen %d/3.", ...
    S.total, S.profile, S.biomass, S.oxygen);
L(end+1) = "";
for k = 1:numel(S.justification), L(end+1) = "* " + S.justification(k); end %#ok<AGROW>
L(end+1) = "";
L(end+1) = "| metric | value | unit / note |";
L(end+1) = "|---|---|---|";
row = @(a,b,c) sprintf("| %s | %s | %s |", a, b, c);
L(end+1) = row("t reached", sprintf("%.4g", M.tEnd), "d");
L(end+1) = row("phi_b peak", sprintf("%.4f", M.phiPeak), sprintf("at z = %.4f m", M.zPeak));
L(end+1) = row("largest rise below peak", sprintf("%.3f %%", 100*M.maxRiseBelow), "of the peak; > 1 % = secondary maximum");
L(end+1) = row("profile change, last window", sprintf("%.3f %%", 100*M.profileChange), sprintf("over %g d", M.steadyWindow));
L(end+1) = row("**biomass, supernatant**", sprintf("%.6f", M.massSup), "kg/m2 (z < -delta)");
L(end+1) = row("**biomass, roughness layer**", sprintf("%.6f", M.massRough), "kg/m2 (-delta <= z < 0)");
L(end+1) = row("**biomass, sand bed**", sprintf("%.6f", M.massBed), "kg/m2 (z >= 0)");
L(end+1) = row("**biomass, TOTAL column**", sprintf("%.6f", M.massTot), "kg/m2, eps-weighted HET+PHO+POM");
L(end+1) = row("supernatant share", sprintf("%.4f %%", 100*M.supFrac), "of the total");
L(end+1) = row("bed integral of eps*phi_b", sprintf("%.5f", M.bedInt), "m — probeChain's `bedInt`, EXPERIMENTS.md column");
L(end+1) = row("supernatant integral of eps*phi_b", sprintf("%.3g", M.supInt), "m");
L(end+1) = row("largest drop in total mass", sprintf("%.3f %%", 100*M.maxDropFrac), "of the maximum");
L(end+1) = row("top-2 cm biomass, dry reading", sprintf("%.1f", M.top2Dry), sprintf("ug C/g, kappa = %.1f, %d cells", M.kappa, M.top2Cells));
L(end+1) = row("top-2 cm biomass, wet reading", sprintf("%.1f", M.top2Wet), sprintf("ug C/g, dry / %g (f_dry unresolved; the factor has no independent provenance)", opts.WetFactor));
L(end+1) = row("**scoring convention**", upper(M.convention), "the Campos biomass point is scored on this reading only; the score is conditional on it");
L(end+1) = row("biomass species summed", "HET + PHO + POM", "dosed PAT is NOT included, so a marker pulse cannot inflate the biomass axis");
L(end+1) = row("phototrophs above sand", sprintf("%.6f", M.phoAbove), "kg/m2 in the roughness layer");
if M.hasPair
    L(end+1) = row("lit/dark, phototrophs above sand", sprintf("%.3f", M.ratioPhoAbove), "x");
    L(end+1) = row("lit/dark, top-2 cm", sprintf("%.3f", M.ratioTop2), "x");
    L(end+1) = row("lit/dark, bed biomass", sprintf("%.3f", M.ratioBed), "x (HET+PHO+POM mass)");
    L(end+1) = row("lit/dark, bed eps*phi_b", sprintf("%.3f", M.ratioBedInt), "x (EXPERIMENTS.md convention)");
end
L(end+1) = row("H/H0", sprintf("%.3f", M.HH0), "Kozeny-Carman, z >= delta");
L(end+1) = row("effluent O2 mean / min / max", sprintf("%.3f / %.3f / %.3f", M.o2Mean, M.o2Min, M.o2Max), sprintf("mg/L, last %.2g d", M.o2Window));
L(end+1) = row("O2 consumed", sprintf("%.3f", M.o2Consumed), sprintf("mg/L (influent %.2f)", M.o2In));
L(end+1) = "";
L(end+1) = sprintf("_Generated by analysis/scoreRun.m, %s._", string(datetime("now")));
L(end+1) = "";
s = strjoin(L, newline);
end
