function [T, P] = collectCover(stage, opts)
% COLLECTCOVER  Gather cover-sweep cells into tables and paired ratios.
%
%   [T, P] = collectCover("A")        % long table T, pairs/ladder table P
%   [T, P] = collectCover("all", WriteCsv=false)
%
% Writes, under analysis/results/cover_sweep/:
%   summary_<stage>.csv     one row per cell, every QoI
%   pairs_<stage>.csv       one row per pairKey member vs its uncovered reference
%   limitation_<stage>.csv  binding substrate and min-Monod per region/phase
%
% GROUPING ------------------------------------------------------------------
% Cells are grouped by pairKey (the tag with the cover field removed), never by
% parsing filenames. A group with exactly 2 members is a PAIR; a group with more
% is a LADDER (Stage A sweeps five cover levels against one reference). Both are
% reduced the same way: every member is divided by the group's cover == 1 member,
% so a ladder is just a pair with extra rungs and the same code handles both.
%
% INTEGRITY -----------------------------------------------------------------
% Each .mat carries the cfg it was produced from. Before using a cell, that cfg
% is compared field-by-field against the current registry over
% coverConfigFields() -- the physics fields only, since `stage` is bookkeeping
% and the same configuration may be requested by two stages. A mismatch is an
% ERROR, not a silent reuse: it means the registry was edited after the run.
%
% MASS GUARD ----------------------------------------------------------------
% Trusting a ratio rests on the two arms sharing the transport operator so that
% the open mass-conservation defect is common-mode. That is checked, not
% assumed: a pair whose members' mass_ratio differ by more than MassTol (5% by
% default) is flagged massOK = false.
%
% See also COVERCASES, COVERTAG, PROBECOVER.

arguments
    stage (1,1) string {mustBeMember(stage, ["A","B","Bprime","C","smoke","all"])} = "all"
    opts.WriteCsv (1,1) logical = true
    opts.MassTol  (1,1) double  = 0.05
end

here = fileparts(mfilename("fullpath"));            % .../analysis
W    = fileparts(here);                             % repo root
addpath(fullfile(here, "probes")); addpath(genpath(fullfile(W,"src")));

dataDir = fullfile(here, "probes", "data");
outDir  = fullfile(here, "results", "cover_sweep");
if ~isfolder(outDir), mkdir(outDir); end

cases = coverCases(stage);
fields = coverConfigFields();

rows = {}; limRows = {}; missing = strings(0);
for i = 1:numel(cases)
    cfg = cases(i);
    f = fullfile(dataDir, "cover_" + cfg.tag + ".mat");
    if ~isfile(f), missing(end+1) = cfg.tag; continue, end %#ok<AGROW>
    S = load(f);

    % --- integrity: stored cfg must still match the registry ---------------
    for fn = fields
        a = S.cfg.(fn); b = cfg.(fn);
        assert(isequal(a, b), ...
            "STALE ARTEFACT %s: field %s is %s on disk but %s in the registry. " + ...
            "Delete the .mat and rerun, or revert the registry.", ...
            cfg.tag, fn, string(a), string(b));
    end

    r = S.res;
    rows{end+1,1} = packRow(cfg, r); %#ok<AGROW>
    limRows{end+1,1} = packLim(cfg, r); %#ok<AGROW>
end

if isempty(rows)
    warning("collectCover:empty", "no completed cells found for stage %s (%d missing)", stage, numel(missing));
    T = table(); P = table(); return
end

T = vertcat(rows{:});
L = vertcat(limRows{:});

%% ---- group into pairs / ladders ---------------------------------------
P = table();
keys = unique(T.pairKey, "stable");
for k = 1:numel(keys)
    g = T(T.pairKey == keys(k), :);
    ref = g(g.cover == 1, :);
    if height(ref) ~= 1
        warning("collectCover:noref", "pairKey %s has %d uncovered references; skipped", keys(k), height(ref));
        continue
    end
    for j = 1:height(g)
        if g.cover(j) == 1, continue, end       % the reference is not a row of its own
        P = [P; ratioRow(g(j,:), ref, opts.MassTol)]; %#ok<AGROW>
    end
end

%% ---- report ------------------------------------------------------------
fprintf("\n=== cover sweep: stage %s ===\n", stage);
fprintf("cells found %d / %d", height(T), numel(cases));
if ~isempty(missing), fprintf("  (missing: %s)", strjoin(missing, ", ")); end
fprintf("\n");
bad = T.flag ~= "OK";
if any(bad), fprintf("NON-OK FLAGS: %s\n", strjoin(T.tag(bad), ", ")); end

if ~isempty(P)
    fprintf("\n--- covered / uncovered ratios (targets: R_B_0_2 0.12-0.25, |dTOC| < 3pp, R_mat < 0.1) ---\n");
    for j = 1:height(P)
        fprintf("  %-46s s=%-7.4g R_B02=%6.3f R_B010=%6.3f R_mat=%6.3f dDOC=%+6.2fpp eps_L=%+6.3f %s%s\n", ...
            P.pairKey(j), P.cover(j), P.R_B_0_2(j), P.R_B_0_10(j), P.R_mat(j), ...
            100*P.dDOC(j), P.eps_L(j), ...
            ternaryStr(P.jointPass(j), "PASS", "    "), ...
            ternaryStr(~P.massOK(j), "  [MASS GUARD FAILED]", ""));
    end
end

fprintf("\n--- light and limitation (final-window) ---\n");
for j = 1:height(L)
    fprintf("  %-46s Ihat surf=%.3g z0=%.3g bed1=%.3g | PHO sup bio: %s %.3g | PHO top2 bio: %s %.3g\n", ...
        L.tag(j), L.Ihat_surface(j), L.Ihat_z0(j), L.Ihat_bed1(j), ...
        L.PHO_sup_bio_species(j), L.PHO_sup_bio_monod(j), ...
        L.PHO_top2_bio_species(j), L.PHO_top2_bio_monod(j));
end

%% ---- write -------------------------------------------------------------
if opts.WriteCsv
    writetable(T, fullfile(outDir, "summary_" + stage + ".csv"));
    writetable(L, fullfile(outDir, "limitation_" + stage + ".csv"));
    if ~isempty(P), writetable(P, fullfile(outDir, "pairs_" + stage + ".csv")); end
    fprintf("\nwrote %s/{summary,limitation%s}_%s.csv\n", outDir, ternaryStr(~isempty(P), ",pairs", ""), stage);
end
end


% =========================================================================
function row = packRow(cfg, r)
row = table();
row.tag = string(cfg.tag);  row.pairKey = string(cfg.pairKey);  row.stage = string(cfg.stage);
row.cover = cfg.cover;  row.phoIn = cfg.phoIn;  row.nh4In = cfg.nh4In;
row.hpo4In = cfg.hpo4In;  row.icIn = cfg.icIn;  row.tsim = cfg.tsim;
row.ncells = cfg.ncells;  row.etaSand = cfg.etaSand;  row.kinetics = string(cfg.kinetics);
row.flag = string(r.flag);  row.tFinal = r.tFinal;  row.wallMin = r.wallMin;
for fn = ["B_0_2","B_0_10","B_0_2_pom","B_0_10_pom","Bpho_0_2","Bpho_0_10","Bhet_0_2", ...
          "Bpeak","zPeak","phib_0_2","phib_0_10","phib_max","phib_sup", ...
          "mat_mass","mat_depth","mat_present","doc_removal","toc_removal", ...
          "o2_out","o2_bed_min","psFrac_sup","psFrac_0_2","psFrac_deep", ...
          "stock_raw","stock_eps","supplied","mass_ratio"]
    row.(fn) = r.(fn);
end
% Added after Stage A; absent from earlier artefacts.
for fn = ["outflux","closure_residual"]
    if isfield(r, fn), row.(fn) = r.(fn); else, row.(fn) = NaN; end
end
end

% =========================================================================
function row = packLim(cfg, r)
row = table();
row.tag = string(cfg.tag);  row.cover = cfg.cover;
L = r.lim;
row.Ihat_surface = L.Ihat_surface;  row.Ihat_z0 = L.Ihat_z0;  row.Ihat_bed1 = L.Ihat_bed1;
for org = ["Heterotroph","Phototroph"]
    short = extractBefore(org, 4);        % "Het" / "Pho"
    for regn = ["sup","top2"]
        for phase = ["biofilm","flowing"]
            src = org + "_" + regn + "_" + phase;
            dst = upper(short) + "_" + regn + "_" + extractBefore(phase, 4);
            if isfield(L, src + "_species")
                row.(dst + "_species") = string(L.(src + "_species"));
                row.(dst + "_frac")    = L.(src + "_frac");
                row.(dst + "_monod")   = L.(src + "_monod");
            end
        end
    end
end
end

% =========================================================================
function row = ratioRow(g, ref, massTol)
% Covered / uncovered. Every reported quantity is a ratio or a difference, so
% the constants that are unresolved at model level (beta, rho_P, the attachment
% prefactor, the absolute influent scale) divide out.
% NB: do NOT write this as arithmetic masking -- NaN*0 is NaN in MATLAB, so
% `a/b*(b>tol) + NaN*(b<=tol)` returns NaN for EVERY input. Caught by the smoke
% test, which reported R_B_0_2 = NaN for a pair whose true ratio was 1.0.
saferatio = @(a, b) safediv(a, b);
row = table();
row.pairKey = g.pairKey;  row.tag = g.tag;  row.stage = g.stage;
row.cover = g.cover;
row.R_B_0_2   = saferatio(g.B_0_2,   ref.B_0_2);
row.R_B_0_10  = saferatio(g.B_0_10,  ref.B_0_10);
row.R_Bpho_0_2= saferatio(g.Bpho_0_2,ref.Bpho_0_2);
row.R_phib_0_2= saferatio(g.phib_0_2,ref.phib_0_2);
row.R_mat     = saferatio(g.mat_mass, ref.mat_mass);
row.dDOC      = g.doc_removal - ref.doc_removal;      % fraction, printed as pp
row.dTOC      = g.toc_removal - ref.toc_removal;
row.dO2       = g.o2_out - ref.o2_out;

% Elasticity of the biomass response to the cover multiplier. eps_L = 1 would be
% a perfectly proportional response; ~0 means the cover does nothing.
if g.cover > 0 && g.cover < 1 && isfinite(row.R_B_0_2) && row.R_B_0_2 > 0
    row.eps_L = log(row.R_B_0_2) / log(g.cover);
else
    row.eps_L = NaN;                                   % s = 0 has no log
end

% Campos joint criterion: the biomass ratio must be won WITHOUT destroying
% removal and WITH the schmutzdecke suppressed.
row.jointPass = row.R_B_0_2 >= 0.12 && row.R_B_0_2 <= 0.25 ...
             && abs(row.dTOC) < 0.03 && row.R_mat < 0.1;

% Common-mode argument, checked rather than assumed.
%
% Compares the CLOSURE RESIDUAL (in - out - dStock)/in, not stock/supplied. The
% latter was the first attempt and is not a numerical diagnostic at all: two arms
% that grow differently hold different stock for physical reasons, so it flagged
% every Stage A pair while nothing was wrong. The residual is arm-independent for
% a conservative scheme even when the arms carry very different biomass.
% Falls back to the old test for .mat files written before the field existed.
if ismember("closure_residual", g.Properties.VariableNames) && ~isnan(g.closure_residual)
    row.massOK = abs(g.closure_residual - ref.closure_residual) <= massTol;
    row.dResidual = g.closure_residual - ref.closure_residual;
else
    row.massOK = NaN;                 % undetermined, not "passed"
    row.dResidual = NaN;
end
end

% =========================================================================
function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end

% =========================================================================
function q = safediv(a, b)
if abs(b) > realmin, q = a/b; else, q = NaN; end
end
