function cases = coverCases(stage)
% COVERCASES  Registry for the covered/uncovered driver sweep.
%
% Single source of truth: probeCover runs cases(idx), collectCover pairs them.
% Nothing downstream parses filenames -- tags come from coverTag() and pairing
% comes from the pairKey field, so an edit here propagates everywhere.
%
% Design: .claude/plans/2026-08-22-cover-sweep-and-referee-experiments.md and its
% 2026-08-22 amendment. Motivation: Campos2002 measured a 4x biomass difference
% between two full-scale Thames beds (uncovered Bed 9 / covered Bed 10, 104 d,
% same pretreatment) with NO schmutzdecke on the covered bed, yet IDENTICAL
% TOC/DOC removal (25/23% vs 23/23%). Our E2 gives 4%, and 0.97 on the corrected
% stack. Four candidate blockers, tested by the stages below.
%
% Stage A  light ladder + metric audit    5 cells, 30 d
% Stage B  nutrient x influent, paired   14 cells, 30 d   (supernatant limb)
% Stage Bp euphotic-depth counterfactual  4 cells, 10 d   (bed limb; SCHEDULED,
%                                                          not contingent -- see
%                                                          the amendment)
% Stage C  confirmation at Campos length  4 cells, 104 d
%
%   cases = coverCases("A" | "B" | "Bprime" | "C" | "smoke" | "all")
%
% Fields per case:
%   stage    string   which stage it belongs to
%   tag      string   unique, encodes every varied field (coverTag)
%   pairKey  string   tag minus the cover field -> collectCover groups on this
%   cover    double   multiplier on LightIrradiation (1 = uncovered, 0 = dark)
%   phoIn    double   influent PHO   [kg/m3]
%   nh4In    double   influent NH4   [kg/m3]
%   hpo4In   double   influent HPO4  [kg/m3]
%   icIn     double   influent IC    [kg/m3]
%   tsim     double   simulated days
%   ncells   double   grid points
%   maxdt    double   AdaptiveMaxDt
%   kinetics string   "corrected" | "manuscript"
%   etaSand  double   LightAttenuationCoeffSand  (Bprime counterfactual)
%
% See also PROBECOVER, COLLECTCOVER, COVERTAG.

arguments
    stage (1,1) string {mustBeMember(stage, ["A","B","Bprime","C","smoke","all"])} = "all"
end

% ---- nominal levels -----------------------------------------------------
% Table B.1 (manuscriptExperiments/baseInfluent) unless a factor overrides.
PHO_TABLE = 1.00e-2;    % Table B.1
PHO_FIELD = 5.00e-4;    % Campos2006b Fig.1(a)
NH4_NOM   = 2.00e-5;    % Table B.1 (Chan2018)
HPO4_NOM  = 0.0;        % Table B.1 -- exactly zero, the recycling-loop driver
IC_NOM    = 6.23e-3;    % Table B.1 (Campeau2017)

% Add-back levels. Chosen from the limitation diagnostic
% (.claude/decisions/2026-08-22-liebig-limitation-diagnostic.md), NOT from the
% original plan: on the corrected stack the binding species are IC (supernatant
% biofilm, 99%) and HPO4 (flowing, 100%) -- NOT NH4, which the first draft of the
% plan would have relieved to no effect.
HPO4_ADD  = 5.0e-4;     % ~25x the corrected K_HPO4 (2e-5); Thames-era SRP scale
IC_ADD    = 6.23e-2;    % 10x Table B.1; ~50x the corrected K_IC (1.2e-3)
NH4_ADD   = 1.0e-2;     % retained for the manuscript-stack continuity cell only

ETA_SAND_NOM = 1500;    % SandFilter default
ETA_SAND_CF  = 150;     % Bprime counterfactual: euphotic depth 5 mm -> 5 cm

base = struct( ...
    "stage", "A", "tag", "", "pairKey", "", "cover", 1.0, ...
    "phoIn", PHO_TABLE, "nh4In", NH4_NOM, "hpo4In", HPO4_NOM, "icIn", IC_NOM, ...
    "tsim", 30.0, "ncells", 100, "maxdt", 3e-6, ...
    "kinetics", "corrected", "etaSand", ETA_SAND_NOM);

cases = base([]);   % 0x1 struct array with the right fields
add = @(c) c;       % readability only

% ---- Stage A: light ladder ----------------------------------------------
% Uniformly spaced in eta: -ln(s) = 0, 2.30, 4.61, 6.91, Inf.
% s = 0.01 is the code of record (manuscriptExperiments.m:287); s = 0.001 is the
% manuscript text (PARAMETERS.md:124); s = 0 is the physical bound of a
% light-excluding membrane and costs one run to learn the ceiling.
if any(stage == ["A","all"])
    for s = [1.0, 0.1, 0.01, 0.001, 0.0]
        c = base; c.stage = "A"; c.cover = s;
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
end

% ---- Stage B: nutrient add-back x influent, paired on cover -------------
% 16-cell grid; the two (cover in {1,0.001}, N = none, I = TABLE) cells are
% already in Stage A, so 14 are new here. Kept in the registry anyway -- the
% collector deduplicates by tag, and a self-contained Stage B is easier to
% resubmit.
if any(stage == ["B","all"])
    for phoIn = [PHO_TABLE, PHO_FIELD]
        for nutrient = ["none", "P", "IC", "ALL"]
            for s = [1.0, 0.001]
                c = base; c.stage = "B"; c.cover = s; c.phoIn = phoIn;
                switch nutrient
                    case "P",   c.hpo4In = HPO4_ADD;
                    case "IC",  c.icIn   = IC_ADD;
                    case "ALL", c.hpo4In = HPO4_ADD; c.icIn = IC_ADD; c.nh4In = NH4_ADD;
                end
                cases(end+1,1) = add(c); %#ok<AGROW>
            end
        end
    end
end

% ---- Stage Bprime: euphotic-depth counterfactuals -----------------------
% SCHEDULED, not contingent. Diagnostic measurement: Ihat = 0.0608 at z = 0 and
% 1.39e-27 one cell in, so the bed is optically dark and Campos's 0-2 cm layer is
% two dark cells. Neither arm is a proposed model change -- they exist to show
% whether euphotic depth is what gates the bed limb.
%   eta150 : LightAttenuationCoeffSand 1500 -> 150  (euphotic depth 5 mm -> 5 cm)
%   n300   : NCells 100 -> 300 (dz 9.95 mm -> 3.32 mm), resolves the euphotic layer
if any(stage == ["Bprime","all"])
    for s = [1.0, 0.0]
        c = base; c.stage = "Bprime"; c.cover = s; c.tsim = 10.0;
        c.etaSand = ETA_SAND_CF;
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
    for s = [1.0, 0.0]
        c = base; c.stage = "Bprime"; c.cover = s; c.tsim = 10.0;
        c.ncells = 300; c.maxdt = 1e-6;   % finer grid tightens the CFL bound
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
end

% ---- Stage C: confirmation at Campos's own run length -------------------
% 104 d is Bed 9's actual run (Campos2002 Table 1), removing a duration confound
% for ~14 min per run. Baseline pair = the null; field+ALL pair = the best case
% Stage B can offer. One manuscript-stack pair is the continuity control that
% must recover the known ~4% / 0.97.
if any(stage == ["C","all"])
    for s = [1.0, 0.001]
        c = base; c.stage = "C"; c.cover = s; c.tsim = 104.0;
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
    for s = [1.0, 0.001]
        c = base; c.stage = "C"; c.cover = s; c.tsim = 104.0;
        c.phoIn = PHO_FIELD; c.hpo4In = HPO4_ADD; c.icIn = IC_ADD; c.nh4In = NH4_ADD;
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
    for s = [1.0, 0.01]
        c = base; c.stage = "C"; c.cover = s; c.tsim = 104.0;
        c.kinetics = "manuscript";
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
end

% ---- smoke: tiny end-to-end check of the whole chain --------------------
if any(stage == "smoke")
    for s = [1.0, 0.0]
        c = base; c.stage = "smoke"; c.cover = s;
        c.tsim = 0.2; c.ncells = 30;
        cases(end+1,1) = add(c); %#ok<AGROW>
    end
end

% ---- tags, pair keys, dedup ---------------------------------------------
for i = 1:numel(cases)
    cases(i).tag = coverTag(cases(i));
    cases(i).pairKey = coverTag(cases(i), OmitCover=true);
end
[~, keep] = unique([cases.tag], "stable");
cases = cases(keep);
end
