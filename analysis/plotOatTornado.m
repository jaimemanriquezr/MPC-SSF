function T = plotOatTornado(resultsDir, opts)
% PLOTOATTORNADO  Manuscript tornado figure [FIG:oat-tornado] from an OAT campaign.
%
% Reads measures.csv from RESULTSDIR. Two layouts are accepted:
%   - one campaign dir holding a single measures.csv (e.g. log_oat_pulse_v2), or
%   - a parent dir of per-parameter dirs log_oat_pulse_<param>/measures.csv
%     (the slurm/oat_pulse.sbatch layout), concatenated here.
% Horizontal bars of the primary ranking I_rms (log-units of removal per e-fold),
% sorted, coloured by parameter block; a marker flags strongly asymmetric
% (nonlinear) responses, and clog-driver parameters are annotated. Fonts are
% sized for a single manuscript column (reviewer R2.4).
%
%   plotOatTornado("analysis/results/log_oat_pulse_v2")   % pre-audit anchor
%   plotOatTornado("analysis/results/oat")                % working-set campaign
arguments
    resultsDir (1,1) string = "analysis/results/oat"
    opts.OutDir (1,1) string = "analysis/results/figures/oat"
    opts.AsymThreshold (1,1) double = 1.0   % asym >= this flags nonlinearity
    opts.Name (1,1) string = "fig_oat_tornado"
end
T = readOatMeasures(resultsDir);
T = sortrows(T, "I_rms", "ascend");   % ascend: largest ends up at the TOP of barh
if ~isfolder(opts.OutDir), mkdir(opts.OutDir); end

blocks = unique(T.block, "stable");
% One colour per block, colourblind-safe (Wong 2011)
palette = [0.00 0.45 0.70; 0.90 0.60 0.00; 0.00 0.62 0.45; 0.80 0.40 0.70; ...
           0.35 0.70 0.90; 0.95 0.90 0.25; 0.55 0.55 0.55];

fig = figure(Visible="off", Units="centimeters", Position=[1 1 14 0.55*height(T)+3]);
ax = axes(fig); hold(ax, "on");
% one barh with per-row colours: a per-block barh call misplaces bars when a
% block's row indices are sparse (bar thickness scales with subset spacing)
hb = barh(ax, 1:height(T), T.I_rms, 0.65, EdgeColor="none", FaceColor="flat");
blockColor = @(b) palette(1 + mod(find(blocks == b, 1)-1, size(palette,1)), :);
for i = 1:height(T)
    hb.CData(i, :) = blockColor(T.block(i));
end
hBlocks = gobjects(numel(blocks), 1);
for b = 1:numel(blocks)   % invisible stubs so the legend has one handle per block
    hBlocks(b) = patch(ax, NaN, NaN, blockColor(blocks(b)), EdgeColor="none");
end
for i = 1:height(T)
    lab = "";
    if T.asymmetry(i) >= opts.AsymThreshold, lab = lab + " \ast"; end
    if T.clog_driver(i) ~= 0, lab = lab + " (clog)"; end
    if strlength(lab) > 0
        text(ax, T.I_rms(i), i, lab, FontSize=8, VerticalAlignment="middle");
    end
end
yticks(ax, 1:height(T));
yticklabels(ax, prettyParamLabels(T.param));
xlabel(ax, "I_{rms}  [log_{10}-units of removal per e-fold]", FontSize=10);
legend(ax, hBlocks, blocks, Location="southeast", FontSize=9, Box="off");
set(ax, FontSize=9, TickLabelInterpreter="tex", Box="on");
grid(ax, "on"); ax.XGrid = "on"; ax.YGrid = "off";
title(ax, "");   % captions live in the manuscript

exportgraphics(fig, fullfile(opts.OutDir, opts.Name + ".pdf"), ContentType="vector");
exportgraphics(fig, fullfile(opts.OutDir, opts.Name + ".png"), Resolution=300);
close(fig);
fprintf("wrote %s.{pdf,png} (%d parameters)\n", fullfile(opts.OutDir, opts.Name), height(T));
end

function T = readOatMeasures(resultsDir)
single = fullfile(resultsDir, "measures.csv");
if isfile(single)
    T = readtable(single, TextType="string");
else
    d = dir(fullfile(resultsDir, "log_oat_*", "measures.csv"));
    assert(~isempty(d), "plotOatTornado:noMeasures", ...
        "no measures.csv under %s (neither direct nor per-parameter)", resultsDir);
    parts = cell(numel(d), 1);
    for i = 1:numel(d)
        parts{i} = readtable(fullfile(d(i).folder, d(i).name), TextType="string");
    end
    T = vertcat(parts{:});
    % A per-parameter dir layout can hold a stray duplicate from a re-run task
    [~, keep] = unique(T.param, "stable");
    T = T(keep, :);
end
T.block = string(T.block); T.param = string(T.param);
end

function labs = prettyParamLabels(names)
% TeX labels for the 29-parameter workingset design (analysis/PARAMETERS.md).
% A name with no entry falls back to itself verbatim.
map = dictionary( ...
    "temperature",    "T", ...
    "influent_PAT",   "PAT_{in}", ...
    "dispersivity",   "\alpha", ...
    "transport_P",    "b^{tr}_P", ...
    "attach_sand",    "b^{att}_{sand}", ...
    "sand_pathogen",  "b^{att,PAT}_{sand}", ...
    "mu_HET",         "\mu_{HET}", ...
    "mu_PHO",         "\mu_{PHO}", ...
    "d_HET",          "d_{HET}", ...
    "d_PHO",          "d_{PHO}", ...
    "hydrolysis",     "k_{hyd}", ...
    "theta_growth",   "\theta_{20,\mu}", ...
    "theta_death",    "\theta_{20,d}", ...
    "K_O2_HET",       "K^{O_2}_{HET}", ...
    "K_DOM_HET",      "K^{DOM}_{HET}", ...
    "K_HPO4_HET",     "K^{HPO_4}_{HET}", ...
    "K_O2_PAT",       "K^{O_2}_{PAT}", ...
    "K_pred",         "K_{pred}", ...
    "marker_growth",  "\mu_{PAT}", ...
    "inactivation",   "d_{PAT}", ...
    "bacterivory",    "p_{PAT}", ...
    "zeta_0",         "\zeta_0", ...
    "kappa",          "\kappa", ...
    "zeta_1",         "\zeta_1", ...
    "detach_scale",   "k_{det}", ...
    "light_att_water","\eta_{water}", ...
    "light_att_sand", "\eta_{sand}", ...
    "attenuation_P",  "\eta_P", ...
    "beta_porosity",  "\beta");
labs = strings(numel(names), 1);
for i = 1:numel(names)
    if isKey(map, names(i)), labs(i) = map(names(i)); else, labs(i) = names(i); end
end
end
