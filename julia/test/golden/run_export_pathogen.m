% Driver: export the pathogen golden-master reference from the authoritative
% slow-sand-filtration code. Run from anywhere:
%
%   matlab -batch "run('/abs/path/to/SSF.jl/test/golden/run_export_pathogen.m')"
%
% golden_dir is self-locating; adjust SSF_DIR if the checkout lives elsewhere.

golden_dir = fileparts(mfilename("fullpath"));   % this file's own directory
SSF_DIR    = "/Users/jaime/Research/SSF/code/1d/slow-sand-filtration";

cd(SSF_DIR);                      % thesis_model.mat + classes live here
restoredefaultpath; initpath;     % put SDfilter/SDmodel/run_pathogen on the path
addpath(golden_dir);              % so export_pathogen_reference is callable

% Baseline: dark_respiration = 0, dark diel light -> light floor inactive.
export_pathogen_reference(fullfile(golden_dir, "reference_pathogen"));

% Light variant: constant light + dark_respiration > 0, with a phototroph-heavy
% seed, so the dark-respiration light floor max(fdark, I·e^{1-I}) is exercised
% (both branches straddle fdark across depth). Validates the light-model
% reconciliation against the authoritative floor form.
export_pathogen_reference(fullfile(golden_dir, "reference_pathogen_light"), ...
    dark_respiration = 0.1, light_const = 0.8, ...
    seed_matrix = [0.05, 0.1, 0.02, 0.01], seed_enclosed_liquids = 0.05);
