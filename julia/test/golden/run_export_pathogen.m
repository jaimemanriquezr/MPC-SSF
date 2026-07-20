% Driver: export the pathogen golden-master reference from the authoritative
% slow-sand-filtration code. Run from anywhere:
%
%   matlab -batch "run('/abs/path/to/MPC-SSF/julia/test/golden/run_export_pathogen.m')"
%
% Adjust the two roots below if the checkouts live elsewhere.

golden_dir = "/Users/jaime/Research/SSF/MPC-SSF/julia/test/golden";
SSF_DIR    = "/Users/jaime/Research/SSF/slow-sand-filtration";
OUTDIR     = fullfile(golden_dir, "reference_pathogen");

cd(SSF_DIR);                      % thesis_model.mat + classes live here
restoredefaultpath; initpath;     % put SDfilter/SDmodel/run_pathogen on the path
addpath(golden_dir);              % so export_pathogen_reference is callable

export_pathogen_reference(OUTDIR);
