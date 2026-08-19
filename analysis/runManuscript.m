% RUNMANUSCRIPT  Full manuscript-experiment suite on the corrected model.
% Stage order and dependencies live in manuscriptExperiments.m ("all").
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, "..", "src")));
addpath(here);
fprintf("manuscript suite start (%s)\n", string(datetime));
try
    tic; manuscriptExperiments("all");
    fprintf("manuscript suite complete in %.1f h (%s)\n", toc/3600, string(datetime));
catch err
    fprintf(2, "manuscript suite FAILED: %s\n%s\n", err.message, ...
        getReport(err, "extended", "hyperlinks", "off"));
end
