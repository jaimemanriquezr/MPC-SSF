% RUNMATLABCAMPAIGN  Sequential driver for the MATLAB anchor campaign.
% Order: PHO bloom study first (new science), then pulse (the decided
% methodology's main scenario), then startup, then flowstep.
here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, "..", "src")));
addpath(here);

stages = { ...
    @() phoBloomStudy(), ...
    @() logOatCampaign("pulse"), ...
    @() logOatCampaign("startup"), ...
    @() logOatCampaign("flowstep")};
names = ["phoBloom", "pulse", "startup", "flowstep"];

for k = 1:length(stages)
    fprintf("\n===== stage %d/%d: %s  (%s) =====\n", k, length(stages), names(k), string(datetime));
    try
        tic; stages{k}();
        fprintf("===== stage %s done in %.1f min =====\n", names(k), toc/60);
    catch err
        fprintf(2, "===== stage %s FAILED: %s =====\n", names(k), err.message);
        fprintf(2, "%s\n", getReport(err, "extended", "hyperlinks", "off"));
    end
end
fprintf("\ncampaign complete (%s)\n", string(datetime));
