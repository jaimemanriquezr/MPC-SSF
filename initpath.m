function initpath()
% INITPATH  Put the mpcssf toolbox on the MATLAB path (from any folder).
    addpath(fullfile(fileparts(mfilename("fullpath")), "src"))
end
