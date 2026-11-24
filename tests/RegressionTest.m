classdef (TestTags = {'Regression'}) RegressionTest < matlab.unittest.TestCase
% REGRESSIONTEST  The gate: each configuration re-run and compared bit for bit
% (isequaln) with its reference in tests/regression/data/.
%
%   runtests("tests/", Tag="Regression")
%   runtests("tests/RegressionTest.m", ParameterName="theta_*")   % a subset
%
% One test per configuration of regressionConfigs. checkBaseline prints where a
% failing configuration first differs (component, phase, frame, cell).

    properties (TestParameter)
        config = RegressionTest.configNames()
    end

    methods (TestClassSetup)
        function addPaths(tc)
            root = fileparts(fileparts(mfilename("fullpath")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "src")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "tests", "regression")));
        end
    end

    methods (Test)
        function matchesBaseline(tc, config)
            tc.verifyTrue(checkBaseline(string(config)), ...
                "output differs from tests/regression/data/baseline_" + config + ".mat");
        end
    end

    methods (Static)
        function names = configNames()
            here = fileparts(mfilename("fullpath"));
            oldPath = addpath(fullfile(here, "regression"));
            restore = onCleanup(@() path(oldPath));
            names = cellstr([regressionConfigs().name]);
        end
    end
end
