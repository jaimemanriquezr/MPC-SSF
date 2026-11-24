classdef (TestTags = {'Unit'}) UnitTest < matlab.unittest.TestCase
% UNITTEST  Seconds-scale property and option checks.
%
%   runtests("tests/", Tag="Unit")
%
% Each test runs one assertion script in tests/unit/; an assert that fails
% inside the script fails the test with its message.

    methods (TestClassSetup)
        function addPaths(tc)
            root = fileparts(fileparts(mfilename("fullpath")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "src")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "tests", "unit")));
        end
    end

    methods (Test)
        function cardinalTemperatureResponse(~)
            testCardinal
        end

        function phototrophRespiration(~)
            testRespiration
        end

        function temperatureAsymmetryOptions(~)
            testAsymmetryOptions
        end
    end
end
