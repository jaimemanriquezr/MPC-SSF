function ok = runTests(suite)
% RUNTESTS  Run the matlab.unittest suites by tag.
%
%   ok = runTests("unit")          % seconds   -- property/option checks
%   ok = runTests("regression")    % ~40 s     -- 19 configurations, bit-identical
%   ok = runTests("diagnostics")   % ~30 min   -- correctness checks that simulate
%   ok = runTests("all")           % all three, in that order
%
% Equivalent to runtests("tests/", Tag="Unit") etc. Returns TRUE only if every
% test that ran passed, and prints the result table.
%
% LAYOUT
%   tests/UnitTest.m, RegressionTest.m, DiagnosticTest.m   the TestCase classes
%   tests/unit/         assertion scripts run by UnitTest
%   tests/regression/   regressionConfigs, runConfig, seedState, captureBaseline,
%                      checkBaseline, compareBaseline, selectConfigs,
%                      checkGateVision, and data/ (the 19 references)
%   tests/diagnostics/  the probes run by DiagnosticTest and their helpers
%
% NOT run from here: checkGateVision(), which asks whether the gate can SEE, not
% whether it passes. Run it by hand after any captureBaseline and whenever a
% config is added.

arguments
    suite (1,1) string {mustBeMember(suite, ["unit", "regression", "diagnostics", "all"])} = "all"
end

here = fileparts(mfilename("fullpath"));
tags = dictionary("unit", "Unit", "regression", "Regression", "diagnostics", "Diagnostics");
if suite == "all"
    selected = values(tags);
else
    selected = tags(suite);
end

tests = matlab.unittest.Test.empty;
for tag = selected(:)'
    tests = [tests, matlab.unittest.TestSuite.fromFolder(here, Tag=tag)]; %#ok<AGROW>
end
results = matlab.unittest.TestRunner.withTextOutput().run(tests);
disp(table(results));
ok = all([results.Passed]);
fprintf("OVERALL: %s\n", string(ok).replace(["true", "false"], ["PASS", "FAIL"]));
end
