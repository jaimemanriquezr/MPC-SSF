classdef (TestTags = {'Diagnostics'}) DiagnosticTest < matlab.unittest.TestCase
% DIAGNOSTICTEST  Correctness checks that must simulate (minutes).
%
%   runtests("tests/", Tag="Diagnostics")
%
% The probes in tests/diagnostics/ compute; the pass criteria live here. Probe
% output goes to tests/diagnostics/data/ (ignored by git).

    methods (TestClassSetup)
        function addPaths(tc)
            root = fileparts(fileparts(mfilename("fullpath")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "src")));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, "tests", "diagnostics")));
        end
    end

    methods (Test)
        function cahnHilliardSecondOrder(tc)
            % Manufactured solution: the kappa operator converges at order 2.
            res = probeMMSCahnHilliard();
            tc.verifyLessThan(abs(res.orderInterior(end) - 2), 0.15, ...
                "interior order of the Cahn-Hilliard operator is not 2");
        end

        function zeroBiologyMassClosure(tc)
            % With every reaction off, each particle species' eps-weighted stock
            % changes by supply - export to 1e-4 of supply. Days=5: at the probe's
            % default 10 d the zero-biology bed clogs at t = 7.84 d.
            out = probeMassClosure(Days=5);
            for n = string(out.species(:))'
                tc.verifyLessThan(abs(out.(n).residualOverSupply), 1e-4, ...
                    n + " does not close to 1e-4 of supply");
            end
        end

        function restartChainingMatchesContinuous(tc)
            tc.verifyTrue(probeRestart().pass, "a chained run differs from a continuous one");
        end

        function cohesionBoundaryConditions(tc)
            tc.verifyTrue(probeBCSelfTest().pass, "Cahn-Hilliard boundary-condition gate failed");
        end

        function pathogenPreset(~)
            % The audited pathogenModel preset, field by field plus a PAT budget.
            testPathogen
        end
    end
end
