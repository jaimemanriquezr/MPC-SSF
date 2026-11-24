function plan = buildfile
% BUILDFILE  Build tasks for MPC-SSF.
%
%   buildtool               % check + test
%   buildtool check         % Code Analyzer on src, tests and examples: fail on errors
%   buildtool test          % Unit and Regression suites (~1 min)
%   buildtool diagnostics   % Diagnostics suite (~30 min)
%   buildtool package       % check + test, then release/MPCSSF-<version>.mltbx
%
% Name, version and identifier come from resources/mpackage.json (the metadata
% mpminstall also reads); bump the version there.

import matlab.buildtool.tasks.CodeIssuesTask
import matlab.buildtool.tasks.TestTask

plan = buildplan(localfunctions);

plan("check") = CodeIssuesTask(["src", "tests", "examples"]);
plan("check").Description = "Code Analyzer on src, tests and examples";

plan("test") = TestTask("tests/", Tag=["Unit", "Regression"]);
plan("test").Description = "Unit and Regression suites";

plan("diagnostics") = TestTask("tests/", Tag="Diagnostics");
plan("diagnostics").Description = "Diagnostics suite (simulates; ~30 min)";

plan("package").Dependencies = ["check", "test"];

plan.DefaultTasks = ["check", "test"];
end

function packageTask(context)
% Build release/MPCSSF-<version>.mltbx: src/ on the path, plus examples, license, readme.
root = context.Plan.RootFolder;
meta = jsondecode(fileread(fullfile(root, "resources", "mpackage.json")));

opts = matlab.addons.toolbox.ToolboxOptions(root, meta.id);
opts.ToolboxFiles = fullfile(root, ["src", "examples", "LICENSE", "README.md"]);
opts.ToolboxMatlabPath = fullfile(root, "src");
release = fullfile(root, "release");
if ~isfolder(release), mkdir(release); end
opts.OutputFile = fullfile(release, "MPCSSF-" + meta.version + ".mltbx");

matlab.addons.toolbox.packageToolbox(opts);
fprintf("Wrote %s\n", opts.OutputFile);
end
