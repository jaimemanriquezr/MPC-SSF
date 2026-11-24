# MPC-SSF

Multiphase continuum simulations of slow sand filtration (one-dimensional), in MATLAB.

## Requirements

- MATLAB R2026b or later.
- Optional: Symbolic Math Toolbox, used only by `isequal` on `mpcssf.Model` and
  `mpcssf.SandFilter` to compare function handles. Simulations run without it.

MATLAB and the toolbox can be installed from the command line with
[MATLAB Package Manager (`mpm`)](https://www.mathworks.com/products/mpm.html), e.g.

```sh
mpm install --release=R2026b --destination=$HOME/matlab --products MATLAB Symbolic_Math_Toolbox
```

## Installation

Three ways; each puts `src/` (the `mpcssf` namespace) on the MATLAB path.

**1. MATLAB Package Manager.** Install from File Exchange by name:

```matlab
mpminstall("mpcssf")       % install
mpmuninstall("mpcssf")     % remove
```

**2. From source.** Clone the repository and run `initpath`, which adds `src/` to the
path for the current session:

```sh
git clone https://github.com/jaimemanriquezr/MPC-SSF.git
```

```matlab
run("MPC-SSF/initpath.m")
```

For development, open the MATLAB project instead (`MPCSSF.prj`, or
`openProject("MPC-SSF/MPCSSF.prj")`): it puts `src/` on the path while it is open and
registers `tests/` as its test folder.

**3. Toolbox file.** Download `MPCSSF-<version>.mltbx` from the repository's Releases
page and double-click it, or run
`matlab.addons.toolbox.installToolbox("MPCSSF-<version>.mltbx")`. Uninstall from the
Add-On Manager.

## Quick start

After any of the three installs:

```matlab
filter = mpcssf.SandFilter(Temperature=19);
filter = filter.addGridPoints(100);

model = mpcssf.presets.pathogenModel();
% influent, in the order of model.Components: HET PHO POM PAT O2 IC NH4 HPO4 DOM
inflow = [2.68e-3, 1.0e-2, 0, 0, 9.1e-3, 6.23e-3, 2.0e-5, 0, 1.75e-4];

results = mpcssf.State(filter, model).simulate(InflowConcentrations=inflow, SimulationTime=1.0);
results.plotBiofilm();
```

Everything lives in the `mpcssf` namespace; `import mpcssf.*` lets you drop the prefix.
`mpcssf.presets.modelLund()` and `mpcssf.presets.pathogenModel()` build the
parameterisation of Manríquez et al. (2026); read the header of
`src/+mpcssf/+presets/modelLund.m` before reusing its values.

## Building an ecological model

```matlab
import mpcssf.*

model = Model();
model.WaterDensity = 1.0;

P  = Particle(Name="Microorganism");  P.Density = 1.1;
P0 = Particle(Name="POM");            P0.Density = 1.1;
L  = Liquid(Name="Nutrient");         L.Density = 1.0;

growth = Reaction(Name="Growth", Order=dictionary(P, 1), ...
                  HalfSaturationConstants=dictionary(L, 1e-3), ...
                  StoichiometricCoefficients=dictionary(P, 1.0, L, -1.0));
death  = Reaction(Name="Inactivation", Order=dictionary(P, 1), ...
                  StoichiometricCoefficients=dictionary(P, -1.0, P0, 1.0));

model.Components = [P, P0, L];
model.Reactions  = [growth, death];
```

More in `examples/`.

## Tests

The tests are `matlab.unittest` classes in `tests/` (the project's test folder),
tagged `Unit`, `Regression` and `Diagnostics`. From the project root:

```matlab
runtests("tests/", Tag="Unit")          % seconds
runtests("tests/", Tag="Regression")    % ~30 s, 19 configurations checked bit-for-bit
runtests("tests/", Tag="Diagnostics")   % ~30 min
runtests("tests/")                      % everything
```

`addpath("tests"); runTests("unit")` does the same and returns a single pass/fail.

## Building

With the project open, from its root:

```matlab
buildtool                 % check (Code Analyzer) + test (Unit, Regression)
buildtool diagnostics     % the ~30 min Diagnostics suite
buildtool package         % check + test, then release/MPCSSF-<version>.mltbx
```

Name, version and identifier live in `resources/mpackage.json`, which both
`buildtool package` and `mpminstall` read; bump the version there.

## Contributing

Work with the project open. Add new files to the project (Project tab, or
`addFile(currentProject, "path/to/file.m")`) so they are tracked with it, and keep
`buildtool` green before committing; the Regression suite compares results bit for
bit, so a change that moves any number must update the references deliberately
(`tests/regression/captureBaseline.m`).

## Layout

| path | contents |
|---|---|
| `src/+mpcssf/` | the model: `SandFilter`, `Model`, `State` (`simulate`), `Results`, `SolverOptions`, `Particle`, `Liquid`, `Reaction`, `CahnHilliardModel` |
| `src/+mpcssf/+presets/` | `modelLund`, `pathogenModel` |
| `examples/` | short scripts |
| `tests/` | unit, regression and diagnostic suites |
| `MPCSSF.prj`, `resources/project/` | MATLAB project definition (path: `src/`; test folder: `tests/`) |
| `resources/mpackage.json` | package metadata (`mpminstall` and the `.mltbx`) |
| `buildfile.m` | `buildtool` tasks |
| `initpath.m` | adds `src/` to the path when not using the project |

## License

MIT, see `LICENSE`.
