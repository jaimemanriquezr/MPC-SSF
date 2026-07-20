# Golden-master validation (Julia ↔ MATLAB)

Proves the Julia solver reproduces the MATLAB reference numerically, not just
structurally. Both sides run an **identical** simulation of `SimpleModel`
(`data/SimpleModel.mat`), the minimal model the repo's own example uses — chosen
because it runs to completion (`flag = OK`) on a clean filter, unlike the full
Lund preset which trips a negativity guard immediately.

## Shared run spec

Defined identically in `export_reference.m` and `compare.jl`:

- `SandFilter()` defaults, `addGridPoints(20)`
- model: `SimpleModel.mat` with the example's overrides — cohesion `Kappa=1e-2`,
  `Zeta0=1.0`, detachment `@(v) sqrt(abs(v))` — mirrored by Julia `simpleModel()`
- inflow: `Microorganism = 1e-2`, `Nutrient = 1.0`
- `SimulationTime = 1e-3`, `TimeStep = 1e-5`, `FrameNumber = 5`

## Running it

1. Export the MATLAB reference (needs MATLAB on PATH), from the repo root:

   ```
   cp julia/test/golden/export_reference.m .
   matlab -batch "initpath; export_reference('/abs/path/to/refdir')"
   rm export_reference.m
   ```

   (`export_reference.m` must sit where `initpath` puts the classes on the path;
   copying it to the repo root is the simplest way.)

2. Compare, either standalone (detailed diff table):

   ```
   julia --project=julia julia/test/golden/compare.jl /abs/path/to/refdir
   ```

   or through the test suite (the `golden-master vs MATLAB` testset runs when the
   env var is set, and is skipped otherwise):

   ```
   MPCSSF_GOLDEN_REF=/abs/path/to/refdir julia --project=julia -e 'using Pkg; Pkg.test()'
   ```

## Current status

Flags agree (`OK`), and every field matches to machine precision **except** the
two flowing-phase concentrations (`Microorganism/flowing`, `Nutrient/flowing`),
which differ by a max abs of ~9e-9 (rel ~7e-4). The identical error magnitude in
both fields points to a single localized cause in the flowing-phase
convection/dispersion update (SOLVER B) — under investigation. Tolerance is
`rtol=1e-8` (a field passes on abs OR rel within tol).
