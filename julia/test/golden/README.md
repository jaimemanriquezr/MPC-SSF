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

The comparison runs automatically in `Pkg.test` against the **committed**
reference in `reference/` (no MATLAB needed) — see the `golden-master vs MATLAB`
testset. To regenerate the reference or compare against a fresh MATLAB run:

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

## Status: validated

Flags agree (`OK`). Every field except the flowing phase matches to ~1e-15 or
better. The flowing-phase concentrations agree to **~9e-9 absolute (~6.5e-8
relative)**, which is accumulated floating-point difference over 100 nonlinear
time steps, **not** a port bug:

- At step 1 every flowing reaction term is identically zero in both codes (the
  fields are still clean; inflow enters via convection during the step), so the
  two solvers are structurally identical — the difference only accumulates once
  concentrations become nonzero.
- The `Microorganism/flowing` and `Nutrient/flowing` errors are **equal and
  opposite to 8 significant figures** (micro low by 8.9456e-9, nutrient high by
  8.9456e-9). So `micro + nutrient` is conserved between the codes and only the
  *split* differs — the signature of the Growth reaction (Nutrient→Micro) rate
  differing in its last bits, fed by the sparse solve `lhsCH \ rhsCH` and
  reduction orderings. MATLAB (UMFPACK / its BLAS) and Julia (OpenBLAS) legitimately
  differ below ~1e-8; bit-identical agreement there is not achievable.

Tolerance is therefore `atol=1e-7`, `rtol=1e-6` (a field passes on abs OR rel
within tol) — ~5 orders below the microorganism inflow (1e-2) and 7 below the
nutrient inflow (1.0).
