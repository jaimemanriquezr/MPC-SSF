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
   SSF_GOLDEN_REF=/abs/path/to/refdir julia --project=julia -e 'using Pkg; Pkg.test()'
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

## Adaptive-CFL golden-master

`export_adaptive_reference.m` / `compare_adaptive.jl` validate the
`time_step=:adaptive` path on **modelLund** (which, unlike SimpleModel, has
nonzero dispersivity and the reaction structure the CFL bound assumes). It
compares the per-step **dt trajectory** (`step_times`) as well as the frame
concentrations/velocities. Reference: `reference_adaptive/`.

Because the adaptive MATLAB `simulate` lives on branch `matlab-claude`, generate
the reference by overlaying those two files, then restoring:

```
git checkout matlab-claude -- src/@State/simulate.m src/@Model/Model.m
cp julia/test/golden/export_adaptive_reference.m .
matlab -batch "initpath; export_adaptive_reference('/abs/path/to/refdir')"
rm export_adaptive_reference.m
git checkout HEAD -- src/@State/simulate.m src/@Model/Model.m
```

Result: **exact match**. Both codes take the same 228 steps and every field
(dt trajectory included) agrees to ~1e-15 — the adaptive stepping is bit-for-bit
deterministic between MATLAB and Julia over this horizon (the run spec keeps a
short physical time so concentrations stay small and the CFL bound is dominated
by advection). The testset also asserts identical step counts.

## Pathogen golden-master

`export_pathogen_reference.m` / `compare_pathogen.jl` validate `modelPathogen()`
against the **authoritative** pathogen code — `@SDfilter/run_pathogen.m` +
`thesis_model.mat` in the sibling `slow-sand-filtration` checkout (this is the
only golden master that drives that repo rather than MPC-SSF's own `src/`, since
the pathogen model exists only there). Reference: `reference_pathogen/`.

The run spec removes every cross-codebase confounder so the two solvers are
numerically comparable:

- **grid**: `SDfilter.add_cells(20)` and Julia `addgridpoints(20)` build the
  identical grid (both `dz = 2/(2n+1)`, same centers).
- **temperature**: 293 K (MATLAB) / 20 °C (Julia). The two codebases use
  *different* temperature-correction formulas — MATLAB `θ^(293−T_K)`, Julia
  `θ^(T_K/293 − 1)` — but both collapse to `μ = nominal` at 293 K / 20 °C, so the
  correction is neutralized.
- **light**: the baseline uses `dark_respiration = 0` and the default diel
  forcing (which is dark, `light = 0`, over the short run), so the light factor
  is trivially equal on both sides. The `dark_respiration > 0` case is now
  reconciled and covered by its own variant — see **Light-model reconciliation**
  below.
- **detachment**: set to `sqrt(|v|/7.2)` on both sides.

To exercise the pathogen physics the run **seeds a uniform mature biofilm** (a
clean filter stays clean, since growth is order-1 in existing biomass). That
makes the biofilm-phase reactions — death, hydrolysis, **pathogen inactivation**,
and **bacterivory** — fire, and combined with `sand_pathogen = 0.1` and
`water_factor = 1e-3` it covers the two pathogen-specific solver knobs
(`Particle.sand_attachment_factor`, `Reaction.efficiency_flowing`).

Regenerate the reference (needs MATLAB and the `slow-sand-filtration` sibling
checkout):

```
matlab -batch "run('julia/test/golden/run_export_pathogen.m')"
```

Result at the committed **10-step** spec (`TimeStep=1e-6`, `SimulationTime=1e-5`):
**match** within the standard `atol=1e-7` / `rtol=1e-6` (worst ~2.5e-8 abs). The
residual is accumulated cross-implementation floating point, confirmed by two
observations: it is exactly **0 at t = 0**, and it scales **super-linearly** with
step count (≈1.1e-8 abs over 10 steps → ≈1.6e-6 over 100), i.e. last-bit
differences in the stiff Cahn-Hilliard + osmosis sparse solve (`zeta_0 = 1`,
`τ = 1e-3`) compounding through the dynamics — not a systematic algorithmic
difference. Unlike the adaptive master (Julia vs MPC-SSF, same code lineage, so
~1e-15), this master crosses to the *independent* publication implementation, so
some per-step divergence is expected; the short horizon keeps it within tol.

### Light-model reconciliation

The authoritative code applies a **dark-respiration floor** to the light factor:
`I = max(fdark, I_eff·e^{1−I_eff})` (`@SDfilter/run_biofilm.m:272`,
`run_pathogen.m`), with a global `fdark = dark_respiration`. MPC-SSF's
`src/@State/simulate.m:237` and (mirroring it) the original Julia port instead
used an **additive** form `max(0, minimumLight + I_eff·e^{1−I_eff})`, which
double-counts the baseline at high light. Since slow-sand-filtration takes
precedence, the port was switched to the floor form (`_light_factor_floor` in
`simulate.jl`); the per-reaction `minimum_light_factor` stands in for the global
`fdark` (only the light-dependent reaction carries it). The two forms coincide
whenever `fdark = 0` **or** light `= 0`, which is why every prior golden master
(all run in darkness) was unaffected and still matches.

Two things validate the reconciliation:

- **Unit test** `light-factor floor (dark respiration)` (runtests.jl) checks
  `_light_factor_floor` reproduces `max(fdark, ·)` on both branches and differs
  from the additive form. This is the *discriminating* test — see next point for
  why a golden master cannot be.
- **Light-active golden** `reference_pathogen_light/` runs `run_pathogen.m` with
  constant light (`0.8`), `dark_respiration = 0.1`, and a phototroph-heavy seed
  (both `max` branches straddle `fdark` across depth), exported by the same
  `run_export_pathogen.m`. Julia matches it within `atol = 1e-6` (worst ~3e-7).
  **This golden does *not* isolate the light form**: the ~3e-7 residual is
  reaction-driven cross-implementation float from the stiff CH+osmosis solve, and
  it is *larger* than the light-form difference itself (the additive form scores
  the same worst-abs), so it masks it. The golden is therefore an end-to-end
  "light-active run reproduces the authoritative code" check; the unit test is
  what proves the floor form is the correct one.
