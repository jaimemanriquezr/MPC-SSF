# TODO

## MATLAB project: reach parity with publication-ready version (`.old-src/`)

The current `src/` is a clean-API rollback of `.old-src/`. Goal: incrementally
re-introduce the publication features on the new class-based design.
(Details: `.claude/old-src-notes.md`.)

### Feature parity backlog
- [ ] **Adaptive time-stepping (CFL).** Port `@SDmodel/CFL_bound.m` logic onto
      `Model`; add `TimeStep="adaptive"` path to `@State/simulate.m` with options
      `CFLFactor`, `AdaptiveVelocityFactor`, `AdaptiveTimeTolerance`,
      `AdaptiveInitialDt`, `AdaptiveMaxDt`. (Currently fixed `TimeStep=1e-5`.)
- [ ] **Pathogen model.** Add pathogen particle component + dispatch (biofilm vs.
      pathogen run), sand-vs-biofilm attachment (`sand_pathogen`), bacterivory /
      predation (`kPred`), and flowing-phase `water_factor`.
- [x] **Richer results & plotting.** DONE: `get_*` accessors on `Results`
      (`reaction_rates`/`reaction_names`, `get_volume_fractions`, `times`,
      `depths`, `concentration`) plus six Makie plots via a package extension
      (`ext/MPCSSFMakieExt.jl`): `plot_concentration`, `plot_concentration_heatmap`,
      `plot_volume_fractions`, `plot_velocity`, `plot_cfl`, `plot_reaction_rates`.
      Core stays Makie-free; plots load on `using CairoMakie`. Accessors and the
      no-backend fallback path are covered in the test suite.
- [x] **Run chaining.** `concatenate`/`plus` to resume/append simulations.
- [ ] **HPC/batch tooling.** Decide whether to re-add `slurm/` + `simulate_filter`
      entry point.
- [ ] **Diagnostics.** `print_*` and `*_error` helpers as needed.

### Verification
- [ ] Cross-check refactored `src/` against `.old-src/` outputs on shared presets
      before considering a feature "done".

## Julia port (branch `julia-port`, local only)

See `.claude/julia-port-outline.md` for the full plan. Sequencing note: the port
tracks the MATLAB parity work — port the **current `src/` feature set first**
(fixed-step biofilm model), then adopt each publication feature in Julia **after**
it lands in MATLAB, so both codebases advance in lockstep and every Julia stage
has a MATLAB reference to validate against.

- [x] Scaffold package (`Project.toml`, `src/`, `test/`, `examples/`).
- Port (fixed-step parity with current `src/`), leaf-first:
  - [x] Ecological layer — Component, Particle, Liquid, Reaction
        (`compute_rate`, `lookup_{order,half_saturation_constants,
        stoichiometric_coefficients,quotients}`) + unit tests. NOT yet run in a
        Julia runtime (none in sandbox); user to run `Pkg.test`.
  - [x] Cohesion (CahnHilliardModel) — kappa/zeta_0/zeta_1 + mobility &
        potential-gradient callables, with tests.
  - [x] SandFilter — grid (`addgridpoints`/`gridsize`/`gridzero`), porosity,
        light attenuation (water/sand), and Cahn-Hilliard matrix assembly
        (`get_cahn_hilliard_matrices` → `Triplets`), with tests. Not yet run in
        a Julia runtime.
  - [x] Model — component views (particles/liquids), reaction-rate vector,
        stoichiometric/half-saturation/order matrices, quotients, global params
        (β, τ, water density, detachment), with tests.
  - [x] State — container (zero-allocated regions, velocities, enclosed water)
        + dependent accessors (biofilm/flowing concentrations, volume
        fractions), with tests.
  - [x] simulate (core fixed-step solver) — full block-by-block port: reaction
        kernel (listK/order), light, attachment/detachment/transfer, eco
        reactions, SOLVER A (Cahn-Hilliard implicit solve → biofilm velocity),
        SOLVER B (upwind convection + dispersion FV update), clogging/negativity
        guards, frame capture. Smoke test (zero-rate model) passes structurally;
        NOT yet validated numerically vs MATLAB (needs golden-master data).
  - [x] Results — frame container + accessors (`concentration`,
        `get_volume_fractions`, `times`, `depths`); plotting deferred. With tests.
- [x] Golden-master tests vs. MATLAB reference outputs. DONE (branch
      `matlab-claude`): `julia/test/golden/` runs an identical SimpleModel
      simulation in MATLAB (`export_reference.m`) and Julia (`compare.jl`) and
      diffs every field. Ran in a real Julia runtime — flags agree (OK) and all
      fields match (non-flowing ~1e-15, flowing ~9e-9 abs / ~6.5e-8 rel, the
      latter accumulated cross-language float, not a bug). Committed reference
      under `test/golden/reference/`; runs by default in `Pkg.test`. Also ported
      `modelLund` and `simpleModel` presets. (SimpleModel is the validation
      vehicle because clean modelLund trips a negativity guard at step 1.)
- [x] Add adaptive time-stepping. DONE: `time_step=:adaptive` in simulate.jl
      ports the MATLAB CFL path (options cfl_factor, adaptive_velocity_factor,
      adaptive_time_tolerance, adaptive_initial_dt, adaptive_max_dt). Validated:
      capped adaptive == fixed step bit-for-bit; Julia CFL plateau dt 9.3887e-8
      matches MATLAB 9.389e-8. Assumes Lund model structure, as the reference
      does. (MATLAB side is on branch `matlab-claude`.)
- [x] Adaptive-CFL golden-master. DONE: `export_adaptive_reference.m` /
      `compare_adaptive.jl` diff the per-step dt trajectory + frames of adaptive
      modelLund vs MATLAB. Exact match — same 228 steps, all fields ~1e-15.
      Reference under `test/golden/reference_adaptive/`; runs in `Pkg.test`.
- [x] Run chaining (Julia). DONE: `final_state(results)` rebuilds a resumable
      State from a frame (port of `@SDresults/get_filter_state`), and
      `concatenate(r1,r2)` / `r1 + r2` join runs. Validated: resuming reproduces
      a continuous run bit-for-bit. (`src/chaining.jl`.)
- [ ] Add pathogen model (after MATLAB).
- [x] Port richer results/plotting. DONE: Makie plots via package extension
      (`ext/MPCSSFMakieExt.jl`) + `Results` accessors; see the plotting item above.
