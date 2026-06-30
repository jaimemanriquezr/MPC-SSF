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
- [ ] **Richer results & plotting.** Port useful `@SDresults` plots
      (`plot_CFL`, `plot_biofilm_*`, `plot_concentration_3D/line`,
      `plot_subphases`, `plot_velocity`) and `get_*` accessors into `@Results`.
- [ ] **Run chaining.** `concatenate`/`plus` to resume/append simulations.
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
- [ ] Golden-master tests vs. MATLAB reference outputs (the remaining gap:
      everything is structurally tested but NOT yet validated numerically
      against MATLAB — needs reference data exported from a MATLAB run).
- [ ] Add adaptive time-stepping (after MATLAB).
- [ ] Add pathogen model (after MATLAB).
- [ ] Port richer results/plotting.
