# TODO

## MATLAB project: reach parity with publication-ready version (`old-src/`)

The current `src/` is a clean-API rollback of `old-src/`. Goal: incrementally
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
- [ ] Cross-check refactored `src/` against `old-src/` outputs on shared presets
      before considering a feature "done".

## Julia port (branch `julia-port`, local only)

See `.claude/julia-port-outline.md` for the full plan. Sequencing note: the port
tracks the MATLAB parity work — port the **current `src/` feature set first**
(fixed-step biofilm model), then adopt each publication feature in Julia **after**
it lands in MATLAB, so both codebases advance in lockstep and every Julia stage
has a MATLAB reference to validate against.

- [ ] Scaffold package (`Project.toml`, `src/`, `test/`, `examples/`).
- [ ] Port ecological layer → cohesion → SandFilter → Model → State → simulate →
      Results (fixed-step parity with current `src/`).
- [ ] Golden-master tests vs. MATLAB reference outputs.
- [ ] Add adaptive time-stepping (after MATLAB).
- [ ] Add pathogen model (after MATLAB).
- [ ] Port richer results/plotting.
