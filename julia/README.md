# SSF.jl

Julia port of the MATLAB MPC-SSF code — multiphase continuum simulation of slow
sand filtration (Diehl et al.). Lives on the local `julia-port` branch only.

## Status

Scaffold. Structs and the public API surface are stubbed; numerical methods are
not yet implemented. **Phase 1** targets parity with the current MATLAB `src/`
(fixed-step biofilm model); adaptive time-stepping and the pathogen model follow
in later phases, tracking the MATLAB parity work.

Planning docs (untracked, in `../.claude/`): `julia-port-outline.md`,
`old-src-notes.md`. Task tracking: `../TODO.md`.

## Layout

    src/
      SSF.jl              top module (include order = port order)
      ecological/            Component, Particle, Liquid, Reaction
      cohesion/              CahnHilliardModel
      SandFilter.jl          geometry, grid, porosity, light, CH matrices
      Model.jl               assembled ecological model
      State.jl               state container + initial conditions
      simulate.jl            core solver (block-by-block port of simulate.m)
      Results.jl             post-processing
    test/runtests.jl
    examples/simulation_example.jl

## Develop

```julia
julia> using Pkg; Pkg.activate("julia"); Pkg.instantiate()
julia> Pkg.test()          # tests the active project (SSF)
```

## Porting workflow

Leaf-first: ecological → cohesion → SandFilter → Model → State → simulate →
Results. Validate each unit against MATLAB reference output before moving on.
