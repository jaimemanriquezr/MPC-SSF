"""
    MPCSSF

Julia port of the MATLAB multiphase continuum model of slow sand filtration
(MPC-SSF, Diehl et al.).

Phase 1 goal: numerical parity with the current MATLAB `src/` — a fixed-step
biofilm model. Adaptive time-stepping and the pathogen model are later phases
that track the MATLAB parity work (see `../TODO.md` and
`../.claude/julia-port-outline.md`).

Files are included in dependency order (leaf-first), mirroring the MATLAB
classes:

    ecological/  -> SandFilter -> Model -> State -> simulate -> Results

Each unit is a stub to be filled in and validated against MATLAB reference
output before moving to the next.
"""
module MPCSSF

using LinearAlgebra
using SparseArrays

# ---- ecological layer (leaf) ------------------------------------------------
include("ecological/Component.jl")
include("ecological/Particle.jl")
include("ecological/Liquid.jl")
include("ecological/Reaction.jl")

# ---- cohesion submodel ------------------------------------------------------
include("cohesion/CahnHilliardModel.jl")

# ---- filter geometry / mesh -------------------------------------------------
include("SandFilter.jl")

# ---- assembled model --------------------------------------------------------
include("Model.jl")

# ---- state container + solver ----------------------------------------------
include("State.jl")
include("simulate.jl")

# ---- post-processing --------------------------------------------------------
include("Results.jl")

export Component, Particle, Liquid, Reaction
export CahnHilliardModel
export SandFilter, addgridpoints, computeporosity
export Model
export State, simulate
export Results

end # module MPCSSF
