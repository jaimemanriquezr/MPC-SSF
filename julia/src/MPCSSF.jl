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
export compute_rate, lookup_order, lookup_half_saturation_constants,
       lookup_stoichiometric_coefficients, lookup_quotients
export CahnHilliardModel
export SandFilter, Grid, Triplets, addgridpoints, computeporosity,
       gridsize, gridzero, light_attenuation_eta_water,
       light_attenuation_eta_sand, get_cahn_hilliard_matrices
export Model, particles, liquids, compute_reaction_rates,
       stoichiometric_coefficients, stoichiometric_matrix_particles,
       stoichiometric_matrix_liquids, half_saturation_constants,
       reaction_orders, quotients
export State, simulate
export Results

end # module MPCSSF
