# Port of src/@State/simulate.m  (the core solver)
#
# This is the heart of the model and the largest single unit. It will be ported
# block-by-block, preserving the MATLAB section structure so it can be diffed
# against the source and validated frame-by-frame:
#
#   I.   Load parameters (filter geometry, model params, inflow handling)
#   II.  Model parameters (rates, stoichiometry, cohesion handles)
#   III. Initial conditions (volume fractions, initial velocities)
#   IV.  Pre-allocation (frame buffers, sparse CH matrices)
#   V.   Output Results scaffolding
#   VI.  Time integration loop:
#          - global/local concentrations & volume fractions
#          - clogging guard
#          - reaction terms (light, attachment/detachment, transfer, ecological)
#          - SOLVER A: Cahn-Hilliard implicit solve -> biofilm velocity vb
#          - SOLVER B: finite-volume update of concentrations (upwind fluxes,
#            dispersion, porosity-weighted faces)
#          - negative/NaN guards, frame capture
#
# Phase 1 implements the FIXED-step path only (parity with current src/).
# Adaptive (CFL) time-stepping and the pathogen variant are later phases
# (see ../TODO.md and ../.claude/old-src-notes.md).

"""
    simulate(state::State; inflow_concentrations=nothing, simulation_time=1.0,
             time_step=1e-5, n_frames=200, clogging_fraction=0.99,
             upwinded=false, quiet=false) -> Results

Run the slow-sand-filtration simulation from `state`. Keyword names mirror the
MATLAB `simulate` options. Returns a [`Results`](@ref).

NOT YET IMPLEMENTED — scaffold only.
"""
function simulate(state::State;
                  inflow_concentrations=nothing,
                  simulation_time::Real=1.0,
                  time_step=1e-5,
                  n_frames::Integer=200,
                  clogging_fraction::Real=0.99,
                  upwinded::Bool=false,
                  quiet::Bool=false)
    error("simulate not yet ported — port src/@State/simulate.m block-by-block " *
          "(see header in this file)")
end
