#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J refine
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH --array=0-7
#SBATCH -o slurm/logs/refine_%a.out
#SBATCH -e slurm/logs/refine_%a.err
#
# Grid refinement x boundary condition. Tasks 0-3 = Neumann (current, manuscript
# Appendix A.1), tasks 4-7 = Dirichlet (phi_b continuous at z = 0, ghost taken from
# the first bed cell, which follows its own ODE -- the Diehl2025 condition).
#
# N = 100 is cap-bound (capFrac 0.96) and cheap; N >= 500 is CFL-bound with dt
# ~ dz^2, so cost ~ N^3 from there. N=1000 is the long pole; the array runs the
# arms concurrently so wall clock is that arm alone.
#
# Read zeroCellRatio down each BC column. Flat => no surface source, refinement is
# meaningful. Growing like 1/dz => phi_b(0) has no continuum limit and only Sigma,
# L and the deposition flux may be reported.
module load matlab/2025b
NS=(100 200 500 1000 100 200 500 1000)
BCS=(neumann neumann neumann neumann dirichlet dirichlet dirichlet dirichlet)
N=${NS[$SLURM_ARRAY_TASK_ID]}
BC=${BCS[$SLURM_ARRAY_TASK_ID]}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeRefinement($N, \"$BC\", Kappa=1e-6, Days=3)"
