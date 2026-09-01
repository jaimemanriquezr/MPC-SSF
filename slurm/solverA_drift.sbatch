#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J drftSch
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH --array=0-2
#SBATCH -o slurm/logs/drift_scheme_%a.out
#SBATCH -e slurm/logs/drift_scheme_%a.err
#
# Free-running Solver A drift for all THREE cohesion schemes, one per array task
# so wall clock is the slowest scheme rather than their sum.
#
# The question: does "matched" reduce the Solver A/B drift relative to "shin",
# and has "bailo" already absorbed that advantage? Bailo's u-row flux is
# M(u_i,u_j) = zeta_0 (u_up)^+ (1-u_down)^+ times v, and since the 2026-08-24
# v_b fix Solver B applies exactly u_up * zeta_0 (1-u_down) * grad(mu) -- the
# same operator. If bailo's drift is <= matched's, matched has no unique virtue
# and the scheme set collapses to two.
#
# NOTE: this arm is only meaningful with the 2026-08-24 fix to simulate.m's
# free-running diagnostic. Before it, "matched" fell through to `lhsPar`, which
# is built from the SHIN stencil, so a matched arm measured shin.
module load matlab/2025b
SCHEMES=(shin matched bailo)
S=${SCHEMES[$SLURM_ARRAY_TASK_ID]}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeSolverADrift(\"$S\", NCells=100, Kappa=1e-6, Days=20)"
