#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J impdisp
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/impdisp-%A_%a.out
#SBATCH -e slurm/logs/impdisp-%A_%a.err
#
# Implicit dispersion, measured in the SAME regime Phase 0 used (3 d), not a
# short startup run.
#
# A 0.02 d run at N=500 gave a 1.00x step reduction and reported enclosed L
# binding with w_s = 0.81 -- but that is a startup artefact, not the production
# regime: w_s in the enclosed-liquid region carries
#   L_e = -sum(sigma .* mu .* X ./ (localEnclosedS + K))
# which blows up while the enclosed liquids are still near zero. At 3 d Phase 0
# measured flowing P/L binding at 0.96 with dispersion at 94% of that sum. The
# claim that "dispersion is regime-independent so a short run is legitimate" was
# wrong: dispersion is, but the terms it competes against are not.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) N=200 ;;
  2) N=500 ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeImplicitDispersion(NCells=$N, Kappa=1e-6, Days=3)"
