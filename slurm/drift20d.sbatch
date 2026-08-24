#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J drift20
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/drift20-%j.out
#SBATCH -e slurm/logs/drift20-%j.err
#
# 20 d free-running Solver A drift, both schemes. At 3 d the drift PEAKS creep
# upward rather than saturating (+2.18e-03/d shin, +3.20e-03/d bailo), and a
# linear extrapolation to 20 d gives 37% / 52% of max phi_b -- which would sink
# the component-elimination design. But that projects 6.7x beyond the measured
# span from four peaks, and the filter is still developing at 3 d (phi_b = 0.145
# and rising). This measures it instead of extrapolating.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  plotSolverADrift(NCells=100, Kappa=1e-6, Days=20)"
