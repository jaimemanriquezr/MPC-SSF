#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J val500
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/val500-%j.out
#SBATCH -e slurm/logs/val500-%j.err
# Three-way VALUE comparison at N=500 with the optimised Bailo. A and B's cost is
# already known here (2158-2289 s, 460-467 s) but their FIELDS were never stored,
# so all three are re-run under one harness to make the profiles comparable.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  compareSchemeValues(NCells=500, Kappa=1e-6, Days=3, WithBailo=true)"
