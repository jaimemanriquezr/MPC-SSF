#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J bcgate
#SBATCH -t 06:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/bc_selftest.out
#SBATCH -e slurm/logs/bc_selftest.err
#
# GATE for slurm/refinement.sbatch. Do not submit the refinement array until this
# reports PASS: a wrong sign on the Dirichlet terms would make every dirichlet arm
# quietly wrong, and a 24 h array is an expensive way to find that out.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeBCSelfTest(NCells=100, Kappa=1e-6, Days=3)"
