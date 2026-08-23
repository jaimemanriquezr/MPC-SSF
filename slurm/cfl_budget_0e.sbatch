#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J cfl0e
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/cfl0e-%j.out
#SBATCH -e slurm/logs/cfl0e-%j.err
# Phase 0 measurement 0e: does Solver A's phi_b match the diagnostic component
# sum? Gates the component-elimination design. N=100 is enough to answer a
# consistency question; 3 d so the regime is active.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis/probes'); \
  probeCflBudget(NCells=100, Kappa=1e-6, MatureDays=3)"
