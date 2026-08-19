#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J msE3
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/msE3-%j.out
#SBATCH -e slurm/logs/msE3-%j.err
# Stage 1 of the manuscript suite at article-grade resolution: the 90 d summer
# run (E3/E1-summer/X1) which also builds the mature30_summer cache that every
# other experiment consumes. Submit stage 2 with --dependency=afterok:<this>.
# Deploy: rsync this tree to ~/MPC-SSF-manuscript over the ControlMaster socket.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  manuscriptExperiments('E3', NCells=500, MaxDt=1e-5)"
