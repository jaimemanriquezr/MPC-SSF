#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J k6E3
#SBATCH -t 72:00:00
#SBATCH -c 1
#SBATCH -o slurm/logs/k6E3-%j.out
#SBATCH -e slurm/logs/k6E3-%j.err
#
# E3 (+ E1 winter) at 20 d, kappa = 1e-6, N = 500 -- the first run of these
# experiments with the Cahn-Hilliard interfacial term actually live. Until the
# 2026-08-23 fix kappa was inert (block (1,1), discarded), so every previously
# reported result was produced with no cohesion at all.
#
# Why kappa=1e-6 at N=500: the interface width is ell = sqrt(kappa/Psi''), not
# sqrt(kappa); with Psi'' = 3*phi*(phi - zeta_1) ~ 0.26 at phi = 0.3 that gives
# ell = 1.96 mm, and dz = 1/500.5 = 2.00e-3 m, i.e. dz/ell = 1.02 -- resolved.
# 1e-6 is also what modelLund and modelPathogen already carry. The old 1e-7
# default rested on an instability that kappa demonstrably could not have
# caused (see the note in manuscriptExperiments.m).
#
# Runs 'E1', which is runSummer (all E3 outputs) followed by e1_winter. Summer
# writes an E3_longterm/done.txt sentinel, so RESUBMITTING THIS SCRIPT skips a
# completed summer and goes straight to winter -- that sentinel is the only
# checkpoint, hence the paired submission in submit_k1e6_20d.sh.
#
# Deploy: rsync this tree to $HOME/MPC-SSF-manuscript over the ControlMaster socket.
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  manuscriptExperiments('E1', Days=20, NCells=500, Kappa=1e-6, \
                        OutRoot='results/manuscript_k1e6_20d')"
