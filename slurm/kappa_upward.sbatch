#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J kupward
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH -a 1-2
#SBATCH -o slurm/logs/kupward-%A_%a.out
#SBATCH -e slurm/logs/kupward-%A_%a.err
#
# Is the upward biofilm migration caused by the kappa fix, or pre-existing?
#
# E3 at 20 d shows bed biomass PEAKING at day 10 then falling 40%, while the
# supernatant triples; by day 20, 57% of the biofilm is above the sand, reaching
# 34 cm up with phi_b = 0.14 at 30 cm. Not a schmutzdecke.
#
# kappa = 0 reproduces the PRE-FIX scheme exactly (the fix restored only the
# interfacial term), so this isolates it. Note the Psi' part of the cohesive flux
# was always active and acts as POSITIVE diffusion wherever Psi'' = 3u(u-zeta_1)
# > 0, i.e. everywhere above u = zeta_1 -- so the prior expectation is that the
# spread is pre-existing and kappa, which sharpens interfaces, opposes it.
module load matlab/2025b
case $SLURM_ARRAY_TASK_ID in
  1) KAP=1e-6 ;;
  2) KAP=0    ;;
esac
matlab -batch "addpath(genpath('src')); addpath('analysis'); \
  manuscriptExperiments('E3', Days=20, NCells=200, Kappa=$KAP, \
                        OutRoot='results/kupward_k$SLURM_ARRAY_TASK_ID')"
