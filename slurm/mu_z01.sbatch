#!/bin/bash -l
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -J muz01
#SBATCH -t 12:00:00
#SBATCH -c 1
#SBATCH --array=0-1
#SBATCH -o slurm/logs/muz01_%a.out
#SBATCH -e slurm/logs/muz01_%a.err
#
# mu_PHO = 20 /d (x10 the preset) on the WORKING SET, zeta0 = 1 + linear detachment.
# Re-runs E5 (which predated the zeta0=1 decision and ran at zeta0=5 with sqrt detachment)
# so the lit/dark ratio is comparable with E4. x20 omitted: E5 showed 40/d identical to 20/d.
module load matlab/2025b
LIGHT=(1 0); TAG=(mu_z01_x10_lit mu_z01_x10_dark)
i=$SLURM_ARRAY_TASK_ID
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeChain(\"${TAG[$i]}\", Zeta0=1, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10, \
             FramesPerLeg=241, LightScale=${LIGHT[$i]}, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
             DetachForm=\"linear\", MuPHO=20.0, \
             Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"
