#!/bin/bash -l
# Manufactured-solution verification of the restored Cahn-Hilliard kappa operator.
#
#   sbatch slurm/mms_ch.sbatch        (single task; assembly only, no time stepping)
#
# Verifies block (2,1) of the mixed system converges at second order against a
# manufactured u with zero normal derivative at both faces of the CH domain -- the
# homogeneous Neumann condition the index clamping imposes. Reports interior and
# full-domain error norms separately so a degraded boundary treatment would show.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:30:00
#SBATCH -J mmsCH
#SBATCH -o slurm/logs/mmsch_%j.out
#SBATCH -e slurm/logs/mmsch_%j.err
set -euo pipefail
mkdir -p slurm/logs
echo "=== mms_ch start=$(date -Is) host=$(hostname) ==="
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
  probeMMSCahnHilliard(NList=[20 40 80 160 320 640 1280 2560], Kappa=1e-7); \
  fprintf('\n--- repeat at the golden kappa ---\n'); \
  probeMMSCahnHilliard(NList=[20 40 80 160 320 640], Kappa=1e-2, Save=false);"
echo "=== done $(date -Is) ==="
