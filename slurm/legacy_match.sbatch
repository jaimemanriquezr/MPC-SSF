#!/bin/bash -l
# Does MPC-SSF reproduce the PUBLISHED seasons figure?
#
# Three variants x two seasons, decomposing the divergence between
# probeSeason90("*","manuscript") and the legacy run that made Fig. 7:
#   published          5 reactions, dark floor 0.01, I_opt = 1.08   <- the published config
#   plusResp           + PG-excess phototroph respiration 0.55/d    <- isolates the added reaction
#   publishedIopt1814  5 reactions, I_opt = 1.814e-2 as printed     <- isolates the table value
#
#   sbatch --array=1-6 slurm/legacy_match.sbatch
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 06:00:00
#SBATCH -J legmatch
#SBATCH -o slurm/logs/legmatch_%A_%a.out
#SBATCH -e slurm/logs/legmatch_%A_%a.err
set -euo pipefail
mkdir -p slurm/logs
SEASONS=(x summer winter summer winter summer winter)
VARIANTS=(x published published plusResp plusResp publishedIopt1814 publishedIopt1814)
S=${SEASONS[$SLURM_ARRAY_TASK_ID]}
V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "=== legacy_match task=$SLURM_ARRAY_TASK_ID season=$S variant=$V start=$(date -Is) ==="
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
               probeLegacyMatch(\"$S\", Variant=\"$V\")"
echo "=== done $(date -Is) ==="
