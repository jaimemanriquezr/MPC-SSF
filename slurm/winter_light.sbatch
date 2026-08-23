#!/bin/bash -l
# Realistic Scandinavian winter photoperiod, published configuration.
#
# The manuscript winter curve gives 13.54 h of daylight and 0.63x summer's daily
# light integral. Lund (55.7 N) on 21 December has 6 h 55 m and, by solar
# elevation, ~0.22x the peak irradiance -- so the manuscript understates the
# seasonal light contrast ~6x. Offset 0.81 matches BOTH measures independently.
#
# Only the WINTER arm is new: the summer curve is already realistic (16.92 h
# against an actual 17.5 h) and legmatch_summer_published_n100 is already done,
# so it serves as the paired reference.
#
#   sbatch --array=1-2 slurm/winter_light.sbatch
#     1  winter, published, realistic curve   <- the new arm
#     2  winter, plusResp,  realistic curve   <- checks the interaction
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 08:00:00
#SBATCH -J wlight
#SBATCH -o slurm/logs/wlight_%A_%a.out
#SBATCH -e slurm/logs/wlight_%A_%a.err
set -euo pipefail
mkdir -p slurm/logs
VARIANTS=(x published plusResp)
V=${VARIANTS[$SLURM_ARRAY_TASK_ID]}
echo "=== winter_light task=$SLURM_ARRAY_TASK_ID variant=$V start=$(date -Is) ==="
module load matlab/2025b
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
               probeLegacyMatch(\"winter\", Variant=\"$V\", WinterLight=\"realistic\")"
echo "=== done $(date -Is) ==="
