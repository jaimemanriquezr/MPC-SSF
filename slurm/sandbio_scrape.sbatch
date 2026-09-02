#!/bin/bash -l
# Stage 2: scrape and challenge the stage-1 chains. Submit with
#   sbatch --dependency=afterok:<stage1 jobid> slurm/sandbio_scrape.sbatch
#
# Each arm builds its OWN uniquely named scraped snapshot, so the arms cannot race.
# The model must carry the SAME SandBiomass that grew the chain -- probePulse rebuilds
# the model, so a mismatch would be a mongrel.
#
# dL = L(intact) - L(scraped) per factor, against Schijven2013's 0.6-1.6 log.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 2:00:00
#SBATCH -J sandbio2
#SBATCH -o slurm/logs/sandbio_scrape_%A_%a.out
#SBATCH -e slurm/logs/sandbio_scrape_%A_%a.err
#SBATCH --array=0-5

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=sb01; SB=0.1; MODE=intact  ;;
    1) TAG=sb03; SB=0.3; MODE=intact  ;;
    2) TAG=sb10; SB=1.0; MODE=intact  ;;
    3) TAG=sb01; SB=0.1; MODE=scraped ;;
    4) TAG=sb03; SB=0.3; MODE=scraped ;;
    5) TAG=sb10; SB=1.0; MODE=scraped ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== sandbio_scrape tag=$TAG SandBiomass=$SB mode=$MODE host=$(hostname) start=$(date -Is) ==="
module load matlab/2025b

if [ "$MODE" = "scraped" ]; then
    SNAP="scrape2cm_${TAG}.mat"
    matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
        reproduceManriquez2026(\"scrape\", Depth=0.02, Days=0, FromLeg=2, NCells=500, \
            LitTag=\"$TAG\", SnapshotTag=\"scrape2cm_${TAG}\")"
else
    SNAP="chain_${TAG}_leg2.mat"
fi

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probePulse(\"sb_${MODE}_${TAG}\", Snapshot=\"$SNAP\", NCells=500, SandBiomass=$SB, \
        PulseFactor=1, PulseT0=0.0, PulseT1=3.0, TPost=3.0)"

echo "=== done $(date -Is) ==="
