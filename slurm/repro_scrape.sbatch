#!/bin/bash -l
# Scraping runs for the Manriquez2026 reproduction (fig:month-scrape-0/4/8/12
# and fig:month-scrape-mass), on the AUDITED WORKING SET.
#
#   sbatch slurm/repro_scrape.sbatch        # array 0..3 = 0, 4, 8, 12 cm
#
# Each arm takes the end state of leg 3 of chain_fld2x_lit (t = 30 d, the E4
# working set), removes the biofilm down to the scraping depth with the legacy
# refill-sand rule (port of @SDresults/scrape.m, see reproduceManriquez2026.m),
# and runs the regrowth. The manuscript figure shows 20 d of regrowth; this runs
# 10 d -- half the published window, so the figure is a PARTIAL reproduction and
# is labelled as such in analysis/results/figures/repro/README.md.
#
# Prerequisite on cosmos: analysis/probes/data/chain/chain_fld2x_lit_leg3.mat
# and analysis/reproduceManriquez2026.m must be present (rsync them first).
#
# Wall estimate: ~190 s per simulated day at N = 500 on cosmos (E6), so ~32 min
# per arm for 10 d; 4 h is a wide margin. Output:
#   analysis/results/figures/repro/data/scrape_GP{0,4,8,12}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH -J reproscrape
#SBATCH -o slurm/logs/repro_scrape_%A_%a.out
#SBATCH -e slurm/logs/repro_scrape_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) DEPTH=0.00 ;;
    1) DEPTH=0.04 ;;
    2) DEPTH=0.08 ;;
    3) DEPTH=0.12 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== repro_scrape depth=${DEPTH} m task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
module load matlab/2025b

# DAYS is overridable: sbatch --export=ALL,DAYS=20 ... (default 10, the 2026-08-27 partial).
DAYS=${DAYS:-10}
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    reproduceManriquez2026(\"scrape\", Depth=$DEPTH, Days=$DAYS, FromLeg=3, NCells=500)"

echo "=== done $(date -Is) ==="
