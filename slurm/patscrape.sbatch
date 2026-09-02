#!/bin/bash -l
# Pathogen removal after scraping, for fig:pat-scraping (Manriquez2026 revision).
# Replaces the April figure `results_flow_Scraped_GP-pathogen_filtration.pdf`,
# which predates the pathogen-model audit (decision 2026-08-27).
#
#   sbatch slurm/patscrape.sbatch        # array 0..3 = 0, 4, 8, 12 cm
#
# Each arm, in ONE MATLAB process:
#   1. scrape the end state of leg 3 of chain_fld2x_lit (t = 30 d, E4 working
#      set) to DEPTH, then run DAYS = 20 d of regrowth, saving a probeChain-
#      format restart snapshot (reproduceManriquez2026 SnapshotTag option);
#   2. probePulse on that snapshot with a CONSTANT feed of marker at the
#      Table B.1 nominal CRef (PulseFactor = 1 over the whole window) -- this is
#      what the fig:pat-scraping caption describes ("a constant feed of marker");
#   3. probePulse again with the 10x spike on [0.1, 0.3) d, the same disturbance
#      logOatCampaign's `pulse` scenario uses;
#   4. probePulse a third time at 100x (the BigPulse factor) on the same window.
#
# THREE feeds per arm, per Jaime 2026-09-02: the two pulse magnitudes compare how
# challenge magnitude affects removal, the constant feed shows the effect of feed
# LENGTH. The caption of fig:pat-scraping currently describes the constant feed
# only. The 20 d regrowth is the expensive part and is shared across all three,
# so each extra feed costs only ~25 min.
#
# ORDER MATTERS. The feeds run cheapest-risk first. checkRunFlag errors on a
# non-OK flag and `set -e` would kill the task, so the 100x arm is LAST: if it
# clogs, the constant-feed and 10x results are already on disk. (probePulse also
# writes an ABORTED file in that case -- the x100 arm of job 3549870 was lost
# before that existed.)
#
# Prerequisite on cosmos: analysis/probes/data/chain/chain_fld2x_lit_leg3.mat.
#
# Wall estimate: ~190 s per simulated day at N = 500 on cosmos (E6) -> ~65 min
# for the 20 d regrowth, plus ~25 min per 3 d feed x3; ~2.5 h per arm. 6 h is margin
# for an arm that clogs and drives the adaptive step down.
#
# Output:
#   analysis/probes/data/chain/scrape_GP{0,4,8,12}_20d.mat      (snapshots)
#   analysis/probes/data/pulse/pulse_scrapeGP{...}_feed.mat     (constant feed)
#   analysis/probes/data/pulse/pulse_scrapeGP{...}_pulse10.mat  (10x spike)
#   analysis/probes/data/pulse/pulse_scrapeGP{...}_pulse100.mat (100x spike)
#   analysis/results/figures/repro/data/scrape_GP{...}.mat       (regrowth, 20 d)
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -J patscrape
#SBATCH -o slurm/logs/patscrape_%A_%a.out
#SBATCH -e slurm/logs/patscrape_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) DEPTH=0.00; CM=0  ;;
    1) DEPTH=0.04; CM=4  ;;
    2) DEPTH=0.08; CM=8  ;;
    3) DEPTH=0.12; CM=12 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

DAYS=${DAYS:-20}
SNAP="scrape_GP${CM}_${DAYS}d"

echo "=== patscrape depth=${DEPTH} m days=${DAYS} task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_fld2x_lit_leg3.mat

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    reproduceManriquez2026(\"scrape\", Depth=$DEPTH, Days=$DAYS, FromLeg=3, NCells=500, \
        SnapshotTag=\"$SNAP\"); \
    probePulse(\"scrapeGP${CM}_feed\", Snapshot=\"$SNAP.mat\", NCells=500, \
        PulseFactor=1, PulseT0=0.0, PulseT1=3.0, TPost=3.0); \
    probePulse(\"scrapeGP${CM}_pulse10\", Snapshot=\"$SNAP.mat\", NCells=500, \
        PulseFactor=10, PulseT0=0.1, PulseT1=0.3, TPost=3.0); \
    probePulse(\"scrapeGP${CM}_pulse100\", Snapshot=\"$SNAP.mat\", NCells=500, \
        PulseFactor=100, PulseT0=0.1, PulseT1=0.3, TPost=3.0)"

echo "=== done $(date -Is) ==="
