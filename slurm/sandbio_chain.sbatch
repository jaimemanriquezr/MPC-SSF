#!/bin/bash -l
# Stage 1: 20 d chains at three BIOMASS sand-attachment factors.
#
#   sbatch slurm/sandbio_chain.sbatch          # array 0-2
#
# E13 showed the model cannot reproduce Schijven2013's 0.6-1.6 log scraping loss at
# any MARKER sand factor -- dL is 0.02 log at best, because removal is spread over
# 2-25 cm and the top 2 cm carries only 18 % of it. The marker cannot fix that: it
# does not build the filter. HET and PHO do, and both currently attach to bare sand
# at full strength.
#
# HYPOTHESIS. Lowering the biomass sand factor makes biomass accrete preferentially
# on existing biofilm rather than on clean grains -- a positive feedback that should
# thicken the Schmutzdecke, thin the deep bed, and make a 2 cm scrape cost more.
#
# The marker stays at its published s = 0 so the two axes are not confounded.
#
# Wall: ~3.2 min per simulated day at N = 500 on cosmos -> ~65 min for 20 d. 3 h
# covers an arm that clogs; a lower sand factor concentrates biomass, which raises
# the clogging risk, so this is not a formality.
#
# Output: analysis/probes/data/chain/chain_sb{01,03,10}_leg{1,2}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 3:00:00
#SBATCH -J sandbio1
#SBATCH -o slurm/logs/sandbio_chain_%A_%a.out
#SBATCH -e slurm/logs/sandbio_chain_%A_%a.err
#SBATCH --array=0-2

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=sb01; SB=0.1 ;;
    1) TAG=sb03; SB=0.3 ;;
    2) TAG=sb10; SB=1.0 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== sandbio_chain tag=$TAG SandBiomass=$SB host=$(hostname) start=$(date -Is) ==="
module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=20, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
               DetachForm=\"linear\", SandBiomass=$SB, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
