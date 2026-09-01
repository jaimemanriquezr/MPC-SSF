#!/bin/bash -l
# Missing cases for the Manriquez2026 figure RECREATION (2026-08-28):
#   task 0  fld2x_winter  90 d at 3 C with the manuscript winter light curve
#           -> the winter arm of fig:seasons-results (summer arm = E7 lit at 90 d)
#   task 1  fld2x_cov01   30 d at LightScale = 0.01 (the manuscript's "covered")
#           -> the covered arm of fig:roofed-results (uncovered = E4 lit leg 3;
#              the existing fld2x_dark is 0 % light, a stronger contrast than published)
#
# Both otherwise verbatim the E4 working set. Requires the probeChain of
# 2026-08-28 (options Temperature, LightForm). Submit FROM THE REPO ROOT.
#
# Wall: 90 d ~ 4.8 h at the 193 s/d cosmos rate, and 3 C runs faster (P1 cold
# arms were 2.2x faster); 30 d ~ 1.6 h. 12 h covers everything.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 12:00:00
#SBATCH -J winterCov
#SBATCH -o slurm/logs/winter_cov_%A_%a.out
#SBATCH -e slurm/logs/winter_cov_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=fld2x_winter ; DAYS=90 ; TEMP=3  ; LSCALE=1    ; LFORM=winter ;;
    1) TAG=fld2x_cov01  ; DAYS=30 ; TEMP=19 ; LSCALE=0.01 ; LFORM=chain  ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== winter_cov tag=$TAG task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=$DAYS, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=$LSCALE, LightForm=\"$LFORM\", Temperature=$TEMP, \
               KDOM=3e-4, TransferScale=10, KHPO4=1e-6, DetachForm=\"linear\", \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
