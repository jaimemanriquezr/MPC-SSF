#!/bin/bash -l
# 2x2 attribution: {mu_HET 2.0, 4.8} x {nu_P 0.094, 52}, lit, 30 d, N = 500.
#
#   sbatch slurm/tenore_2x2.sbatch          # array 0..3
#
# WHY. Commit b1fcd03 changed two presets at once (Tenore2021: mu_HET 2.0 -> 4.8,
# nu_P 0.094 -> 52). The obvious baseline, chain_fld2x_lit (grown 2026-08-26), is
# NOT usable for attribution: ee0e4d3 rewrote 523 lines of simulate.m plus both
# presets and was committed 2026-09-01, but the work was uncommitted when the
# chain ran, so git cannot establish which code produced it. This array instead
# runs all four corners under TODAY'S code, so the only differences are the two
# parameters and the design is self-contained.
#
#   arm 0  mu 2.0  nu 0.094   old / old   <- the control the fld2x chain cannot be
#   arm 1  mu 4.8  nu 0.094   mu_HET only
#   arm 2  mu 2.0  nu 52      nu_P only
#   arm 3  mu 4.8  nu 52      both (= the shipped presets)
#
# arm3 - arm1 and arm2 - arm0 both isolate nu_P; arm1 - arm0 and arm3 - arm2 both
# isolate mu_HET. Disagreement between the two estimates of either effect IS the
# interaction, which is the point of running the full square rather than a fork.
#
# LIT ONLY. nu_P is inert without light, so a dark arm cannot see it; mu_HET is
# recovered from the lit contrasts above.
#
# Supersedes cosmos 3564903 (slurm/tenore_chain.sbatch), whose lit arm duplicates
# arm 3 and which was submitted before the confound was found.
#
# Wall: ~105 s per simulated day at N = 500 locally, cosmos ~1.8x (E6) -> ~95 min
# per arm for 30 d. 4 h covers an arm that clogs. Each leg writes its own file.
#
# Output: analysis/probes/data/chain/chain_t22_{m20n009,m48n009,m20n52,m48n52}_leg{1,2,3}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -J tenore2x2
#SBATCH -o slurm/logs/tenore_2x2_%A_%a.out
#SBATCH -e slurm/logs/tenore_2x2_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=t22_m20n009; MU=2.0; NU=0.094 ;;
    1) TAG=t22_m48n009; MU=4.8; NU=0.094 ;;
    2) TAG=t22_m20n52;  MU=2.0; NU=52    ;;
    3) TAG=t22_m48n52;  MU=4.8; NU=52    ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== tenore_2x2 tag=$TAG mu_HET=$MU nu_P=$NU task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_${TAG}_leg*.mat 2>/dev/null || echo "no legs for $TAG -- starting from t = 0 (expected)"

module load matlab/2025b

# E4 call (slurm/fld2x_104d.sbatch) with Days 104 -> 30, plus the two swept options.
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=30, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=1, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
               DetachForm=\"linear\", MuHET=$MU, Attenuation=$NU, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
