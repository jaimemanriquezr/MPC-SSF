#!/bin/bash -l
# 104 d lit/dark pair at the working set (E4 parameters), N = 500 -- Campos2002
# Bed 9 run length. Settles the head-loss trajectory and the d97 top-2 cm biomass
# against Campos like-for-like, which the 60 d E4 run could only extrapolate.
#
#   sbatch slurm/fld2x_104d.sbatch          # array 0 = lit, 1 = dark
#
# Submit FROM THE REPO ROOT on cosmos (~/MPC-SSF-manuscript); every addpath is
# relative to it. `#!/bin/bash -l` is required or `module` is not found.
#
# RESUME, NOT A FRESH RUN. probeChain skips legs already on disk and guards every
# option against the leg-1 record, so this EXTENDS E4's chain_fld2x_{lit,dark}
# legs 1-6 (0-60 d) with legs 7-11 (60-104 d; leg 11 is the 4 d remainder). The
# call below must therefore stay IDENTICAL to the E4 call in EXPERIMENTS.md apart
# from Days. Before submitting, make sure legs 1-6 are on cosmos:
#
#   rsync -av -e "ssh -o ControlPath=$HOME/.ssh/cm-cosmos" \
#         analysis/probes/data/chain/chain_fld2x_{lit,dark}_leg{1..6}.mat \
#         cosmos:MPC-SSF-manuscript/analysis/probes/data/chain/
#
# If they are absent the job silently starts from t = 0 and recomputes all 104 d.
# If the f_dry decision changes the influent, change the TAG too (a new chain),
# drop the rsync step, and expect the full-length wall time.
#
# Wall estimate (E4: ~105 s per simulated day locally at N = 500; cosmos ~1.8x
# slower, E6): 44 d resumed ~ 2.3 h per arm; 104 d from scratch ~ 5.5 h. The
# 12 h limit covers either. Each leg writes its own file, so a wall-clock kill
# loses at most one leg; resubmit the same script to continue.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 12:00:00
#SBATCH -J fld2x104
#SBATCH -o slurm/logs/fld2x104_%A_%a.out
#SBATCH -e slurm/logs/fld2x104_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=fld2x_lit;  LIGHT=1 ;;
    1) TAG=fld2x_dark; LIGHT=0 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== fld2x_104d tag=$TAG task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_${TAG}_leg*.mat 2>/dev/null || echo "no legs on disk -- running from t = 0"

module load matlab/2025b

# E4 call (EXPERIMENTS.md, 2026-08-26) with Days 60 -> 104. Do not edit other options.
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=104, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=$LIGHT, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
               DetachForm=\"linear\", \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
