#!/bin/bash -l
# LOAD LADDER at the working set (E4 parameters), N = 500, 60 d, lit + dark.
#
#   sbatch slurm/ladder_n500.sbatch     # array 0..3
#
# Four arms, identical to E4 (EXPERIMENTS.md) in every option EXCEPT the
# influent heterotroph/phototroph concentrations:
#
#   0  fld1x_lit    HET 1.5e-4  PHO 5.0e-4   (1x field, Campos2006b Fig 1a)   light 1
#   1  fld1x_dark   idem                                                      light 0
#   2  fld4x_lit    HET 6.0e-4  PHO 2.0e-3   (4x field)                       light 1
#   3  fld4x_dark   idem                                                      light 0
#
# The 2x rung of the ladder is E4/E7 itself (tags fld2x_lit / fld2x_dark), so it
# is deliberately NOT repeated here.
#
# RESUME semantics: probeChain skips legs already on disk and guards every option
# against the leg-1 record. `chain_fld1x_{lit,dark}_leg1.mat` (0-10 d) already
# exist from E1c/E1d with exactly these options, so the fld1x arms resume at
# t = 10 d and compute legs 2-6; the fld4x arms are new chains and run 1-6.
# If any option below drifts from the E1 record, the arm dies immediately with
# probeChain:parameterMismatch -- check the .out before assuming it is running.
#
# Wall estimate: E4 cost ~105 s per simulated day locally at N = 500 and cosmos
# is ~1.8x slower (E6), i.e. ~190 s/d => ~32 min per 10 d leg, ~3.2 h for the
# fld4x arms (6 legs) and ~2.6 h for the fld1x arms (5 legs). Higher load costs
# more per day, hence the 12 h limit.
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 12:00:00
#SBATCH -J ladder
#SBATCH -o slurm/logs/ladder_%A_%a.out
#SBATCH -e slurm/logs/ladder_%A_%a.err
#SBATCH --array=0-3

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=fld1x_lit;  LIGHT=1; HET=1.5e-4; PHO=5.0e-4 ;;
    1) TAG=fld1x_dark; LIGHT=0; HET=1.5e-4; PHO=5.0e-4 ;;
    2) TAG=fld4x_lit;  LIGHT=1; HET=6.0e-4; PHO=2.0e-3 ;;
    3) TAG=fld4x_dark; LIGHT=0; HET=6.0e-4; PHO=2.0e-3 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== ladder_n500 tag=$TAG light=$LIGHT HET=$HET PHO=$PHO task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_${TAG}_leg*.mat 2>/dev/null || echo "no legs on disk -- running from t = 0"

module load matlab/2025b

# E4 call verbatim apart from Influent(1:2) and LightScale.
matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=60, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=$LIGHT, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
               DetachForm=\"linear\", \
               Influent=[$HET, $PHO, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
