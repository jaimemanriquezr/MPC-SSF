#!/bin/bash -l
# 30 d lit/dark chain under the TENORE PRESETS (commit b1fcd03: mu_HET 2.0 -> 4.8,
# nu_P 0.094 -> 52), to see how the two changes alter the evolution of the filter.
#
#   sbatch slurm/tenore_chain.sbatch        # array 0 = lit, 1 = dark
#
# WHY TWO ARMS. The dark arm has LightScale = 0, so no light reaches the column and
# nu_P is inert: any difference from the old dark chain is mu_HET alone. The lit arm
# carries both changes. Differencing the pair attributes the response.
#
# NEW TAG ON PURPOSE. probeChain RESUMES from legs already on disk, so reusing
# fld2x_lit/fld2x_dark would silently continue the old-preset chain and produce a
# mongrel. `tenore_*` has no legs, so both arms start from t = 0.
#
# COMPARE AGAINST: chain_fld2x_{lit,dark}_leg1..3 (same call, old presets, 2026-08-26).
# The call below is E4's verbatim (EXPERIMENTS.md / slurm/fld2x_104d.sbatch) except
# Days 104 -> 30 and the tag, so the presets are the ONLY intended difference.
#
# CAVEAT, to state in any report: the old chain was grown on 2026-08-26. Any src/
# change since then that is not a preset is a confound. The only preset change since
# is b1fcd03; a diff of src/ over that window should be checked before attributing.
#
# Wall estimate: ~105 s per simulated day at N = 500 locally, cosmos ~1.8x slower
# (E6) -> ~95 min for 30 d per arm. 4 h covers an arm that clogs and drives the
# adaptive step down. Each leg writes its own file, so a kill loses at most one leg.
#
# Output: analysis/probes/data/chain/chain_tenore_{lit,dark}_leg{1,2,3}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -J tenorechain
#SBATCH -o slurm/logs/tenore_chain_%A_%a.out
#SBATCH -e slurm/logs/tenore_chain_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=tenore_lit;  LIGHT=1 ;;
    1) TAG=tenore_dark; LIGHT=0 ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== tenore_chain tag=$TAG task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
echo "--- presets in effect ---"
grep -n "attenuationParticle = \|NominalRate=4.80" src/presets/modelLund.m || true
ls -l analysis/probes/data/chain/chain_${TAG}_leg*.mat 2>/dev/null || echo "no legs for $TAG -- starting from t = 0 (expected)"

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=30, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=$LIGHT, KDOM=3e-4, TransferScale=10, KHPO4=1e-6, \
               DetachForm=\"linear\", \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, 5.0e-6, 1.0e-3])"

echo "=== done $(date -Is) ==="
