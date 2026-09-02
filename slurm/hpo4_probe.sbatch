#!/bin/bash -l
# Is the winter sub-surface biofilm bulge driven by PHOSPHATE limitation?
#
#   sbatch slurm/hpo4_probe.sbatch          # array 0 = HPO4 influent 0, 1 = 5.0e-6 control
#
# WHY. At 90 d the winter profile is non-monotonic: phi_b dips just under the
# surface and peaks at ~8 cm. Diagnosis (2026-09-02, from chain_fld2x_winter_leg9):
# HPO4 is the LIMITING substrate for heterotroph growth through the top 8 cm, and
# its Monod factor RISES with depth (0.241 at 2 mm -> 0.288 at 8 cm) because
# phototrophs, which deplete it, thin out with depth. Below 8 cm DOM takes over as
# limiting and falls. The net specific rate therefore peaks at the crossover, 8 cm,
# which is exactly where the bulge sits. Setting the influent HPO4 to zero removes
# the external phosphate supply (recycling from decay remains) and should move or
# destroy the bulge if that diagnosis is right.
#
# 60 DAYS, NOT 30. The bulge is ABSENT at 30 d (max at z = 4 mm) and present at
# 60 d (0.1808 at z = 88 mm) and 90 d (0.1843 at 84 mm). A 30 d probe would show
# nothing and prove nothing.
#
# OLD PRESETS PINNED. chain_fld2x_winter was grown with mu_HET = 2.0 and
# nu_P = 0.094; today's presets are 4.8 and 52 (b1fcd03). Both are pinned back via
# MuHET= and the new Attenuation= option so that HPO4 is the ONLY difference from
# the reference chain. Arm 1 repeats the reference influent as a self-contained
# control, so the comparison does not depend on trusting the Aug-26 run.
#
# WINTER SETTINGS, read off chain_fld2x_winter_leg9 rather than assumed:
# T = 3 C, LightForm = "winter" (light is 0 at t = 0.75 d where summer gives 0.37),
# LightScale = 1, zeta0 = 1, zeta1 = 0.27, KDOM = 3e-4, KHPO4 = 1e-6,
# TransferScale = 10, linear detachment, N = 500, MaxDt = 5e-5.
#
# Wall: ~190 s per simulated day at N = 500 on cosmos -> ~3.2 h for 60 d. 6 h
# covers an arm that clogs. Each leg writes its own file.
#
# Output: analysis/probes/data/chain/chain_hpo4{zero,base}_winter_leg{1..6}.mat
#
#SBATCH -A lu2026-2-100
#SBATCH -p lu48
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -J hpo4probe
#SBATCH -o slurm/logs/hpo4_probe_%A_%a.out
#SBATCH -e slurm/logs/hpo4_probe_%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
mkdir -p slurm/logs

case "$SLURM_ARRAY_TASK_ID" in
    0) TAG=hpo4zero_winter; HPO4=0.0     ;;
    1) TAG=hpo4base_winter; HPO4=5.0e-6  ;;
    *) echo "unexpected array index $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "=== hpo4_probe tag=$TAG HPO4_in=$HPO4 task=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is) ==="
ls -l analysis/probes/data/chain/chain_${TAG}_leg*.mat 2>/dev/null || echo "no legs for $TAG -- starting from t = 0 (expected)"

module load matlab/2025b

matlab -batch "addpath(genpath('src')); addpath('analysis'); addpath('analysis/probes'); \
    probeChain(\"$TAG\", Zeta0=1, Zeta1=0.27, NCells=500, Days=60, MaxDt=5e-5, LegDays=10, \
               FramesPerLeg=241, LightScale=1, LightForm=\"winter\", Temperature=3, \
               KDOM=3e-4, TransferScale=10, KHPO4=1e-6, DetachForm=\"linear\", \
               MuHET=2.0, Attenuation=0.094, \
               Influent=[3.0e-4, 1.0e-3, 0, 0, 9.10e-3, 6.23e-3, 2.0e-5, $HPO4, 1.0e-3])"

echo "=== done $(date -Is) ==="
